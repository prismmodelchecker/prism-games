package strat;

import java.util.Map;
import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;
import java.math.BigInteger;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Map.Entry;

import prism.PrismException;
import explicit.STPG;
import explicit.Distribution;
import explicit.Model;
import explicit.STPGExplicit;
import explicit.MDP;
import explicit.SMG;
import explicit.MDPSimple;
import explicit.SMGModelChecker;
import explicit.PPLSupport;
import explicit.rewards.STPGRewards;
import parser.State;

import org.apache.commons.math3.optim.linear.SimplexSolver;
import org.apache.commons.math3.optim.linear.Relationship;
import org.apache.commons.math3.optim.linear.LinearConstraintSet;
import org.apache.commons.math3.optim.linear.LinearConstraint;
import org.apache.commons.math3.optim.linear.LinearObjectiveFunction;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.MaxIter;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.fraction.BigFraction;
import org.apache.commons.math3.optim.linear.NoFeasibleSolutionException;
import parma_polyhedra_library.*;

public class MultiObjectiveStrategy implements Strategy
{

    // number of states
    //protected int n;

    String info = "No information available.";

    // last state in history
    protected int lastState;
    // last memory element in history (paired with lastState)
    protected int lastCorner;

    // MEMORY UPDATE FUNCTION
    // first dimension is source state: t
    // second dimension is target state: w
    // first key is action: u
    // second key is corner point for source: p
    // value is a distribution over corner poins of target: choose q_i with probability beta_i
    protected Map<Integer,Map<Integer,Map<Integer,Double>>>[][] memoryUpdateFunction;
    // memory size
    int memorySize = 0;


    // INITIAL DISTRIBUTION FUNCTION
    // first dimension is state
    // key is corner point
    // value is probability
    protected Map<Integer,Double>[] initialDistributionFunction;


    // NEXT STATE FUNCTION
    // first dimension is state: t
    // Map is from corner points p to actions u
    protected Map<Integer,Integer>[] nextMoveFunction;


    protected List<double[]>[] LIST_gsX;


    protected Integer sample(Map<Integer,Double> distribution) throws PrismException
    {
	double r = Math.random();
	for(Map.Entry<Integer,Double> kv : distribution.entrySet()) {
	    r -= kv.getValue();
	    if(r <= 0)
		return kv.getKey();
	}
	throw new PrismException("Distribution invalid.");
    }

    @Override
    public void init(int state) throws InvalidStrategyStateException
    {
	lastState = state;
	try {
	    lastCorner = sample(initialDistributionFunction[state]);
	} catch (PrismException e) {
	    throw new InvalidStrategyStateException("Not possible to sample from corrupt distribution.");
	}
    }

    @Override
    public Distribution getNextMove(int state) throws InvalidStrategyStateException
    {
	if (state != lastState) {
	    throw new InvalidStrategyStateException("Stored state doesn't agree with observed state.");
	}
	Distribution result = new Distribution();
	result.add(nextMoveFunction[lastState].get(lastCorner), 1.0);
	return result;
    }

    @Override
    public void reset()
    {
	lastCorner = -1;
	lastState = -1;
    }

    @Override
    public void exportToFile(String file)
    {
	System.out.println("Exporting to file not supported yet.");
    }

    @Override
    public Model buildProduct(Model model) throws PrismException
    {
	if(model instanceof SMG) {
	    return buildProductSMG((SMG) model);
	}
	throw new UnsupportedOperationException("The product building is not supported for this class of strategies and models.");
    }

    private MDP buildProductSMG(SMG G) throws PrismException
    {
	System.out.println("-------------- BUILDING STATE SPACE ----------------");
	List<State> S = G.getStatesList();
	List<State> newS = new ArrayList<State>();

	// initial state
	State s_init = new State(2);
	s_init.setValue(0, "init");
	s_init.setValue(1, 0);
	newS.add(s_init);

	Map<Integer,Map<Integer,State>> oldSandCornerToNewS = new HashMap<Integer,Map<Integer,State>>();
	// create new states
	// need to precompute, as the indices are later needed
	for(int s = 0; s < S.size(); s++) {
	    State state = S.get(s);
	    oldSandCornerToNewS.put(s, new HashMap<Integer,State>());
	    Map<Integer,Map<Integer,Double>> corners = memoryUpdateFunction[s][0].entrySet().iterator().next().getValue();
	    for(Integer corner : corners.keySet()) {
		State new_state = new State(2);
		new_state.setValue(0, state);
		new_state.setValue(1, corner);
		newS.add(new_state);
		oldSandCornerToNewS.get(s).put(corner, new_state);
	    }
	}

	// printing state space
	for(int s = 0; s < newS.size(); s++) {
	    if(s!=0) System.out.printf(", ");
	    System.out.printf("%d: %s", s, newS.get(s));
	}
	System.out.println();

	// the initial distribution
	Distribution d_init = new Distribution();
	for(int p = 0; p < LIST_gsX[0].size(); p++) {
	    State new_s_init = oldSandCornerToNewS.get(0).get(p);
	    if(initialDistributionFunction[0].containsKey(p)) { // otherwise is zero anyway
		d_init.add(newS.indexOf(new_s_init), initialDistributionFunction[0].get(p));
	    }
	}

	// first create a new MDP
	MDPSimple mdp = new MDPSimple(newS.size());
	mdp.setStatesList(newS);
	mdp.addChoice(0, d_init);
	mdp.addInitialState(0);

	System.out.println("------------- BUILDING TRANSITIONS ----------------");

	for (int s = 0; s < S.size(); s++) { // go through all states in S
	    // for each action need to build a new distribution
	    if(G.getPlayer(s) == STPGExplicit.PLAYER_1) { // Player 1 states become stochastic states
		for(int p = 0; p < LIST_gsX[s].size(); p++) {
		    //System.out.printf("getting p %d from s %d\n", p, s);
		    int u = nextMoveFunction[s].get(p); // u is a stochastic state
		    State origin = oldSandCornerToNewS.get(s).get(p);
		    Distribution d = new Distribution();
		    
		    int nsu = G.getNumTransitions(s,u);
		    Iterator<Entry<Integer,Double>> dsu = G.getTransitionsIterator(s,u);
		    for(int w = 0; w < nsu; w++) { // for each successor w of u
			Entry<Integer,Double> e_w = dsu.next();
			int key_w = e_w.getKey();
			double val_w = e_w.getValue();
			Map<Integer,Double> cd = memoryUpdateFunction[s][u].get(key_w).get(p);
			// build a distribution
			for(Entry<Integer,Double> ecd : cd.entrySet()) {
			    Integer nextCorner = ecd.getKey();
			    State destination = oldSandCornerToNewS.get(key_w).get(nextCorner);
			    //System.out.printf("w:%d, corner:%d: %s\n", key_w, nextCorner, destination);
			    Double nextProb = ecd.getValue();
			    d.add(newS.indexOf(destination), nextProb*val_w);
			}
		    }

		    mdp.addChoice(newS.indexOf(origin),d);
		    System.out.printf("%d(%d)--%d-->: %s\n", s, p, u, d.toString());
		}
	    } else if (G.getPlayer(s) == STPGExplicit.PLAYER_2) { // Player 2 retains decisions
		for(int u = 0; u < G.getNumChoices(s); u++) { // for each action (i.e. stochastic state)
		    for(int p = 0; p < LIST_gsX[s].size(); p++) {
			State origin = oldSandCornerToNewS.get(s).get(p);
			Distribution d = new Distribution();
			
			int nsu = G.getNumTransitions(s,u);
			Iterator<Entry<Integer,Double>> dsu = G.getTransitionsIterator(s,u);
			for(int w = 0; w < nsu; w++) { // for each successor w of u
			    Entry<Integer,Double> e_w = dsu.next();
			    int key_w = e_w.getKey();
			    double val_w = e_w.getValue();
			    Map<Integer,Double> cd = memoryUpdateFunction[s][u].get(key_w).get(p);
			    // build a distribution
			    for(Entry<Integer,Double> ecd : cd.entrySet()) {
				Integer nextCorner = ecd.getKey();
				State destination = oldSandCornerToNewS.get(key_w).get(nextCorner);
				//System.out.printf("w:%d, corner:%d: %s\n", key_w, nextCorner, destination);
				Double nextProb = ecd.getValue();
			    d.add(newS.indexOf(destination), nextProb*val_w);
			    }
			}
			
			mdp.addChoice(newS.indexOf(origin),d);
			System.out.printf("%d(%d)--%d-->: %s\n", s, p, u, d.toString());
		    }
		}
	    }
 	}

	System.out.println("------------- RESULTING MDP ----------------");

	System.out.println(mdp.toString());

	return mdp;

    }

    @Override
    public void setInfo(String info)
    {
	this.info = info;
    }

    @Override
    public String getInfo()
    {
	return info;
    }
    
    @Override
    public int getMemorySize()
    {
	return memorySize;
    }
    
    @Override
    public String getType()
    {
	return "rCMQ HR strategy.";
    }
    
    @Override
    public Object getCurrentMemoryElement()
    {
	return new Object[] { lastState, lastCorner };
    }
    
    @Override
    public void setMemory(Object memory) throws InvalidStrategyStateException
    {
	System.out.println("Setting memory not supported.");
    }

    @Override
    public String getStateDescription()
    {
	return null;
    }

    @Override
    public int getInitialStateOfTheProduct(int s)
    {
	return 0;
    }


    private <T> List<List<T>> selectGenerators(List<T> gs, List<T> tuple, int l, List<T> not)
    {
	if(not == null) {
	    not = new ArrayList<T>();
	}

	List<List<T>> output = new ArrayList<List<T>>();
	if(l==0){
	    if(tuple != null) {
		output.add(tuple);
	    }
	} else {
	    outer:
	    for(T g : gs) {
		if(not != null) {
		    for(T n : not) {
			if(g == n) { // ref comparison should be fine here
			    continue outer;
			}
		    }
		}
		List<T> new_tuple = new ArrayList<T>();
		if(tuple != null) {
		    new_tuple.addAll(tuple);
		}
		new_tuple.add(g);
		not.add(g);
		List<T> new_not = new ArrayList<T>();
		new_not.addAll(not);
		output.addAll(selectGenerators(gs, new_tuple, l-1, new_not));
	    }
	}
	
	return output;
    }




    // select all combinations of q^u_i from a list of length l of tuples
    private <T> List<List<List<T>>> selectMultiGenerators(List<List<List<T>>> tuples, int successor_number, List<List<T>> multiTuple)
    {
	List<List<List<T>>> output = new ArrayList<List<List<T>>>();

	if(successor_number >= tuples.size()) { // base case
	    // generate new list of multituples
	    // put in currently computed tuple - if it exists
	    if(multiTuple != null) {
		output.add(multiTuple);
	    }
	} else {
	    for(int i = 0; i < tuples.get(successor_number).size(); i++) {
		List<List<T>> new_multiTuple = new ArrayList<List<T>>();
		if(multiTuple != null) {
		    new_multiTuple.addAll(multiTuple);
		}
		new_multiTuple.add(tuples.get(successor_number).get(i));
		output.addAll(selectMultiGenerators(tuples, successor_number+1, new_multiTuple));
	    }
	}

	return output;

    }


    private List<double[]> gsToList(Generator_System gs, int dim)
    {
	List<double[]> result = new ArrayList<double[]>(gs.size());
	for(Generator g : gs) {
	    double[] p = new double[dim];
	    Linear_Expression le = g.linear_expression();
	    BigInteger d = g.divisor().getBigInteger();
	    Map<Variable, BigInteger> map = new HashMap<Variable, BigInteger>();
	    PPLSupport.getCoefficientsFromLinearExpression(le, false, BigInteger.ONE, map);
	    for(Variable k : map.keySet()) {
		if(k!=null) {
		    BigFraction value = new BigFraction(map.get(k), d);
		    p[k.id()] = value.doubleValue();
		}
	    }
	    result.add(p);
	}
	return result;
    }

    // G is the stochastic game
    // X is a polyhedron for each good or bad state
    // Y are the polyhedra of stochastic states:
    //     the first index is the good or bad state
    //     the second index is the action taken
    public MultiObjectiveStrategy(STPG G, int initial_state, double[] v, Map<Integer,Polyhedron> X, List<List<Polyhedron>> Y, List<STPGRewards> stpgRewards)
    {
	// memory is the list of tuples (t, p), where p is in X(t)

	int N = G.getNumStates(); // number of states in game (excluding stochastic states)
	int L = ((int) X.get(0).space_dimension()); // total number of goals
	int M = L - stpgRewards.size(); // total number of probabilistic goals

	SimplexSolver solver = new SimplexSolver(1.0e-2, 10);
	List<List<double[]>> LIST_tuples;

	System.out.println("-------- CANONICAL ORDER -------------");
	// establish canonical order of corners

	/*
	gsX = new Generator_System[N];
	for(int t = 0; t < N; t++) {
	    gsX[t] = X.get(t).minimized_generators();
	    Polyhedron gsXpoly = new C_Polyhedron(gsX[t]);
	    gsXpoly.add_space_dimensions_and_project(L-gsXpoly.space_dimension());
	    SMGModelChecker.printReachabilityPolyhedron(gsXpoly, ((int)gsXpoly.space_dimension()), t);
	}
	*/
	LIST_gsX = new List[N];
	for(int t = 0; t < N; t++) {
	    LIST_gsX[t] = gsToList(X.get(t).minimized_generators(), L);
	}

	

	System.out.println("--------- INITIAL DISTRIBUTION -------------");

	initialDistributionFunction = new Map[N];
	
	initialDistributionFunction[initial_state] = new HashMap<Integer,Double>();
	// find q_i^u and beta_i in X(t)
	search_for_distribution:
	for(int l = 1; l < L+1; l++) { // first find l
	    // compute all possible combinations of q_i^u
	    LIST_tuples = selectGenerators(LIST_gsX[initial_state], null, l, null);
	    
	    // preparation for LP
	    double[] coeffs_beta = new double[l];
	    for(int i = 0; i < l; i++) {
		coeffs_beta[i] = 1;
	    }
	    // check for each such tuple
	    iteration_through_tuples:
	    //for(List<Generator> tuple : tuples) {
	    for(List<double[]> tuple : LIST_tuples) {
		// now formulate an LP for beta_i
		
		// max_{beta_i} sum_i beta
		// s.t. sum_i beta_i q_i^u >= v - rewards(t)
		//      sum_i beta_i <= 1
		
		// describe the optimization problem
		// optimization function - maximize betas
		LinearObjectiveFunction f = new LinearObjectiveFunction(coeffs_beta, 0);
		
		// constraints
		List<LinearConstraint> constraints = new ArrayList<LinearConstraint>();
		double[][] coeffs_q = new double[L][l];
		for(int i = 0; i < l; i++) {
		    // get coefficients from tuple.get(i)
		    for(int k = 0; k < L; k++) {
			coeffs_q[k][i] = tuple.get(i)[k];
		    }
		    
		}
		
		// TODO: put value of v - reward(t) into x
		double[] bounds = new double[L];
		for(int i = 0; i < L; i++) {
		    // lower bound on sum of betas
		    bounds[i] = v[i];
		    constraints.add(new LinearConstraint(coeffs_q[i], Relationship.GEQ, bounds[i]));
		}
		for(int i = 0; i < l; i++) {
		    // lower bound on beta_i
		    double[] onlyone = new double[l];
		    onlyone[i] = 1.0;
		    constraints.add(new LinearConstraint(onlyone, Relationship.GEQ, 0.0));
		}
		// upper bound on sum of beta
		constraints.add(new LinearConstraint(coeffs_beta, Relationship.LEQ, 1));
		
		PointValuePair solution;
		try{
		    solution = solver.optimize(f,
					       new LinearConstraintSet(constraints),
					       GoalType.MAXIMIZE,
					       new MaxIter(10000));
		} catch ( NoFeasibleSolutionException e) {
		    // tuple not feasible, try a different one
		    continue iteration_through_tuples;
		}
		
		// there has been no exception, so the problem was fasible
		// can extract the distribution now from the solution
		for(int i = 0; i < l; i++) {
		    initialDistributionFunction[initial_state].put(LIST_gsX[initial_state].indexOf(tuple.get(i)),solution.getPoint()[i]);
		    //System.out.printf("dim%d: %.6f\n", i, solution.getPoint()[i]);
		}
		break search_for_distribution;
	    }
	}
    
	System.out.printf("State %d: %s\n", initial_state, initialDistributionFunction[initial_state].toString());
	

	// compute memory update function
	System.out.println("--------- MEMORY UPDATE AND NEXT MOVE -------------");

	memoryUpdateFunction = new Map[N][]; // for N states
	nextMoveFunction = new Map[N];

	for(int t = 0; t < N; t++) { // for each state (good or bad)
	    // do some garbage collection of the previous iteration
	    System.gc();

	    //System.out.printf("t: %d\n", t);
	    memoryUpdateFunction[t] = new Map[G.getNumChoices(t)];
	    // compute memory size
	    memorySize += LIST_gsX[t].size();

	    if(G.getPlayer(t) == STPGExplicit.PLAYER_1) { // next move only defined for good guy
		nextMoveFunction[t] = new HashMap<Integer,Integer>(LIST_gsX[t].size());
	    }

	    for (int u = 0; u < G.getNumChoices(t); u++) { // for each stochastic successor (i.e. action u)
		//System.out.printf("u: %d\n", u);
		//tuGenerator_System gsYtu = Y.get(t).get(u).minimized_generators();
		List<double[]> gsYtu = gsToList(Y.get(t).get(u).minimized_generators(), L);
		//Polyhedron gsYpoly = new C_Polyhedron(gsYtu);
		//gsYpoly.add_space_dimensions_and_project(L-gsYpoly.space_dimension());
		//System.out.printf("%d-%d->: ", t, u);
		//SMGModelChecker.printReachabilityPolyhedron(gsYpoly, ((int)gsYpoly.space_dimension()), t);
		int ntu = G.getNumTransitions(t,u);
		Iterator<Entry<Integer,Double>> dtu = G.getTransitionsIterator(t,u);
		memoryUpdateFunction[t][u] = new HashMap<Integer,Map<Integer,Map<Integer,Double>>>();
		for(int w = 0; w < ntu; w++) {
		    int key_w = dtu.next().getKey();
		    memoryUpdateFunction[t][u].put(key_w, new HashMap<Integer,Map<Integer,Double>>());
		}

		// temporarily store results for stochastic states
		// specific to u
		Map<Integer,Map<Integer, Double>>[] stochastic = new Map[ntu];


		///// STOCHASTIC //////

		// interpret u as a stochastic state, and look at all its successors w
		// first get tuples for each successor

		//System.out.printf("Stochastic %d,%d corners: %d\n", t, u, gsYtu.size());

		List<List<List<double[]>>> LIST_succ_tuples = new ArrayList<List<List<double[]>>>();
		dtu = G.getTransitionsIterator(t,u);
		for(int w = 0; w < ntu; w++) {
		    int key_w = dtu.next().getKey();
		    stochastic[w] = new HashMap<Integer, Map<Integer,Double>>(); // initialize for each successor w of u
		    //System.out.printf("Generate tuples out of %s\n", gsXws[w].toString());
		    List<List<double[]>> LIST_succ_tuple;
		    look_for_nonempty_tuple:
		    for(int l = L; l >= 1; l--) {
			LIST_succ_tuple = selectGenerators(LIST_gsX[key_w], null, l, null);
			//System.out.printf("l: %d, selected st: %s\n", l, LIST_succ_tuple);
			if(LIST_succ_tuple.size()!=0) {
			    if(l < L) {
				for(List<double[]> point : LIST_succ_tuple) {
				    for(int ll = point.size(); ll < L; ll++) {
					point.add(point.get(0));
				    }
				}
			    }
			    
			    //System.out.printf("L:%d, succ_tuple: %s\n", L, succ_tuple);
			    LIST_succ_tuples.add(LIST_succ_tuple);
			    break look_for_nonempty_tuple;
			}

			
		    }
		}

		//System.out.printf("SUCCTUPLES: %s\n", succ_tuples.toString());

		List<List<List<double[]>>> LIST_multiTuples = selectMultiGenerators(LIST_succ_tuples, 0, null);

		// here choose the distributions of the stochastic states
		for (int p = 0; p < gsYtu.size(); p++) { // for each corner point in u
	    
		    double[] bounds = new double[L];
		    /*
		    Linear_Expression le_p = gsYtu.get(p).linear_expression();
		    Coefficient d_p = gsYtu.get(p).divisor();
		    Map<Variable, BigInteger> map_p = new HashMap<Variable, BigInteger>();
		    PPLSupport.getCoefficientsFromLinearExpression(le_p, false, BigInteger.ONE, map_p);
		    for(Variable k : map_p.keySet()) {
			if(k !=null) {
			    BigFraction c_p = new BigFraction(map_p.get(k), d_p.getBigInteger());
			    bounds[k.id()] = c_p.doubleValue();
			}
		    }
		    */
		    for(int k = 0; k < L; k++) {
			bounds[k] = gsYtu.get(p)[k];
		    }

		    // do p - reward
		    for(int k = M; k < L; k++) {
			bounds[k] -= stpgRewards.get(k-M).getStateReward(t);
		    }
		    //System.out.printf("t:%d, u:%d, p:%d, bounds: %s\n", t, u, p, Arrays.toString(bounds));
		    //System.out.printf("%d multituples for %s\n", multiTuples.size(), gsYtu.get(p));

		    double[] coeffs_beta = new double[ntu*L];
		    double[][] coeffs_beta_indiv = new double[ntu][L*ntu];
		    for(int i = 0; i < L; i++) {
			for(int w = 0; w < ntu; w++){
			    coeffs_beta_indiv[w][w*L+i] = 1;
			    coeffs_beta[w*L+i] = 1;
			}
		    }
		    List<LinearConstraint> constraints = new ArrayList<LinearConstraint>();

		    boolean nothingfound = true;
		    iteration_through_multi_tuples:
		    for(List<List<double[]>> LIST_multiTuple : LIST_multiTuples) { // for each combination of tuples
			// formulate an LP that contains the following constraints
			// sum_{w} /\(u,w) sum_i beta^w_i q^w_i
			
			// the objective function is sum_w sum_i beta^w_i
			
			// The LP then optimizes for all successors simultaneously,
			// so get q^w_i and beta^w_i for each successor w
			
			// build objective function
			// note: For now take L points in successor and don't try to optimize yet.
			//       Would get a lot of optimization problems to actually calculate the
			//       smallest number of points necessary.

			LinearObjectiveFunction f = new LinearObjectiveFunction(coeffs_beta, 0);
		    
			// constraints
			constraints.clear();
			// now that all combinations of tuples are computed, can build the constraints
			// first dimension is constraint
			// second dimension is beta^w_i index
			double[][] coeffs_q = new double[L][ntu*L];
			dtu = G.getTransitionsIterator(t,u);
			for(int w = 0; w < ntu; w++) { // for each successor w
			    double delta_uw = dtu.next().getValue();
			    for(int i = 0; i < L; i++) { // for each component
				// get coefficients from tuple.get(i)
				//System.out.printf("w: %d, i: %d, mt: %s\n", w, i, multiTuple);
				for(int k = 0; k < L; k++) {
				    coeffs_q[k][w*L+i] = delta_uw * LIST_multiTuple.get(w).get(i)[k];
				}
			    }
			}

			/*
			System.out.printf("t: %d, u:%d, p:%d, m-tuple: %d ... coefficients\n", t, u, p, LIST_multiTuples.indexOf(LIST_multiTuple));
			System.out.println(Arrays.deepToString(coeffs_q));
			System.out.println("...");
			*/

			for(int i = 0; i < L; i++) {
			    // lower bound on sum of betas
			    constraints.add(new LinearConstraint(coeffs_q[i], Relationship.GEQ, bounds[i]));
			}
			for(int i = 0; i < L*ntu; i++) {
			    // lower bound on beta^w_i
			    double[] onlyone = new double[L*ntu];
			    onlyone[i] = 1.0;
			    constraints.add(new LinearConstraint(onlyone, Relationship.GEQ, 0.0));
			}
			// upper bound on sums of betas
			for(int w = 0; w < ntu; w++) { // for each successor w
			    constraints.add(new LinearConstraint(coeffs_beta_indiv[w], Relationship.LEQ, 1));
			}
			

			PointValuePair solution;
			try{
			    solution = solver.optimize(f,
						       new LinearConstraintSet(constraints),
						       GoalType.MAXIMIZE,
						       new MaxIter(10000));
			} catch ( NoFeasibleSolutionException e) {
			    // tuple not feasible, try a different one
			    continue iteration_through_multi_tuples;
			}
			nothingfound = false;
			
			//System.out.printf("t:%d, u:%d, bounds:%s\n",t, u, Arrays.toString(bounds));
			//System.out.printf("t:%d, u:%d, coeffs_q:%s\n", t, u, Arrays.deepToString(coeffs_q));
			//System.out.printf("multituple: %s\n", LIST_multiTuple);

			
			//System.out.printf("state: %d, L: %d, ntu: %d, solnsize: %d\n", t, L, ntu, solution.getPoint().length);

			// there has been no exception, so the problem wasa feasible
			// can extract the distribution now from the solution
			dtu = G.getTransitionsIterator(t,u);
			for(int w = 0; w < ntu; w++) { // for each successor
			    stochastic[w].put(p, new HashMap<Integer, Double>()); // initialize for each corner p of u
			    int key_w = dtu.next().getKey();
			    for(int i = 0; i < L; i++) { // for each dimension
				Integer index = LIST_gsX[key_w].indexOf(LIST_multiTuple.get(w).get(i));
				//System.out.printf("key_w: %d, w: %d, i:%d: ", key_w, w, i);
				//System.out.println(LIST_multiTuple.get(w).get(i));
				//System.out.println(LIST_multiTuple.get(w));
				//System.out.println(LIST_gsX[key_w]);
				double prob = solution.getPoint()[L*w+i];
				prob = prob > 1.0 ? 1.0 : prob;
				stochastic[w].get(p).put(index, prob);
				//System.out.printf("%d:%d: %.6f\n", u, index, solution.getPoint()[L*w+i]);
			    }
			}
			break iteration_through_multi_tuples;
		    }

		    //System.out.printf("p:%d: %s\n", p, nothingfound);
		}
		


		/*
		System.out.printf("STOCHASTIC t:%d, u:%d\n", t,u);
		for(int w = 0; w < ntu; w++) {
	            System.out.printf("w:%d: %s\n", w, stochastic[w]);
		}
		*/

		///// STOCHASTIC END ////


		// now for each corner point p for t, need to find a distribution
		// that is, find l, and l coefficients beta_i summing to one such that
		// for good and bad states: sum_i beta_i q_i^u >= p - rewards(t)
		// and for stochastic states: ...
		// for q^i_u in Y(t,u) - need to actually find these
		for (int p = 0; p < LIST_gsX[t].size(); p++) { // for each corner point in t
		    //System.out.printf("p: %d\n", p);
		    dtu = G.getTransitionsIterator(t,u);
		    for(int w = 0; w < ntu; w++) {
			int key_w = dtu.next().getKey();
			memoryUpdateFunction[t][u].get(key_w).put(p, new HashMap<Integer, Double>());
		    }

		    // put value of p - reward(t) into bounds
		    double[] bounds = new double[L];
		    for(int k = 0; k < L; k++) {
			bounds[k] = LIST_gsX[t].get(p)[k];
		    }
		    // do p - reward
		    for(int k = M; k < L; k++) {
			bounds[k] -= stpgRewards.get(k-M).getStateReward(t);
		    }
		    //System.out.printf("t:%d, u:%d, p:%d, bounds: %s\n", t, u, p, Arrays.toString(bounds));

		    
		    // find q_i^u and beta_i in Y(t,u)
		    search_for_distribution:
		    for(int l = 1; l < L+1; l++) { // first find l
			//System.out.printf("l: %d\tout of %d\n", l, L);
			// compute all possible combinations of q_i^u
			LIST_tuples = selectGenerators(gsYtu, null, l, null);
			
			// preparation for LP
			double[] coeffs_beta = new double[l];
			for(int i = 0; i < l; i++) {
			    coeffs_beta[i] = 1;
			}
			// check for each such tuple
			iteration_through_tuples:
			for(List<double[]> LIST_tuple : LIST_tuples) {
			    //System.out.printf("tuple: %s\n", tuple.toString());
			    // now formulate an LP for beta_i
			    
			    // max_{beta_i} sum_i beta
			    // s.t. sum_i beta_i q_i^u >= p - rewards(t)
			    //      sum_i beta_i <= 1
			    
			    // describe the optimization problem
			    // optimization function - maximize betas
			    LinearObjectiveFunction f = new LinearObjectiveFunction(coeffs_beta, 0);
			    
			    // constraints
			    List<LinearConstraint> constraints = new ArrayList<LinearConstraint>();
			    double[][] coeffs_q = new double[L][l];
			    for(int i = 0; i < l; i++) {
				// get coefficients from tuple.get(i)
				for(int k = 0; k < L; k++) {
				    coeffs_q[k][i] = LIST_tuple.get(i)[k];
				}
			    }
			    //System.out.println(Arrays.deepToString(coeffs_q));
			    
			    for(int i = 0; i < L; i++) {
				// lower bound on sum of betas
				constraints.add(new LinearConstraint(coeffs_q[i], Relationship.GEQ, bounds[i]));
			    }
			    for(int i = 0; i < l; i++) {
				// lower bound on beta_i
				double[] onlyone = new double[l];
				onlyone[i] = 1.0;
				constraints.add(new LinearConstraint(onlyone, Relationship.GEQ, 0.0));
			    }
			    // upper bound on sum of beta
			    constraints.add(new LinearConstraint(coeffs_beta, Relationship.LEQ, 1));
			    
			    PointValuePair solution;
			    try{
				solution = solver.optimize(f,
							   new LinearConstraintSet(constraints),
							   GoalType.MAXIMIZE,
							   new MaxIter(10000));
			    } catch ( NoFeasibleSolutionException e) {
				// tuple not feasible, try a different one
				continue iteration_through_tuples;
			    }
			    
			    // there has been no exception, so the problem was fasible
			    // can extract the distribution now from the solution
			    dtu = G.getTransitionsIterator(t,u);
			    for(int w = 0; w < ntu; w++) {
				int key_w = dtu.next().getKey();
				for(int i = 0; i < l; i++) {
				    Integer q_index = gsYtu.indexOf(LIST_tuple.get(i));
				    Map<Integer,Double> action = stochastic[w].get(q_index);
				    //System.out.printf("PLAYER %d at state %d, action %d, target %d(%dth)\n", G.getPlayer(t), t, u, key_w, w);
				    //System.out.printf("action: %s\n", action);
				    //System.out.printf("q_index: %d\n", q_index);
				    if(G.getPlayer(t) != STPGExplicit.PLAYER_1 || action != null) {
					double beta = solution.getPoint()[i]; // beta^u_i
					for(Integer j : action.keySet()) {
					    double prob = beta*action.get(j);
					    prob = prob > 1.0 ? 1.0 : prob;
					    if(prob!=0.0) {
						memoryUpdateFunction[t][u].get(key_w).get(p).put(j, prob);
					    }
					}
				    }
				}
			    }
			    // moreover, if good guy state, know that now u can be picked in t for p
			    if(G.getPlayer(t) == STPGExplicit.PLAYER_1) { // next move only defined for good guy
				nextMoveFunction[t].put(p, u);
			    }

			    break search_for_distribution;
			    
			}
		    }
		}
		// print results for state t and action u

		Iterator<Entry<Integer,Double>> pdtu = G.getTransitionsIterator(t,u);
		for(int w = 0; w < G.getNumTransitions(t, u); w++) {
		    int key_w = pdtu.next().getKey();
		    System.out.printf("Transition %d-%d->%d: %s\n", t, u, key_w, memoryUpdateFunction[t][u].get(key_w).toString());
		}
	    }
	}

	// print results
	for(int t = 0; t < N; t++) {
	    for(int u = 0; u < G.getNumChoices(t); u++) {
		
	    }
	}

	System.out.println("--------- NEXT MOVE -------------");
	for(int t = 0; t < nextMoveFunction.length; t++){
	    if(nextMoveFunction[t] != null) {
		System.out.printf("State %d: %s\n", t, nextMoveFunction[t].toString());
	    }
	}

	System.out.println("--------- STRATEGY DONE -------------");
    }


    @Override
    public void updateMemory(int action, int state) throws InvalidStrategyStateException
    {
	// distribution between next corner is: 
	Map<Integer, Double> distribution = memoryUpdateFunction[lastState][action].get(state).get(lastCorner);
	try {
	    lastCorner = sample(distribution);
	} catch (PrismException e) {
	    throw new InvalidStrategyStateException("Not possible to sample from corrupt distribution.");
	}
	    
    }

}
