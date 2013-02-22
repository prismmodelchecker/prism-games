package strat;

import java.util.Map;
import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;
import java.math.BigInteger;
import java.util.Arrays;

import prism.PrismException;
import explicit.STPG;
import explicit.Distribution;
import explicit.Model;
import explicit.STPGExplicit;
import explicit.SMGModelChecker;
import explicit.PPLSupport;

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
    protected int n;

    // MEMORY UPDATE FUNCTION
    // first dimension is source state: t
    // second dimension is target state: u
    // key is corner point for source: p
    // value is a distribution over corner poins of target: choose q_i with probability beta_i
    Map<Integer,Map<Integer,Double>>[][] memoryUpdateFunction;


    // INITIAL DISTRIBUTION FUNCTION
    // first dimension is state
    // value is probability
    protected double[] initialDistributionFunction;


    // NEXT STATE FUNCTION
    // first dimension is state: t
    // Map is from corner points p to successors u
    protected Map<Integer,Integer>[] nextMoveFunction;

    @Override
    public void init(int state) throws InvalidStrategyStateException
    {
    }

    @Override
    public Distribution getNextMove(int state) throws InvalidStrategyStateException
    {
	return null;
    }

    @Override
    public void reset()
    {
    }

    @Override
    public void exportToFile(String file)
    {
    }

    @Override
    public Model buildProduct(Model model) throws PrismException
    {
	return null;
    }

    @Override
    public void setInfo(String info)
    {
    }

    @Override
    public String getInfo()
    {
	return null;
    }
    
    @Override
    public int getMemorySize()
    {
	return 0;
    }
    
    @Override
    public String getType()
    {
	return "rCMQ HR strategy.";
    }
    
    @Override
    public Object getCurrentMemoryElement()
    {
	return null;
    }
    
    @Override
    public void setMemory(Object memory) throws InvalidStrategyStateException
    {
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

    private List<List<Generator>> selectGenerators(Generator_System gs, List<Generator> tuple, int l, List<Generator> not)
    {
	if(not == null) {
	    not = new ArrayList<Generator>();
	}

	List<List<Generator>> output = new ArrayList<List<Generator>>();
	if(l==0){
	    if(tuple != null) {
		output.add(tuple);
	    }
	    return output;
	} else {
	    outer:
	    for(Generator g : gs) {
		if(not != null) {
		    for(Generator n : not) {
			if(g == n) { // ref comparison should be fine here
			    continue outer;
			}
		    }
		}
		List<Generator> new_tuple = new ArrayList<Generator>();
		if(tuple != null) {
		    new_tuple.addAll(tuple);
		}
		new_tuple.add(g);
		not.add(g);
		List<Generator> new_not = new ArrayList<Generator>();
		new_not.addAll(not);
		output.addAll(selectGenerators(gs, new_tuple, l-1, new_not));
	    }
	    return output;
	}
    }

    // G is the stochastic game
    // X is a polyhedron for each good or bad state
    // Y are the polyhedra of stochastic states:
    //     the first index is the good or bad state
    //     the second index is the action taken
    public MultiObjectiveStrategy(STPG G, List<Polyhedron> X, List<List<Polyhedron>> Y)
    {

	// memory is the list of tuples (t, p), where p is in X(t)

	int N = G.getNumStates();
	int L = ((int) X.get(0).space_dimension());

	Polyhedron extend;
	Generator_System extand_gs = new Generator_System();
	for(int i = 0; i < L; i++) {
	    
	}

	
	
	List<List<Generator>> tuples = selectGenerators(X.get(0).generators(), null, 2, null);
	int ti = 0;
	for(List<Generator> tuple : tuples) {
	    ti++;
	    Generator_System gs = new Generator_System();
	    for(Generator g : tuple) {
		gs.add(g);
	    }
	    Polyhedron P = new C_Polyhedron(gs);
	    SMGModelChecker.printReachabilityPolyhedron(P, ((int)P.space_dimension()), ti);
	}


	// compute memory update function
	System.out.println("--------- MEMORY UPDATE -------------");

	memoryUpdateFunction = new Map[N][];
	for(int t = 0; t < N; t++) { // for each state (good or bad)
	    System.out.printf("t: %d\n", t);
	    memoryUpdateFunction[t] = new Map[G.getNumChoices(t)];
	    Generator_System gsXt = X.get(t).minimized_generators();
	    for (int u = 0; u < G.getNumChoices(t); u++) { // for each stochastic successor
		System.out.printf("u: %d\n", u);
		memoryUpdateFunction[t][u] = new HashMap<Integer,Map<Integer,Double>>();
		Generator_System gsXu = X.get(u).minimized_generators();
		// now for each corner point p for t, need to find a distribution
		// that is, find l, and l coefficients beta_i summing to one such that
		// for good and bad states: sum_i beta_i q_i^u >= p - rewards(t)
		// and for stochastic states: ...
		// for q^i_u in Y(t,u) - need to actually find these
		for (int p = 0; p < gsXt.size(); p++) {
		    System.out.printf("p: %d\n", p);
		    memoryUpdateFunction[t][u].put(p, new HashMap<Integer,Double>());
		    if(G.getPlayer(t) == STPGExplicit.PLAYER_1) {

			// find q_i^u and beta_i in Y(t,u)
			for(int l = 1; l < L+1; l++) { // first find l
			    System.out.printf("l: %d\tout of %d\n", l, L);
			    // compute all possible combinations of q_i^u
			    tuples = selectGenerators(gsXt, null, l, null);

			    // preparation for LP
			    double[] coeffs_beta = new double[l];
			    for(int i = 0; i < l; i++) {
				coeffs_beta[i] = 1;
			    }
			    // check for each such tuple
			    iteration_through_tuples:
			    for(List<Generator> tuple : tuples) {
				System.out.printf("tuple: %s\n", tuple.toString());
				// now formulate an LP feasibility check for beta_i

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
				    Linear_Expression le = tuple.get(i).linear_expression();
				    Coefficient d = tuple.get(i).divisor();
				    Map<Variable, BigInteger> map = new HashMap<Variable, BigInteger>();
				    PPLSupport.getCoefficientsFromLinearExpression(le, false, BigInteger.ONE, map);
				    for(Variable k : map.keySet()) {
					if(k != null) {
					    BigFraction c = new BigFraction(map.get(k), d.getBigInteger());
					    coeffs_q[k.id()][i] = c.doubleValue();
					}
				    }
				}
				System.out.println(Arrays.deepToString(coeffs_q));

				// put value of p - reward(t) into x
				double[] bounds = new double[L];
				Linear_Expression le = gsXt.get(p).linear_expression();
				Coefficient d = gsXt.get(p).divisor();
				Map<Variable, BigInteger> map = new HashMap<Variable, BigInteger>();
				PPLSupport.getCoefficientsFromLinearExpression(le, false, BigInteger.ONE, map);
				for(Variable k : map.keySet()) {
				    if(k !=null) {
					BigFraction c = new BigFraction(map.get(k), d.getBigInteger());
					bounds[k.id()] = c.doubleValue();
				    }
				}
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
				// TODO: can coeffs_beta be reused?
				// upper bound on sum of beta
				constraints.add(new LinearConstraint(coeffs_beta, Relationship.LEQ, 1));

				SimplexSolver solver = new SimplexSolver(1.0e-5, 10);
				PointValuePair solution;
				try{
				    solution = solver.optimize(f,
									      new LinearConstraintSet(constraints),
									      GoalType.MAXIMIZE,
									      new MaxIter(10000));
				} catch ( NoFeasibleSolutionException e) {
				    // tuple not feasible
				    continue iteration_through_tuples;
				}
									  


				// there now should be l values in the solution:
				for(int i = 0; i < l; i++) {
				    System.out.printf("dim%d: %.6f\n", i, solution.getPoint()[i]);
				}
				

				/* EXAMPLE:
				// describe the optimization problem
				LinearObjectiveFunction f = new LinearObjectiveFunction(new double[] { -2, 1 }, -5);
				List<LinearConstraint> constraints = new ArrayList();
				constraints.add(new LinearConstraint(new double[] { 1, 2 }, Relationship.LEQ, 6));
				constraints.add(new LinearConstraint(new double[] { 3, 2 }, Relationship.LEQ, 12));
				constraints.add(new LinearConstraint(new double[] { 0, 1 }, Relationship.GEQ, 0));


				SimplexSolver solver = new SimplexSolver();
				PointValuePair solution = solver.optimize(f, new LinearConstraintSet(constraints), GoalType.MINIMIZE, new MaxIter(10000));
				double[] values = solution.getPointRef();

				// get the solution
				double x = solution.getPoint()[0];
				double y = solution.getPoint()[1];
				double min = solution.getValue();

				System.out.printf("x: %f, y: %f, min: %f", x, y, min);

				*/

			    }
			}
		    }
		}
		
	    }
	}

	// compute next move function
	System.out.println("--------- NEXT MOVE -------------");

	nextMoveFunction = new Map[N];
	for (int t = 0; t < N; t++) { // for each state (good or bad, but want good only)
	    //System.out.printf("Player %d: %d =?= %d\n", t, G.getPlayer(t), STPGExplicit.PLAYER_1);
	    if(G.getPlayer(t) == STPGExplicit.PLAYER_1) { // only defined for good guy
		// need to find l, and l coefficients beta_i summing to one such that
		// sum_i beta_i q_i^u >= p - rewards(t)
		// where q^i_u must be in Y(t,u), but only require existence here, no need to store
		
		// do this by iterating through all corner points p,
		// and for each of them iterate through all successors,
		// and try to find whether p - rewards(t) lies in the polygon Y(t,u)
		// if so, the condition is fulfilled (but don't know which beta_i and q^u_i)
		Generator_System gsXt = X.get(t).minimized_generators();
		nextMoveFunction[t] = new HashMap<Integer,Integer>(gsXt.size());
		for (int p = 0; p < gsXt.size(); p++) { // for each corner point p in X(t)
		    for (int u = 0; u < G.getNumChoices(t); u++) { // for each stochastic successor
			// if Y(t,u) contains p - rewards(t) then
			// construct auxiliary polyhedron containing only X(t)(p) to test containment
			Generator_System test_gs = new Generator_System();
			//System.out.printf("%d, %d -> %d\n", t, p, u);
			test_gs.add(gsXt.get(p)); // TODO: do p - reward
			Polyhedron test_poly = new C_Polyhedron(test_gs);
			// equate dimensions
			test_poly.add_space_dimensions_and_project(Y.get(t).get(u).space_dimension() - test_poly.space_dimension());
			//SMGModelChecker.printReachabilityPolyhedron(test_poly, ((int)test_poly.space_dimension()), t);
			//SMGModelChecker.printReachabilityPolyhedron(Y.get(t).get(u), ((int) Y.get(t).get(u).space_dimension()), u);
			if(Y.get(t).get(u).contains(test_poly)) {
			    nextMoveFunction[t].put(p, u);
			    break;
			}
		    }
		}

	    }
	}


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
	
    }

}
