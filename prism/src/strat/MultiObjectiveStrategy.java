package strat;

import java.util.Map;
import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;

import prism.PrismException;
import explicit.STPG;
import explicit.Model;
import explicit.STPGExplicit;

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

    public void getNextMove(int state) throws InvalidStrategyStateException
    {
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

    // G is the stochastic game
    // X is a polyhedron for each good or bad state
    // Y are the polyhedra of stochastic states:
    //     the first index is the good or bad state
    //     the second index is the action taken
    public MultiObjectiveStrategy(STPG G, List<Polyhedron> X, List<List<Polyhedron>> Y)
    {

	// memory is the list of tuples (t, p), where p is in X(t)

	int N = G.getNumStates();

	// compute memory update function

	memoryUpdateFunction = new Map[N][];
	for(int t = 0; t < n; t++) { // for each state (good or bad)
	    memoryUpdateFunction[t] = new Map[G.getNumChoices(t)];
	    Generator_System gsXt = X.get(t).minimized_generators();
	    for (int u = 0; u < G.getNumChoices(t); u++) { // for each stochastic successor
		memoryUpdateFunction[t][u] = new HashMap<Integer,Map<Integer,Double>>();
		Generator_System gsXu = X.get(u).minimized_generators();
		// now for each corner point p for t, need to find a distribution
		// that is, find l, and l coefficients beta_i summing to one such that
		// for good and bad states: sum_i beta_i q_i^u >= p - rewards(t)
		// and for stochastic states: ...
		// for q^i_u in Y(t,u) - need to actually find these
		for (int p = 0; p < gsXt.size(); p++) {
		    memoryUpdateFunction[t][u].put(p, new HashMap<Integer,Double>());
		    if(G.getPlayer(t) == STPGExplicit.PLAYER_1) {
			// find q_i^u and beta_i in Y(t,u)
			
		    }
		}
		
	    }
	}

	// compute next move function
	nextMoveFunction = new Map[N];
	for (int t = 0; t < n; t++) { // for each state (good or bad, but want good only)
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
			Generator_System test_gs = new Generator_System();
			test_gs.add(gsXt.get(p)); // TODO: do p - reward
			Polyhedron test_poly = new C_Polyhedron(test_gs);
			if(Y.get(t).get(p).contains(test_poly)) {
			    nextMoveFunction[t].put(p, u);
			    break; // inner loop
			}
		    }
		}

	    }
	}

    }

    @Override
    public void updateMemory(int action, int state) throws InvalidStrategyStateException
    {
	
    }

}
