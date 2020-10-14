//==============================================================================
//	
//	Copyright (c) 2013-
//	Authors:
//	* Clemens Wiltsche <clemens.wiltsche@cs.ox.ac.uk> (University of Oxford)
//	
//------------------------------------------------------------------------------
//	
//	This file is part of PRISM.
//	
//	PRISM is free software; you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation; either version 2 of the License, or
//	(at your option) any later version.
//	
//	PRISM is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//	
//	You should have received a copy of the GNU General Public License
//	along with PRISM; if not, write to the Free Software Foundation,
//	Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//	
//==============================================================================

package strat;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.math.BigInteger;
import java.text.NumberFormat;
import java.util.AbstractMap.SimpleEntry;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Scanner;
import java.util.StringTokenizer;

import org.apache.commons.math3.fraction.BigFraction;
import org.apache.commons.math3.optim.MaxIter;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.linear.LinearConstraint;
import org.apache.commons.math3.optim.linear.LinearConstraintSet;
import org.apache.commons.math3.optim.linear.LinearObjectiveFunction;
import org.apache.commons.math3.optim.linear.NoFeasibleSolutionException;
import org.apache.commons.math3.optim.linear.Relationship;
import org.apache.commons.math3.optim.linear.SimplexSolver;
import org.apache.commons.math3.optim.linear.UnboundedSolutionException;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;

import explicit.Distribution;
import explicit.Model;
import explicit.PPLSupport;
import explicit.Pareto;
import explicit.SMG;
import explicit.rewards.SMGRewards;
import parma_polyhedra_library.Generator;
import parma_polyhedra_library.Generator_System;
import parma_polyhedra_library.Generator_Type;
import parma_polyhedra_library.Linear_Expression;
import parma_polyhedra_library.Variable;
import prism.Prism.StrategyExportType;
import prism.PrismException;
import prism.PrismLangException;
import prism.PrismLog;
import prism.PrismUtils;

public class StochasticUpdateStrategy implements Strategy
{
	// turn on for specific debugging
        private boolean log_problem = false;
        // logging for user
        private boolean logStrategy = false;
        private PrismLog mainLog = null;

	protected String info = "No information available.";

	// initial state the strategy is tailored for
	protected int initial_state;
	// last state in history
	protected int lastState;
	// last memory element in history (paired with lastState)
	protected int lastCorner;

	/**
	 * INITIAL DISTRIBUTION
	 * first index: corner
	 **/
	protected Distribution alpha;

	/**
	 * MEMORY UPDATE FUNCTION: PLAYER STATES
	 * first index: current state
	 * second index: currend corner
	 * third index: next (stochastic) state = next move
	 * returns: distribution over next corner
	 **/
	protected Map<Integer, Map<Integer, Distribution>>[] pi_t; // pi_u(t, p, u)[q]
	/**
	 * MEMORY UPDATE FUNCTION: STOCHASTIC STATES
	 * first index: current state
	 * second index: current stochastic state = current move
	 * third index: currend corner at stochastic state = current corner at move
	 * fourth index: next state
	 * returns: distribution over next corner
	 **/
	protected Map<Integer, Map<Integer, Map<Integer, Distribution>>>[] pi_u; // pi_u((t,u), q, w)[j]

	/**
	 * NEXT STATE FUNCTION
	 * first index: current state
	 * second index: current corner
	 * returns: distribution over next move (next corner is determined by pi_u)
	 **/
	protected Map<Integer, Distribution>[] pi_n; // pi_n(t, p) = u

	// memory size
	protected int memorySize = -1;

	// accuracy
	protected double varepsilon;

	@Override
	public void init(int state) throws InvalidStrategyStateException
	{
		if (state != initial_state) {
			throw new InvalidStrategyStateException(String.format("Strategy tailored to initial state %d, but asked to start from %d", initial_state, state));
		}
		lastState = initial_state;
		try {
			lastCorner = alpha.sampleFromDistribution();
		} catch (PrismException e) {
			throw new InvalidStrategyStateException("Initial distribution invalid. Recompute.");
		}
	}

	@Override
	public Distribution getNextMove(int state) throws InvalidStrategyStateException
	{
	        if (log_problem) {
			System.out.printf("getting next move: %d (last_state=%d, last_corner=%d)\n", state, lastState, lastCorner);
			System.out.printf("pi_n[state]: %s\n", pi_n[state].toString());
		}

		if (state != lastState)
			throw new InvalidStrategyStateException(String.format("Strategy thinks game is at %d, but you ask to proceed from %d", lastState, state));

		if(pi_n.length <= state || pi_n[state] == null)
		    throw new InvalidStrategyStateException(String.format("No choice for state %i specified", state));
		Distribution result = pi_n[state].get(lastCorner);
		return result == null ? new Distribution() : result;
	}

	@Override
	public void updateMemory(int action, int state) throws InvalidStrategyStateException
	{
		if (log_problem)
			System.out.printf("update mem: current state: %d, current mem: %d, action: %d, next state: %d\n", lastState, lastCorner, action, state);
		try {
		        if(pi_t.length <= lastState
			   || pi_t[lastState] == null
			   || pi_t[lastState].get(lastCorner) == null
			   || pi_t[lastState].get(lastCorner).get(action) == null)
			    throw new InvalidStrategyStateException("Cannot proceed to states not selected by the strategy. No stochastic memory update present");

			// first go to stochastic state, according to the action
			int tempCorner = pi_t[lastState].get(lastCorner).get(action).sampleFromDistribution();
		        if(pi_u.length <= lastState
			   || pi_u[lastState] == null
			   || pi_u[lastState].get(action) == null
			   || pi_u[lastState].get(action).get(tempCorner) == null
			   || pi_u[lastState].get(action).get(tempCorner).get(state) == null)
			    throw new InvalidStrategyStateException("Cannot proceed to states not selected by the strategy. No stochastic memory update present");

			// then go to the next state
			lastCorner = pi_u[lastState].get(action).get(tempCorner).get(state).sampleFromDistribution();
			// finally, update the next state
			lastState = state;

		} catch (PrismException e) {
			throw new InvalidStrategyStateException("Something went wrong when sampling from the memory distribution");
		}
	}

	public Distribution memoryUpdate(int action, int state) throws InvalidStrategyStateException
	{
	        if(pi_t.length <= lastState
		   || pi_t[lastState] == null
		   || pi_t[lastState].get(lastCorner) == null
		   || pi_t[lastState].get(lastCorner).get(action) == null)
		    throw new InvalidStrategyStateException("Cannot proceed to states not selected by the strategy. No stochastic memory update present");

		// first go to stochastic state, according to the action
		Distribution result = new Distribution();
		Distribution state_to_action = pi_t[lastState].get(lastCorner).get(action);
		for (Integer tempCorner : state_to_action.getSupport()) { // for each corner at the stochastic state
			double p_tC = state_to_action.get(tempCorner); // probability to go to tempCorner

			if(pi_u.length <= lastState
			   || pi_u[lastState] == null
			   || pi_u[lastState].get(action) == null
			   || pi_u[lastState].get(action).get(tempCorner) == null
			   || pi_u[lastState].get(action).get(tempCorner).get(state) == null)
			    throw new InvalidStrategyStateException("Cannot proceed to states not selected by the strategy. No stochastic memory update present");

			Distribution action_to_state = pi_u[lastState].get(action).get(tempCorner).get(state);
			if (action_to_state != null) {
				for (Integer nextCorner : action_to_state.getSupport()) { // for each corner at the next state
					double p_nC = action_to_state.get(nextCorner); // probability to go to nextCorner
					result.add(nextCorner, p_tC * p_nC);
				}
			}
		}
		return result;
	}

    public String memoryUpdateString(int state, int choice, int next, NumberFormat df) throws InvalidStrategyStateException
    {
		// display probability and memory update
		Distribution dist = getNextMove(state);
		Distribution mu = memoryUpdate(choice, next);
		String label = df.format(dist.get(choice)) + " mu: {";
		boolean first = true;
		for (Integer m : mu.getSupport()) {
			if (first)
				first = false;
			else
				label += ", ";
			label += String.format("(%d, %d)=", next, m) + df.format(mu.get(m));
		}
		label += "}";
		return label;
    }
    
	@Override
	public void reset()
	{
		lastCorner = -1;
		lastState = -1;
	}

	public Distribution parseDistribution(String line)
	{
		Distribution d = new Distribution();
		StringTokenizer st = new StringTokenizer(line.trim(), ",={}");
		while (st.hasMoreTokens()) {
			d.add(Integer.parseInt(st.nextToken().trim()), Double.parseDouble(st.nextToken().trim()));
		}
		return d;
	}

	// construct strategy from file
	public StochasticUpdateStrategy(Scanner scan)
	{
		int states = 0;
		info = ""; // clear startegy info
		String nextLine;
		nextLine = scan.nextLine();
		while (scan.hasNext()) { // parse portions of strategy
			if (nextLine.startsWith("States:")) { // parse number of states
				nextLine = scan.nextLine();
				inner: while (!nextLine.startsWith("InitState:") && scan.hasNext()) {
					if (nextLine.startsWith("//")) {
						nextLine = scan.nextLine();
						continue inner; // ignore comments
					}
					states = Integer.parseInt(nextLine);
					// initialise arrays to hold strategy
					pi_n = new Map[states];
					pi_t = new Map[states];
					pi_u = new Map[states];
					nextLine = scan.nextLine();
				}
			} else if (nextLine.startsWith("InitState:")) { // initial state
				nextLine = scan.nextLine();
				inner: while (!nextLine.startsWith("Init:") && scan.hasNext()) {
					if (nextLine.startsWith("//")) {
						nextLine = scan.nextLine();
						continue inner; // ignore comments
					}
					initial_state = Integer.parseInt(nextLine);
					nextLine = scan.nextLine();
				}
			} else if (nextLine.startsWith("Init:")) { // parse initial distribution
				nextLine = scan.nextLine();
				inner: while (!nextLine.startsWith("Next:") && scan.hasNext()) {
					if (nextLine.startsWith("//")) {
						nextLine = scan.nextLine();
						continue inner; // ignore comments
					}
					alpha = parseDistribution(nextLine);
					nextLine = scan.nextLine();
				}
			} else if (nextLine.startsWith("Next:")) { // parse next state function
				nextLine = scan.nextLine();
				inner: while (!nextLine.startsWith("MemUpdStates:") && scan.hasNext()) {
					if (nextLine.startsWith("//")) {
						nextLine = scan.nextLine();
						continue inner; // ignore comments
					}
					Scanner local = new Scanner(nextLine);
					int s = local.nextInt(); // state
					if (pi_n[s] == null) {
						pi_n[s] = new HashMap<Integer, Distribution>();
					}
					int p = local.nextInt(); // corner
					Distribution d = parseDistribution(local.nextLine()); // distribution
					pi_n[s].put(p, d);
					nextLine = scan.nextLine();
				}
			} else if (nextLine.startsWith("MemUpdStates:")) { // parse memory update: states
				nextLine = scan.nextLine();
				inner: while (!nextLine.startsWith("MemUpdMoves:") && scan.hasNext()) {
					if (nextLine.startsWith("//")) {
						nextLine = scan.nextLine();
						continue inner; // ignore comments
					}
					Scanner local = new Scanner(nextLine);
					int s = local.nextInt(); // state
					if (pi_t[s] == null) {
						pi_t[s] = new HashMap<Integer, Map<Integer, Distribution>>();
					}
					int p = local.nextInt(); // corner
					if (!pi_t[s].containsKey(p)) {
						pi_t[s].put(p, new HashMap<Integer, Distribution>());
					}
					int u = local.nextInt(); // next state
					Distribution d = parseDistribution(local.nextLine()); // distribution
					pi_t[s].get(p).put(u, d);
					nextLine = scan.nextLine();
				}
			} else if (nextLine.startsWith("MemUpdMoves:")) { // parse memory update: moves
				nextLine = scan.nextLine();
				inner: while (!nextLine.startsWith("Info:") && scan.hasNext()) {
					if (nextLine.startsWith("//")) {
						nextLine = scan.nextLine();
						continue inner; // ignore comments
					}
					Scanner local = new Scanner(nextLine);
					int s = local.nextInt(); // state
					if (pi_u[s] == null) {
						pi_u[s] = new HashMap<Integer, Map<Integer, Map<Integer, Distribution>>>();
					}
					int u = local.nextInt(); // move
					if (!pi_u[s].containsKey(u)) {
						pi_u[s].put(u, new HashMap<Integer, Map<Integer, Distribution>>());
					}
					int q = local.nextInt(); // corner at move
					if (!pi_u[s].get(u).containsKey(q)) {
						pi_u[s].get(u).put(q, new HashMap<Integer, Distribution>());
					}
					int w = local.nextInt(); // next state
					Distribution d = parseDistribution(local.nextLine()); // distribution
					pi_u[s].get(u).get(q).put(w, d);
					nextLine = scan.nextLine();
				}
			} else if (nextLine.startsWith("Info:")) { // parse strategy info
				nextLine = scan.nextLine();
				boolean first = true;
				inner: while (!nextLine.startsWith("endstrategy") && scan.hasNext()) {
					if (!first)
						info += "\n";
					first = false;
					if (nextLine.startsWith("//")) {
						nextLine = scan.nextLine();
						continue inner; // ignore comments
					}
					info += nextLine;
					nextLine = scan.nextLine();
				}
			} else if (nextLine.startsWith("endstrategy")) {
				break;
			} else {
				nextLine = scan.nextLine();
			}
		}
	}

	@Override
	public void exportToFile(String filename)
	{
		FileWriter out = null;
		try {
			out = new FileWriter(filename);
			out.write(this.toString());
			out.flush();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (out != null)
				try {
					out.close();
				} catch (IOException e) {
				}
		}
	}

	@Override
	public String toString()
	{

		ByteArrayOutputStream stream = new ByteArrayOutputStream();
		PrintStream out = new PrintStream(stream);

		// header
		out.print(Strategies.FORMAT_STRING_SU_STRAT_MONO + "\n");
		out.print("// Stochastic Memory Update Strategy\n");
		out.print("start strategy\n");
		out.print("States:\n");
		out.print(String.format("%d\n", pi_n.length));
		out.print("// Initial state\n");
		out.print("InitState:\n");
		out.print(String.format("%d\n", initial_state));

		// initial distribution
		out.print("// initial distribution\n");
		out.print("Init:\n");
		out.print(alpha.toString());
		out.print("\n");

		// next state function
		out.print("// next state function\n");
		out.print("// note: only P1 states\n");
		out.print("Next:\n");
		out.print("// first index: current state\n");
		out.print("// second index: current corner\n");
		for (int s = 0; s < pi_n.length; s++) { // go through states
			if (pi_n[s] != null) {
				for (Integer p : pi_n[s].keySet()) { // go through corners
					out.print(String.format("%d %d %s\n", s, p, pi_n[s].get(p)));
				}
			}
		}

		// memory update function: player states
		out.print("// memory update function: player states\n");
		out.print("MemUpdStates:\n");
		out.print("// first index: current state\n");
		out.print("// second index: current corner\n");
		out.print("// third index: next move\n");
		for (int s = 0; s < pi_t.length; s++) { // go through states
			if (pi_t[s] != null) {
				for (Integer p : pi_t[s].keySet()) { // go through current corner
					if (pi_t[s].get(p) != null) {
						for (Integer u : pi_t[s].get(p).keySet()) { // go through next state
							out.print(String.format("%d %d %d %s\n", s, p, u, pi_t[s].get(p).get(u)));
						}
					}
				}
			}
		}

		// memory update function: moves
		out.print("// memory update function: moves\n");
		out.print("MemUpdMoves:\n");
		out.print("// first index: current state\n");
		out.print("// second index: current move\n");
		out.print("// third index: curent corner (at move)\n");
		out.print("// fourth index: next state\n");
		for (int s = 0; s < pi_u.length; s++) { // go through states
			if (pi_u[s] != null) {
				for (Integer u : pi_u[s].keySet()) { // go through moves
					if (pi_u[s].get(u) != null) {
						for (Integer q : pi_u[s].get(u).keySet()) { // go through corners at move
							if (pi_u[s].get(u).get(q) != null) {
								for (Integer w : pi_u[s].get(u).get(q).keySet()) { // go through next state
									out.print(String.format("%d %d %d %d %s\n", s, u, q, w, pi_u[s].get(u).get(q).get(w)));
								}
							}
						}
					}
				}
			}
		}

		// strategy info
		out.print("Info:\n");
		out.print(info);

		// footer
		out.print("\nendstrategy\n");

		out.flush();
		out.close();
		return stream.toString();
	}

	// Note that stochastic memory update doesn't work with product!
	@Override
	public Model buildProduct(Model model) throws PrismException
	{
		if (!model.getClass().equals(SMG.class)) {
			throw new PrismLangException("Unsupported model type");
		}
		throw new PrismLangException("Product with model not supported for stochastic memory update strategy.");
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
		return Strategies.FORMAT_STRING_SU_STRAT_MONO;
	}

	@Override
	public Object getCurrentMemoryElement()
	{
		Entry mem = new SimpleEntry<Integer, Integer>(lastState, lastCorner);
		return mem;
		// note that this can only be of a non-stochastic state,
		// as strategy steps over these
	}

	@Override
	public void setMemory(Object memory) throws InvalidStrategyStateException
	{
		if (memory instanceof Entry) {
			this.lastState = (Integer) ((Entry) memory).getKey();
			this.lastCorner = (Integer) ((Entry) memory).getValue();
		} else {
			throw new InvalidStrategyStateException("Memory has to be integer for this strategy.");
		}
	}

	@Override
	public String getDescription()
	{
		String desc = "";
		desc += "Stochastic update strategy\n";
		desc += "Size of memory: " + getMemorySize() + "\n";
		return desc;
	}

	@Override
	public int getInitialStateOfTheProduct(int s)
	{
		return -1; // not available for SU strategies
	}

	// produces the list of next multi-generators, i.e. one generator for each successor
	private List<List<double[]>> selectMultiGenerator(List<List<double[]>> gss, List<List<double[]>> previous_tuples, int l) throws PrismException
	{
		List<List<double[]>> result = new ArrayList<List<double[]>>();
		// 0. initialise
		if (previous_tuples == null) {
			for (int w = 0; w < gss.size(); w++) { // for each successor
				result.add(selectGenerator(gss.get(w), null, l));
			}
			return result;
		}
		// 1. determine if there is a next generator at the current index - and handly carrys
		int current_multi_index = 0;
		while (current_multi_index < gss.size()) {
			List<double[]> i_tuple = selectGenerator(gss.get(current_multi_index), previous_tuples.get(current_multi_index), l);
			if (i_tuple == null) { // maximum generator reached at current_multi_index
				// reset current_multi_index
				i_tuple = selectGenerator(gss.get(current_multi_index), null, l);
				// put at position in result
				result.add(i_tuple);
				// carry - look at next index
				current_multi_index++;
			} else { // current_multi_index may be increased
				// put at position in result
				result.add(i_tuple);
				// keep the remaining tuples
				for (int w = current_multi_index + 1; w < gss.size(); w++) {
					result.add(previous_tuples.get(w));
				}
				// done
				break;
			}
		}
		if (current_multi_index == gss.size()) { // exited loop because no more multi-tuple present
			return null;
		}
		return result;
	}

	// select a generator that has l components, even if duplicates
	private List<double[]> selectGenerator(List<double[]> gs, List<double[]> previous_tuple, int l) throws PrismException
	{
		for (int ll = l; ll >= 1; ll--) {
			try {
				List<double[]> tuple = selectGeneratorCandidate(gs, previous_tuple, ll);
				if (tuple == null) {
					return null; // no more tuple to be found
				} else if (ll < l) { // if a non-full-sized tuple is selected
					for (int lll = tuple.size(); lll < l; lll++) { // fill up
						tuple.add(tuple.get(0));
					}
					return tuple;
				} else { // full tuple found
					return tuple;
				}
			} catch (PrismException e) { // thrown if too many generators required
				continue; // try fewer generators
			}
		}
		// fall through if nothing found
		return null;
	}

	// produces the next list of generators
	private List<double[]> selectGeneratorCandidate(List<double[]> gs, List<double[]> previous_tuple, int l) throws PrismException
	{
		List<double[]> result = new ArrayList<double[]>();
		// 0. initialise
		if (previous_tuple == null) {
			if (l <= gs.size()) { // initialise by putting the first l generators in the tuple
				for (int i = 0; i < l; i++) {
					result.add(gs.get(i));
				}
				return result;
			} else { // more generators requested than available in gs
				throw new PrismLangException("Too many generators required");
			}
		}

		// 1. get the list of indices used in the tuple
		List<Integer> indices_used = new ArrayList<Integer>(l);
		for (int i = 0; i < gs.size(); i++) {
			if (previous_tuple.contains(gs.get(i))) { // index used in tuple
				indices_used.add(i);
			}
		}
		// 2. see which one is the highest one that can be moved
		int index_to_move = -1;
		for (int i = indices_used.size() - 1; i >= 0; i--) { // iterate from highest to lowest
			// now test if there is a space to move the index up
			if (indices_used.get(i) == gs.size() - 1) { // already at the top, cannot move ... continue
				continue;
			} else if (i == indices_used.size() - 1) { // i is highest index, and not already at the top ... pick this one
				index_to_move = i;
				break;
			} else { // i is not highest, and not alrady at the top
				// check if the index has space to go up
				if (indices_used.get(i) + 1 < indices_used.get(i + 1)) { // if so, pick this one
					index_to_move = i;
					break;
				}
			}
		}
		// 2.5 if none is found, return null
		if (index_to_move == -1) {
			return null;
		}
		// 3. move that index one up
		indices_used.set(index_to_move, indices_used.get(index_to_move) + 1);
		// the remaining indices need to come in succession afterwards (void if highest index)
		int count = 1;
		for (int i = index_to_move + 1; i < indices_used.size(); i++) {
			indices_used.set(i, indices_used.get(index_to_move) + count);
			count++;
		}
		// 4. return the new tuple
		for (Integer i : indices_used) {
			result.add(gs.get(i));
		}
		return result;
	}

	private static Comparator<double[]> COMPARATOR = new Comparator<double[]>()
	{
		public int compare(double[] g1, double[] g2)
		{
			double l1 = 0.0;
			double l2 = 0.0;
			for (int i = 0; i < g1.length; i++) {
				l1 += g1[i] * g1[i];
				l2 += g2[i] * g2[i];
			}
			return Double.compare(l2, l1);
		}
	};

	/**
	 * Turns a generator system into a list
	 * so that the order of points can be stored canonically
	 **/
	private List<double[]> gsToList(Generator_System gs, int dim)
	{
		List<double[]> result = new ArrayList<double[]>(gs.size());
		for (Generator g : gs) {
			if (g.type() == Generator_Type.POINT) {
				double[] p = new double[dim];
				Linear_Expression le = g.linear_expression();
				BigInteger d = g.divisor().getBigInteger();
				Map<Variable, BigInteger> map = new HashMap<Variable, BigInteger>();
				PPLSupport.getCoefficientsFromLinearExpression(le, false, BigInteger.ONE, map);
				for (Variable k : map.keySet()) {
					if (k != null) {
						BigFraction value = new BigFraction(map.get(k), d);
						p[(int) k.id()] = value.doubleValue();
					}
				}
				result.add(p);
			} // ignore anything that is not a point
		}

		// exclude all non-pareto points
		boolean changed = true;
		while (changed) {
			changed = false;
			List<double[]> to_delete = new ArrayList<double[]>();
			for (double[] p : result) {
				look_for_larger_point: for (double[] q : result) { // if there is a point (q) which is larger than p in all dimensions
					if (p != q) {
						for (int i = 0; i < p.length; i++) {
							if (p[i] > q[i]) {
								continue look_for_larger_point;
							}
						}
						// fall through only if q is larger than p in all dimensions
						to_delete.add(p);
						break look_for_larger_point;
					}
				}
			}
			if (to_delete.size() > 0) {
				result.removeAll(to_delete);
				changed = true;
			}
		}

		// now sort the result to put in generators with decreasing length
		Collections.sort(result, COMPARATOR);

		return result;
	}

	/**
	 * if reachable_only is set, only construct reachable memory and state space.
	 * @param G The game.
	 * @param v The objective at the initial state.
	 * @param X Pareto sets at the player states.
	 * @param Y Pareto sets at the stochastic states (moves).
	 * @param rewards The reward structures for which to construct the strategy.
	 * @param biggest_reward For rounding: the biggest reward in the respective dimensions.
	 * @param baseline_accuracy For rounding: the baseline accuracy.
	 * @param reachable_only Construct strategy only for reachable states and memory elements.
	 * @param rounding Use rounding.
	 * @param varepsilon Tolerance for energy objectives (should be zero otherwise).
	 * @param logStrategy turn logging on or off during strategy construction.
	 * @param mainLog logger during the strategy construction.
	 */
	public StochasticUpdateStrategy(SMG G, double[] v, Pareto[] X, List<Pareto>[] Y, List<SMGRewards> rewards, double[] biggest_reward, long baseline_accuracy,
					boolean reachable_only, boolean rounding, double varepsilon, boolean logStrategy, PrismLog mainLog) throws PrismException
	{
	        // set logger
	        this.logStrategy = logStrategy;
		this.mainLog = mainLog;

		if(logStrategy) mainLog.print("Constructing SU Strategy from the following sets:");
		if(logStrategy) PPLSupport.printReachabilityPolyhedra(X, Y, v.length, mainLog);


		int initial_state = G.getFirstInitialState();

		// tolerance - should be zero if not an enery objective 
		this.varepsilon = varepsilon;

		// make initial state specific to strategy
		this.initial_state = initial_state;
		this.memorySize = 0;

		// store size of game and goal
		int gameSize = G.getNumStates(); // game size
		int n = rewards.size(); // number of rewards

		// accuracy
		long[] accuracy = new long[n];
		if (rounding) {
			for (int i = 0; i < n; i++) {
				accuracy[i] = Math.max(1, ((long) (((double) baseline_accuracy) / biggest_reward[i])));
			}
		}

		// need LP solver to get convex combinations
		SimplexSolver solver = null;
		if (rounding)
			solver = new SimplexSolver(1 / ((double) baseline_accuracy));
		else
			solver = new SimplexSolver();

		//----------------------------------------------------------------------
		// CANONICAL ORDER OF CORNER POINTS
		List<double[]>[] LIST_gsX = new List[gameSize];

		// for constructing reachable memory elements only
		// list to store which corners need to be covered
		BitSet[] c_X = reachable_only ? new BitSet[gameSize] : null;
		Map<Integer, BitSet>[] c_Y = reachable_only ? new Map[gameSize] : null;

		// establish canonical order of corners
		for (int t = 0; t < gameSize; t++) {
			LIST_gsX[t] = gsToList(X[t].get().minimized_generators(), n);
			if (reachable_only)
				c_X[t] = new BitSet(LIST_gsX[t].size());
		}

		//----------------------------------------------------------------------
		// INITIAL DISTRIBUTION
		alpha = new Distribution(); // initialize initial distribution to empty.

		// put value of v - reward(t) into bounds for LP
		double[] bounds = new double[n];
		/*
		for(int i = 0; i < n; i++) {
		    // lower bound on sum of betas
		    if(rounding) {
			bounds[i] = Math.floor((v[i] - rewards.get(i).getStateReward(initial_state))*accuracy[i])/accuracy[i];
		    } else {
			bounds[i] = v[i] - rewards.get(i).getStateReward(initial_state);
		    }
		}
		*/
		for (int i = 0; i < n; i++) {
			// lower bound on sum of betas
			if (rounding) {
				bounds[i] = Math.floor((v[i] - varepsilon) * accuracy[i]) / accuracy[i];
			} else {
				bounds[i] = v[i] - varepsilon;
			}
		}

		boolean nothingfound = true;
		search_for_distribution: for (int l = 1; l < n + 1; l++) { // first find l, the number of corners needed in the convex combination
			List<double[]> tuple = null;

			// preparation for LP
			double[] coeffs_beta = new double[l];
			for (int i = 0; i < l; i++) {
				coeffs_beta[i] = 1.0;
			}
			// check for each such tuple
			nothingfound = true;
			iteration_through_tuples: while (true) {
				tuple = selectGenerator(LIST_gsX[initial_state], tuple, l);
				if (tuple == null)
					break iteration_through_tuples;
				//for(List<double[]> tuple : LIST_tuples) {
				// now formulate an LP for beta_i

				// max_{beta_i} sum_i beta
				// s.t. sum_i beta_i q_i^u >= v - rewards(t)
				//      sum_i beta_i <= 1

				// optimization function - maximize betas
				LinearObjectiveFunction f = new LinearObjectiveFunction(coeffs_beta, 0);

				// constraints
				List<LinearConstraint> constraints = new ArrayList<LinearConstraint>();
				double[][] coeffs_q = new double[n][l];
				for (int i = 0; i < l; i++) {
					// get coefficients from tuple.get(i)
					for (int k = 0; k < n; k++) {
						coeffs_q[k][i] = tuple.get(i)[k];
					}
				}

				// put value of v - reward(t) into bounds
				for (int i = 0; i < n; i++) {
					constraints.add(new LinearConstraint(coeffs_q[i], Relationship.GEQ, bounds[i]));
				}
				double[][] onlyone = new double[l][l];
				for (int i = 0; i < l; i++) {
					// lower bound on beta^w_i
					onlyone[i][i] = 1.0;
					constraints.add(new LinearConstraint(onlyone[i], Relationship.GEQ, 0.0));
				}
				// strict bound on sum of beta ... probability distribution
				constraints.add(new LinearConstraint(coeffs_beta, Relationship.EQ, 1));

				PointValuePair solution = null;
				try {
					solution = solver.optimize(f, new LinearConstraintSet(constraints), GoalType.MAXIMIZE, new MaxIter(10000));
				} catch (NoFeasibleSolutionException e) {
					// tuple not feasible, try a different one
					continue iteration_through_tuples;
				}
				nothingfound = false;

				// there has been no exception, so the problem was fasible
				// can extract the distribution now from the solution
				for (int i = 0; i < l; i++) {
					alpha.add(LIST_gsX[initial_state].indexOf(tuple.get(i)), solution.getPoint()[i]);
				}
				break search_for_distribution;
			}
		}

		if (nothingfound) {
			throw new PrismLangException("Goal not realizable");
		}

		alpha = renormDistribution(alpha);

		// mark all corner points that need to be covered
		if (reachable_only) {
			for (Integer p : alpha.getSupport()) {
				c_X[initial_state].set(p);
			}
		}
		boolean corners_to_cover = true;

		//----------------------------------------------------------------------
		// TRANSITIONS AND MEM UPDATE
		pi_n = new Map[G.getNumStates()];
		pi_t = new Map[G.getNumStates()];
		pi_u = new Map[G.getNumStates()];

		while (corners_to_cover) {
			corners_to_cover = false;
			for (int t = 0; t < G.getNumStates(); t++) { // go through all states in S

				if (log_problem)
					printCorners(c_X);

				if (!reachable_only || c_X[t].cardinality() != 0) { // some corners at state t need to be worked on

					// initialize only if necessary, i.e. if t hasn't been visited yet
					if (pi_n[t] == null)
						pi_n[t] = new HashMap<Integer, Distribution>();
					if (pi_t[t] == null)
						pi_t[t] = new HashMap<Integer, Map<Integer, Distribution>>();
					if (pi_u[t] == null)
						pi_u[t] = new HashMap<Integer, Map<Integer, Map<Integer, Distribution>>>();

					// initialize corner tracking for successors of t
					if (reachable_only) {
						if (c_Y[t] == null) {
							c_Y[t] = new HashMap<Integer, BitSet>();
						}
					}

					if (G.getPlayer(t) == 1) {
						corners_to_cover |= getP1Choices(G, t, Y, LIST_gsX, c_X, c_Y, rewards, accuracy, solver, reachable_only, rounding);
					} else {
						corners_to_cover |= getP2Choices(G, t, Y, LIST_gsX, c_X, c_Y, rewards, accuracy, solver, reachable_only, rounding);
					}

				} // check if c_X cardinality
			}
		}
	}

	private void printCorners(BitSet[] c_X)
	{
		System.out.printf("corners to cover: ");
		boolean comma = false;
		for (int ii = 0; ii < c_X.length; ii++) {
			if (comma) {
				System.out.printf(", ");
				comma = false;
			}
			if (c_X[ii].length() > 0) {
				System.out.printf("%d=%s", ii, c_X[ii]);
				comma = true;
			}
		}
		System.out.printf("\n");
	}

	/**
	 * computes next state and memory update functions for choices of P1
	 *
	 * G ... game to construct strategy for
	 * t ... the P2 state
	 * Yt ... Pareto sets at moves of t
	 * LIST_gsX ... canonical order of points at X
	 * c_X ... corners to cover
	 * c_Y ... corners to cover at moves
	 * rewards ... the reward structures of the game under consideration
	 * accuracy ... accuracy per dimension
	 * reachable_only ... only construct strategy with reachable corners and states
	 * rounding ... whether rounding is active
	 *
	 * return value:
	 * corners_to_cover ... whether new corners need to be covered in the reachable_only search
	 **/
	private boolean getP1Choices(SMG G, int t, List<Pareto>[] Y, List<double[]>[] LIST_gsX, BitSet[] c_X, Map<Integer, BitSet>[] c_Y, List<SMGRewards> rewards,
			long[] accuracy, SimplexSolver solver, boolean reachable_only, boolean rounding) throws PrismException
	{
		int n = rewards.size();
		double[] bounds;
		int nt = G.getNumChoices(t);

		boolean corners_to_cover = false; // return value

		// to handle selfloops
		boolean corners_from_selfloops = false;
		BitSet c_Xt = reachable_only ? new BitSet() : null;
		do { // while corners_from_selfloops
			corners_from_selfloops = false;

			// find next move function and memory update at state
			iterate_through_corners: for (int p = 0; p < LIST_gsX[t].size(); p++) { // for each corner point in t

				if (!reachable_only || c_X[t].get(p) || c_Xt.get(p)) { // if corner p needs to be handled

					if (log_problem) System.out.printf("looking for: t:%d p:%d(%s) ----> \n", t, p, Arrays.toString(LIST_gsX[t].get(p)));
					if (logStrategy) mainLog.print(String.format("looking for: t:%d p:%d(%s) ----> \n", t, p, Arrays.toString(LIST_gsX[t].get(p))));

					if (!pi_t[t].containsKey(p)) {
						pi_t[t].put(p, new HashMap<Integer, Distribution>());
					}

					// put value of p - reward(t) -varepsilon into bounds
					bounds = new double[n];
					for (int k = 0; k < n; k++) {
						if (rounding) {
							bounds[k] = Math.floor((LIST_gsX[t].get(p)[k] - rewards.get(k).getStateReward(t) - varepsilon) * accuracy[k]) / accuracy[k];
						} else {
							bounds[k] = LIST_gsX[t].get(p)[k] - rewards.get(k).getStateReward(t) - varepsilon;
						}
					}

					boolean nothingfound = true;
					for (int l = 1; l < n + 1; l++) { // number of corner points allowed in convex combination

						List<LinearConstraint> constraints = new ArrayList<LinearConstraint>();

						List<List<double[]>> gss = new ArrayList<List<double[]>>();
						int nnt = 0; // keeps track of how many successors have viable corners
						for (int u = 0; u < nt; u++) {
							List<double[]> gsYtu = gsToList(Y[t].get(u).get().minimized_generators(), n);
							// only take those points for which at least one component is larger than the bounds (otherwise doesn't help in convex combination)
							List<double[]> keep = new ArrayList<double[]>();
							look_for_tuples_to_keep:
							for(double[] gYtu : gsYtu) {
								for(int k = 0; k < n; k++) {
									if(gYtu[k] >= bounds[k] - PrismUtils.epsilonDouble) {
										keep.add(gYtu);
										continue look_for_tuples_to_keep;
									}
								}
							}
							if(keep.size()>0)
								nnt++;
							gss.add(keep);
						}
						List<List<double[]>> LIST_multiTuple = null;

						double[] coeffs_gamma = new double[nnt * l];
						for (int i = 0; i < l; i++) {
							for (int u = 0; u < nnt; u++) {
								coeffs_gamma[u * l + i] = 1.0;
							}
						}
						
						int multi_counter = 0;
						iteration_through_multi_tuples: while (true) {
							LIST_multiTuple = selectMultiGenerator(gss, LIST_multiTuple, l);
							if (LIST_multiTuple == null)
								break iteration_through_multi_tuples;
							multi_counter++;
							/*
							mainLog.print(String.format("mTs(%d): [", multi_counter));
							for(List<double[]> mT : LIST_multiTuple) {
							        mainLog.print("[");
								if(mT != null) {
									for(double[] T : mT) {
										mainLog.print(String.format("%s", T!=null ? Arrays.toString(T) : "null"));
									}
								}	
								mainLog.print("].");
							}
							mainLog.print("]\n");
							*/
							// formulate a QP:
							// max f(alpha, beta) = sum_u (alpha^u + sum_i beta^u_i)
							// s.t. sum_u alpha^u sum_i beta^u_i q^u_i >= p - reward(t)
							//      sum_u alpha^u = 1
							//      sum_i beta^u_i = 1 for all u

							// The QP then optimizes for all moves simultaneously,
							// so get alpha^u, q^u_i and beta^u_i for each move u

							// transform the QP into an LP, by letting gamma^u_i = beta^u_i alpha^u
							// max f(gamma) = sum_{i, u} gamma^u_i
							// s.t. sum_{u,i} gamma^u_i q^u_i >= p - reward(t)
							//      sum_{u,i} gamma^u_i = 1

							// then alpha^u = sum_i gamma^u_i
							//  and beta^u_i = gamma^u_i / alpha^u

							// build objective function
							LinearObjectiveFunction f = new LinearObjectiveFunction(coeffs_gamma, 0);

							// constraints
							constraints.clear();

							// now that all combinations of tuples are computed, can build the constraints
							// first dimension is constraint
							// second dimension is gamma^u_i index
							double[][] coeffs_q = new double[n][nnt * l];
							for (int i = 0; i < l; i++) { // for each component
								for (int k = 0; k < n; k++) {
									int uu = 0;
									for (int u = 0; u < nt; u++) { // for each move u
										if(LIST_multiTuple.get(u) != null && LIST_multiTuple.get(u).get(i) != null) {
											coeffs_q[k][uu * l + i] = LIST_multiTuple.get(u).get(i)[k];
											uu++;
										}
										//else
											//coeffs_q[k][u * l + i] = -Double.MAX_VALUE; // avoid in LP
									}
								}
							}

							for (int i = 0; i < n; i++) {
								// lower bound on sum of gammas
								constraints.add(new LinearConstraint(coeffs_q[i], Relationship.GEQ, bounds[i]));
							}
							double[][] onlyone = new double[l * nnt][l * nnt];
							for (int i = 0; i < l * nnt; i++) {
								// lower bound on gamma^u_i
								onlyone[i][i] = 1.0;
								constraints.add(new LinearConstraint(onlyone[i], Relationship.GEQ, 0.0));
							}
							// strict bound on sum of gamma ... probability distribution
							constraints.add(new LinearConstraint(coeffs_gamma, Relationship.EQ, 1.0));

							PointValuePair solution = null;
							try {
							    if(log_problem)
							        mainLog.print(String.format("----\ngoal: %s\n\nbounds: %s\n\ncoeffs: %s\n\ncoeffs_indiv: %s\n\nonlyone: %s\n----\n",
											    Arrays.toString(coeffs_gamma),
											    Arrays.toString(bounds),
											    Arrays.deepToString(coeffs_q),
											    "",
											    Arrays.deepToString(onlyone)));
								solution = solver.optimize(f, new LinearConstraintSet(constraints), GoalType.MAXIMIZE, new MaxIter(10000));	
								/*
								mainLog.print(String.format("solution: %s\n----\n",
											    solution==null?"null":Arrays.toString(solution.getPoint())));
								*/
								
							} catch (NoFeasibleSolutionException e) {
								// tuple not feasible, try a different one
								continue iteration_through_multi_tuples;
							} catch (UnboundedSolutionException e) {
								throw e;
							}
							// there has been no exception, so the problem was feasible
							nothingfound = false;
							if (reachable_only)
								c_X[t].clear(p); // reset corner

							// can extract the distribution now from the solution
							double[] solution_point = new double[l * nt];
							for(int i = 0; i < l; i++) {
								int uu = 0;
								for(int u = 0; u < nt; u++) {
									if(LIST_multiTuple.get(u) != null && LIST_multiTuple.get(u).get(i) != null) {
										solution_point[l * u + i] = solution.getPoint()[l * uu + i];
										uu++;
									}
								}
							}
							Distribution next_move = new Distribution(); // NEXT MOVE - specific to t and p
							go_through_moves: for (int u = 0; u < nt; u++) {
								// calculate alpha^u
								double alpha_u = 0.0;
								for (int i = 0; i < l; i++) {
									alpha_u += solution_point[l * u + i]; // alpha^u = sum_i gamma^u_i
								}
								// alpha^u = 0 implies that this move is not taken
								if (PrismUtils.doublesAreEqual(alpha_u, 0.0))
									continue go_through_moves;
								next_move.add(u, alpha_u);

								Distribution d = new Distribution(); // NEXT MEMORY - specific to t, p and u
								for (int i = 0; i < l; i++) {
									Integer index = gss.get(u).indexOf(LIST_multiTuple.get(u).get(i));
									double prob = solution_point[l * u + i] / alpha_u; // beta^u_i = gamma^u_i / alpha^u
									prob = prob > 1.0 ? 1.0 : prob;
									prob = prob < 0.0 ? 0.0 : prob;
									d.add(index, prob);
								}
								pi_t[t].get(p).put(u, renormDistribution(d));
								memorySize += d.size();
								if (log_problem) System.out.printf("pi_t(%d)\tt:%d, p:%d --u:%d--> %s \n", G.getPlayer(t), t, p, u, d.toString());
								if (logStrategy) mainLog.print(String.format("pi_t(%d)\tt:%d, p:%d --u:%d--> %s \n", G.getPlayer(t), t, p, u, d.toString()));

								// STOCHASTIC
								List<double[]> gsYtu = gsToList(Y[t].get(u).get().minimized_generators(), n);
								// set up corner tracking
								if (reachable_only) {
									if (c_Y[t].get(u) == null) { // make sure not to overwrite
										c_Y[t].put(u, new BitSet(gsYtu.size()));
									}
									// mark corners q at move to be handled, if not already handled
									for (Integer q : d.getSupport()) {
										if (pi_u[t] == null || pi_u[t].get(u) == null || pi_u[t].get(u).get(q) == null) { // need to handle corner point wq at move u of t
											c_Y[t].get(u).set(q);
											corners_to_cover = true;
										}
									}
								}

								if (pi_u[t].get(u) == null) { // make sure to not overwrite
									pi_u[t].put(u, new HashMap<Integer, Map<Integer, Distribution>>());
								}

								int ntu = G.getNumTransitions(t, u);
								Iterator<Entry<Integer, Double>> dtu;

								if (!reachable_only || c_Y[t].get(u).cardinality() != 0) { // some corners at move u of t need to be worked on
									for (int q_index = 0; q_index < gsYtu.size(); q_index++) {
										if (!reachable_only || c_Y[t].get(u).get(q_index)) { // corner q at move u of t needs to be handled
											if (reachable_only)
												c_Y[t].get(u).clear(q_index); // reset

											if (pi_u[t].get(u).get(q_index) == null) { // make sure to not replace
												pi_u[t].get(u).put(q_index, new HashMap<Integer, Distribution>());
											}

											// ll is number of objectives, i.e. number of required corners
											// start with 1, continue up to n, i.e. total number of goals
											Distribution[] action = null;
											search_for_stochastic_corners: for (int ll = 1; ll < n + 1; ll++) {
												try {
													action = getActions(G, ntu, u, t, q_index, ll, rewards, accuracy, solver, gsYtu, LIST_gsX, rounding);
													break search_for_stochastic_corners;
												} catch (PrismException e) {
													//System.out.printf("Nothing found for ll=%d\n", ll);
													continue search_for_stochastic_corners;
												}
											}
											if (action != null) {
												dtu = G.getTransitionsIterator(t, u); // iterates through w
												for (int w = 0; w < ntu; w++) {
													Entry<Integer, Double> e_w = dtu.next();
													int key_w = e_w.getKey();
													double val_w = e_w.getValue();
													Distribution mem_d = renormDistribution(action[w]);
													pi_u[t].get(u).get(q_index).put(key_w, mem_d);
													memorySize += action[w].size();
													if (log_problem) System.out.printf("pi_u(%d)\tt:%d, u:%d, q:%d --w:%d--> %s \n", G.getPlayer(t), t, u, q_index, key_w, action[w].toString());
													if (logStrategy) mainLog.print(String.format("pi_u(%d)\tt:%d, u:%d, q:%d --w:%d--> %s \n", G.getPlayer(t), t, u, q_index, key_w, action[w].toString()));
													// mark corners at the successors to be handled, if not already handled
													if (reachable_only) {
														for (Integer r : mem_d.getSupport()) {
															if (pi_t[key_w] == null || pi_t[key_w].get(r) == null) { // need to still handle corner point r at state w
																if (key_w == t) { // self loop
																	if (log_problem)
																		System.out.printf("selfloop: key_w=%d, r=%d\n", key_w, r);
																	c_Xt.set(r);
																	corners_from_selfloops = true;
																} else {
																	if (log_problem)
																		System.out.printf("explore: key_w=%d, r=%d\n", key_w, r);
																	c_X[key_w].set(r);
																	corners_to_cover = true;
																}
															}
														}
													}
												}
											} else {
												if (log_problem) System.out.printf("PROBLEM(P1):t:%d, u:%d, q:%d --w:?--> ? \n", t, u, q_index);
												if (logStrategy) mainLog.print(String.format("PROBLEM(P1):t:%d, u:%d, q:%d --w:?--> ? \n", t, u, q_index));

												throw new PrismLangException("No distributions for move found");
											}
										}
									}
								}
								// STOCHASTIC END
							}
							// set next move choice at t with corner p
							pi_n[t].put(p, renormDistribution(next_move));

							continue iterate_through_corners;
						}
						// nothing found for l - loop, increase l
					}
					if (nothingfound) {
						throw new PrismLangException("No distribution found for corner");
					}
				}
			}

		} while (corners_from_selfloops);

		// reset c_X at t
		if (log_problem)
			System.out.printf("resetting c_X[%d]\n", t);
		if (reachable_only)
			c_X[t].clear();

		// new corners to cover?
		return corners_to_cover;
	}

	/**
	 * computes memory update functions for choices of P2
	 *
	 * G ... game to construct strategy for
	 * t ... the P2 state
	 * Yt ... Pareto sets at moves of t
	 * LIST_gsX ... canonical order of points at X
	 * c_X ... corners to cover
	 * c_Y ... corners to cover at moves
	 * rewards ... the reward structures of the game under consideration
	 * accuracy ... accuracy per dimension
	 * solver ... the Simplex solver for solving the LPs
	 * reachable_only ... only construct strategy with reachable corners and states
	 * rounding ... whether rounding is active
	 *
	 * return value:
	 * corners_to_cover ... whether new corners need to be covered in the reachable_only search
	 **/
	private boolean getP2Choices(SMG G, int t, List<Pareto>[] Y, List<double[]>[] LIST_gsX, BitSet[] c_X, Map<Integer, BitSet>[] c_Y, List<SMGRewards> rewards,
			long[] accuracy, SimplexSolver solver, boolean reachable_only, boolean rounding) throws PrismException
	{
		int n = rewards.size();
		double[] bounds;

		boolean corners_to_cover = false; // return value

		// to handle selfloops
		boolean corners_from_selfloops = false;
		BitSet c_Xt = reachable_only ? new BitSet() : null;
		do { // while corners_from_selfloops
			corners_from_selfloops = false;
			// for each action need to build a new distribution
			for (int u = 0; u < G.getNumChoices(t); u++) { // for each stochastic successor (i.e. action u)
				if (pi_u[t].get(u) == null) { // make sure to not overwrite
					pi_u[t].put(u, new HashMap<Integer, Map<Integer, Distribution>>());
				}
				if (log_problem) System.out.printf("looking for: t:%d --u:%d--> \n", t, u);
				if (logStrategy) mainLog.print(String.format("looking for: t:%d --u:%d--> \n", t, u));

				List<double[]> gsYtu = gsToList(Y[t].get(u).get().minimized_generators(), n);
				if (reachable_only) {
					if (c_Y[t].get(u) == null) { // make sure not to overwrite
						c_Y[t].put(u, new BitSet(gsYtu.size()));
					}
				}

				int ntu = G.getNumTransitions(t, u);
				Iterator<Entry<Integer, Double>> dtu;

				// now for each corner point p for t, need to find a distribution
				// that is, find l, and l coefficients beta_i summing to one such that
				// for good and bad states: sum_i beta_i q_i^u >= p - rewards(t)
				// and for stochastic states: ...
				// for q^i_u in Y(t,u) - need to actually find these
				for (int p = 0; p < LIST_gsX[t].size(); p++) { // for each corner point in t

					if (!reachable_only || c_X[t].get(p) || c_Xt.get(p)) { // if corner p needs to be handled

						if (log_problem) System.out.printf("looking for: t:%d p:%d(%s) --u:%d--> \n", t, p, Arrays.toString(LIST_gsX[t].get(p)), u);
						if (logStrategy) mainLog.print(String.format("looking for: t:%d p:%d(%s) --u:%d--> \n", t, p, Arrays.toString(LIST_gsX[t].get(p)), u));

						if (!pi_t[t].containsKey(p)) {
							pi_t[t].put(p, new HashMap<Integer, Distribution>());
						}
						// put value of p - reward(t) - varepsilon into bounds
						bounds = new double[n];
						for (int k = 0; k < n; k++) {
							if (rounding) {
								bounds[k] = Math.floor((LIST_gsX[t].get(p)[k] - rewards.get(k).getStateReward(t) - varepsilon) * accuracy[k]) / accuracy[k];
							} else {
								bounds[k] = LIST_gsX[t].get(p)[k] - rewards.get(k).getStateReward(t) - varepsilon;
							}
						}
						// find q_i^u and beta_i in Y(t,u)
						search_for_distribution: for (int l = 1; l < n + 1; l++) { // first find l
							List<double[]> LIST_tuple = null;

							// preparation for LP
							double[] coeffs_beta = new double[l];
							for (int i = 0; i < l; i++) {
								coeffs_beta[i] = 1.0;
							}
							// check for each such tuple
							iteration_through_tuples: while (true) {
								LIST_tuple = selectGenerator(gsYtu, LIST_tuple, l);
								if (LIST_tuple == null)
									break iteration_through_tuples;
								// now formulate an LP for beta_i

								// max_{beta_i} sum_i beta
								// s.t. sum_i beta_i q_i^u >= p - rewards(t)
								//      sum_i beta_i = 1

								// describe the optimization problem
								// optimization function - maximize betas
								LinearObjectiveFunction f = new LinearObjectiveFunction(coeffs_beta, 0);

								// constraints
								List<LinearConstraint> constraints = new ArrayList<LinearConstraint>();
								double[][] coeffs_q = new double[n][l];
								for (int i = 0; i < l; i++) {
									// get coefficients from tuple.get(i)
									for (int k = 0; k < n; k++) {
										coeffs_q[k][i] = LIST_tuple.get(i)[k];
									}
								}

								for (int i = 0; i < n; i++) {
									// lower bound on sum of betas
									constraints.add(new LinearConstraint(coeffs_q[i], Relationship.GEQ, bounds[i]));
								}
								double[][] onlyone = new double[l][l];
								for (int i = 0; i < l; i++) {
									// lower bound on beta^w_i
									onlyone[i][i] = 1.0;
									constraints.add(new LinearConstraint(onlyone[i], Relationship.GEQ, 0.0));
								}
								// strict bound on sum of beta ... probability distribution
								constraints.add(new LinearConstraint(coeffs_beta, Relationship.EQ, 1));

								PointValuePair solution = null;
								// TODO: MaxIter?
								try {
								    if(log_problem)
									mainLog.print(String.format("----\ngoal: %s\n\nbounds: %s\n\ncoeffs: %s\n\ncoeffs_indiv: %s\n\nonlyone: %s\n----\n",
												    Arrays.toString(coeffs_beta),
												    Arrays.toString(bounds),
												    Arrays.deepToString(coeffs_q),
												    "",
												    Arrays.deepToString(onlyone)));
								        solution = solver.optimize(f, new LinearConstraintSet(constraints), GoalType.MAXIMIZE, new MaxIter(10000));
								} catch (NoFeasibleSolutionException e) {
									// tuple not feasible, try a different one
									continue iteration_through_tuples;
								} catch (Exception e2) { // any other exception
									throw new PrismLangException(e2.toString()); // TODO: let's see what happens
									//continue iteration_through_tuples; // TODO: what's happening here? shouldn't continue, this might be a bug!
								}

								// memory-update distribution for non-stochastic state part
								Distribution d = new Distribution(); // specific to t, p, and u
								for (int i = 0; i < l; i++) {
									Integer q_index = gsYtu.indexOf(LIST_tuple.get(i));
									double beta = solution.getPoint()[i]; // beta^u_i
									beta = beta > 1.0 ? 1.0 : beta;
									beta = beta < 0.0 ? 0.0 : beta;
									d.add(q_index, beta);
								}
								pi_t[t].get(p).put(u, renormDistribution(d));
								memorySize += d.size();
								if (log_problem) System.out.printf("pi_t(%d)\tt:%d, p:%d --u:%d--> %s \n", G.getPlayer(t), t, p, u, d.toString());
								if (logStrategy) mainLog.print(String.format("pi_t(%d)\tt:%d, p:%d --u:%d--> %s \n", G.getPlayer(t), t, p, u, d.toString()));

								// mark corners q at move to be handled, if not already handled
								if (reachable_only) {
									for (Integer q : d.getSupport()) {
										if (pi_u[t] == null || pi_u[t].get(u) == null || pi_u[t].get(u).get(q) == null) { // need to handle corner point q at move u of t
											c_Y[t].get(u).set(q);
											corners_to_cover = true;
										}
									}
								}

								break search_for_distribution; // a distribution for this l has been found
							}
						}
					} // check if corner p needs handling
				} // loop through corners p at t

				// STOCHASTIC
				if (!reachable_only || c_Y[t].get(u).cardinality() != 0) { // some corners at move u of t need to be worked on
					for (int q_index = 0; q_index < gsYtu.size(); q_index++) {
						if (!reachable_only || c_Y[t].get(u).get(q_index)) { // corner q at move u of t needs to be handled
							if (reachable_only)
								c_Y[t].get(u).clear(q_index); // reset

							if (pi_u[t].get(u).get(q_index) == null) { // make sure to not replace
								pi_u[t].get(u).put(q_index, new HashMap<Integer, Distribution>());
							}

							// ll is number of objectives, i.e. number of required corners
							// start with 1, continue up to n, i.e. total number of goals
							Distribution[] action = null;
							search_for_stochastic_corners: for (int ll = 1; ll < n + 1; ll++) {
								try {
									action = getActions(G, ntu, u, t, q_index, ll, rewards, accuracy, solver, gsYtu, LIST_gsX, rounding);
									break search_for_stochastic_corners;
								} catch (PrismException e) {
									if (log_problem)
										System.out.printf("Nothing found for ll=%d\n", ll);
									continue search_for_stochastic_corners;
								}
							}
							if (action != null) {
								dtu = G.getTransitionsIterator(t, u); // iterates through w
								for (int w = 0; w < ntu; w++) {
									Entry<Integer, Double> e_w = dtu.next();
									int key_w = e_w.getKey();
									double val_w = e_w.getValue();
									Distribution d = renormDistribution(action[w]);
									pi_u[t].get(u).get(q_index).put(key_w, d);
									memorySize += action[w].size();
									if (log_problem) System.out.printf("pi_u(%d)\tt:%d, u:%d, q:%d --w:%d--> %s \n", G.getPlayer(t), t, u, q_index, key_w, action[w].toString());
									if (logStrategy) mainLog.print(String.format("pi_u(%d)\tt:%d, u:%d, q:%d --w:%d--> %s \n", G.getPlayer(t), t, u, q_index, key_w, action[w].toString()));

									// mark corners at the successors to be handled, if not already handled
									if (reachable_only) {
										for (Integer r : d.getSupport()) {
											if (pi_t[key_w] == null || pi_t[key_w].get(r) == null) { // need to still handle corner point r at state w
												if (key_w == t) { // self loop
													if (log_problem)
														System.out.printf("selfloop: key_w=%d, r=%d\n", key_w, r);
													corners_from_selfloops = true;
													c_Xt.set(r);
												} else {
													if (log_problem)
														System.out.printf("explore: key_w=%d, r=%d\n", key_w, r);
													c_X[key_w].set(r);
													corners_to_cover = true;
												}
											}
										}
									}
								}
							} else {
								if (log_problem) System.out.printf("PROBLEM(P2):t:%d, u:%d, q:%d --w:?--> ? \n", t, u, q_index);
								if (logStrategy) mainLog.print(String.format("PROBLEM(P2):t:%d, u:%d, q:%d --w:?--> ? \n", t, u, q_index));

								throw new PrismLangException("No distributions for move found");
							}
						}
					}
				}
				// STOCHASTIC END
			} // go through all moves

		} while (corners_from_selfloops);

		// reset c_X at t
		if (log_problem)
			System.out.printf("resetting c_X[%d]\n", t);
		if (reachable_only)
			c_X[t].clear();

		// new corners to cover?
		return corners_to_cover;
	}

	/**
	 * computes the memory update functions for moves (stochastic states)
	 *
	 * G ... game to construct strategy for
	 * ntu ... number of successor states following the move u
	 * u ... move considered (a move following t)
	 * t ... state considered
	 * q ... index of corner at move u to be looked at
	 * n ... number of points in convex combination of corners allowed
	 * rewards ... the reward structures of the game under consideration
	 * accuracy ... accuracy per dimension
	 * solver ... the Simplex solver for solving the LPs
	 * reachable_only ... only construct strategy with reachable corners and states
	 * rounding ... whethe rounding is active
	 **/
	private Distribution[] getActions(SMG G, int ntu, int u, int t, int q, int n, List<SMGRewards> rewards, long[] accuracy, SimplexSolver solver,
			List<double[]> gsYtu, List<double[]>[] LIST_gsX, boolean rounding) throws PrismException
	{
		// interpret u as a stochastic state, and look at all its successors w

		// the result:
		Distribution[] stochastic = new Distribution[ntu];

		int dim = rewards.size();

		// bounds in stochastic states are independent of state rewards, but need transition rewards
		double[] bounds = new double[dim];
		for (int k = 0; k < dim; k++) {
			if (rounding) {
				bounds[k] = Math.floor((gsYtu.get(q)[k] - rewards.get(k).getTransitionReward(t, u) - varepsilon) * accuracy[k]) / accuracy[k];
			} else {
				bounds[k] = gsYtu.get(q)[k] - rewards.get(k).getTransitionReward(t, u) - varepsilon;
			}
		}

		double[] coeffs_beta = new double[ntu * n];
		double[][] coeffs_beta_indiv = new double[ntu][n * ntu];
		for (int i = 0; i < n; i++) {
			for (int w = 0; w < ntu; w++) {
				coeffs_beta_indiv[w][w * n + i] = 1.0;
				coeffs_beta[w * n + i] = 1.0;
			}
		}
		List<LinearConstraint> constraints = new ArrayList<LinearConstraint>();

		List<List<double[]>> gss = new ArrayList<List<double[]>>();
		Iterator<Entry<Integer, Double>> dtu = G.getTransitionsIterator(t, u);
		for (int w = 0; w < ntu; w++) {
			int key_w = dtu.next().getKey();
			gss.add(LIST_gsX[key_w]);
		}
		List<List<double[]>> LIST_multiTuple = null;

		int multi_counter = 0;

		boolean nothingfound = true;
		iteration_through_multi_tuples: while (true) {
			LIST_multiTuple = selectMultiGenerator(gss, LIST_multiTuple, n);
			if (LIST_multiTuple == null)
				break iteration_through_multi_tuples;
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
			double[][] coeffs_q = new double[dim][ntu * n];
			dtu = G.getTransitionsIterator(t, u);
			for (int w = 0; w < ntu; w++) { // for each successor w
				double delta_uw = dtu.next().getValue();
				for (int i = 0; i < n; i++) { // for each component
					// get coefficients from tuple.get(i)
					for (int k = 0; k < dim; k++) {
						coeffs_q[k][w * n + i] = delta_uw * LIST_multiTuple.get(w).get(i)[k];
					}
				}
			}

			for (int i = 0; i < dim; i++) {
				// lower bound on sum of betas
				constraints.add(new LinearConstraint(coeffs_q[i], Relationship.GEQ, bounds[i]));
			}
			double[][] onlyone = new double[n * ntu][n * ntu];
			for (int i = 0; i < n * ntu; i++) {
				// lower bound on beta^w_i
				onlyone[i][i] = 1.0;
				constraints.add(new LinearConstraint(onlyone[i], Relationship.GEQ, 0.0));
			}
			// strict bound on sum of beta ... probability distribution
			for (int w = 0; w < ntu; w++) { // for each successor w
				constraints.add(new LinearConstraint(coeffs_beta_indiv[w], Relationship.EQ, 1.0));
			}

			PointValuePair solution = null;
			try {
			    if(log_problem)
				System.out.printf("----\ngoal: %s\n\nbounds: %s\n\ncoeffs: %s\n\ncoeffs_indiv: %s\n\nonlyone: %s\n\nsolution: %s\n----\n",
						  Arrays.toString(coeffs_beta),
						  Arrays.toString(bounds),
						  Arrays.deepToString(coeffs_q),
						  Arrays.deepToString(coeffs_beta_indiv),
						  Arrays.deepToString(onlyone),
						  solution==null?"null":Arrays.toString(solution.getPoint()));
				solution = solver.optimize(f, new LinearConstraintSet(constraints), GoalType.MAXIMIZE, new MaxIter(10000));
			} catch (NoFeasibleSolutionException e) {
				// tuple not feasible, try a different one
				continue iteration_through_multi_tuples;
			} catch (UnboundedSolutionException e) {
				throw e;
			}
			nothingfound = false;
			// there has been no exception, so the problem was feasible
			// can extract the distribution now from the solution
			dtu = G.getTransitionsIterator(t, u);
			for (int w = 0; w < ntu; w++) { // for each successor
				stochastic[w] = new Distribution(); // initialize for corner q of u

				int key_w = dtu.next().getKey();
				for (int i = 0; i < n; i++) { // for each dimension
					Integer index = LIST_gsX[key_w].indexOf(LIST_multiTuple.get(w).get(i));
					double prob = solution.getPoint()[n * w + i];
					prob = prob > 1.0 ? 1.0 : prob;
					prob = prob < 0.0 ? 0.0 : prob;
					stochastic[w].add(index, prob); // adds if index already exists ... don't know which dimension the LP-solver has assigned the probability mass to
				}
			}
			break iteration_through_multi_tuples;
		}
		//System.out.println(nothingfound);
		if (nothingfound) {
			throw new PrismLangException("Nothing found for L");
		}
		return stochastic;
	}

	/**
	 * scales probabilities in a distribution to sum to one
	 **/
	private Distribution renormDistribution(Distribution d)
	{
		Distribution result = new Distribution();
		double total = 0.0;
		for (Entry<Integer, Double> e : d) {
			total += e.getValue();
		}
		for (Entry<Integer, Double> e : d) {
			result.add(e.getKey(), e.getValue() / total);
		}
		return result;
	}

	@Override
	public void exportActions(PrismLog out)
	{
		out.print(this.toString());
	}

	@Override
	public void exportIndices(PrismLog out)
	{
		// TODO
	}

	@Override
	public void exportInducedModel(PrismLog out)
	{
		// TODO
	}

	@Override
	public void clear()
	{
		// TODO Auto-generated method stub

	}

	@Override
	public void exportDotFile(PrismLog out)
	{
		// TODO Auto-generated method stub

	}

	@Override
	public void restrictStrategyToReachableStates() throws PrismException
	{
		// TODO Auto-generated method stub
		throw new PrismException("Reach option is not supported for this strategy type");
	}

	@Override
	public void exportStratToFile(File file, StrategyExportType exportType)
	{
		// TODO Auto-generated method stub
		
	}

	@Override
	public HashMap<String, Double> getNextAction(int state) throws InvalidStrategyStateException
	{
		// TODO Auto-generated method stub
		return null;
	}
}
