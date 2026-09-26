//==============================================================================
//
//	Copyright (c) 2026-
//	Authors:
//	* Gabriel Santos <gabriel.santos@cs.ox.ac.uk> (University of Oxford)
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

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import explicit.CSG;
import explicit.CSGMultiObjectives;
import explicit.Distribution;
import explicit.DistributionOver;
import explicit.MDPSimple;
import parser.State;
import prism.JointAction;
import prism.PrismException;
import prism.PrismLog;
import prism.PrismNotSupportedException;

/**
 * Equilibrium strategy profile for a CSG with more than two coalitions (Nash or correlated), as synthesised by
 * CSGModelCheckerEquilibria. It has memory: the coalitions that are done (D) or have failed (E) and, while some
 * bounded objective is undecided, the number of steps taken (t). The memory is updated as in the model checking
 * (see {@link CSGMultiObjectives}): on entering a state, the step count advances and the coalitions whose objectives
 * become decided there are added to D or E. The choice in (state, memory) is the local equilibrium strategy of that
 * subgame: for Nash equilibria, one distribution per coalition over its actions (independent); for correlated
 * equilibria, one distribution over joint actions. Either way, the decision returned by getChoiceAction is a
 * distribution over the joint actions ({@link JointAction}) of the state's choices.
 * <br>
 * Memory values are indices encoding (D, E, t): a base-3 digit per coalition (0 undecided, 1 done, 2 failed),
 * times (horizon + 2), plus t + 1 (t = -1 once all bounded objectives are decided).
 */
public class CSGEquilibriumStrategy extends StrategyWithStates<Double>
{
	protected final CSG<Double> model;
	protected final CSGMultiObjectives obj;
	protected final int n;
	protected final boolean correlated;
	/** Local strategies, indexed by memory and state (null if undefined) */
	protected final List<Map<BitSet, Double>>[][] local;
	protected final int memorySize;

	/**
	 * @param model The CSG
	 * @param obj The coalitions' objectives
	 * @param correlated Whether the local strategies are correlated (joint) or Nash (one per coalition)
	 * @param localStrategies Local strategies keyed by D (bits 0..n-1), E (bits n..2n-1) and, while bounded objectives
	 *        are undecided, the step t (bit 2n + 1 + t); per state, as lists of maps from (player) action indices to
	 *        probabilities (one per coalition for Nash, a single joint one for correlated)
	 */
	@SuppressWarnings("unchecked")
	public CSGEquilibriumStrategy(CSG<Double> model, CSGMultiObjectives obj, boolean correlated, Map<BitSet, List<Map<BitSet, Double>>[]> localStrategies)
	{
		this.model = model;
		this.obj = obj;
		this.n = obj.n;
		this.correlated = correlated;
		int pow3 = 1;
		for (int c = 0; c < n; c++)
			pow3 *= 3;
		memorySize = pow3 * (obj.horizon + 2);
		local = (List<Map<BitSet, Double>>[][]) new List[memorySize][];
		for (Map.Entry<BitSet, List<Map<BitSet, Double>>[]> e : localStrategies.entrySet()) {
			BitSet k = e.getKey();
			int tb = k.nextSetBit(2 * n + 1);
			local[encode(k.get(0, n), k.get(n, 2 * n), tb < 0 ? -1 : tb - 2 * n - 1)] = e.getValue();
		}
		// state look-up (for the StrategyGenerator interface, e.g. the simulator)
		Map<State, Integer> index = new HashMap<>();
		List<State> states = model.getStatesList();
		if (states != null) {
			for (int s = 0; s < states.size(); s++)
				index.put(states.get(s), s);
		}
		setStateLookUp(state -> index.getOrDefault(state, -1));
	}

	// Memory encoding

	protected int encode(BitSet D, BitSet E, int t)
	{
		int code = 0;
		for (int c = n - 1; c >= 0; c--)
			code = 3 * code + (D.get(c) ? 1 : E.get(c) ? 2 : 0);
		return code * (obj.horizon + 2) + (t + 1);
	}

	/** Decodes memory m into D and E (which are cleared first); returns t */
	protected int decode(int m, BitSet D, BitSet E)
	{
		D.clear();
		E.clear();
		int t = m % (obj.horizon + 2) - 1;
		int code = m / (obj.horizon + 2);
		for (int c = 0; c < n; c++) {
			int d = code % 3;
			if (d == 1)
				D.set(c);
			else if (d == 2)
				E.set(c);
			code /= 3;
		}
		return t;
	}

	/** Memory after entering state s with D, E and step t (t = -1: all bounded objectives already decided) */
	protected int enter(BitSet D, BitSet E, int t, int s)
	{
		obj.classify(s, t, D, E);
		if (t >= 0 && obj.boundedDecided(D, E))
			t = -1;
		return encode(D, E, t);
	}

	// Strategy interface

	@Override
	public Memory memory()
	{
		return Memory.FINITE;
	}

	@Override
	public int getMemorySize()
	{
		return memorySize;
	}

	@Override
	public int getInitialMemory(int sInit)
	{
		return enter(new BitSet(), new BitSet(), obj.bounded.isEmpty() ? -1 : 0, sInit);
	}

	@Override
	public int getUpdatedMemory(int m, Object action, int sNext)
	{
		BitSet D = new BitSet(), E = new BitSet();
		int t = decode(m, D, E);
		return enter(D, E, t < 0 ? -1 : t + 1, sNext);
	}

	@Override
	public String getMemoryString(int m)
	{
		if (m < 0 || m >= memorySize)
			return "?";
		BitSet D = new BitSet(), E = new BitSet();
		int t = decode(m, D, E);
		return "D=" + D + " E=" + E + (t >= 0 ? " t=" + t : "");
	}

	@Override
	public boolean isRandomised()
	{
		return true;
	}

	/** Actions (player action indices, idle actions included) of choice i of state s */
	protected BitSet playerActions(int s, int i)
	{
		BitSet b = new BitSet();
		int[] indexes = model.getIndexes(s, i);
		for (int p = 0; p < indexes.length; p++)
			b.set(indexes[p] > 0 ? indexes[p] : model.getIdles()[p]);
		return b;
	}

	/** Probability of choice i of state s under local strategy strat */
	protected double choiceProbability(List<Map<BitSet, Double>> strat, int s, int i)
	{
		BitSet b = playerActions(s, i);
		if (correlated)
			return strat.get(0).getOrDefault(b, 0.0);
		double p = 1.0;
		for (Map<BitSet, Double> sc : strat) {
			double q = 0.0;
			for (Map.Entry<BitSet, Double> e : sc.entrySet()) {
				BitSet k = (BitSet) e.getKey().clone();
				k.andNot(b);
				if (k.isEmpty()) {
					q = e.getValue();
					break;
				}
			}
			p *= q;
			if (p == 0.0)
				break;
		}
		return p;
	}

	/** Distribution over the choices (indices) of state s in memory m, or null if undefined */
	protected Distribution<Double> choiceDistribution(int s, int m)
	{
		if (m < 0 || m >= memorySize || local[m] == null || local[m][s] == null)
			return null;
		Distribution<Double> d = new Distribution<>();
		for (int i = 0; i < model.getNumChoices(s); i++) {
			double p = choiceProbability(local[m][s], s, i);
			if (p > 0.0)
				d.add(i, p);
		}
		return d.isEmpty() ? null : d;
	}

	@Override
	public Object getChoiceAction(int s, int m)
	{
		Distribution<Double> d = choiceDistribution(s, m);
		if (d == null)
			return UNDEFINED;
		return DistributionOver.create(d, i -> model.getAction(s, i));
	}

	@Override
	public int getChoiceIndex(int s, int m)
	{
		Distribution<Double> d = choiceDistribution(s, m);
		if (d == null || d.size() != 1)
			return -1;
		return d.getSupport().iterator().next();
	}

	/**
	 * Joint actions may also be given as arrays of (1-indexed) player action indices (-1 for idle), which is how
	 * model generators (e.g. the simulator's) represent the actions of CSG choices.
	 */
	protected Object toJointAction(Object act)
	{
		return act instanceof int[] ? new JointAction((int[]) act, model.getActions()) : act;
	}

	@Override
	public Double getChoiceActionProbability(Object decision, Object act)
	{
		return super.getChoiceActionProbability(decision, toJointAction(act));
	}

	@Override
	public boolean isActionChosen(Object decision, Object act)
	{
		return super.isActionChosen(decision, toJointAction(act));
	}

	@Override
	public CSG<Double> getModel()
	{
		return model;
	}

	// Export

	@Override
	public prism.Model<Double> constructInducedModel(StrategyExportOptions options) throws PrismException
	{
		throw new PrismNotSupportedException("CSG strategy product not yet supported");
	}

	@Override
	public void exportActions(PrismLog out, StrategyExportOptions options) throws PrismException
	{
		throw new PrismNotSupportedException("CSG strategy export in this format not yet supported");
	}

	@Override
	public void exportIndices(PrismLog out, StrategyExportOptions options) throws PrismException
	{
		throw new PrismNotSupportedException("CSG strategy export in this format not yet supported");
	}

	@Override
	public void exportInducedModel(PrismLog out, StrategyExportOptions options) throws PrismException
	{
		throw new PrismNotSupportedException("CSG strategy export in this format not yet supported");
	}

	/**
	 * Replays the strategy from the initial state and exports the resulting graph, one node per reachable
	 * (state, memory) pair, as in the two-coalition case: choices are labelled with the local strategy
	 * ("CSG: p: [a] + ... -- p: [b] + ..." for Nash, per coalition; "CSG: p: [a][b][c] + ..." for correlated,
	 * over joint actions) and the memory; once every objective is decided, a self-loop labelled Sat(i)/Unsat(i).
	 */
	@Override
	public void exportDotFile(PrismLog out, StrategyExportOptions options) throws PrismException
	{
		MDPSimple<Double> mdp = new MDPSimple<>();
		List<State> statelist = new ArrayList<>();
		Map<Long, Integer> nodes = new HashMap<>();
		ArrayDeque<long[]> todo = new ArrayDeque<>();
		mdp.setVarList(model.getVarList());
		int s0 = model.getFirstInitialState();
		int m0 = getInitialMemory(s0);
		int n0 = node(mdp, statelist, nodes, todo, s0, m0);
		mdp.addInitialState(n0);
		BitSet D = new BitSet(), E = new BitSet();
		while (!todo.isEmpty()) {
			long[] sm = todo.poll();
			int s = (int) sm[0], m = (int) sm[1];
			int nd = nodes.get(key(s, m));
			int t = decode(m, D, E);
			if (D.cardinality() + E.cardinality() == n) {
				Distribution<Double> d = new Distribution<>();
				d.add(nd, 1.0);
				mdp.addActionLabelledChoice(nd, d, decidedLabel(D));
				continue;
			}
			Distribution<Double> strat = choiceDistribution(s, m);
			if (strat == null) {
				Distribution<Double> d = new Distribution<>();
				d.add(nd, 1.0);
				mdp.addActionLabelledChoice(nd, d, "undefined {" + getMemoryString(m) + "}");
				continue;
			}
			Distribution<Double> d = new Distribution<>();
			for (int i : strat.getSupport()) {
				double p = strat.get(i);
				for (int u : model.getChoice(s, i).getSupport()) {
					int mu = getUpdatedMemory(m, null, u);
					int nu = node(mdp, statelist, nodes, todo, u, mu);
					d.add(nu, p * model.getChoice(s, i).get(u));
				}
			}
			mdp.addActionLabelledChoice(nd, d, label(local[m][s]) + " {" + getMemoryString(m) + "}");
		}
		mdp.setStatesList(statelist);
		mdp.exportToDotFile(out, null, true);
		out.print("\n/*");
		out.print("\n -- Transitions --  \n");
		mdp.exportToPrismExplicitTra(out);
		out.print("\n -- States --  \n");
		mdp.exportStates(0, mdp.getVarList(), out);
		out.print("*/\n");
	}

	private static long key(int s, int m)
	{
		return ((long) s << 32) | (m & 0xffffffffL);
	}

	private int node(MDPSimple<Double> mdp, List<State> statelist, Map<Long, Integer> nodes, ArrayDeque<long[]> todo, int s, int m)
	{
		Integer nd = nodes.get(key(s, m));
		if (nd == null) {
			nd = mdp.addState();
			nodes.put(key(s, m), nd);
			statelist.add(model.getStatesList().get(s));
			todo.add(new long[] { s, m });
		}
		return nd;
	}

	private String actionsString(BitSet acts)
	{
		StringBuilder sb = new StringBuilder();
		for (int i = acts.nextSetBit(0); i >= 0; i = acts.nextSetBit(i + 1))
			sb.append("[").append(model.getActions().get(i - 1)).append("]");
		return sb.toString();
	}

	/** Label of a local strategy: per coalition (Nash) or over joint actions (correlated) */
	private String label(List<Map<BitSet, Double>> strat)
	{
		StringBuilder sb = new StringBuilder("CSG: ");
		for (int c = 0; c < strat.size(); c++) {
			if (c > 0)
				sb.append(" -- ");
			int k = 0;
			for (Map.Entry<BitSet, Double> e : strat.get(c).entrySet()) {
				if (k++ > 0)
					sb.append(" + ");
				sb.append(e.getValue()).append(": ").append(actionsString(e.getKey()));
			}
		}
		return sb.toString();
	}

	private String decidedLabel(BitSet D)
	{
		StringBuilder sb = new StringBuilder();
		for (int c = 0; c < n; c++) {
			if (c > 0)
				sb.append(" -- ");
			sb.append(D.get(c) ? "Sat(" : "Unsat(").append(c).append(")");
		}
		return sb.toString();
	}

	@Override
	public void clear()
	{
	}
}
