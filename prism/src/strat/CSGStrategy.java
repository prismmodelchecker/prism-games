//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford)
//  * Gabriel Santos <gabriel.santos@cs.ox.ac.uk> (University of Oxford)
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
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import explicit.CSG;
import explicit.Distribution;
import explicit.DistributionOver;
import explicit.MDPModelChecker;
import explicit.MDPSimple;
import explicit.rewards.CSGRewards;
import explicit.rewards.MDPRewardsSimple;
import parser.State;
import prism.JointAction;
import prism.PrismException;
import prism.PrismLog;
import prism.PrismNotSupportedException;

/**
 * Strategy of a coalition C in a zero-sum CSG (C against the other players N\C), memoryless (unbounded properties,
 * and next). Its decisions are distributions over C's part of the joint actions (partial joint actions, with the
 * other players' entries undefined): a joint action is chosen whenever its part for C is. The strategy induces an MDP,
 * in which the other players' choices remain nondeterministic.
 */
public class CSGStrategy extends StrategyWithStates<Double> {

	protected CSG<Double> model;
	protected List<List<List<Map<BitSet, Double>>>> csgchoices; // player -> iteration -> state -> indexes -> value
	protected BitSet no;
	protected BitSet yes;
	protected BitSet inf;
	
	public CSGStrategy(CSG<Double> model, List<List<List<Map<BitSet, Double>>>> csgchoices, BitSet no, BitSet yes, BitSet inf) {
		this.model = model;
		this.csgchoices = csgchoices;
		this.no = no;
		this.yes = yes;
		this.inf = inf;
		// state look-up (for the StrategyGenerator interface, e.g. the simulator)
		Map<State, Integer> index = new HashMap<>();
		List<State> states = model.getStatesList();
		if (states != null) {
			for (int s = 0; s < states.size(); s++)
				index.put(states.get(s), s);
		}
		setStateLookUp(state -> index.getOrDefault(state, -1));
	}

	@Override
	public CSG<Double> getModel()
	{
		return model;
	}

	@Override
	public int getNumStates()
	{
		return model.getNumStates();
	}
	
	@Override
	public Memory memory()
	{
		return Memory.NONE;
	}

	@Override
	public boolean isRandomised()
	{
		return true;
	}

	/** The players in C (derived from C's actions, set by setCoalitionActions) */
	protected BitSet coalitionPlayers()
	{
		BitSet players = new BitSet();
		if (coalitionActions != null) {
			for (int p = 0; p < model.getNumPlayers(); p++)
				if (coalitionActions.get(model.getIdleForPlayer(p)))
					players.set(p);
		}
		return players;
	}

	/** C's part of a joint action (the other players' entries undefined) */
	protected JointAction partOfC(JointAction joint)
	{
		BitSet players = coalitionPlayers();
		JointAction part = new JointAction(joint.size());
		for (int p = players.nextSetBit(0); p >= 0 && p < joint.size(); p = players.nextSetBit(p + 1))
			part.set(p, joint.get(p));
		return part;
	}

	/** C's action given as player action indices (idle actions included), as a partial joint action */
	protected JointAction partOfC(BitSet actions)
	{
		int numPlayers = model.getNumPlayers();
		JointAction part = new JointAction(numPlayers);
		BitSet players = coalitionPlayers();
		for (int p = players.nextSetBit(0); p >= 0; p = players.nextSetBit(p + 1)) {
			int idle = model.getIdleForPlayer(p);
			if (actions.get(idle)) {
				part.set(p, JointAction.IDLE_ACTION);
				continue;
			}
			BitSet own = (BitSet) model.getIndexes()[p].clone();
			own.and(actions);
			int i = own.nextSetBit(0);
			part.set(p, i > 0 ? model.getActions().get(i - 1) : null);
		}
		return part;
	}

	/** Joint actions may also be given as arrays of (1-indexed) player action indices (-1 for idle), as by model generators */
	protected JointAction toJointAction(Object act)
	{
		if (act instanceof int[])
			return new JointAction((int[]) act, model.getActions());
		return act instanceof JointAction ? (JointAction) act : null;
	}

	/**
	 * The decision in state s: a distribution over C's (partial) joint actions, or UNDEFINED (states already decided,
	 * or where no strategy was computed). The memory m is ignored (memoryless).
	 */
	@Override
	public Object getChoiceAction(int s, int m)
	{
		if (s < 0 || yes.get(s) || no.get(s) || inf.get(s) || coalitionActions == null)
			return UNDEFINED;
		Map<BitSet, Double> strat = csgchoices.get(0).get(0).get(s);
		if (strat == null)
			return UNDEFINED;
		Distribution<Double> d = new Distribution<>();
		List<JointAction> parts = new ArrayList<>();
		for (Map.Entry<BitSet, Double> e : strat.entrySet()) {
			if (e.getValue() > 0.0) {
				d.add(parts.size(), e.getValue());
				parts.add(partOfC(e.getKey()));
			}
		}
		return d.isEmpty() ? UNDEFINED : DistributionOver.create(d, parts::get);
	}

	/** The probability of C's part of a joint action (or of a partial joint action for C) */
	@SuppressWarnings("unchecked")
	@Override
	public Double getChoiceActionProbability(Object decision, Object act)
	{
		JointAction joint = toJointAction(act);
		if (!(decision instanceof DistributionOver) || joint == null)
			return 0.0;
		return ((DistributionOver<Double, Object>) decision).getProbability(partOfC(joint));
	}

	@Override
	public boolean isActionChosen(Object decision, Object act)
	{
		return getChoiceActionProbability(decision, act) > 0.0;
	}

	/**
	 * The choice picked in state s, if determined: C's decision deterministic and a single joint action matching it
	 * (i.e. the other players have a single action); otherwise -1.
	 */
	@Override
	public int getChoiceIndex(int s, int m)
	{
		Object decision = getChoiceAction(s, m);
		if (decision == UNDEFINED)
			return -1;
		int found = -1;
		for (int i = 0; i < model.getNumChoices(s); i++) {
			double p = getChoiceActionProbability(decision, model.getAction(s, i));
			if (p > 0.0) {
				if (p < 1.0 || found >= 0)
					return -1;
				found = i;
			}
		}
		return found;
	}

	@Override
	public prism.Model<Double> constructInducedModel(StrategyExportOptions options) throws PrismException
	{
		// the other players' choices remain nondeterministic: an MDP (in either mode)
		return buildInducedMDP(0).mdp;
	}

	@Override
	public void exportActions(PrismLog out, StrategyExportOptions options) throws PrismException
	{
		List<State> states = model.getStatesList();
		boolean showStates = options.getShowStates() && states != null;
		for (int s = 0; s < model.getNumStates(); s++) {
			Object decision = getChoiceAction(s, -1);
			if (decision != UNDEFINED)
				out.println((showStates ? states.get(s) : s) + "=" + decision);
		}
	}

	@Override
	public void exportIndices(PrismLog out, StrategyExportOptions options) throws PrismException
	{
		throw new PrismNotSupportedException("Zero-sum CSG strategies fix only the coalition's part of each choice, so cannot be exported as choice indices");
	}

	@Override
	public void exportInducedModel(PrismLog out, StrategyExportOptions options) throws PrismException
	{
		buildInducedMDP(0).mdp.exportToPrismExplicitTra(out, options.getModelPrecision());
	}

	@Override
	public void exportDotFile(PrismLog out, StrategyExportOptions options) throws PrismException
	{
		exportZeroSumStrategy(out);
	}
	
	/** Whether the strategy was computed for the complement of the path formula asked for (e.g. F !a for G a) */
	protected boolean complemented = false;

	/**
	 * Sets whether the strategy was computed for the complement of the path formula asked for (e.g. F !a for G a,
	 * whose probabilities are subtracted from 1): the states labelled Sat/Unsat in the export are then swapped.
	 */
	public void setComplemented(boolean complemented)
	{
		this.complemented = complemented;
	}

	/** C's actions (player action indices, idle actions included): the coalition whose strategy this is */
	protected BitSet coalitionActions = null;

	/**
	 * Sets the actions (player action indices, idle actions included) of the coalition C whose strategy this is;
	 * required to build the MDP induced by the strategy (where the other players' choices remain nondeterministic).
	 */
	public void setCoalitionActions(BitSet coalitionActions)
	{
		this.coalitionActions = coalitionActions;
	}

	/** The MDP induced by the strategy of C (reachable part), with, per state and choice, the joint choices of the CSG and their weights */
	protected static class InducedMDP
	{
		MDPSimple<Double> mdp = new MDPSimple<>();
		List<Integer> stateOf = new ArrayList<>();
		/** Per MDP state and choice: CSG choice index -> probability of C's part */
		List<List<Map<Integer, Double>>> weights = new ArrayList<>();
	}

	/**
	 * Builds the MDP induced by the strategy of C: in each state, one choice per action b of the other players, whose
	 * distribution is the mixture, under C's local strategy, of the transitions of the joint actions (a, b). States
	 * already decided (yes, no, inf) are absorbing.
	 * @param k Iteration of the strategy (0 for unbounded properties)
	 */
	protected InducedMDP buildInducedMDP(int k) throws PrismException
	{
		if (coalitionActions == null)
			throw new PrismException("The coalition's actions are needed to build the MDP induced by its strategy");
		InducedMDP ind = new InducedMDP();
		ind.mdp.setVarList(model.getVarList());
		List<State> statelist = new ArrayList<>();
		Map<Integer, Integer> index = new HashMap<>();
		ArrayDeque<Integer> todo = new ArrayDeque<>();
		int s0 = model.getFirstInitialState();
		ind.mdp.addInitialState(inducedNode(ind, statelist, index, todo, s0));
		while (!todo.isEmpty()) {
			int x = todo.poll();
			int s = ind.stateOf.get(x);
			List<Map<Integer, Double>> ws = new ArrayList<>();
			ind.weights.add(x, ws);
			String end = yes.get(s) ? (complemented ? "Unsat" : "Sat") : no.get(s) ? (complemented ? "Sat" : "Unsat") : inf.get(s) ? "Infinity" : null;
			Map<BitSet, Double> strat = csgchoices.get(0).get(k).get(s);
			if (end != null || strat == null) {
				Distribution<Double> d = new Distribution<>();
				d.add(x, 1.0);
				ind.mdp.addActionLabelledChoice(x, d, end != null ? end : "undefined");
				ws.add(new HashMap<>());
				continue;
			}
			// group the joint choices by the other players' part b
			Map<BitSet, Map<Integer, Double>> byOthers = new LinkedHashMap<>();
			for (int t = 0; t < model.getNumChoices(s); t++) {
				BitSet joint = new BitSet();
				int[] indexes = model.getIndexes(s, t);
				for (int q = 0; q < indexes.length; q++)
					joint.set(indexes[q] > 0 ? indexes[q] : model.getIdles()[q]);
				BitSet a = (BitSet) joint.clone();
				a.and(coalitionActions);
				BitSet b = (BitSet) joint.clone();
				b.andNot(coalitionActions);
				Double p = strat.get(a);
				byOthers.computeIfAbsent(b, __ -> new LinkedHashMap<>());
				if (p != null && p > 0.0)
					byOthers.get(b).put(t, p);
			}
			String cpart = mixtureString(strat);
			for (Map.Entry<BitSet, Map<Integer, Double>> e : byOthers.entrySet()) {
				Distribution<Double> d = new Distribution<>();
				for (Map.Entry<Integer, Double> w : e.getValue().entrySet()) {
					for (Iterator<Map.Entry<Integer, Double>> it = model.getTransitionsIterator(s, w.getKey()); it.hasNext();) {
						Map.Entry<Integer, Double> tr = it.next();
						d.add(inducedNode(ind, statelist, index, todo, tr.getKey()), w.getValue() * tr.getValue());
					}
				}
				if (d.isEmpty())
					continue;
				ind.mdp.addActionLabelledChoice(x, d, cpart + " -- " + actionsString(e.getKey()));
				ws.add(e.getValue());
			}
		}
		ind.mdp.setStatesList(statelist);
		return ind;
	}

	private int inducedNode(InducedMDP ind, List<State> statelist, Map<Integer, Integer> index, ArrayDeque<Integer> todo, int s)
	{
		Integer x = index.get(s);
		if (x == null) {
			x = ind.mdp.addState();
			index.put(s, x);
			ind.stateOf.add(s);
			statelist.add(model.getStatesList().get(s));
			todo.add(x);
		}
		return x;
	}

	private String actionsString(BitSet acts)
	{
		StringBuilder sb = new StringBuilder();
		for (int i = acts.nextSetBit(0); i >= 0; i = acts.nextSetBit(i + 1))
			sb.append("[").append(model.getActions().get(i - 1)).append("]");
		return sb.toString();
	}

	/** C's local strategy, e.g. "0.5: [a1] + 0.5: [b1]" */
	private String mixtureString(Map<BitSet, Double> strat)
	{
		StringBuilder sb = new StringBuilder();
		for (Map.Entry<BitSet, Double> e : strat.entrySet()) {
			if (e.getValue() <= 0.0)
				continue;
			if (sb.length() > 0)
				sb.append(" + ");
			sb.append(e.getValue()).append(": ").append(actionsString(e.getKey()));
		}
		return sb.toString();
	}

	/** Exports the MDP induced by the strategy of C (choices labelled "C's mixture -- others' action") */
	public void exportZeroSumStrategy(PrismLog out) throws PrismException
	{
		MDPSimple<Double> mdp = buildInducedMDP(0).mdp;
		mdp.exportToDotFile(out, null, true);
		out.print("\n/*");
		out.print("\n -- Transitions --  \n");
		mdp.exportToPrismExplicitTra(out);
		out.print("\n -- States --  \n");
		mdp.exportStates(0, mdp.getVarList(), out);
		out.print("*/\n");
	}

	/**
	 * The value the strategy of C guarantees in the initial state (unbounded objectives): the optimal value of the
	 * other players in the induced MDP (minimising if C maximises, and vice versa), for reaching the yes states or,
	 * with rewards, the expected reward to reach them.
	 * @param rewards C's rewards (null for probabilities)
	 * @param min Whether C minimises
	 */
	public double guaranteedValue(MDPModelChecker mc, CSGRewards<Double> rewards, boolean min) throws PrismException
	{
		InducedMDP ind = buildInducedMDP(0);
		int num = ind.mdp.getNumStates();
		BitSet target = new BitSet();
		for (int x = 0; x < num; x++)
			if (yes.get(ind.stateOf.get(x)))
				target.set(x);
		if (rewards == null)
			return mc.computeReachProbs(ind.mdp, target, !min).soln[0];
		MDPRewardsSimple<Double> rew = new MDPRewardsSimple<>(num);
		for (int x = 0; x < num; x++) {
			int s = ind.stateOf.get(x);
			rew.setStateReward(x, rewards.getStateReward(s));
			List<Map<Integer, Double>> ws = ind.weights.get(x);
			for (int i = 0; i < ws.size(); i++) {
				double r = 0.0;
				for (Map.Entry<Integer, Double> w : ws.get(i).entrySet())
					r += w.getValue() * rewards.getTransitionReward(s, w.getKey());
				rew.setTransitionReward(x, i, r);
			}
		}
		return mc.computeReachRewards(ind.mdp, rew, target, !min).soln[0];
	}

	@Override
	public void clear() {
	}
	
	@Override
	public String toString() {
		String[] action;
		String joint, label;
		int c, i, k, p, s;
		boolean chck;
		label = "CSG:";
		k = 0;
		if (csgchoices != null) {
			/*
			for (List<List<Map<BitSet, Double>>> entry : csgchoices) {
				s += entry.toString();
				s += "\n";
			}
			*/
			for (s = 0; s < model.getNumStates(); s++) {
				action = new String[csgchoices.size()];
				chck = true;
				for (p = 0; p < csgchoices.size(); p++) {
					chck = chck && csgchoices.get(p).get(k).get(s) != null;
					action[p] = "";
				}
				label += " $ s:" + s + " -> ";
				for (p = 0; p < csgchoices.size(); p++) {
					if (csgchoices.get(p).get(k).get(s) != null) {
						c = csgchoices.get(p).get(k).get(s).keySet().size();
						for (BitSet act : csgchoices.get(p).get(k).get(s).keySet()) {
							joint = "";
							for (i = act.nextSetBit(0); i >= 0; i = act.nextSetBit(i + 1)) {
								joint += "[" + model.getActions().get(i - 1) + "]";
							}
							c--;
							action[p] += csgchoices.get(p).get(k).get(s).get(act) +": " + joint + ((c > 0)? " + " : ""); 
						}
						label += (p + 1 < csgchoices.size())? action[p] + " -- " : action[p];
					}	
				}
			}
			return label;
		}
		else 
			return label;
	}
}
