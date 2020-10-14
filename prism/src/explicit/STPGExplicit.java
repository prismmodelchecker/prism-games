// ==============================================================================
//	
// Copyright (c) 2002-
// Authors:
// * Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford)
//	
// ------------------------------------------------------------------------------
//	
// This file is part of PRISM.
//	
// PRISM is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//	
// PRISM is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//	
// You should have received a copy of the GNU General Public License
// along with PRISM; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//	
// ==============================================================================

package explicit;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import explicit.rewards.MDPRewards;
import explicit.rewards.STPGRewards;
import prism.ModelType;
import prism.PrismException;
import prism.PrismLog;
import prism.PrismUtils;

/**
 * Simple explicit-state representation of a (turn-based) stochastic two-player game (STPG).
 */
public class STPGExplicit extends MDPSimple implements STPG
{
	/** Which player owns each state,fl i.e. stateOwners[i] is owned by player i (1 or 2) */
	protected List<Integer> stateOwners;

	// Constructors

	/**
	 * Constructor: empty STPG.
	 */
	public STPGExplicit()
	{
		super();
		stateOwners = new ArrayList<Integer>(0);
	}

	/**
	 * Constructor: new STPG with fixed number of states.
	 */
	public STPGExplicit(int numStates)
	{
		super(numStates);
		stateOwners = new ArrayList<Integer>(numStates);
	}

	/**
	 * Construct an STPG from an existing one and a state index permutation,
	 * i.e. in which state index i becomes index permut[i].
	 */
	public STPGExplicit(STPGExplicit stpg, int permut[])
	{
		super(stpg, permut);
		stateOwners = new ArrayList<Integer>(numStates);
		// Create blank array of correct size
		for (int i = 0; i < numStates; i++) {
			stateOwners.add(0);
		}
		// Copy permuted player info
		for (int i = 0; i < numStates; i++) {
			stateOwners.set(permut[i], stpg.stateOwners.get(i));
		}
	}

	/**
	 * Copy constructor
	 */
	public STPGExplicit(STPGExplicit stpg)
	{
		super(stpg);
		stateOwners = new ArrayList<Integer>(stpg.stateOwners);
	}

	// Mutators (for ModelSimple)

	/**
	 * Add a new (player 1) state and return its index.
	 */
	@Override
	public int addState()
	{
		return addState(1);
	}

	/**
	 * Add multiple new (player 1) states.
	 */
	@Override
	public void addStates(int numToAdd)
	{
		super.addStates(numToAdd);
		for (int i = 0; i < numToAdd; i++)
			stateOwners.add(1);
	}

	/**
	 * Add a new (player {@code p}) state and return its index. For an STPG, {@code p} should be 1 or 2. 
	 * @param p Player who owns the new state.
	 */
	public int addState(int p)
	{
		super.addStates(1);
		stateOwners.add(p);
		return numStates - 1;
	}

	/**
	 * Add multiple new states, with owners as given in the list {@code p}
	 * (the number of states to add is dictated by the length of the list).
	 * For an STPG, player indices should be 1 or 2.
	 * @param p List of players owning each new state
	 */
	public void addStates(List<Integer> p)
	{
		super.addStates(p.size());
		stateOwners.addAll(p);
	}

	/**
	 * Set player {@code p} to own state {@code s}. For an STPG, {@code} should be 1 or 2.
	 * It is not checked whether {@code s} or {@code p} are in the correct range.
	 */
	public void setPlayer(int s, int p)
	{
		stateOwners.set(s, p);
	}

	// Accessors (for Model)

	@Override
	public ModelType getModelType()
	{
		// Resolve conflict: STPG interface does not (currently) extend MDP  
		return STPG.super.getModelType();
	}

	@Override
	public void exportToPrismExplicitTra(PrismLog out)
	{
		// Resolve conflict: STPG interface does not (currently) extend MDP  
		STPG.super.exportToPrismExplicitTra(out);
	}

	@Override
	public void exportToPrismLanguage(final String filename) throws PrismException
	{
		// Resolve conflict: STPG interface does not (currently) extend MDP  
		STPG.super.exportToPrismLanguage(filename);
	}

	@Override
	public String infoString()
	{
		// Resolve conflict: STPG interface does not (currently) extend MDP  
		return STPG.super.infoString();
	}

	@Override
	public String infoStringTable()
	{
		// Resolve conflict: STPG interface does not (currently) extend MDP  
		return STPG.super.infoStringTable();
	}

	// Accessors (for STPG)

	@Override
	public int getPlayer(int s)
	{
		return stateOwners.get(s);
	}

	@Override
	public int getNumTransitions(int s, int i)
	{
		return super.getNumTransitions(s, i);
	}

	@Override
	public Iterator<Entry<Integer, Double>> getTransitionsIterator(int s, int i)
	{
		return super.getTransitionsIterator(s, i);
	}

	@Override
	public boolean isChoiceNested(int s, int i)
	{
		// No nested choices
		return false;
	}

	@Override
	public int getNumNestedChoices(int s, int i)
	{
		// No nested choices
		return 0;
	}

	@Override
	public Object getNestedAction(int s, int i, int j)
	{
		// No nested choices
		return null;
	}

	@Override
	public int getNumNestedTransitions(int s, int i, int j)
	{
		// No nested choices
		return 0;
	}

	@Override
	public Iterator<Entry<Integer, Double>> getNestedTransitionsIterator(int s, int i, int j)
	{
		// No nested choices
		return null;
	}

	@Override
	public void prob0step(BitSet subset, BitSet u, boolean forall1, boolean forall2, BitSet result)
	{
		int i, c;
		boolean b1, b2;
		boolean forall = false;

		for (i = 0; i < numStates; i++) {
			if (subset.get(i)) {
				forall = (getPlayer(i) == 1) ? forall1 : forall2;
				c = 0;
				b1 = forall; // there exists or for all
				for (Distribution distr : trans.get(i)) {
					b2 = distr.containsOneOf(u);
					if (forall) {
						if (!b2) {
							b1 = false;
							continue;
						}
					} else {
						if (b2) {
							b1 = true;
							continue;
						}
					}
				}
				result.set(i, b1);
			}
		}
	}
	
	/**
	 * Returns true is state {@code s} has a selfloop with probability 1
	 * @param s The state to be tested.
	 * @return True if probability 1 selfloop.
	 */
	public boolean hasProb1Selfloop(int s)
	{
		int c = 0;
		
		for (Distribution distr : trans.get(s)) {
			if(distr.getSupport().size() == 1 && distr.getSupport().contains(s))
				return true;
		}
		
		return false; // no probability 1 selfloop found
	}
	
	public void reachpositivestep(BitSet u, boolean forall1, boolean forall2, BitSet result)
	{
		int i, c;
		boolean forall = false;
		boolean first;
		Set<Integer> u1;
		
		for (i = 0; i < numStates; i++) {
			if (u.get(i)) {
				forall = (getPlayer(i) == 1) ? forall1 : forall2;
				c = 0;
				u1 = null; // reach in one step
				first = true;
				for (Distribution distr : trans.get(i)) {
					if(first) { 
						u1 = new HashSet<Integer>(distr.getSupport()); // put all successors in reachable states
						first = false;
					} else if (!first && forall) {
						u1.retainAll(distr.getSupport()); // intersect
					} else if (!first & !forall) {
						u1.addAll(distr.getSupport()); // union
					}
				}
				for(int r : u1)
					result.set(r, true);
			}
		}
	}

    /**
     * @param u The subtree so far
     * @param closedPlayer Player for which subtree is closed
     * @param result The subtree after extending
     **/
        public void subtreeStep(BitSet u, int closedPlayer, BitSet result)
	{
		for (int i = 0; i < numStates; i++) {
		        // go only through states in subtree so far,
		        // and only extend subtree if closed for that player,
		        // or if the state has only one choice that is enabled
		        boolean jump = (getNumChoices(i) == 1) && getPlayer(i) != closedPlayer;
		        if (u.get(i) && (getPlayer(i) == closedPlayer || jump)) {
				int c = 0;
				for (Distribution distr : trans.get(i)) {
					for(int r : distr.getSupport()) { // add all successors (no matter which player)
					    result.set(r, true);
					}
				}
			}
		}
	}

	@Override
	public void prob1step(BitSet subset, BitSet u, BitSet v, boolean forall1, boolean forall2, BitSet result)
	{
		int i, c;
		boolean b1, b2;
		boolean forall = false;

		for (i = 0; i < numStates; i++) {
			if (subset.get(i)) {
				forall = (getPlayer(i) == 1) ? forall1 : forall2;
				c = 0;
				b1 = forall; // there exists or for all
				for (Distribution distr : trans.get(i)) {
					b2 = distr.containsOneOf(v) && distr.isSubsetOf(u);
					if (forall) {
						if (!b2) {
							b1 = false;
							continue;
						}
					} else {
						if (b2) {
							b1 = true;
							continue;
						}
					}
				}
				result.set(i, b1);
			}
		}
	}

	@Override
	public void mvMultMinMax(double vect[], boolean min1, boolean min2, double result[], BitSet subset, boolean complement, int adv[])
	{
		int s;
		boolean min = false;
		// Loop depends on subset/complement arguments
		if (subset == null) {
			for (s = 0; s < numStates; s++) {
				min = (getPlayer(s) == 1) ? min1 : min2;
				result[s] = mvMultMinMaxSingle(s, vect, min, adv);
			}
		} else if (complement) {
			for (s = subset.nextClearBit(0); s < numStates; s = subset.nextClearBit(s + 1)) {
				min = (getPlayer(s) == 1) ? min1 : min2;
				result[s] = mvMultMinMaxSingle(s, vect, min, adv);
			}
		} else {
			for (s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s + 1)) {
				min = (getPlayer(s) == 1) ? min1 : min2;
				result[s] = mvMultMinMaxSingle(s, vect, min, adv);
			}
		}
	}

	@Override
	public double mvMultMinMaxSingle(int s, double vect[], boolean min1, boolean min2)
	{
		boolean min = (getPlayer(s) == 1) ? min1 : min2;
		return mvMultMinMaxSingle(s, vect, min, null);
	}

	@Override
	public List<Integer> mvMultMinMaxSingleChoices(int s, double vect[], boolean min1, boolean min2, double val)
	{
		boolean min = (getPlayer(s) == 1) ? min1 : min2;
		return mvMultMinMaxSingleChoices(s, vect, min, val);
	}

	@Override
	public double mvMultGSMinMax(double vect[], boolean min1, boolean min2, BitSet subset, boolean complement, boolean absolute)
	{
		int s;
		double d, diff, maxDiff = 0.0;
		// Loop depends on subset/complement arguments
		if (subset == null) {
			for (s = 0; s < numStates; s++) {
				d = mvMultJacMinMaxSingle(s, vect, min1, min2);
				diff = absolute ? (Math.abs(d - vect[s])) : (Math.abs(d - vect[s]) / d);
				maxDiff = diff > maxDiff ? diff : maxDiff;
				vect[s] = d;
			}
		} else if (complement) {
			for (s = subset.nextClearBit(0); s < numStates; s = subset.nextClearBit(s + 1)) {
				d = mvMultJacMinMaxSingle(s, vect, min1, min2);
				diff = absolute ? (Math.abs(d - vect[s])) : (Math.abs(d - vect[s]) / d);
				maxDiff = diff > maxDiff ? diff : maxDiff;
				vect[s] = d;
			}
		} else {
			for (s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s + 1)) {
				d = mvMultJacMinMaxSingle(s, vect, min1, min2);
				diff = absolute ? (Math.abs(d - vect[s])) : (Math.abs(d - vect[s]) / d);
				maxDiff = diff > maxDiff ? diff : maxDiff;
				vect[s] = d;
			}
		}
		return maxDiff;
	}

	@Override
	public double mvMultJacMinMaxSingle(int s, double vect[], boolean min1, boolean min2)
	{
		boolean min = (getPlayer(s) == 1) ? min1 : min2;
		return mvMultJacMinMaxSingle(s, vect, min, null);
	}

	@Override
	public void mvMultRewMinMax(double vect[], STPGRewards rewards, boolean min1, boolean min2, double result[], BitSet subset, boolean complement, int adv[])
	{
		int s;
		boolean min = false;
		MDPRewards mdpRewards = rewards.buildMDPRewards();
		// Loop depends on subset/complement arguments
		if (subset == null) {
			for (s = 0; s < numStates; s++) {
				min = (getPlayer(s) == 1) ? min1 : min2;
				result[s] = mvMultRewMinMaxSingle(s, vect, mdpRewards, min, adv, 1.0);
			}
		} else if (complement) {
			for (s = subset.nextClearBit(0); s < numStates; s = subset.nextClearBit(s + 1)) {
				min = (getPlayer(s) == 1) ? min1 : min2;
				result[s] = mvMultRewMinMaxSingle(s, vect, mdpRewards, min, adv, 1.0);
			}
		} else {
			for (s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s + 1)) {
				min = (getPlayer(s) == 1) ? min1 : min2;
				result[s] = mvMultRewMinMaxSingle(s, vect, mdpRewards, min, adv, 1.0);
			}
		}
	}

	@Override
	public void mvMultRewMinMax(double vect[], STPGRewards rewards, boolean min1, boolean min2, double result[], BitSet subset, boolean complement, int adv[],
			double disc)
	{
		int s;
		boolean min = false;
		MDPRewards mdpRewards = rewards.buildMDPRewards();
		// Loop depends on subset/complement arguments
		if (subset == null) {
			for (s = 0; s < numStates; s++) {
				min = (getPlayer(s) == 1) ? min1 : min2;
				result[s] = mvMultRewMinMaxSingle(s, vect, mdpRewards, min, adv, disc);
			}
		} else if (complement) {
			for (s = subset.nextClearBit(0); s < numStates; s = subset.nextClearBit(s + 1)) {
				min = (getPlayer(s) == 1) ? min1 : min2;
				result[s] = mvMultRewMinMaxSingle(s, vect, mdpRewards, min, adv, disc);
			}
		} else {
			for (s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s + 1)) {
				min = (getPlayer(s) == 1) ? min1 : min2;
				//System.out.printf("s: %s, min1: %s, min2: %s, min: %s, player: %d\n", s, min1, min2, min, getPlayer(s));
				result[s] = mvMultRewMinMaxSingle(s, vect, mdpRewards, min, adv, disc);
			}
		}
	}

	@Override
	public double mvMultRewMinMaxSingle(int s, double vect[], STPGRewards rewards, boolean min1, boolean min2, int adv[])
	{
		MDPRewards mdpRewards = rewards.buildMDPRewards();
		boolean min = (getPlayer(s) == 1) ? min1 : min2;
		return mvMultRewMinMaxSingle(s, vect, mdpRewards, min, adv);
	}

	@Override
	public List<Integer> mvMultRewMinMaxSingleChoices(int s, double vect[], STPGRewards rewards, boolean min1, boolean min2, double val)
	{
		MDPRewards mdpRewards = rewards.buildMDPRewards();
		boolean min = (getPlayer(s) == 1) ? min1 : min2;
		return mvMultRewMinMaxSingleChoices(s, vect, mdpRewards, min, val);
	}

	// Standard methods

	@Override
	public String toString()
	{
		int i, j, n;
		Object o;
		String s = "";
		s = "[ ";
		for (i = 0; i < numStates; i++) {
			if (i > 0)
				s += ", ";
			s += i + "(PP-" + getPlayer(i) + "): ";
			s += "[";
			n = getNumChoices(i);
			for (j = 0; j < n; j++) {
				if (j > 0)
					s += ",";
				o = getAction(i, j);
				if (o != null)
					s += o + ":";
				s += trans.get(i).get(j);
			}
			s += "]";
		}
		s += " ]\n";
		return s;
	}

	/**
	 * Allows discounting
	 * @param s
	 * @param vect
	 * @param mdpRewards
	 * @param min
	 * @param adv
	 * @param disc
	 * @return
	 */
	public double mvMultRewMinMaxSingle(int s, double vect[], MDPRewards mdpRewards, boolean min, int adv[], double disc)
	{
		int j, k, advCh = -1, c;
		double d, prob, minmax;
		boolean first;
		List<Distribution> step;

		c = 0;
		minmax = 0;
		first = true;
		j = -1;
		step = trans.get(s);
		for (Distribution distr : step) {
			j++;
			// Compute sum for this distribution
			d = mdpRewards.getTransitionReward(s, j);

			for (Map.Entry<Integer, Double> e : distr) {
				k = (Integer) e.getKey();
				prob = (Double) e.getValue();
				d += prob * vect[k] * disc;
			}

			// Check whether we have exceeded min/max so far
			if (first || (min && d < minmax) || (!min && d > minmax)) {
				minmax = d;
				// If adversary generation is enabled, remember optimal choice
				if (adv != null) {
					advCh = j;
				}
			}
			first = false;
		}
		// If adversary generation is enabled, store optimal choice
		if (adv != null & !first) {
			// Only remember strictly better choices (required for max)
			if (adv[s] == -1 || (min && minmax < vect[s]) || (!min && minmax > vect[s]) || this instanceof STPG) {
				adv[s] = advCh;
			}
		}

		// Add state reward (doesn't affect min/max)
		minmax += mdpRewards.getStateReward(s);

		return minmax;
	}

	@Override
	public void exportToDotFile(PrismLog out, BitSet mark)
	{
		exportToDotFile(out, mark, false);
	}

	@Override
	public void exportToDotFile(PrismLog out, BitSet mark, boolean states)
	{
		exportToDotFile(out, mark, null, states);
	}

	public void exportToDotFile(PrismLog out, BitSet mark, BitSet players, boolean states)
	{
		int i, j, numChoices;
		String nij;
		Object action;
		out.print("digraph " + getModelType() + " {\nsize=\"8,5\"\nnode [shape=box];\n");
		for (i = 0; i < numStates; i++) {
			String state_label = null;
			if (states) {
				state_label = "label=\"" + statesList.get(i).toString() + "\"";
			}
			int player = (players == null) ? getPlayer(i) : (players.get(i) ? 1 : 2);
			// Player 1 states are diamonds
			if (player == 1 && mark != null && mark.get(i)) {
				if (state_label != null)
					state_label += ", ";
				else
					state_label = "";
				out.print(i + " [" + state_label + "shape=diamond, style=filled, fillcolor=\"#cccccc\"]\n");
			} else if (mark != null && mark.get(i)) {
				if (state_label != null)
					state_label += ", ";
				else
					state_label = "";
				out.print(i + " [" + state_label + "style=filled  fillcolor=\"#cccccc\"]\n");
			} else if (player == 1) {
				if (state_label != null)
					state_label += ", ";
				else
					state_label = "";
				out.print(i + " [" + state_label + "shape=diamond]\n");
			} else if (state_label != null) {
				out.print(i + " [" + state_label + "]\n");
			}
			numChoices = getNumChoices(i);
			for (j = 0; j < numChoices; j++) {
				action = getAction(i, j);
				nij = "n" + i + "_" + j;
				out.print(i + " -> " + nij + " [ arrowhead=none,label=\"" + j);
				if (action != null)
					out.print(":" + action);
				out.print("\" ];\n");
				out.print(nij + " [ shape=point,width=0.1,height=0.1,label=\"\" ];\n");
				Iterator<Map.Entry<Integer, Double>> iter = getTransitionsIterator(i, j);
				while (iter.hasNext()) {
					Map.Entry<Integer, Double> e = iter.next();
					if (PrismUtils.doublesAreEqual(e.getValue(), 1.0)) {
						out.print(nij + " -> " + e.getKey() + ";\n");
					} else {
						out.print(nij + " -> " + e.getKey() + " [ label=\"" + e.getValue() + "\" ];\n");
					}
				}
			}
		}
		out.print("}\n");
	}

	@Override
	public void checkForDeadlocks(BitSet except) throws PrismException
	{
		for (int i = 0; i < numStates; i++) {
			if (trans.get(i).isEmpty() && (except == null || !except.get(i)))
			    throw new PrismException("Game has a deadlock in state " + i +
						     (statesList==null ? "" : ": " + statesList.get(i)));
		}
	}
}
