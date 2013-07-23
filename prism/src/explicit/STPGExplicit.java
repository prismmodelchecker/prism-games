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
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import prism.ModelType;
import prism.PrismException;
import explicit.rewards.MDPRewards;
import explicit.rewards.STPGRewards;

/**
 * Simple explicit-state representation of a stochastic two-player game
 * (STPG). States can be labelled arbitrarily with player 1 player 2.
 */
public class STPGExplicit extends MDPSimple implements STPG
{
	// state labels
	public static final int PLAYER_1 = 1;
	public static final int PLAYER_2 = 2;

	protected List<Integer> stateLabels;

	public STPGExplicit()
	{
		super();

		// initialising state labels
		stateLabels = new ArrayList<Integer>(numStates);
	}

	public STPGExplicit(int n)
	{
		super(n);

		// initialising state labels
		stateLabels = new ArrayList<Integer>(n);
	}

	/**
	 * Construct an STPG from an existing one and a state index permutation,
	 * i.e. in which state index i becomes index permut[i].
	 */
	public STPGExplicit(STPGExplicit stpg, int permut[])
	{
		super(stpg, permut);
		stateLabels = new ArrayList<Integer>(numStates);
		// Create blank array of correct size
		for (int i = 0; i < numStates; i++) {
			stateLabels.add(0);
		}
		// Copy permuted player info
		for (int i = 0; i < numStates; i++) {
			stateLabels.set(permut[i], stpg.stateLabels.get(i));
		}
	}

	/**
	 * Copy constructor
	 */
	public STPGExplicit(STPGExplicit stpg)
	{
		super(stpg);
		stateLabels = new ArrayList<Integer>(stpg.stateLabels);
	}

	/**
	 * Adds one state, assigned to player 1
	 */
	@Override
	public int addState()
	{
		return addState(PLAYER_1);
	}

	/**
	 * Adds specified number of states all assigned to player 1
	 */
	@Override
	public void addStates(int numToAdd)
	{
		super.addStates(numToAdd);
		for (int i = 0; i < numToAdd; i++)
			stateLabels.add(PLAYER_1);
	}

	/**
	 * Adds state assigned to the specified player
	 * 
	 * @param player state owner
	 * @return state id
	 */
	public int addState(int player)
	{
		checkPlayer(player);

		super.addStates(1);
		stateLabels.add(player);

		//System.out.println("State " + (numStates - 1) + " player " + player);

		return numStates - 1;
	}

	/**
	 * Adds the number of states the same as number of Integer in the list, each
	 * assigned to the corresponding player
	 * 
	 * @param players list of players (to which corresponding state belongs)
	 */
	public void addStates(List<Integer> players)
	{
		checkPlayers(players);

		super.addStates(players.size());
		stateLabels.addAll(players);
	}

	/**
	 * labels the given state with the given player
	 * 
	 * @param s state
	 * @param player player
	 */
	public void setPlayer(int s, int player)
	{
		checkPlayer(player);
		if (s < stateLabels.size())
			stateLabels.set(s, player);
	}

	// Accessors (for Model)

	@Override
	public ModelType getModelType()
	{
		return ModelType.STPG;
	}

	// Accessors (for STPG)

	@Override
	public int getPlayer(int s)
	{
		return stateLabels.get(s);
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

				if (getPlayer(i) == PLAYER_1)
					forall = forall1;
				else if (getPlayer(i) == PLAYER_2)
					forall = forall2;

				c = 0;
				b1 = forall; // there exists or for all
				for (Distribution distr : trans.get(i)) {
					// ignoring the choice if it is disabled
					if (someChoicesDisabled && disabledChoices.containsKey(i)
							&& disabledChoices.get(i).get(c++) == true)
						continue;
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

	@Override
	public void prob1step(BitSet subset, BitSet u, BitSet v, boolean forall1, boolean forall2, BitSet result)
	{
		int i, c;
		boolean b1, b2;
		boolean forall = false;

		for (i = 0; i < numStates; i++) {
			if (subset.get(i)) {

				if (getPlayer(i) == PLAYER_1)
					forall = forall1;
				else if (getPlayer(i) == PLAYER_2)
					forall = forall2;

				c = 0;
				b1 = forall; // there exists or for all
				for (Distribution distr : trans.get(i)) {
					// ignoring the choice if it is disabled
					if (someChoicesDisabled && disabledChoices.containsKey(i)
							&& disabledChoices.get(i).get(c++) == true)
						continue;
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

	// TODO fix the method
	@Override
	public void mvMultMinMax(double vect[], boolean min1, boolean min2, double result[], BitSet subset,
			boolean complement, int adv[])
	{
		int s;
		boolean min = false;
		// Loop depends on subset/complement arguments
		if (subset == null) {
			for (s = 0; s < numStates; s++) {
				if (getPlayer(s) == PLAYER_1)
					min = min1;
				else if (getPlayer(s) == PLAYER_2)
					min = min2;

				result[s] = mvMultMinMaxSingle(s, vect, min, adv);
			}
		} else if (complement) {
			for (s = subset.nextClearBit(0); s < numStates; s = subset.nextClearBit(s + 1)) {
				if (getPlayer(s) == PLAYER_1)
					min = min1;
				else if (getPlayer(s) == PLAYER_2)
					min = min2;

				result[s] = mvMultMinMaxSingle(s, vect, min, adv);
			}
		} else {
			for (s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s + 1)) {
				if (getPlayer(s) == PLAYER_1)
					min = min1;
				else if (getPlayer(s) == PLAYER_2)
					min = min2;
				result[s] = mvMultMinMaxSingle(s, vect, min, adv);
			}
		}
	}

	@Override
	public double mvMultMinMaxSingle(int s, double vect[], boolean min1, boolean min2)
	{
		boolean min = getPlayer(s) == PLAYER_1 ? min1 : getPlayer(s) == PLAYER_2 ? min2 : false;
		return mvMultMinMaxSingle(s, vect, min, null);
	}

	@Override
	public List<Integer> mvMultMinMaxSingleChoices(int s, double vect[], boolean min1, boolean min2, double val)
	{
		boolean min = getPlayer(s) == PLAYER_1 ? min1 : getPlayer(s) == PLAYER_2 ? min2 : false;
		return mvMultMinMaxSingleChoices(s, vect, min, val);
	}

	@Override
	public double mvMultGSMinMax(double vect[], boolean min1, boolean min2, BitSet subset, boolean complement,
			boolean absolute)
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
		boolean min = getPlayer(s) == PLAYER_1 ? min1 : getPlayer(s) == PLAYER_2 ? min2 : false;
		return mvMultJacMinMaxSingle(s, vect, min, null);
	}

	@Override
	public void mvMultRewMinMax(double vect[], STPGRewards rewards, boolean min1, boolean min2, double result[],
			BitSet subset, boolean complement, int adv[])
	{
		int s;
		boolean min = false;
		MDPRewards mdpRewards = rewards.buildMDPRewards();
		// Loop depends on subset/complement arguments
		if (subset == null) {
			for (s = 0; s < numStates; s++) {
				if (getPlayer(s) == PLAYER_1)
					min = min1;
				else if (getPlayer(s) == PLAYER_2)
					min = min2;
				result[s] = mvMultRewMinMaxSingle(s, vect, mdpRewards, min, adv, 1.0);
			}
		} else if (complement) {
			for (s = subset.nextClearBit(0); s < numStates; s = subset.nextClearBit(s + 1)) {
				if (getPlayer(s) == PLAYER_1)
					min = min1;
				else if (getPlayer(s) == PLAYER_2)
					min = min2;
				result[s] = mvMultRewMinMaxSingle(s, vect, mdpRewards, min, adv, 1.0);
			}
		} else {
			for (s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s + 1)) {
				if (getPlayer(s) == PLAYER_1)
					min = min1;
				else if (getPlayer(s) == PLAYER_2)
					min = min2;
				result[s] = mvMultRewMinMaxSingle(s, vect, mdpRewards, min, adv, 1.0);
			}
		}
	}
	
	@Override
	public void mvMultRewMinMax(double vect[], STPGRewards rewards, boolean min1, boolean min2, double result[],
			BitSet subset, boolean complement, int adv[], double disc)
	{
		int s;
		boolean min = false;
		MDPRewards mdpRewards = rewards.buildMDPRewards();
		// Loop depends on subset/complement arguments
		if (subset == null) {
			for (s = 0; s < numStates; s++) {
				if (getPlayer(s) == PLAYER_1)
					min = min1;
				else if (getPlayer(s) == PLAYER_2)
					min = min2;
				result[s] = mvMultRewMinMaxSingle(s, vect, mdpRewards, min, adv, disc);
			}
		} else if (complement) {
			for (s = subset.nextClearBit(0); s < numStates; s = subset.nextClearBit(s + 1)) {
				if (getPlayer(s) == PLAYER_1)
					min = min1;
				else if (getPlayer(s) == PLAYER_2)
					min = min2;
				result[s] = mvMultRewMinMaxSingle(s, vect, mdpRewards, min, adv, disc);
			}
		} else {
			for (s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s + 1)) {
				if (getPlayer(s) == PLAYER_1)
					min = min1;
				else if (getPlayer(s) == PLAYER_2)
					min = min2;
				result[s] = mvMultRewMinMaxSingle(s, vect, mdpRewards, min, adv, disc);
			}
		}
	}

	@Override
	public double mvMultRewMinMaxSingle(int s, double vect[], STPGRewards rewards, boolean min1, boolean min2,
			int adv[])
	{
		MDPRewards mdpRewards = rewards.buildMDPRewards();
		boolean min = getPlayer(s) == PLAYER_1 ? min1 : getPlayer(s) == PLAYER_2 ? min2 : false;
		return mvMultRewMinMaxSingle(s, vect, mdpRewards, min, adv);
	}

	@Override
	public List<Integer> mvMultRewMinMaxSingleChoices(int s, double vect[], STPGRewards rewards, boolean min1,
			boolean min2, double val)
	{
		MDPRewards mdpRewards = rewards.buildMDPRewards();
		boolean min = getPlayer(s) == PLAYER_1 ? min1 : getPlayer(s) == PLAYER_2 ? min2 : false;
		return mvMultRewMinMaxSingleChoices(s, vect, mdpRewards, min, val);
	}

	// Accessors (other)

	/** Checks whether the given player is valid and throws exception otherwise **/
	private void checkPlayer(int player)
	{
		switch (player) {
		case PLAYER_1:
			return;
		case PLAYER_2:
			return;
		}
		throw new IllegalArgumentException("Player " + player + " is undefined!");
	}

	/**
	 * Checks whether every player in the list is valid and throws exception
	 * otherwise
	 **/
	private void checkPlayers(List<Integer> players)
	{
		for (Integer p : players)
			checkPlayer(p);
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
			s += i + "(P-" + getPlayer(i) + "): ";
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

			// ignoring the choice if it is disabled
			if (someChoicesDisabled && disabledChoices.containsKey(s) && disabledChoices.get(s).get(c++) == true)
				continue;

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
				if (adv != null)
				{
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
	
	public static void main(String[] args)
	{
		STPGModelChecker mc;
		STPGExplicit stpg;
		DistributionSet set;
		Distribution distr;
		ModelCheckerResult res;
		BitSet target;

		// Simple example: Create and solve the stochastic game from:
		// Mark Kattenbelt, Marta Kwiatkowska, Gethin Norman, David Parker
		// A Game-based Abstraction-Refinement Framework for Markov Decision
		// Processes
		// Formal Methods in System Design 36(3): 246-280, 2010

		try {
			// Build game
			stpg = new STPGExplicit();

			/*
			// state 0
			stpg.addState(PLAYER_1);
			set = stpg.newDistributionSet(null);
			distr = new Distribution();
			distr.set(1, 1.0);
			set.add(distr);
			stpg.addDistributionSet(0, set);

			// state 1
			stpg.addState(PLAYER_2);
			set = stpg.newDistributionSet(null);
			distr = new Distribution();
			distr.set(2, 1.0);
			set.add(distr);
			stpg.addDistributionSet(1, set);

			// state 2
			stpg.addState(PLAYER_1);
			set = stpg.newDistributionSet(null);
			distr = new Distribution();
			distr.set(3, 1.0);
			set.add(distr);
			distr = new Distribution();
			distr.set(4, 1.0);
			set.add(distr);
			stpg.addDistributionSet(2, set);

			// state 3
			stpg.addState(PLAYER_2);
			set = stpg.newDistributionSet(null);
			distr = new Distribution();
			distr.set(5, 1.0);
			set.add(distr);
			distr = new Distribution();
			distr.set(2, 1.0);
			set.add(distr);
			stpg.addDistributionSet(3, set);

			// state 4
			stpg.addState(PLAYER_2);
			set = stpg.newDistributionSet(null);
			distr = new Distribution();
			distr.set(6, 1.0);
			set.add(distr);
			distr = new Distribution();
			distr.set(5, 0.5);
			distr.set(6, 0.5);
			set.add(distr);
			stpg.addDistributionSet(4, set);

			// state 5
			stpg.addState(PLAYER_1);
			set = stpg.newDistributionSet(null);
			distr = new Distribution();
			distr.set(7, 1.0);
			set.add(distr);
			stpg.addDistributionSet(5, set);

			// state 6
			stpg.addState(PLAYER_1);
			set = stpg.newDistributionSet(null);
			distr = new Distribution();
			distr.set(8, 1.0);
			set.add(distr);
			stpg.addDistributionSet(6, set);

			// state 7
			stpg.addState(PLAYER_2);
			set = stpg.newDistributionSet(null);
			distr = new Distribution();
			distr.set(7, 1.0);
			set.add(distr);
			stpg.addDistributionSet(7, set);

			// state 8
			stpg.addState(PLAYER_2);
			set = stpg.newDistributionSet(null);
			distr = new Distribution();
			distr.set(8, 1.0);
			set.add(distr);
			stpg.addDistributionSet(8, set);
			*/

			// Print game
			System.out.println(stpg);

			// Model check
			mc = new STPGModelChecker(null);
			// mc.setVerbosity(2);
			target = new BitSet();
			target.set(8);
			stpg.exportToDotFile("stpg.dot", target);
			System.out.println("min min: " + mc.computeReachProbs(stpg, target, true, true).soln[0]);
			System.out.println("max min: " + mc.computeReachProbs(stpg, target, false, true).soln[0]);
			System.out.println("min max: " + mc.computeReachProbs(stpg, target, true, false).soln[0]);
			System.out.println("max max: " + mc.computeReachProbs(stpg, target, false, false).soln[0]);
		} catch (PrismException e) {
			System.out.println(e);
		}
	}
}
