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
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.Map.Entry;

import explicit.rewards.MDPRewards;
import explicit.rewards.STPGRewards;

import prism.ModelType;

/**
 * Simple explicit-state representation of a stochastic multi-player game
 * (SMG). States can be labelled arbitrarily with player 1..n, player 0
 * has a special purpose of scheduling the moves of other players
 */
public class SMG extends MDPSimple implements STPG
{

	// State labels: states with label i are controlled by player i
	protected List<Integer> stateLabels;
	
	// Set of players which form a coalition
	protected Set<Integer> coalition;

	public SMG()
	{
		super();
		stateLabels = new ArrayList<Integer>(numStates);
	}

	/**
	 * Construct an SMG from an existing one and a state index permutation,
	 * i.e. in which state index i becomes index permut[i].
	 */
	public SMG(SMG smg, int permut[])
	{
		super(smg, permut);
		stateLabels = new ArrayList<Integer>(numStates);
		// Create blank array of correct size
		for (int i = 0; i < numStates; i++) {
			stateLabels.add(0);
		}
		// Copy permuted player info
		for (int i = 0; i < numStates; i++) {
			stateLabels.set(permut[i], smg.stateLabels.get(i));
		}

	}
	
	/**
	 * Returns the list of states that belong to the scheduler
	 * 
	 * @return the list of states that belong to the scheduler
	 */
	public Set<Integer> getSchedulerStates()
	{
		Set<Integer> ret = new HashSet<Integer>();
//		for (int i = 0; i < stateLabels.size(); i++)
//			if (stateLabels.get(i) == 0)
//				ret.add(i);
		return ret;
	}
	
	
	/**
	 * Adds one state, assigned to player 0
	 */
	@Override
	public int addState()
	{
		return addState(0);
	}

	/**
	 * Adds specified number of states all assigned to player 0
	 */
	@Override
	public void addStates(int numToAdd)
	{
		for (int i = 0; i < numToAdd; i++)
			stateLabels.add(0);
	}

	/**
	 * Adds state assigned to the specified player
	 * 
	 * @param player state owner
	 * @return state id
	 */
	public int addState(int player)
	{
		super.addStates(1);
		stateLabels.add(player);

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
		if (s < stateLabels.size())
			stateLabels.set(s, player);
	}
	
	public void setCoalition(Set<Integer> coalition)
	{
		this.coalition = coalition;
	}

	/**
	 * Makes a half-deep (up to one reference level) copy of itself
	 */
	public SMG clone()
	{
		SMG smg = new SMG();
		smg.copyFrom(this);
		smg.actions = new ArrayList<List<Object>>(this.actions);
		smg.allowDupes = this.allowDupes;
		smg.maxNumDistrs = this.maxNumDistrs;
		smg.maxNumDistrsOk = this.maxNumDistrsOk;
		smg.numDistrs = this.numDistrs;
		smg.numTransitions = this.numTransitions;
		smg.stateLabels = new ArrayList<Integer>(this.stateLabels);
		smg.trans = new ArrayList<List<Distribution>>(this.trans);
		return smg;
	}

	// Accessors (for Model)

	@Override
	public ModelType getModelType()
	{
		return ModelType.SMG;
	}

	// Accessors (for STPG)

	@Override
	public int getPlayer(int s)
	{
		return stateLabels.get(s);
	}
	
	@Override
	public Object getAction(int s, int i)
	{
		return super.getAction(s, i);
	}

	@Override
	public int getNumTransitions(int s, int i)
	{
		return super.getNumTransitions(s, i);
	}

	@Override
	public Iterator<Entry<Integer,Double>> getTransitionsIterator(int s, int i)
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
		int i;
		boolean b1, b2;
		boolean forall = false;

		for (i = 0; i < numStates; i++) {
			if (subset.get(i)) {

				if (coalition.contains(stateLabels.get(i)))
					forall = forall1;
				else //if (stateLabels.get(i).equals(2))
					forall = forall2;

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

	@Override
	public void prob1step(BitSet subset, BitSet u, BitSet v, boolean forall1, boolean forall2, BitSet result)
	{
		int i;
		boolean b1, b2;
		boolean forall = false;

		for (i = 0; i < numStates; i++) {
			if (subset.get(i)) {

				if (coalition.contains(stateLabels.get(i)))
					forall = forall1;
				else //if (stateLabels.get(i).equals(2))
					forall = forall2;

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
	public double mvMultGSMinMax(double[] vect, boolean min1, boolean min2, BitSet subset, boolean complement, boolean absolute)
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
	public double mvMultJacMinMaxSingle(int s, double[] vect, boolean min1, boolean min2)
	{
		boolean min = coalition.contains(stateLabels.get(s)) ? min1 : min2;
		return mvMultJacMinMaxSingle(s, vect, min);
	}

	@Override
	public void mvMultMinMax(double[] vect, boolean min1, boolean min2, double[] result, BitSet subset, boolean complement, int[] adv)
	{
		int s;
		boolean min = false;
		// Loop depends on subset/complement arguments
		if (subset == null) {
			for (s = 0; s < numStates; s++) {
				if (coalition.contains(stateLabels.get(s)))
					min = min1;
				else //if (stateLabels.get(s).equals(2))
					min = min2;

				result[s] = mvMultMinMaxSingle(s, vect, min, adv);
			}
		} else if (complement) {
			for (s = subset.nextClearBit(0); s < numStates; s = subset.nextClearBit(s + 1)) {
				if (coalition.contains(stateLabels.get(s)))
					min = min1;
				else //if (stateLabels.get(s).equals(2))
					min = min2;

				result[s] = mvMultMinMaxSingle(s, vect, min, adv);
			}
		} else {
			for (s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s + 1)) {
				if (coalition.contains(stateLabels.get(s)))
					min = min1;
				else //if (stateLabels.get(s).equals(2))
					min = min2;
				result[s] = mvMultMinMaxSingle(s, vect, min, adv);
			}
		}
	}

	@Override
	public double mvMultMinMaxSingle(int s, double[] vect, boolean min1, boolean min2)
	{
		boolean min = coalition.contains(stateLabels.get(s)) ? min1 : min2;
		return mvMultMinMaxSingle(s, vect, min, null);

	}

	@Override
	public List<Integer> mvMultMinMaxSingleChoices(int s, double[] vect, boolean min1, boolean min2, double val)
	{
		boolean min = coalition.contains(stateLabels.get(s)) ? min1 : min2;
		return mvMultMinMaxSingleChoices(s, vect, min, val);
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
				if (coalition.contains(stateLabels.get(s)))
					min = min1;
				else
					min = min2;
				result[s] = mvMultRewMinMaxSingle(s, vect, mdpRewards, min, adv);
			}
		} else if (complement) {
			for (s = subset.nextClearBit(0); s < numStates; s = subset.nextClearBit(s + 1)) {
				if (coalition.contains(stateLabels.get(s)))
					min = min1;
				else
					min = min2;
				result[s] = mvMultRewMinMaxSingle(s, vect, mdpRewards, min, adv);
			}
		} else {
			for (s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s + 1)) {
				if (coalition.contains(stateLabels.get(s)))
					min = min1;
				else
					min = min2;
				result[s] = mvMultRewMinMaxSingle(s, vect, mdpRewards, min, adv);
			}
		}
	}

	@Override
	public double mvMultRewMinMaxSingle(int s, double vect[], STPGRewards rewards, boolean min1, boolean min2, int adv[])
	{
		MDPRewards mdpRewards = rewards.buildMDPRewards();
		boolean min = coalition.contains(stateLabels.get(s)) ? min1 : min2;		
		return mvMultRewMinMaxSingle(s, vect, mdpRewards, min, adv);
	}

	@Override
	public List<Integer> mvMultRewMinMaxSingleChoices(int s, double vect[], STPGRewards rewards, boolean min1, boolean min2, double val)
	{
		MDPRewards mdpRewards = rewards.buildMDPRewards();
		boolean min = coalition.contains(stateLabels.get(s)) ? min1 : min2;
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
			s += i + "(P-" + stateLabels.get(i) + " " + statesList.get(i) + "): ";
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
}
