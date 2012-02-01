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
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Map.Entry;

import explicit.rewards.MDPRewards;
import explicit.rewards.STPGRewards;

import prism.ModelType;
import prism.PrismException;

/**
 * Simple explicit-state representation of a stochastic multi-player game (SMG).
 * States can be labelled arbitrarily with player 1..n, player 0 has a special
 * purpose of scheduling the moves of other players
 */
public class SMG extends STPGExplicit implements STPG {

	// State labels: states with label i are controlled by player i
	protected List<Integer> stateLabels;

	// Set of players which form a coalition
	protected Set<Integer> coalition;

	public SMG() {
		super();
		stateLabels = new ArrayList<Integer>(numStates);
	}

	/**
	 * Construct an SMG from an existing one and a state index permutation, i.e.
	 * in which state index i becomes index permut[i].
	 */
	public SMG(SMG smg, int permut[]) {
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
	public Set<Integer> getSchedulerStates() {
		Set<Integer> ret = new HashSet<Integer>();
		return ret;
	}

	/**
	 * Adds one state, assigned to player 0
	 */
	@Override
	public int addState() {
		return addState(0);
	}

	/**
	 * Adds specified number of states all assigned to player 0
	 */
	@Override
	public void addStates(int numToAdd) {
		for (int i = 0; i < numToAdd; i++)
			stateLabels.add(0);
	}

	/**
	 * Adds state assigned to the specified player
	 * 
	 * @param player
	 *            state owner
	 * @return state id
	 */
	public int addState(int player) {
		super.addStates(1);
		stateLabels.add(player);
		return numStates - 1;
	}

	/**
	 * Adds the number of states the same as number of Integer in the list, each
	 * assigned to the corresponding player
	 * 
	 * @param players
	 *            list of players (to which corresponding state belongs)
	 */
	public void addStates(List<Integer> players) {
		super.addStates(players.size());
		stateLabels.addAll(players);
	}

	/**
	 * labels the given state with the given player
	 * 
	 * @param s
	 *            state
	 * @param player
	 *            player
	 */
	public void setPlayer(int s, int player) {
		if (s < stateLabels.size())
			stateLabels.set(s, player);
	}

	/**
	 * Sets the coalition (representing Player 1)
	 * @param coalition
	 */
	public void setCoalition(Set<Integer> coalition) {
		this.coalition = coalition;
	}

	/**
	 * Makes a half-deep (up to one reference level) copy of itself
	 */
	public SMG clone() {
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
	public ModelType getModelType() {
		return ModelType.SMG;
	}

	// Accessors (for STPG)

	/**
	 * Returns the 1 if the state belong to coalition and 2 otherwise
	 */
	@Override
	public int getPlayer(int s) {
		return coalition.contains(stateLabels.get(s))?1:2;
	}



	// Methods for intervals
	// TODO fix the method
	public List<List<Interval>> mvMultIntervals(boolean min1, boolean min2,
			List<List<Interval>> intervals) throws PrismException {
		int s;
		boolean min = false;
		List<List<Interval>> intervals_ = new ArrayList<List<Interval>>(
				intervals.size());
		for (s = 0; s < numStates; s++) {
			if (getPlayer(s)==1)
				min = min1;
			else
				min = min2;
			// System.out.println("Checking state " + s);
			intervals_.add(mvMultIntervalsSingle(s, intervals, min));
		}
		return intervals_;
	}

	public List<Interval> mvMultIntervalsSingle(int s,
			List<List<Interval>> intervals, boolean min) throws PrismException {
		int k;
		List<Distribution> step;
		List<List<Interval>> distInts;
		List<Interval> ints;
		Interval int1, int2;

		int states[] = new int[2];
		double probs[] = new double[2];

		distInts = new ArrayList<List<Interval>>(2);
		// collect intervals for distributions
		step = trans.get(s);
		if (step.size() == 1)
			step.add(step.get(0));

		// System.out.println("Checking distributions");
		for (Distribution distr : step) {

			k = 0;
			// gets transition probabilities
			for (Map.Entry<Integer, Double> e : distr) {
				states[k] = (Integer) e.getKey();
				probs[k++] = (Double) e.getValue();
			}

			ints = new ArrayList<Interval>(10);
			for (int i = 0; i < intervals.get(states[0]).size(); i++) {
				int1 = intervals.get(states[0]).get(i);
				for (int j = 0; j < intervals.get(states[1]).size(); j++) {
					int2 = intervals.get(states[1]).get(j);
					ints.add(new Interval((probs[0] * int1.lhs)
							+ (probs[1] * int2.lhs), (probs[0] * int1.rhs)
							+ (probs[1] * int2.rhs)));
				}
			}

			// System.out.println("Filtering");
			// System.out.println(ints);
			// removing obsolete intervals
			filterIntervals(ints);
			// System.out.println(ints);

			// add intervals to the distr
			distInts.add(ints);

		}

		// System.out.println("Taking unions and intersections");
		// System.out.println(distInts);
		// generating a new list
		ints = new LinkedList<Interval>();

		// if player maximising player's state take the 'union'
		if (!min) {
			// System.out.println("Union");
			// joining all the lists
			for (List<Interval> l : distInts)
				ints.addAll(l);
		}
		// if player is minimising player's state: take the 'intersection'
		else {
			// System.out.println("Intersection");
			// Adding unions
			for (int i = 0; i < distInts.get(0).size(); i++)
				for (int j = 0; j < distInts.get(1).size(); j++)
					ints.add(Interval.getUnion(distInts.get(0).get(i), distInts
							.get(1).get(j)));

			// Adding intersections
			for (Interval i1 : distInts.get(0))
				for (Interval i21 : distInts.get(1))
					for (Interval i22 : distInts.get(1))
						if (i1.achievableFrom(i21) || i1.achievableFrom(i21)
								|| i1.achievableFrom(i21, i22))
							ints.add(i1);
			for (Interval i2 : distInts.get(1))
				for (Interval i11 : distInts.get(0))
					for (Interval i12 : distInts.get(0))
						if (i2.achievableFrom(i11) || i2.achievableFrom(i11)
								|| i2.achievableFrom(i11, i12))
							ints.add(i2);

		}
		// System.out.println("Filtering");
		// System.out.println(ints);
		// removing obsolete intervals
		filterIntervals(ints);
		// System.out.println(ints);

		// System.out.println("Finished with this state");

		return ints;

	}

	private void filterIntervals(List<Interval> ints) {

		Interval int1;
		boolean remove;

		// remove intervals which are the same
		BitSet toRemove = new BitSet(ints.size());

		for (int i = 0; i < ints.size(); i++) {
			int1 = ints.get(i);

			remove = false;

			// checking for equality
			for (int j = i + 1; j < ints.size(); j++) {
				if (int1.equals(ints.get(j))) {
					toRemove.set(i);
					remove = true;
					break;
				}
			}
			if (remove)
				continue;

//			// checking for individual containment
//			for (int j = 0; j < ints.size(); j++) {
//				if (j == i || toRemove.get(j))
//					continue;
//				if (int1.achievableFrom(ints.get(j))) {
//					toRemove.set(i);
//					remove = true;
//					break;
//				}
//			}
//			if (remove)
//				continue;

			// check for convex containment
			for (int j = 0; j < ints.size(); j++) {
				if (j == i || toRemove.get(j))
					continue;
				for (int k = 0; k < ints.size(); k++) {
					if (k == i || toRemove.get(k))
						continue;
					if (int1.achievableFrom(ints.get(j), ints.get(k))) {
						toRemove.set(i);
						remove = true;
						break;
					}
				}
				if (remove)
					break;
			}
		}

		removeElements(ints, toRemove);

		// if(ints.isEmpty())
		// ints.add(new Interval(0,1));
	}

	private void removeElements(List<Interval> ints, BitSet remove) {

		for (int i = ints.size() - 1; i >= 0; i--)
			if (remove.get(i))
				ints.remove(i);
	}

	// Standard methods

	@Override
	public String toString() {
		int i, j, n;
		Object o;
		String s = "";
		s = "[ ";
		for (i = 0; i < numStates; i++) {
			if (i > 0)
				s += ", ";
			s += i + "(P-" + stateLabels.get(i) + " " + statesList.get(i)
					+ "): ";
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
