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
import java.util.List;

import prism.PrismException;

/**
 * Simple explicit-state representation of a stochastic two-player game
 * (STPG). States can be labeled arbitrarily with player 1 player 2.
 */
public class STPGExplicit extends STPGAbstrSimple
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

	// TODO fix the method
	@Override
	public void mvMultMinMax(double vect[], boolean min1, boolean min2, double result[], BitSet subset,
			boolean complement)
	{
		int s;
		boolean p1min, p2min;
		p1min = p2min = false;

		// Loop depends on subset/complement arguments
		if (subset == null) {
			for (s = 0; s < numStates; s++) {
				if (stateLabels.get(s) == PLAYER_1)
					p1min = p2min = min1;
				else if (stateLabels.get(s) == PLAYER_2)
					p1min = p2min = min2;

				result[s] = mvMultMinMaxSingle(s, vect, p1min, p2min);
			}
		} else if (complement) {
			for (s = subset.nextClearBit(0); s < numStates; s = subset.nextClearBit(s + 1)) {
				if (stateLabels.get(s) == PLAYER_1)
					p1min = p2min = min1;
				else if (stateLabels.get(s) == PLAYER_2)
					p1min = p2min = min2;

				result[s] = mvMultMinMaxSingle(s, vect, p1min, p2min);
			}
		} else {
			for (s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s + 1)) {
				if (stateLabels.get(s) == PLAYER_1)
					p1min = p2min = min1;
				else if (stateLabels.get(s) == PLAYER_2)
					p1min = p2min = min2;
				result[s] = mvMultMinMaxSingle(s, vect, p1min, p2min);
			}
		}
	}

	@Override
	public void prob0step(BitSet subset, BitSet u, boolean forall1, boolean forall2, BitSet result)
	{
		int i;
		boolean b1, b2, b3;
		boolean forall11, forall22;
		forall11 = forall22 = false;

		for (i = 0; i < numStates; i++) {
			if (subset.get(i)) {

				if (stateLabels.get(i) == PLAYER_1)
					forall11 = forall22 = forall1;
				else if (stateLabels.get(i) == PLAYER_2)
					forall11 = forall22 = forall2;

				b1 = forall11; // there exists or for all player 1 choices
				for (DistributionSet distrs : trans.get(i)) {
					b2 = forall22; // there exists or for all player 2 choices
					for (Distribution distr : distrs) {
						b3 = distr.containsOneOf(u);
						if (forall22) {
							if (!b3)
								b2 = false;
						} else {
							if (b3)
								b2 = true;
						}
					}
					if (forall11) {
						if (!b2)
							b1 = false;
					} else {
						if (b2)
							b1 = true;
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
		boolean b1, b2, b3;
		boolean forall11, forall22;
		forall11 = forall22 = false;

		for (i = 0; i < numStates; i++) {
			if (subset.get(i)) {

				if (stateLabels.get(i) == PLAYER_1)
					forall11 = forall22 = forall1;
				else if (stateLabels.get(i) == PLAYER_2)
					forall11 = forall22 = forall2;

				b1 = forall11; // there exists or for all player 1 choices
				for (DistributionSet distrs : trans.get(i)) {
					b2 = forall22; // there exists or for all player 2 choices
					for (Distribution distr : distrs) {
						b3 = distr.containsOneOf(v) && distr.isSubsetOf(u);
						if (forall22) {
							if (!b3)
								b2 = false;
						} else {
							if (b3)
								b2 = true;
						}
					}
					if (forall11) {
						if (!b2)
							b1 = false;
					} else {
						if (b2)
							b1 = true;
					}
				}
				result.set(i, b1);
			}
		}
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

	/**
	 * Adds distribution to the distribution set of the given state
	 * 
	 * @param s state
	 * @param dist distribution
	 */
	public void addDistribution(int s, Distribution dist)
	{
		// initialising distribution set
		DistributionSet distSet;
		if (trans.get(s).size() == 0) {
			distSet = newDistributionSet(null);
			addDistributionSet(s, distSet);
		} else {
			distSet = trans.get(s).get(0);
		}

		// adding distribution to the set
		distSet.add(dist);

	}

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

	/**
	 * Get transition function as string.
	 */
	public String toString()
	{
		int i;
		boolean first;
		String s = "";
		first = true;
		s = "[ ";
		for (i = 0; i < numStates; i++) {
			if (first)
				first = false;
			else
				s += ", ";
			s += i + "(P-" + stateLabels.get(i) + "): " + trans.get(i) + "\n";
		}
		s += " ]\n";
		return s;
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

			// Print game
			System.out.println(stpg);

			// Model check
			mc = new STPGModelChecker();
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
