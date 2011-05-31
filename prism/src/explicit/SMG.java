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
import java.util.List;
import java.util.Map;

import prism.ModelType;

/**
 * Simple explicit-state representation of a stochastic multiplayer game
 * (SMG). States can be labelled arbitrarily with player 1..n, player 0
 * has a special purpose of scheduling the moves of other players
 */
public class SMG extends MDPSimple
{
	public static final int SCHED_RANDOM = 0;
	public static final int SCHED_NONDET = 1;
	public static final int SCHED_ALTERNATING = 2;

	// state labels: state labelled with label i
	// is controlled by player i
	protected List<Integer> stateLabels;

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
	 * Method transforms SMG to STPG using the specified parameters. 
	 * Scheduler (player 0) has to be included  in player1 or player2
	 *
	 * @param player1 list of players to represent player 1
	 * @param player2 list of players to represent player 2
	 * @param schedType type of scheduler
	 * @return the new SMG which has only two players - indexed 1 and 2
	 */
	public SMG constructSTPG(List<Integer> player1, List<Integer> player2, int schedType)
	{
		// resolving scheduling
		SMG stpg = this.clone();

		// relabeling players
		
		return null;
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

	@Override
	public ModelType getModelType()
	{
		return ModelType.SMG;
	}

	/**
	 * Makes a deep copy of itself
	 */
	public SMG clone()
	{
		SMG smg = new SMG();
		smg.copyFrom(this);
		smg.actions = this.actions;
		smg.allowDupes = this.allowDupes;
		smg.maxNumDistrs = this.maxNumDistrs;
		smg.maxNumDistrsOk = this.maxNumDistrsOk;
		smg.numDistrs = this.numDistrs;
		smg.numTransitions = this.numTransitions;
		smg.stateLabels = this.stateLabels;
		smg.trans = this.trans;
		smg.transRewards = this.transRewards;
		smg.transRewardsConstant = this.transRewardsConstant;
		
		return smg;
	}
	
	/**
	 * Get transition function as string.
	 */
	public String toString()
	{
		int i, j, n;
		Object o;
		String s = "";
		s = "[ ";
		for (i = 0; i < numStates; i++) {
			if (i > 0)
				s += ", ";
			s += i + "(P-" + stateLabels.get(i) + "): ";
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
