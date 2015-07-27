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
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import parser.ast.Coalition;
import prism.ModelType;
import prism.PrismException;

/**
 * Simple explicit-state representation of a (turn-based) stochastic multi-player game (SMG).
 * States can be labelled arbitrarily with player 1..n, player 0 has a special
 * purpose of scheduling the moves of other players
 */
public class SMG extends STPGExplicit implements STPG
{
	// NB: We re-use the existing stateOwners list in the superclass to assign states to players

	// A definition of the players in the game, i.e., the (integer) index and name of each one.
	// It is stored as a mapping from indices to names.
	// Indices can be arbitrary and do not need to be contiguous.
	// Names are optional and can be null or "" if undefined,
	// but all normally defined names should be unique.
	// The size of this map defines the number of players in the game.
	protected Map<Integer, String> playerNames;
	
	// Optionally, a mapping from player indices to 1 or 2,
	// as induced by a coalition of players (those in the coalition
	// are mapped to 1, those who are not are mapped to 2).
	// The mapped index of the player with original index {@code p}
	// is stored as {@code coalitionPlayerMap[p]}, regardless of the
	// range of player indices actually used.
	protected int[] coalitionPlayerMap;

	// Constructors

	/**
	 * Constructor: empty SMG.
	 */
	public SMG()
	{
		super();
		playerNames = new HashMap<Integer, String>();
		coalitionPlayerMap = null;
	}

	/**
	 * Constructor: new SMG with fixed number of states.
	 */
	public SMG(int numStates)
	{
		super(numStates);
		playerNames = new HashMap<Integer, String>();
		coalitionPlayerMap = null;
	}

	/**
	 * Construct an SMG from an existing one and a state index permutation,
	 * i.e. in which state index i becomes index permut[i].
	 * Player and coalition info is also copied across.
	 */
	public SMG(SMG smg, int permut[])
	{
		super(smg, permut);
		playerNames = new HashMap<Integer, String>(smg.playerNames);
		coalitionPlayerMap = smg.coalitionPlayerMap == null ? null : smg.coalitionPlayerMap.clone(); 
	}

	/**
	 * Copy constructor
	 */
	public SMG(SMG smg)
	{
		super(smg);
		playerNames = new HashMap<Integer, String>(smg.playerNames);
		coalitionPlayerMap = smg.coalitionPlayerMap == null ? null : smg.coalitionPlayerMap.clone(); 
	}

	// Mutators

	/**
	 * Add a new (player 1) state and return its index.
	 */
	@Override
	public int addState()
	{
		return addState(0);
	}

	/**
	 * Add multiple new (player 1) states.
	 */
	@Override
	public void addStates(int numToAdd)
	{
		for (int i = 0; i < numToAdd; i++) {
			addState();
		}
	}

	/**
	 * Set the info about players, i.e., the (integer) index and name of each one.
	 * This is given as a mapping from indices to names.
	 * Indices can be arbitrary and do not need to be contiguous.
	 * Names are optional and can be null or "" if undefined,
	 * but all normally defined names should be unique.
	 */
	public void setPlayerInfo(Map<Integer, String> playerNames)
	{
		this.playerNames = new HashMap<Integer, String>(playerNames);
	}

	/**
	 * Copy player info from another SMG.
	 */
	public void copyPlayerInfo(SMG smg)
	{
		setPlayerInfo(smg.playerNames);
	}

	/**
	 * Set a coalition of players for this SMG
	 * (which effectively makes it an STPG with player 1 representing the coalition and 2 the rest).
	 * Pass null to remove any coalition info from this SMG.
	 * 
	 * @param coalition Coalition info object 
	 */
	public void setCoalition(Coalition coalition) throws PrismException
	{
		// Clear info if coalition is null
		if (coalition == null) {
			coalitionPlayerMap = null;
			return;
		}
		
		// This is first time we really need the {@code playerNames} info,
		// so if it has not been set, we tried to create it based on {@code stateOwners} 
		if (playerNames.isEmpty()) {
			for (int i = 0; i < numStates; i++) {
				int p = stateOwners.get(i);
				if (!playerNames.containsKey(p)) {
					playerNames.put(p, null);
				}
			}
		}
		
		// Find max player index
		int maxIndex = 0;
		for (int index : playerNames.keySet()) {
			maxIndex = Math.max(maxIndex, index);
		}
		
		// Construct mapping
		coalitionPlayerMap = new int[maxIndex + 1];
		for (int i = 0; i < maxIndex + 1; i++) {
			coalitionPlayerMap[i] = -1;
		}
		for (Entry<Integer, String> entry : playerNames.entrySet()) {
			int playerIndex = entry.getKey();
			boolean inCoalition = coalition.isPlayerIndexInCoalition(playerIndex, playerNames);
			// In coalition => player 1; not in coalition (or undefined) => player 2
			coalitionPlayerMap[playerIndex] = inCoalition ? 1 : 2;
		}
	}

	/**
	 * Copy coalition info from another SMG.
	 */
	public void copyCoalitionInfo(SMG smg)
	{
		coalitionPlayerMap = smg.coalitionPlayerMap == null ? null : smg.coalitionPlayerMap.clone();
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
		smg.stateOwners = new ArrayList<Integer>(this.stateOwners);
		smg.trans = new ArrayList<List<Distribution>>(this.trans);
		return smg;
	}

	// Accessors (for Model)

	@Override
	public ModelType getModelType()
	{
		return ModelType.SMG;
	}

	// Accessors (for STPG/SMG)

	//@Override
	public int getNumPlayers()
	{
		return playerNames.size();
	}
	
	@Override
	public int getPlayer(int s)
	{
		int playerIndex = stateOwners.get(s); 
		if (coalitionPlayerMap == null) {
			// No coalition: just return index
			return playerIndex;
		} else {
			// Coalition defined: look up if player 1 or 2
			// (note: undefined players are mapped to player 2)
			return playerIndex == -1 ? 2 : coalitionPlayerMap[playerIndex];
		}
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
			s += i + "(P-" + stateOwners.get(i) + " " + statesList.get(i) + "): ";
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
