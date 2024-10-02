//==============================================================================
//	
//	Copyright (c) 2020-
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford)
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

package explicit;

import java.util.BitSet;
import java.util.List;

import prism.ModelType;
import prism.PlayerInfoOwner;

/**
 * Interface for classes that provide (read) access to an explicit-state concurrent stochastic game (CSG).
 */
public interface CSG<Value> extends MDP<Value>, PlayerInfoOwner
{
	// Accessors (for Model) - default implementations
	
	@Override
	public default ModelType getModelType()
	{
		return ModelType.CSG;
	}
	
	// Accessors (for CSG)
	
	/**
	 * Get the list of all action labels
	 */
	public List<Object> getActions();
	
	/**
	 * Get the list of player action indices for choice {@code i} of state {@code s},
	 * i.e., as an array of integers, giving the (1-indexed) indices of
	 * the actions for each player attached to the choice.
	 * An index of -1 indicates that a player idles.
	 */
	public int[] getIndexes(int s, int i);
	
	/**
	 * Get the indices of actions taken in some choice of state {@code s} by player {@code p},
	 * including special "idle" actions.
	 * @param s Index of state (0-indexed)
	 * @param p Index of player (0-indexed)
	 */
	public BitSet getIndexesForPlayer(int s, int p);
	
	/**
	 * Get the list of player actions for choice {@code i} of state {@code s},
	 * as an array of strings, and where "&lt;p&gt;" denotes player p idling.
	 */
	public String[] getActions(int s, int i);
	
	/**
	 * Get the indices of all actions owned by each player,
	 * i.e., a BitSet of (1-indexed) action indices for each player.
	 */
	public BitSet[] getIndexes();
	
	/**
	 * Get the indices of the actions representing "idle" for each player.
	 */
	public int[] getIdles();
	
	/**
	 * Get the index of the actions representing "idle" for player {@code p}.
	 * @param p Index of player (0-indexed)
	 */
	public int getIdleForPlayer(int p);
	
	/**
	 * Get a BitSet specifying whether each player has more than 1 action in state {@code s}.
	 */
	public BitSet getConcurrentPlayers(int s);
	
	// Temp:

	public Distribution<Value> getChoice(int s, int i);
}
