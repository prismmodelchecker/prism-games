//==============================================================================
//	
//	Copyright (c) 2020-
//	Authors:
//	* Dave Parker <d.a.parker@cs.bham.ac.uk> (University of Birmingham)
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

package prism;

import java.util.List;

import parser.ast.Coalition;

/**
 * Interface implemented by game model classes that have info
 * about players, stored as a {@link PlayerInfo} object.
 * 
 * Also provides some convenience methods to access this info
 * in the form of default method implementations.
 */
public interface PlayerInfoOwner
{
	/**
	 * Get the {@link PlayerInfo} object.
	 */
	public PlayerInfo getPlayerInfo();
	
	// Default implementations of methods to modify player info
	
	/**
	 * Add a player, specifying its name.
	 * Names are optional and can be null or "" if undefined,
	 * but all normally defined names should be unique.
	 */
	public default void addPlayer(String playerName)
	{
		getPlayerInfo().addPlayer(playerName);
	}

	/**
	 * Set the names of all players, provided as a list.
	 * Names are optional and can be null or "" if undefined,
	 * but all normally defined names should be unique.
	 */
	public default void setPlayerNames(List<String> playerNames)
	{
		getPlayerInfo().setPlayerNames(playerNames);
	}
	
	/**
	 * Set a coalition of players for the model
	 * (which effectively makes it a 2-player model,
	 * with player 0 representing the coalition and 1 the rest).
	 * Pass null to remove any coalition info from the model.
	 * @param coalition Coalition info object
	 */
	public default void setCoalition(Coalition coalition) throws PrismException
	{
		getPlayerInfo().setCoalition(coalition);
	}

	// Default implementations of methods to access player info
	
	/**
	 * Get the number of players in the game.
	 */
	public default int getNumPlayers()
	{
		return getPlayerInfo().getNumPlayers();
	}
	
	/**
	 * Get the name of player {@code i}. Returns "" if unnamed.
	 * @param i Index of player (0-indexed)
	 */
	public default String getPlayerName(int i)
	{
		return getPlayerInfo().getPlayerName(i);
	}
	
	/**
	 * Get the list of player names. A name is "" if the player is unnamed.
	 */
	public default List<String> getPlayerNames()
	{
		return getPlayerInfo().getPlayerNames();
	}
}
