//==============================================================================
//	
//	Copyright (c) 2016-
//	Authors:
//	* Gabriel Santos <gabriel.santos@cs.ox.ac.uk> (University of Oxford)
//	* Dave Parker <david.parker@cs.ox.ac.uk> (University of Oxford)
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

package symbolic.model;

import java.util.Vector;

import jdd.JDD;
import jdd.JDDNode;
import jdd.JDDVars;
import parser.Values;
import parser.VarList;
import prism.ModelType;
import prism.PlayerInfo;
import prism.PlayerInfoOwner;
import prism.PrismException;

/**
 * Class for symbolic (BDD-based) representation of an SMG.
 */
public class GamesModel extends NondetModel implements PlayerInfoOwner
{
	/** All DD vars for players */
	private JDDVars allDDPlayerVars;
	
	/** Cube over DD vars for each player (i.e. one-hot encoding) */
	private JDDNode[] ddPlayerCubes;

	/** Player + coalition information */
	protected PlayerInfo playerInfo;

	// Constructor

	public GamesModel(JDDNode trans, JDDNode start, JDDVars allDDRowVars, JDDVars allDDColVars, JDDVars allDDNondetVars, ModelVariablesDD modelVariables,
					   VarList varList, JDDVars[] varDDRowVars, JDDVars[] varDDColVars, JDDVars allDDPlayerVars, JDDNode[] ddPlayerCubes, PlayerInfo playerInfo)
	{
		super(trans, start, allDDRowVars, allDDColVars, allDDNondetVars, modelVariables, varList, varDDRowVars, varDDColVars);
		this.allDDPlayerVars = allDDPlayerVars;
		this.ddPlayerCubes = ddPlayerCubes;
		this.playerInfo = playerInfo;
	}

	// Accessors (for Model)

	@Override
	public ModelType getModelType()
	{
		return ModelType.SMG;
	}

	@Override
	public void clear()
	{
		super.clear();
		allDDPlayerVars.derefAll();
		JDD.DerefArray(ddPlayerCubes, getNumPlayers());
	}

	// Accessors (for PlayerInfoOwner)

	@Override
	public PlayerInfo getPlayerInfo()
	{
		return playerInfo;
	}

	// Accessors (for SMG)

	/**
	 * Get the DD variables for encoding players
	 */
	public JDDVars getAllDDPlayerVars()
	{
		return allDDPlayerVars;
	}

	/**
	 * Get a cube over DD vars for each player (i.e. one-hot encoding).
	 */
	public JDDNode getDdPlayerCube(int p)
	{
		return ddPlayerCubes[p];
	}

	/**
	 * Convert to an MDP.
	 */
	public NondetModel toMDP() throws PrismException
	{
		JDDNode transMDP = JDD.SumAbstract(trans.copy(), allDDPlayerVars);
		JDDNode transRewardsMDP[] = new JDDNode[transRewards.length];

		for(int i = 0; i < transRewards.length; i++) {
			transRewardsMDP[i] = JDD.SumAbstract(transRewards[i], allDDPlayerVars);
		}
		
		JDDVars allDDNondetVarsMDP = allDDNondetVars.copy();
		allDDNondetVarsMDP.removeVars(allDDPlayerVars);
		
		NondetModel equivMDP =  new NondetModel(transMDP, start, allDDRowVars, allDDColVars, allDDNondetVarsMDP, modelVariables,
												varList, varDDRowVars, varDDColVars);
		equivMDP.setRewards(stateRewards, transRewardsMDP, rewardStructNames);

		equivMDP.setSynchs((Vector<String>)synchs);
	
		equivMDP.setTransInd(transInd);
		equivMDP.setTransSynch(transSynch);
		
		equivMDP.setTransActions(transActions);

		equivMDP.doReachability();
		equivMDP.filterReachableStates();
		
		return equivMDP;
	}
}