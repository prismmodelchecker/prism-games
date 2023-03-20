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

package prism;

import java.util.Vector;

import jdd.JDD;
import jdd.JDDNode;
import jdd.JDDVars;
import parser.Values;
import parser.VarList;

/*
 * Class for MTBDD-based storage of a PRISM model that is an SMG.
 */
public class GamesModel extends NondetModel implements PlayerInfoOwner
{
	/** All DD vars for players */
	private JDDVars allDDPlayerVars;
	
	/** Cube over DD vars for each player (i.e. one-hot encoding) */
	private JDDNode[] ddPlayerCubes;

	/** Player + coalition information */
	protected PlayerInfo playerInfo;
	
	public GamesModel(JDDNode tr, JDDNode s, JDDNode[] sr, JDDNode[] trr, String[] rsn, JDDVars arv, JDDVars acv, JDDVars asyv, JDDVars asv,
					  JDDVars achv, JDDVars andv, ModelVariablesDD mvdd, int nm, String[] mn, JDDVars[] mrv, JDDVars[] mcv, int nv, 
					  VarList vl, JDDVars[] vrv, JDDVars[] vcv, Values cv, JDDVars apvs, JDDNode[] ddPlayerCubes, PlayerInfo playerInfo) {
		super(tr, s, sr, trr, rsn, arv, acv, asyv, asv, achv, andv, mvdd, nm, mn, mrv,
				mcv, nv, vl, vrv, vcv, cv);
		this.allDDPlayerVars = apvs;
		this.ddPlayerCubes = ddPlayerCubes;
		this.playerInfo = playerInfo;
	}

	public ModelType getModelType()
	{
		return ModelType.SMG;
	}
	
	public JDDVars getAllDDPlayerVars()
	{
		return allDDPlayerVars;
	}

	public JDDNode getDdPlayerCube(int p)
	{
		return ddPlayerCubes[p];
	}

	@Override
	public PlayerInfo getPlayerInfo()
	{
		return playerInfo;
	}
	
	public NondetModel toMDP() throws PrismException
	{
		JDDNode transMDP = JDD.SumAbstract(trans.copy(), allDDPlayerVars);
		JDDNode transRewardsMDP[] = new JDDNode[transRewards.length];

		for(int i = 0; i < transRewards.length; i++) {
			transRewardsMDP[i] = JDD.SumAbstract(transRewards[i], allDDPlayerVars);
		}
		
		JDDVars allDDNondetVarsMDP = allDDNondetVars.copy();
		allDDNondetVarsMDP.removeVars(allDDPlayerVars);
		
		NondetModel equivMDP =  new NondetModel(transMDP, start, stateRewards, transRewardsMDP, rewardStructNames, allDDRowVars, allDDColVars,
			     						allDDSynchVars, allDDSchedVars, allDDChoiceVars, allDDNondetVarsMDP, modelVariables,
			     						numModules, moduleNames, moduleDDRowVars, moduleDDColVars,
			     						numVars, varList, varDDRowVars, varDDColVars, constantValues);
		
		equivMDP.setSynchs((Vector<String>)synchs);
	
		equivMDP.setTransInd(transInd);
		equivMDP.setTransSynch(transSynch);
		
		equivMDP.setTransActions(transActions);
		equivMDP.setTransPerAction(transPerAction);
		
		equivMDP.doReachability();
		equivMDP.filterReachableStates();
		
		return equivMDP;
	}
	
	@Override
	public void clear()
	{
		super.clear();
		allDDPlayerVars.derefAll();
		JDD.DerefArray(ddPlayerCubes, getNumPlayers());
	}
}