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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Vector;

import prism.ModelType;
import prism.PlayerInfo;
import prism.PlayerInfoOwner;

public class CSGSimple extends MDPSimple implements CSG {

	protected Map<Integer, Map<String, Map<String, Map<String, Double>>>> rewards; 
	
	protected List<List<int[]>> transIndexes;
	
	protected BitSet[] indexes; // indexes of actions for each player
	protected Vector<String> actions; // names of actions
 	
	protected int[] idles;
	
	/**
	 * Player information (coalition part not used for now)
	 */
	protected PlayerInfo playerInfo;
	
	public CSGSimple() {
		super();
		transIndexes = new ArrayList<List<int[]>>();
		playerInfo = new PlayerInfo();
	}
	
	public CSGSimple(CSGSimple csg, int[] permut) {
		super(csg, permut);
		//int idle;
		//System.out.println("-- permut " + Arrays.toString(permut));
		transIndexes = new ArrayList<List<int[]>>();
		for (int s = 0; s < csg.getNumStates(); s++) {
			transIndexes.add(s, null);
		}
		for (int s = 0; s < csg.getNumStates(); s++) {
			//System.out.println("-- state " + s);
			//csg.getTransIndexes(s).stream().map(t -> Arrays.toString(t)).forEach(System.out::println);
			transIndexes.set(permut[s], csg.getTransIndexes(s));
		}
		indexes = csg.getIndexes();
		playerInfo = new PlayerInfo(csg.playerInfo);
		actions = csg.getActions();
		idles = csg.getIdles();
		/*
		idle = actions.size() + 1;
		for (int p = 0; p < indexes.length; p++) {
			indexes[p].set(idle);
			idle++;
		}
		*/
	}
	
	@Override
	public void setPlayerNames(List<String> playerNames)
	{
		// Override setPlayerNames because we need to set up indexing
		CSG.super.setPlayerNames(playerNames);
		initIndexes();
	}
	
	/**
	 * Copy the player info from another model
	 */
	public void copyPlayerInfo(PlayerInfoOwner model)
	{
		playerInfo = new PlayerInfo(model.getPlayerInfo());
		initIndexes();
	}
	
	protected void initIndexes() {
		int numPlayers = getNumPlayers();
		idles = new int[numPlayers];
		indexes = new BitSet[numPlayers];
		for (int j = 0; j < numPlayers; j++) {
			indexes[j] = new BitSet();
		}
	}
	
	// Accessors (for PlayerInfoOwner)

	@Override
	public PlayerInfo getPlayerInfo()
	{
		return playerInfo;
	}
	
	// Accessors (for CSG)
	
	public BitSet[] getIndexes() {
		return indexes;
	}
	
	public Vector<String> getActions() {
		return actions;
	}

	public List<List<int[]>> getTransIndexes() {
		return transIndexes;
	}
	
	public List<int[]> getTransIndexes(int s) {
		return transIndexes.get(s);
	}
	
	public int[] getIndexes(int s, int i) {
		return transIndexes.get(s).get(i);
	}
	
	public int[] getIdles() {
		return idles;
	}
	
	public int getIdleForPlayer(int p) {
		return idles[p];
	}
	
	public String[] getActions(int s, int i) {
		int[] indexes = getIndexes(s, i);
		String[] result = new String[indexes.length];
		for (int a = 0; a < indexes.length; a++) {
			result[a] = (indexes[a] > 0)? actions.get(indexes[a] - 1) : "<" + a + ">";
		}
		return result;
	}
	
	public List<String> getActionsForPlayer(int p) {
		List<String> result = new ArrayList<String>();
		for (int i = indexes[p].nextSetBit(0); i >= 0; i = indexes[p].nextSetBit(i + 1)) {	
			result.add(actions.get(i));
		}
		return result;
	}
	
	public Set<String> getActionsForPlayer(int s, int p) {
		Set<String> result = new HashSet<String>();
		String[] actions;
		for (int t = 0; t < getNumChoices(s); t++) {
			actions = getActions(s, t);
			if (!actions[p].equals("-"))
				result.add(actions[p]);
		}
		return result;
	}
	
	public BitSet getIndexesForPlayer(int s, int p) {
		BitSet result = new BitSet();
		int[] indexes;
		for (int t = 0; t < getNumChoices(s); t++) {
			indexes = getIndexes(s, t);
			if (indexes[p] > 0)
				result.set(indexes[p]);
			else 
				result.set(idles[p]);
		}
		return result;
	}
	
	public int[] getNumActions(int s) {
		int numPlayers = getNumPlayers();
		BitSet[] acc = new BitSet[numPlayers];
		int[] tmp;
		int[] result = new int[numPlayers];
		Arrays.fill(result, 0);
		for (int p = 0; p < numPlayers; p++) {
			acc[p] = new BitSet();
			for (int t = 0; t < getNumChoices(s); t++) {
				tmp = transIndexes.get(s).get(t);
				if (tmp[p] > 0) {
					if (!acc[p].get(tmp[p])) {
						acc[p].set(tmp[p]);
						result[p] += 1;
					}
				}
			}
		}
		return result;
	}
	
	public BitSet getActivePlayers(int s) {
		BitSet result = new BitSet();
		int[] numActions = getNumActions(s);
		int numPlayers = getNumPlayers();
		for (int p = 0; p < numPlayers; p++) {
			if (numActions[p] > 0)
				result.set(p);
		}
		return result;
	}
	
	public BitSet getConcurrentPlayers(int s) {
		BitSet result = new BitSet();
		int[] numActions = getNumActions(s);
		int numPlayers = getNumPlayers();
		for (int p = 0; p < numPlayers; p++) {
			if (numActions[p] >= 2)
				result.set(p);
		}
		return result;
	}
	
	public void setIndexes(int s, int t, int[] indexes) {
		this.transIndexes.get(s).add(t, indexes);
	}
	
	public void setIndexes(BitSet[] indexes) {
		this.indexes = indexes;
	}
	
	public void setIdles(int[] idles) {
		this.idles = idles;
	}
	
	public void setActions(List<Object> actions) {
		this.actions = new Vector<String>();
		for (Object action : actions) {
			this.actions.add(action.toString());
		}
	}
	
	public void setActions(Vector<String> actions) {
		this.actions = new Vector<String>(actions);
	}
	
	@Override
	public int addState() {
		addStates(1);
		transIndexes.add(numStates - 1, new ArrayList<int[]>());
		return numStates - 1;
	}
	
	public int addChoice(int s, Distribution distr, int[] indexes) {
		int i = super.addChoice(s, distr);
		for (int j = 0; j < indexes.length; j++) {
			if (indexes[j] >= 0)
				this.indexes[j].set(indexes[j]);
		}
		transIndexes.get(s).add(i, indexes);
		return i;
	}
	
	public void addIdleIndexes() {
		int max = actions.size() + 1;
		for (int p = 0; p < idles.length; p++) {
			actions.add(max - 1, "<" + p + ">");
			idles[p] = max++;
		}
	}
	
	public int addActionLabelledChoice(int s, Distribution distr, int[] indexes) {
		//System.out.println("\n-- Adding choice for state " + s);
		//System.out.println("-- distr " + distr);
		//System.out.println("-- indexes " + Arrays.toString(indexes));		
		int i, j;
		String label = "";
		for (j = 0; j < indexes.length; j++) {
			if(indexes[j] >= 0)
				label += "[" +  actions.get(indexes[j] - 1) + "]";
			else 
				label += "<" + j + ">";
		} 
		i = super.addActionLabelledChoice(s, distr, label);
		//System.out.println("-- label" + label);
		for (j = 0; j < indexes.length; j++) {
			if (indexes[j] >= 0)
				this.indexes[j].set(indexes[j]);
			//else 
				//this.indexes[j].set(idles[j]);
		}
		transIndexes.get(s).add(i, indexes);
		return i;
	}
	
	public void fixDeadlock(int s) {
		int numPlayers = getNumPlayers();
		Distribution distr = new Distribution();
		distr.add(s, 1.0);
		int[] indexes = new int[numPlayers];
		for (int p = 0; p < numPlayers; p++) {
			indexes[p] = -1;
		}
		addActionLabelledChoice(s, distr, indexes);
	}
	
	public void printModelInfo() {
		System.out.println(actions);
		for (int s = 0; s < getNumStates(); s++) {
			System.out.print("\n## state " + s + " : " + Arrays.toString(getNumActions(s)));
			for (int p = 0; p < getNumPlayers(); p++) {
				System.out.print(" : " + getActionsForPlayer(s, p) + " : " + getIndexesForPlayer(s, p));
			}
			System.out.println();
			for (int t = 0; t < getNumChoices(s); t++) {
				System.out.println(Arrays.toString(getIndexes(s, t)) + " : " + Arrays.toString(getActions(s, t)));
			}
		}
	}
}
