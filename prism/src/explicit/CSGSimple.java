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
import java.util.Set;
import java.util.Vector;

import prism.JointAction;
import prism.PlayerInfo;
import prism.PlayerInfoOwner;

/**
 * Simple explicit-state representation of a (multi-player) concurrent stochastic game (CSG).
 */
public class CSGSimple<Value> extends MDPSimple<Value> implements CSG<Value>
{
	/** List of all action labels */
	protected List<Object> actions;
	
	/** List of player action indices for each state/choice,
	 * stored as an array giving the (1-indexed) index for the action
	 * performed by each player in this transition, and -1 indicates that the player idles. */
	protected List<List<int[]>> transIndexes;
	
	/** Indices of actions owned by each player,
	 * i.e., a BitSet of (1-indexed) action indices for each player. */
	protected BitSet[] indexes;

	/** Indices of the actions representing "idle" for each player. */
	protected int[] idles;

	/**
	 * Player information (coalition part not used for now)
	 */
	protected PlayerInfo playerInfo;

	// Constructors

	/**
	 * Constructor: empty CSG.
	 */
	public CSGSimple()
	{
		super();
		transIndexes = new ArrayList<List<int[]>>();
		playerInfo = new PlayerInfo();
	}

	/**
	 * Construct a CSG from an existing one and a state index permutation,
	 * i.e. in which state index i becomes index permut[i].
	 * Player and coalition info is also copied across.
	 */
	public CSGSimple(CSGSimple<Value> csg, int[] permut)
	{
		super(csg, permut);
		transIndexes = new ArrayList<List<int[]>>();
		for (int s = 0; s < csg.getNumStates(); s++) {
			transIndexes.add(s, null);
		}
		for (int s = 0; s < csg.getNumStates(); s++) {
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

	// Mutators

	@Override
	public void clearState(int s)
	{
		super.clearState(s);
		transIndexes.set(s, new ArrayList<int[]>());
	}

	@Override
	public int addState()
	{
		int i = super.addState();
		transIndexes.add(new ArrayList<int[]>());
		return i;
	}

	@Override
	public void addStates(int numToAdd)
	{
		super.addStates(numToAdd);
		for (int i = 0; i < numToAdd; i++) {
			transIndexes.add(new ArrayList<int[]>());
		}
	}

	/**
	 * Add a choice (distribution {@code distr}) to state {@code s} (which must exist).
	 * Behaves the same as {@link MDPSimple#addActionLabelledChoice(int, Distribution, Object)},
	 * but the action must be a {@link JointAction}, from which indexing info is generated
	 * (it can then be accessed via {@link #getIndexes(int, int)}).
	 */
	@Override
	public int addActionLabelledChoice(int s, Distribution<Value> distr, Object action)
	{
		if (!(action instanceof JointAction)) {
			throw new RuntimeException("CSG action labels must be JointAction objects");
		}
		JointAction jointAction = (JointAction) action;
		// Store choice in underlying MDP
		int i = super.addActionLabelledChoice(s, distr, jointAction);
		// Generate indexing info
		int numPlayers = getNumPlayers();
		int[] indexes = new int[numPlayers];
		for (int p = 0; p < numPlayers; p++) {
			if (jointAction.get(p) == JointAction.IDLE_ACTION) {
				indexes[p] = -1;
			} else {
				int j = getActions().indexOf(jointAction.get(p));
				if (j == -1) {
					throw new RuntimeException("CSG action labels must be valid JointAction objects");
				}
				indexes[p] = j + 1;
				this.indexes[p].set(indexes[p]);

			}
		}
		transIndexes.get(s).add(i, indexes);
		return i;
	}

	/**
	 * Add a choice (distribution {@code distr}) to state {@code s} (which must exist).
	 * Behaves the same as {@link MDPSimple#addActionLabelledChoice(int, Distribution, Object)},
	 * but {@code indexes} is an array storing the (1-indexed) index for the action
	 * performed by each player in this transition, and -1 indicates that the player idles.
	 * A representation of this is stored as a {@link JointAction} (accessible via e.g.
	 * {@link #getAction(int, int)}), whereas the array of indices can be accessed via
	 * {@link #getIndexes(int, int)}.
	 */
	public int addActionLabelledChoice(int s, Distribution<Value> distr, int[] indexes)
	{
		int numPlayers = getNumPlayers();
		// Create joint action and store choice in underlying MDP
		JointAction jointAction = new JointAction(indexes, actions);
		int i = super.addActionLabelledChoice(s, distr, jointAction);
		// Store indexing info
		for (int p = 0; p < numPlayers; p++) {
			if (indexes[p] >= 0) {
				this.indexes[p].set(indexes[p]);
			}
		}
		transIndexes.get(s).add(i, indexes);
		return i;
	}

	/**
	 * Set the list of player action indices for choice {@code i} of state {@code s}:
	 * {@code indexes} is an array storing the (1-indexed) index for the action
	 * performed by each player in this transition, and -1 indicates that the player idles.
	 * A string representation of this is stored as an action label (accessible via e.g.
	 * {@link #getAction(int, int)}), whereas the array of indices can be accessed via
	 * {@link #getIndexes(int, int)}.
	 */
	public void setIndexes(int s, int i, int[] indexes)
	{
		// TODO: call setAction too?
		this.transIndexes.get(s).add(i, indexes);
	}

	/**
	 * Set the list of all action labels
	 */
	public void setActions(List<Object> actions)
	{
		this.actions = new ArrayList<>(actions);
	}

	public void setIndexes(BitSet[] indexes)
	{
		this.indexes = indexes;
	}

	public void setIdles(int[] idles)
	{
		this.idles = idles;
	}

	public void addIdleIndexes()
	{
		int max = actions.size() + 1;
		for (int p = 0; p < idles.length; p++) {
			actions.add(max - 1, "<" + p + ">");
			idles[p] = max++;
		}
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

	protected void initIndexes()
	{
		int numPlayers = getNumPlayers();
		idles = new int[numPlayers];
		indexes = new BitSet[numPlayers];
		for (int j = 0; j < numPlayers; j++) {
			indexes[j] = new BitSet();
		}
	}

	public void fixDeadlock(int s)
	{
		int numPlayers = getNumPlayers();
		Distribution<Value> distr = new Distribution<>();
		distr.add(s, getEvaluator().one());
		int[] indexes = new int[numPlayers];
		for (int p = 0; p < numPlayers; p++) {
			indexes[p] = -1;
		}
		addActionLabelledChoice(s, distr, indexes);
	}

	// Accessors (for PlayerInfoOwner)

	@Override
	public PlayerInfo getPlayerInfo()
	{
		return playerInfo;
	}

	// Accessors (for CSG)

	@Override
	public List<Object> getActions()
	{
		return actions;
	}

	@Override
	public int[] getIndexes(int s, int i)
	{
		return transIndexes.get(s).get(i);
	}

	@Override
	public BitSet getIndexesForPlayer(int s, int p)
	{
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

	@Override
	public String[] getActions(int s, int i)
	{
		int[] indexes = getIndexes(s, i);
		String[] result = new String[indexes.length];
		for (int a = 0; a < indexes.length; a++) {
			result[a] = (indexes[a] > 0) ? actions.get(indexes[a] - 1).toString() : "<" + a + ">";
		}
		return result;
	}

	@Override
	public BitSet[] getIndexes()
	{
		return indexes;
	}

	@Override
	public int[] getIdles()
	{
		return idles;
	}

	@Override
	public int getIdleForPlayer(int p)
	{
		return idles[p];
	}

	@Override
	public BitSet getConcurrentPlayers(int s)
	{
		BitSet result = new BitSet();
		int[] numActions = getNumActions(s);
		int numPlayers = getNumPlayers();
		for (int p = 0; p < numPlayers; p++) {
			if (numActions[p] >= 2)
				result.set(p);
		}
		return result;
	}

	// Local accessors / utility methods
	
	/**
	 * Get the list of player action indices for all states/choices
	 */
	public List<List<int[]>> getTransIndexes()
	{
		return transIndexes;
	}

	/**
	 * Get the list of player action indices for the choices of state {@code s}
	 */
	public List<int[]> getTransIndexes(int s)
	{
		return transIndexes.get(s);
	}


	public List<Object> getActionsForPlayer(int p)
	{
		List<Object> result = new ArrayList<>();
		for (int i = indexes[p].nextSetBit(0); i >= 0; i = indexes[p].nextSetBit(i + 1)) {
			result.add(actions.get(i));
		}
		return result;
	}

	public Set<String> getActionsForPlayer(int s, int p)
	{
		Set<String> result = new HashSet<String>();
		String[] actions;
		for (int t = 0; t < getNumChoices(s); t++) {
			actions = getActions(s, t);
			if (!actions[p].equals("-"))
				result.add(actions[p]);
		}
		return result;
	}

	/**
	 * Get the number of distinct (non-idle) actions taken by each player in state {@code s}.
	 */
	public int[] getNumActions(int s)
	{
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

	public BitSet getActivePlayers(int s)
	{
		BitSet result = new BitSet();
		int[] numActions = getNumActions(s);
		int numPlayers = getNumPlayers();
		for (int p = 0; p < numPlayers; p++) {
			if (numActions[p] > 0)
				result.set(p);
		}
		return result;
	}

	public void printModelInfo()
	{
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
