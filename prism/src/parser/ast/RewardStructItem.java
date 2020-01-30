//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford, formerly University of Birmingham)
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

package parser.ast;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import parser.visitor.*;
import prism.PrismLangException;

public class RewardStructItem extends ASTElement
{
	// Synchronising action(s):
	// * null = none (i.e. state reward)
	// * otherwise, action-label(s) for a transition reward,
	//   usually exactly 1, but can be more for concurrent games: 
	//   -  "" = empty/unlabelled/asynchronous action
	//   -  "act" = "act"-labelled action
	private List<String> synchs;
	// Index of action label in model's list of all actions ("synchs")
	// This is 1-indexed, with 0 denoting independent (unlabelled).
	// -1 denotes either none (i.e. state reward, synch==null) or not (yet) known.
	private List<Integer> synchIndices;
	// Guard expression
	private Expression states;
	// Reward expression
	private Expression reward;
	
	// Constructors
	
	public RewardStructItem(String synch, Expression states, Expression reward)
	{
		setSynch(synch);
		this.states = states;
		this.reward = reward;
	}
	
	public RewardStructItem(List<String> synchs, Expression states, Expression reward)
	{
		setSynchs(synchs);
		this.states = states;
		this.reward = reward;
	}
	
	// Set methods
	
	/**
	 * Set the synchronising action label for this reward struct item
	 * (as a string; if it is unlabelled, then this is "").
	 * If it is null, this means there is no action label, i.e. a state reward.
	 * This is the method normally used for setting this;
	 * {@link #setSynchs(List)} is for the multi-action case. 
	 */
	public void setSynch(String synch)
	{
		// Null case
		if (synch == null) {
			synchs = null;
			synchIndices = null;
		}
		// Otherwise reset to a size 1 list containing just s
		synchs = new ArrayList<String>(1);
		synchs.add(synch);
		// Also update synchIndices to matching size list (value -1)
		synchIndices = new ArrayList<Integer>(1);
		synchIndices.add(-1);
	}
	
	/**
	 * Set multiple synchronising action labels for this reward struct item
	 * (each as a string; if it is unlabelled, then this is "").
	 * These are passed in as a list, which is stored directly, not copied.
	 * If it is null, this means there is no action label, i.e. a state reward.
	 * An empty list is treated as a singleton list containing "".
	 */
	public void setSynchs(List<String> synchs)
	{
		// Null case
		if (synchs == null) {
			this.synchs = null;
			synchIndices = null;
		}
		// Otherwise store (checking for empty list)
		// NB: Also update synchIndices to matching size list (value -1)
		else {
			if (synchs.isEmpty()) {
				this.synchs = new ArrayList<String>(1);
				this.synchs.add("");
				synchIndices = new ArrayList<Integer>(1);
				synchIndices.add(-1);
			} else {
				this.synchs = synchs;
				synchIndices = new ArrayList<>(Collections.nCopies(synchs.size(), -1));
			}
		}
	}
	
	/**
	 * Find and cache the index of any action labels,
	 * using a passed in list of all synchronising actions for the index.
	 * The cached index starts from 1; 0 means unlabelled, -1 means unknown.
	 * Throws an exception if an action cannot be found.
	 */
	public void setSynchIndices(List<String> allSynchs) throws PrismLangException
	{
		// State rewards - nothing to do
		if (synchs == null) {
			return;
		}
		// Transition rewards:
		int numSynchs = synchs.size();
		for (int i = 0; i < numSynchs; i++) {
			String synch = synchs.get(i);
			// For independent actions, the index is 0
			if (synch.equals("")) {
				synchIndices.set(i, 0);
				continue;
			}
			// Otherwise, see if action name exists
			int j = allSynchs.indexOf(synch);
			if (j != -1) {
				// If so, set the index (starts from 1)
				synchIndices.set(i, j + 1);
				continue;
			}
			// Otherwise, there is a problem.
			throw new PrismLangException("Unknown action name " + synch + " in reward structure item", this);
		}
	}
	
	public void setStates(Expression states)
	{
		this.states = states;
	}
	
	public void setReward(Expression reward)
	{
		this.reward = reward;
	}
	
	// Get methods
	
	/**
	 * Get the action label for this reward struct item, if present.
	 * If none (i.e. a state reward), it is null.
	 * For a transition reward, it is a string;
	 * "" denotes the independent (unlabelled) case.
	 */
	public String getSynch()
	{
		// If non-null, the list should never be non-empty
		return synchs == null ? null : synchs.get(0);
	}
	
	/**
	 * Get the list of action labels for this reward struct item, if present.
	 * If none (i.e. a state reward), it is null.
	 * For a transition reward, it is a list of strings;
	 * "" denotes the independent (unlabelled) case.
	 * Usually (apart from concurrent games), there is a single
	 * action and you can just use {@link #getSynch()}.
	 */
	public List<String> getSynchs() 
	{
		return synchs;
	}
		
	/**
	 * Get the index of the action label for this reward struct item (in the model's list of actions).
	 * This is 1-indexed, with 0 denoting the independent (unlabelled) case.
	 * -1 denotes either none (i.e. state reward, synch==null) or not (yet) known.
	 */
	public int getSynchIndex()
	{
		// If synchs is non-null, the list should never be non-empty
		return synchs == null ? -1 : synchIndices.get(0);
	}
	
	/**
	 * Get the list of indices for all action labels for this reward struct item (in the model's list of actions).
	 * Each is 1-indexed, with 0 denoting the independent (unlabelled) case.
	 * -1 denotes either none (i.e. state reward, synch==null) or not (yet) known.
	 */
	public List<Integer> getSynchIndices() 
	{
		return synchIndices;
	}
	
	public Expression getStates()
	{
		return states;
	}
	
	public Expression getReward()
	{
		return reward;
	}
	
	/**
	 *	Returns whether this reward is a state (false) or transition (true) reward
	 */
	public boolean isTransitionReward() 
	{
		return (synchs != null);
	}

	// Methods required for ASTElement:
	
	/**
	 * Visitor method.
	 */
	public Object accept(ASTVisitor v) throws PrismLangException 
	{
		return v.visit(this);
	}
	
	/**
	 * Convert to string.
	 */
	public String toString()
	{
		String s = "";
		if (synchs != null) {
			s += "[" + String.join(",", synchs) + "] ";
		}
		s += states + " : " + reward + ";";	
		return s;
	}
	
	/**
	 * Perform a deep copy.
	 */
	public ASTElement deepCopy() 
	{
		List<String> synchsCopy = synchs == null ? null : new ArrayList<String>(synchs);
		RewardStructItem ret = new RewardStructItem(synchsCopy, states.deepCopy(), reward.deepCopy());
		ret.synchIndices = synchIndices == null ? null : new ArrayList<Integer>(synchIndices);
		ret.setPosition(this);
		return ret;
	}
}

//------------------------------------------------------------------------------
