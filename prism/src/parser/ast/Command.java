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

public class Command extends ASTElement
{
	// Action label(s)
	// Usually just 1 label, but can be more for concurrent games
	// There is always at least one string in the list. We use "" for the unlabelled case.  
	private List<String> synchs;
	// Cached indices of each action label in the model's list of all actions ("synchs")
	// This is 1-indexed, with 0 denoting independent (unlabelled).
	// -1 denotes not (yet) known.
	private List<Integer> synchIndices;
	// Guard
	private Expression guard;
	// List of updates
	private Updates updates;
	// Parent module
	private Module parent;
		
	// Constructor
	
	public Command()
	{
		setSynch("");
		guard = null;
		updates = null;
	}
	
	// Set methods
	
	/**
	 * Set the synchronising action label for this command
	 * (as a string; if it is unlabelled, then this is "", not null).
	 * This is the method normally used for setting this;
	 * {@link #setSynchs(List)} is for the multi-action case. 
	 */
	public void setSynch(String synch)
	{
		// Reset to a size 1 list containing just s
		synchs = new ArrayList<String>(1);
		synchs.add(synch);
		// Also update synchIndices to matching size list (value -1)
		synchIndices = new ArrayList<Integer>(1);
		synchIndices.add(-1);
	}
	
	/**
	 * Set multiple synchronising action labels for this command
	 * (each as a string; if it is unlabelled, then this is "").
	 * These are passed in as a list, which is stored directly, not copied.
	 * An empty list is treated as a singleton list containing "".
	 */
	public void setSynchs(List<String> synchs)
	{
		// NB: Also update synchIndices to matching size list (value -1)
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
	
	/**
	 * Find and cache the index of any action labels,
	 * using a passed in list of all synchronising actions for the index.
	 * The cached index starts from 1; 0 means unlabelled, -1 means unknown.
	 * Throws an exception if an action cannot be found.
	 */
	public void setSynchIndices(List<String> allSynchs) throws PrismLangException
	{
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
			throw new PrismLangException("Unknown action name " + synch + " in command", this);
		}
	}
	
	public void setGuard(Expression g)
	{
		guard = g;
	}
	
	public void setUpdates(Updates u)
	{
		updates = u;
		u.setParent(this);
	}
	
	public void setParent(Module m)
	{
		parent = m;
	}

	// Get methods
	
	/**
	 * Get the action label for this command. For independent (unlabelled) commands,
	 * this is the empty string "" (it should never be null).
	 */
	public String getSynch()
	{
		// The list should never be non-empty
		return synchs.get(0);
	}
	
	/**
	 * Get the list of action labels for this command, if present.
	 * It is a list of strings; "" denotes the independent (unlabelled) case.
	 * Usually (apart from concurrent games), there is a single
	 * action and you can just use {@link #getSynch()}.
	 */
	public List<String> getSynchs() 
	{
		return synchs;
	}
	
	/**
	 * Get the index of the action label for this command, if present
	 * (in the model's list of actions). It is 1-indexed, with 0 denoting
	 * the independent (unlabelled) case. -1 denotes not (yet) known.
	 */
	public int getSynchIndex()
	{
		// The list should never be non-empty
		return synchIndices.get(0);
	}
	
	/**
	 * Get the list of indices for the action labels for this command, if present
	 * (in the model's list of actions). Each is 1-indexed, with 0 denoting
	 * the independent (unlabelled) case. -1 denotes not (yet) known.
	 * Usually (apart from concurrent games), there is a single
	 * action and you can just use {@link #getSynchIndex()}.
	 */
	public ArrayList<Integer> getSynchIndices() 
	{
		return (ArrayList<Integer>) synchIndices;
	}
	
	/**
	 * Returns true is this is an unlabelled command ("[] ..."),
	 * i.e. there is no synchronous action label attached to it.
	 */
	public boolean isUnlabelled()
	{
		// NB: this works even for concurrent games,
		// because the modelling language requires either a non-empty list
		// of action labels or a single empty action ("")  
		return "".equals(getSynch());
	}
	
	public Expression getGuard()
	{
		return guard;
	}
	
	public Updates getUpdates()
	{
		return updates;
	}
	
	public Module getParent()
	{
		return parent;
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
		return "[" + String.join(",", synchs) + "] " + guard + " -> " + updates;
	}
	
	/**
	 * Perform a deep copy.
	 */
	public ASTElement deepCopy()
	{
		Command ret = new Command();
		ret.synchs = new ArrayList<String>(getSynchs());
		ret.synchIndices = new ArrayList<Integer>(getSynchIndices());
		ret.setGuard(getGuard().deepCopy());
		ret.setUpdates((Updates)getUpdates().deepCopy());
		ret.setPosition(this);
		return ret;
	}
}

//------------------------------------------------------------------------------
