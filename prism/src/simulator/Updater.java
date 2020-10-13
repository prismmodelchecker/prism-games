//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <d.a.parker@cs.bham.ac.uk> (University of Birmingham/Oxford)
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

package simulator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Vector;

import parser.State;
import parser.VarList;
import parser.ast.Command;
import parser.ast.Module;
import parser.ast.ModulesFile;
import parser.ast.Update;
import parser.ast.Updates;
import prism.ModelType;
import prism.PrismComponent;
import prism.PrismException;
import prism.PrismLangException;
import prism.PrismSettings;

public class Updater extends PrismComponent
{
	// Settings:
	// Do we check that probabilities sum to 1?
	protected boolean doProbChecks = true;
	// The precision to which we check probabilities sum to 1
	protected double sumRoundOff = 1e-5;
	
	// Info on model being explored
	protected ModulesFile modulesFile;
	protected ModelType modelType;
	protected int numModules;
	protected VarList varList;
	// Synchronising action info
	protected Vector<String> synchs;
	protected int numSynchs;
	protected int synchModuleCounts[];
	// Model info/stats
	protected int numRewardStructs;

	// Temporary storage:

	// Element i,j of updateLists is a list of the updates from module i labelled with action j
	// (where j=0 denotes independent, otherwise 1-indexed action label)
	protected List<List<List<Updates>>> updateLists;
	// Bit j of enabledSynchs is set iff action j is currently enabled
	// (where j=0 denotes independent, otherwise 1-indexed action label)
	protected BitSet enabledSynchs;
	// Element j of enabledModules is a BitSet showing modules which enable action j
	// (where j=0 denotes independent, otherwise 1-indexed action label)
	protected BitSet enabledModules[];
	// Number of players
	protected int numPlayers; 
	
	// Element i of playersIndexes is the index of the player for module i
	protected int[] playersIndexes;
	// Element i of playersActionsIndexes contains the indexes of the actions for player i
	protected BitSet[] playersActionsIndexes;
	// Entry i of actionIndexPlayerMap contains a mapping from indexes to players
	protected Map<Integer,Integer> actionIndexPlayerMap;
	
	protected ArrayList<ArrayList<Set<BitSet>>> expansions;
	
	public Updater(ModulesFile modulesFile, VarList varList)
	{
		this(modulesFile, varList, null);
	}
	
	public Updater(ModulesFile modulesFile, VarList varList, PrismComponent parent)
	{
		// Store some settings
		doProbChecks = parent.getSettings().getBoolean(PrismSettings.PRISM_DO_PROB_CHECKS);
		sumRoundOff = parent.getSettings().getDouble(PrismSettings.PRISM_SUM_ROUND_OFF);
		
		// Get info from model
		this.modulesFile = modulesFile;
		modelType = modulesFile.getModelType();
		numModules = modulesFile.getNumModules();
		synchs = modulesFile.getSynchs();
		numSynchs = synchs.size();
		numRewardStructs = modulesFile.getNumRewardStructs();
		this.varList = varList;
		
		// Compute count of number of modules using each synch action
		// First, compute and cache the synch actions for each of the modules
		List<HashSet<String>> synchsPerModule = new ArrayList<HashSet<String>>(numModules);
		for (int i = 0; i < numModules; i++) {
			synchsPerModule.add(new HashSet<String>(modulesFile.getModule(i).getAllSynchs()));
		}
		// Second, do the counting
		synchModuleCounts = new int[numSynchs];
		for (int j = 0; j < numSynchs; j++) {
			synchModuleCounts[j] = 0;
			String s = synchs.get(j);
			for (int i = 0; i < numModules; i++) {
				if (synchsPerModule.get(i).contains(s))
					synchModuleCounts[j]++;
			}
		}
		// Build lists/bitsets for later use
		updateLists = new ArrayList<List<List<Updates>>>(numModules);
		for (int i = 0; i < numModules; i++) {
			updateLists.add(new ArrayList<List<Updates>>(numSynchs + 1));
			for (int j = 0; j < numSynchs + 1; j++) {
				updateLists.get(i).add(new ArrayList<Updates>());
			}
		}
		enabledSynchs = new BitSet(numSynchs + 1);
		enabledModules = new BitSet[numSynchs + 1];
		for (int j = 0; j < numSynchs + 1; j++) {
			enabledModules[j] = new BitSet(numModules);
		}
		numPlayers = modulesFile.getNumPlayers();
	}

	/**
	 * Initialise auxiliary variables for building CSGs.
	 * @throws PrismLangException
	 */
	public void initialiseCSG() throws PrismLangException {
		playersActionsIndexes = new BitSet[numPlayers];
		actionIndexPlayerMap = new HashMap<Integer, Integer>();
		expansions = new ArrayList<ArrayList<Set<BitSet>>>();
		if (numPlayers > 0) {	
			int index;
			BitSet seen = new BitSet();
			playersIndexes = new int[numModules];
			Arrays.fill(playersIndexes, -1);
			for (int p = 0; p < numPlayers; p++) {
				playersActionsIndexes[p] = new BitSet();
			}
			for (int m = 0; m < numModules; m++) {
				playersIndexes[m] = modulesFile.getPlayerForModule(modulesFile.getModuleName(m));
				expansions.add(m, new ArrayList<Set<BitSet>>());
				if (playersIndexes[m] != -1) {
					for (int c = 0; c < modulesFile.getModule(m).getNumCommands(); c++) {
						if (modulesFile.getModule(m).getCommand(c).isUnlabelled()) {
							throw new PrismLangException("Commands in a player-owned module cannot be unlabelled", modulesFile.getModule(m).getCommand(c));
						}
						expansions.get(m).add(c, new HashSet<BitSet>());
						index = modulesFile.getModule(m).getCommand(c).getSynchIndices().get(0);
						if (!seen.get(index) || playersActionsIndexes[playersIndexes[m]].get(index)) {
							seen.set(index);
							playersActionsIndexes[playersIndexes[m]].set(index);
							actionIndexPlayerMap.put(index, playersIndexes[m]);
						}
						else
							throw new PrismLangException("Action " + modulesFile.getModule(m).getCommand(c).getSynchs().get(0) + " of module "
																   + modulesFile.getModule(m).getName() 
																   + " had already been associated to a different module. Action sets must be disjoint");
					}
				}
				else {
					for (int c = 0; c < modulesFile.getModule(m).getNumCommands(); c++) {
						expansions.get(m).add(c, new HashSet<BitSet>());
					}
				}
			}
			//System.out.println("-- actionIndexPlayerMap");
			//System.out.println(actionIndexPlayerMap);
			//System.out.println("-- playersActionsIndexes");
			//System.out.println(Arrays.toString(playersActionsIndexes));
			for (int m = 0; m < numModules; m++) {
				if (playersIndexes[m] == -1) {
					for (int c = 0; c < modulesFile.getModule(m).getNumCommands(); c++) {
						for (int i : modulesFile.getModule(m).getCommand(c).getSynchIndices()) {
							if (!actionIndexPlayerMap.keySet().contains(i) && i != 0)
								throw new PrismLangException("Label \"" + modulesFile.getSynch(i - 1) + "\" of commandd " + c  
																+ " of module " + modulesFile.getModuleName(m) + " is not associated to any player." 
																+ " Independent modules can only synchronize on players' actions.");
						}
					}
				}
			}
		}
	}
	
	/**
	 * Set the precision to which we check that probabilities sum to 1.
	 */
	public void setSumRoundOff(double sumRoundOff)
	{
		this.sumRoundOff = sumRoundOff;
	}

	/**
	 * Get the precision to which we check that probabilities sum to 1.
	 */
	public double getSumRoundOff()
	{
		return sumRoundOff;
	}
	
	/**
	 * Determine the set of outgoing transitions from state 'state' and store in 'transitionList'.
	 * @param state State from which to explore
	 * @param transitionList TransitionList object in which to store result
	 */
	public void calculateTransitions(State state, TransitionList transitionList) throws PrismException
	{
		List<ChoiceListFlexi> chs;
		int i, j, k, l, n, count;
		
		// Clear lists/bitsets
		transitionList.clear();
		for (i = 0; i < numModules; i++) {
			for (j = 0; j < numSynchs + 1; j++) {
				updateLists.get(i).get(j).clear();
			}
		}
		enabledSynchs.clear();
		for (i = 0; i < numSynchs + 1; i++) {
			enabledModules[i].clear();
		}

		// Calculate the available updates for each module/action
		// (update information in updateLists, enabledSynchs and enabledModules)
		for (i = 0; i < numModules; i++) {
			calculateUpdatesForModule(i, state);
		}

		// Add independent transitions for each (enabled) module to list
		for (i = enabledModules[0].nextSetBit(0); i >= 0; i = enabledModules[0].nextSetBit(i + 1)) {
			for (Updates ups : updateLists.get(i).get(0)) {
				ChoiceListFlexi ch = processUpdatesAndCreateNewChoice(-(i + 1), ups, state);
				if (ch.size() > 0)
					transitionList.add(ch);
			}
		}
		
		// Add synchronous transitions to list
		chs = new ArrayList<ChoiceListFlexi>();
		for (i = enabledSynchs.nextSetBit(1); i >= 0; i = enabledSynchs.nextSetBit(i + 1)) {
			chs.clear();
			// Check counts to see if this action is blocked by some module
			if (enabledModules[i].cardinality() < synchModuleCounts[i - 1])
				continue;
			// If not, proceed...
			for (j = enabledModules[i].nextSetBit(0); j >= 0; j = enabledModules[i].nextSetBit(j + 1)) {
				count = updateLists.get(j).get(i).size();
				// Case where there is only 1 Updates for this module
				if (count == 1) {
					Updates ups = updateLists.get(j).get(i).get(0);
					// Case where this is the first Choice created
					if (chs.size() == 0) {
						ChoiceListFlexi ch = processUpdatesAndCreateNewChoice(i, ups, state);
						if (ch.size() > 0) {
							chs.add(ch);
						}
					}
					// Case where there are existing Choices
					else {
						// Product with all existing choices
						for (ChoiceListFlexi ch : chs) {
							processUpdatesAndAddToProduct(ups, state, ch);
						}
					}
				}
				// Case where there are multiple Updates (i.e. local nondeterminism)
				else {
					// Case where there are no existing choices
					if (chs.size() == 0) {
						for (Updates ups : updateLists.get(j).get(i)) {
							ChoiceListFlexi ch = processUpdatesAndCreateNewChoice(i, ups, state);
							if (ch.size() > 0) {
								chs.add(ch);
							}
						}
					}
					// Case where there are existing Choices
					else {
						// Duplicate (count-1 copies of) current Choice list
						n = chs.size();
						for (k = 0; k < count - 1; k++)
							for (l = 0; l < n; l++) {
								chs.add(new ChoiceListFlexi(chs.get(l)));
							}
						// Products with existing choices
						for (k = 0; k < count; k++) {
							Updates ups = updateLists.get(j).get(i).get(k);
							for (l = 0; l < n; l++) {
								processUpdatesAndAddToProduct(ups, state, chs.get(k * n + l));
							}
						}
					}
				}
			}
			// Add all new choices to transition list
			for (ChoiceListFlexi ch : chs) {
				transitionList.add(ch);
			}
		}
		
		// For a DTMC, we need to normalise across all transitions
		// This is partly to handle "local nondeterminism"
		// and also to handle any dubious trickery done by disabling probability checks
		if (modelType == ModelType.DTMC) {
			double probSum = transitionList.getProbabilitySum();
			transitionList.scaleProbabilitiesBy(1.0 / probSum);
		}
	
		// Check validity of the computed transitions
		// (not needed currently)
		//transitionList.checkValid(modelType);
		
		// Check for errors (e.g. overflows) in the computed transitions
		//transitionList.checkForErrors(state, varList);
		
		//System.out.println(transitionList);
	}
	
	public void calculateEnabledCommands(int m, BitSet active, State state) throws PrismLangException 
	{
		//System.out.println("\n## state " + state);
		//System.out.println("\n## calculateEnabledCommands, module " + m);
		Module module;
		Command command;
		BitSet tmp = new BitSet();
		BitSet indexes = new BitSet();
		HashSet<BitSet> seen = new HashSet<BitSet>();
		int i, e, n, p;
		
		module = modulesFile.getModule(m);
		n = module.getNumCommands();
		e = -1;

		for (i = 0; i < n; i++) {
			command = module.getCommand(i);
			expansions.get(m).get(i).clear();
			if (command.getSynchIndices().get(0) == 0) {
				p = playersIndexes[m];
				if (p != -1) {
					throw new PrismLangException("Module " + modulesFile.getModuleName(m) 
				   										   + " from to player " + p
				   										   + " has an unlabelled command");
				}
				else if (command.getGuard().evaluateBoolean(state)) {
					if (e == -1) {
						active.set(i);
						e = i;
					}
					else 
						throw new PrismLangException("Module " + modulesFile.getModuleName(m) 
															   + " has multiple unlabeled active commands in state " + state);
				}
			}
			else {
				if (command.getGuard().evaluateBoolean(state)) {
					indexes.clear();
					for(int j : command.getSynchIndices()) {
						indexes.set(j);
					}
					for (p = 0; p < numPlayers; p++) {
						tmp.clear();
						tmp.or(indexes);
						tmp.andNot(playersActionsIndexes[p]);
						if (tmp.cardinality() < indexes.cardinality() - 1)
							throw new PrismLangException("Module " + modulesFile.getModuleName(m) + " has multiple actions associated to player " 
															   	   + modulesFile.getPlayerName(p) + " in command " + i);
					}
					if (!seen.contains(indexes)) {
						seen.add((BitSet) indexes.clone());
						active.set(i);
					}
					else 
						throw new PrismLangException("Module " + modulesFile.getModuleName(m) + " has multiple active commands labeled "
															   + command.getSynchs() + " in state " + state);
				}
			}
		}
	}
		
	/**
	 * Makes the product of all indexes in indexes[p] (for player p) and stores the result in products.   
	 * @param products Set where the products are stored. 
	 * @param indexes Array of bitsets with the indexes of actions for each player. 
	 * @param product A single product, passed recursively.  
	 * @param p Index of player p
	 */
	public void indexProduct(Set<BitSet> products, BitSet[] indexes, BitSet product, int p) 
	{
		BitSet newproduct;
		int i;
		if (p < numPlayers - 1) {
			if (!indexes[p].isEmpty()) {
				for (i = indexes[p].nextSetBit(0); i >= 0; i = indexes[p].nextSetBit(i + 1)) {	
					newproduct = new BitSet();
					newproduct.or(product);
					newproduct.set(i);
					indexProduct(products, indexes, newproduct, p + 1);
				}
			}
			else {
				newproduct = new BitSet();
				newproduct.or(product);
				indexProduct(products, indexes, newproduct, p + 1);	
			}
		}
		else {
			if (!indexes[p].isEmpty()) {
				for (i = indexes[p].nextSetBit(0); i >= 0; i = indexes[p].nextSetBit(i + 1)) {	
					newproduct = new BitSet();
					newproduct.or(product);
					newproduct.set(i);
					products.add(newproduct);
				}
			}
			else {
				newproduct = new BitSet();
				newproduct.or(product);
				products.add(newproduct);
			}
		}
	}	
	
	public void expandCommand(Set<BitSet> products, int m, int c)
	{
		BitSet bidx = new BitSet();
		BitSet tmp;
		for (int i : modulesFile.getModule(m).getCommand(c).getSynchIndices())
			if (i != 0)
				bidx.set(i);
		for (BitSet prod : products) {
			tmp = new BitSet();
			tmp.or(bidx);
			tmp.andNot(prod);
			if (tmp.isEmpty()) {
				expansions.get(m).get(c).add(prod);
			}
		}
	}
	
	public void calculateTransitionsCSG(State state, TransitionList transitionList) throws PrismLangException 
	{
		//System.out.println("\n## state " + state);
		//System.out.println("## playersActionsIndexes " + Arrays.toString(playersActionsIndexes));
		//System.out.println("## actionIndexPlayerMap " + actionIndexPlayerMap);
		List<ChoiceListFlexi> chs = new ArrayList<ChoiceListFlexi>(); // used to store the final choices from a given state
		List<Integer> indexes; // used to store the indexes of commands at different points
		Set<BitSet> products = new HashSet<BitSet>(); // used to store store products of indexes at different points
		Set<BitSet> synchs;
		BitSet[] active = new BitSet[numModules];
		BitSet[] moves = new BitSet[numPlayers]; // indexes of possible actions taken by each player
		BitSet[] enabled  = new BitSet[numPlayers]; // number of actions which will be actually finally available for each player, used only for checking
		BitSet tmp = new BitSet();
		ChoiceListFlexi chfl;
		Command cmd1, cmd2;
		Updates ups; // (temp) used to store updates at different points
		String warning, missing;
		int[] actions; // (temp) used to store ordered indexes at different points
		int cidx, i, id, j, m, n, nchs, p;
		boolean[] msynch = new boolean[numModules];
		transitionList.clear();				
		Arrays.fill(moves, null);
		for (m = 0; m < numModules; m++) {
			active[m] = new BitSet();
			calculateEnabledCommands(m, active[m], state);
			p = playersIndexes[m];
			if (p != -1) {
				if (moves[p] == null)
					moves[p] = new BitSet();
				for (cidx = active[m].nextSetBit(0); cidx >= 0; cidx = active[m].nextSetBit(cidx + 1)) {	
					indexes = modulesFile.getModule(m).getCommand(cidx).getSynchIndices();
					moves[p].set(indexes.get(0));			
				}
			}
		}
		//System.out.println("-- active " + Arrays.toString(active));
		//System.out.println("-- moves " + Arrays.toString(moves));
		indexProduct(products, moves, new BitSet(), 0);
		for (m = 0; m < numModules; m++) {
			for (i = active[m].nextSetBit(0); i >= 0; i = active[m].nextSetBit(i + 1)) {	
				//System.out.println(modulesFile.getModule(m).getCommand(i));
				expandCommand(products, m, i);
			}
		}
		for (m = 0; m < numModules; m++) {
			for (i = active[m].nextSetBit(0); i >= 0; i = active[m].nextSetBit(i + 1)) {	
				cmd1 = modulesFile.getModule(m).getCommand(i);
				if (!(cmd1.getSynchIndices().get(0) == 0)) {
					for (j = active[m].nextSetBit(0); j >= 0; j = active[m].nextSetBit(j + 1)) {	
						cmd2 = modulesFile.getModule(m).getCommand(j);
						if (i != j) {
							for (BitSet prod : expansions.get(m).get(i)) {
								if (!(cmd2.getSynchIndices().get(0) == 0)) {
									if (cmd2.getSynchIndices().get(0) == cmd1.getSynchIndices().get(0)) {
										if (cmd2.getSynchIndices().size() < cmd1.getSynchIndices().size()) {
											expansions.get(m).get(j).remove(prod);
										}
									}
								}
								else { 
									expansions.get(m).get(j).remove(prod);
								}
							}
						}
					}
				}
			}
		}
		//System.out.println("-- products " + products);
		//System.out.println("-- expansions " + expansions);						
		synchs = new HashSet<BitSet>(products);
		for (BitSet prod : products) {
			tmp.clear();
			tmp.or(prod);
			//System.out.println("-- product " + prod);
			for (m = 0; m < numModules; m++) {
				if (playersIndexes[m] != -1) {
					if (active[m].isEmpty())
						continue;
					for (i = active[m].nextSetBit(0); i >= 0; i = active[m].nextSetBit(i + 1)) {	
						if (expansions.get(m).get(i).contains(prod)) {
							tmp.clear(modulesFile.getModule(m).getCommand(i).getSynchIndices().get(0));
						}
					}
				}
			}
			//System.out.println("-- tmp " + tmp);
			if (!tmp.isEmpty()) {
				warning = "Missing specification for action product [";
				missing = "";
				for (id = prod.nextSetBit(0); id >= 0; id = prod.nextSetBit(id + 1)) {	
					missing += modulesFile.getSynch(id - 1) + ((prod.nextSetBit(id + 1) > 0)? "," : "]");
				}
				warning += missing;
				warning += " in state " + state + ".";
				mainLog.printWarning(warning);
				synchs.remove(prod);
			}
		}
		nchs = 0;
		//System.out.println("-- active " + Arrays.toString(active));
		for (BitSet prod : synchs) {
			actions = new int[numPlayers];
			Arrays.fill(actions, -1);
			chfl = null;
			//System.out.println("\n-- prod");
			//System.out.println(prod);
			for (m = 0; m < numModules; m++) {
				msynch[m] = false;
				//System.out.println("--\n module " + m + " " + modulesFile.getModuleName(m));
				for (i = active[m].nextSetBit(0); i >= 0; i = active[m].nextSetBit(i + 1)) {	
					ups = modulesFile.getModule(m).getCommand(i).getUpdates();
					//System.out.println(modulesFile.getModule(m).getCommand(i));
					//System.out.println("-- expansions " + i);
					//System.out.println(expansions.get(m).get(i));
					//System.out.println(modulesFile.getModule(m).getCommand(i).getSynchIndices());
					if (expansions.get(m).get(i).contains(prod)) {
						if (!msynch[m]) {
							//System.out.println("## selected " + i);
							if (chfl == null) {
								chfl = processUpdatesAndCreateNewChoice(nchs, ups, state);
								nchs++;
							}
							else
								processUpdatesAndAddToProduct(ups, state, chfl);
							msynch[m] = true;
						}
						else {
							//System.out.println();
							//System.out.println(prod);
							throw new PrismLangException("Module " + modulesFile.getModuleName(m) + " has multiple active commands for action "
																								  + "\'" + modulesFile.getModule(m).getCommand(i).getSynch() + "\'" + " in state " + state);
						}
					}
				}
			}
			for (i = prod.nextSetBit(0); i >= 0; i = prod.nextSetBit(i + 1)) {
				actions[actionIndexPlayerMap.get(i)] = i;
			}
			if (chfl != null) {
				//System.out.println("-- chfl " + chfl);
				chfl.setActions(actions);
				chs.add(chfl);
			}
		}	
		n = 1;
		for (p = 0; p < numPlayers; p++) {
			enabled[p] = new BitSet();
			for (ChoiceListFlexi ch : chs) { 
				// Setting to size + 1 in the case of idle actions as they are given index -1 which cannot be set
				enabled[p].set((ch.getActions()[p] > 0)? ch.getActions()[p] : modulesFile.getSynchs().size() + 1); 
			}
			n *= enabled[p].cardinality();
		}				
		products.clear();
		/*
		System.out.println("\n-- chs " + chs);
		System.out.println(modulesFile.getSynchs());
		for (ChoiceListFlexi ch : chs) {
			System.out.println(ch);
			System.out.println(Arrays.toString(ch.getActions()));
		}	
		*/
		if (n != chs.size()) {
			indexProduct(products, enabled, new BitSet(), 0);
			BitSet tmp1 = new BitSet();
			for (ChoiceListFlexi ch : chs) { 
				tmp1.clear();
				for (int idx : ch.getActions()) {
					tmp1.set(idx);
				}
				products.remove(tmp1);
			}
			for (BitSet tmp2 : products) {
				warning = "Missing specification for action product [";
				missing = "";
				for (id = tmp2.nextSetBit(0); id >= 0; id = tmp2.nextSetBit(id + 1)) {	
					missing += modulesFile.getSynch(id - 1) + ((tmp2.nextSetBit(id + 1) > 0)? "," : "]");
				}
				warning += missing;
				warning += " in state " + state + ".";
				mainLog.printWarning(warning);
			}
			throw new PrismLangException("Error in model specification.");
		}
		for (ChoiceListFlexi ch : chs) {
			transitionList.add(ch, ch.getActions());
		}	
	}
		
	// Private helpers
	
	/**
	 * Determine the enabled updates for the 'm'th module from (global) state 'state'.
	 * Update information in updateLists, enabledSynchs and enabledModules.
	 * @param m The module index
	 * @param state State from which to explore
	 */
	protected void calculateUpdatesForModule(int m, State state) throws PrismLangException
	{
		Module module;
		Command command;
		int i, j, n;

		module = modulesFile.getModule(m);
		n = module.getNumCommands();
		for (i = 0; i < n; i++) {
			command = module.getCommand(i);
			if (command.getGuard().evaluateBoolean(state)) {
				j = command.getSynchIndex();
				updateLists.get(m).get(j).add(command.getUpdates());
				enabledSynchs.set(j);
				enabledModules[j].set(m);
			}
		}
	}

	/**
	 * Create a new Choice object (currently ChoiceListFlexi) based on an Updates object
	 * and a (global) state. Check for negative probabilities/rates and, if appropriate,
	 * check probabilities sum to 1 too.
	 * @param moduleOrActionIndex Module/action for the choice, encoded as an integer (see Choice)
	 * @param ups The Updates object 
	 * @param state Global state
	 */
	private ChoiceListFlexi processUpdatesAndCreateNewChoice(int moduleOrActionIndex, Updates ups, State state) throws PrismLangException
	{
		ChoiceListFlexi ch;
		List<Update> list;
		int i, n;
		double p, sum;

		// Create choice and add all info
		ch = new ChoiceListFlexi();
		ch.setModuleOrActionIndex(moduleOrActionIndex);
		n = ups.getNumUpdates();
		sum = 0;
		for (i = 0; i < n; i++) {
			// Compute probability/rate
			p = ups.getProbabilityInState(i, state);
			// Check for non-finite/NaN probabilities/rates
			if (!Double.isFinite(p) || p < 0) {
				String s = modelType.choicesSumToOne() ? "Probability" : "Rate";
				s += " is invalid (" + p + ") in state " + state.toString(modulesFile);
				// Note: we indicate error in whole Updates object because the offending
				// probability expression has probably been simplified from original form.
				throw new PrismLangException(s, ups);
			}
			// Skip transitions with zero probability/rate
			if (p == 0)
				continue;
			sum += p;
			list = new ArrayList<Update>();
			list.add(ups.getUpdate(i));
			ch.add(p, list);
		}
		// For now, PRISM treats empty (all zero probs/rates) distributions as an error.
		// Later, when errors in symbolic model construction are improved, this might be relaxed.
		if (ch.size() == 0) {
			String msg = modelType.probabilityOrRate();
			msg += (ups.getNumUpdates() > 1) ? " values sum to " : " is ";
			msg += "zero for updates in state " + state.toString(modulesFile);
			throw new PrismLangException(msg, ups);
		}
		// Check distribution sums to 1 (if required, and if is non-empty)
		if (doProbChecks && ch.size() > 0 && modelType.choicesSumToOne() && Math.abs(sum - 1) > sumRoundOff) {
			throw new PrismLangException("Probabilities sum to " + sum + " in state " + state.toString(modulesFile), ups);
		}
		return ch;
	}

	/**
	 * Create a new Choice object (currently ChoiceListFlexi) based on the product
	 * of an existing ChoiceListFlexi and an Updates object, for some (global) state.
	 * If appropriate, check probabilities sum to 1 too.
	 * @param ups The Updates object 
	 * @param state Global state
	 * @param ch The existing Choices object
	 */
	private void processUpdatesAndAddToProduct(Updates ups, State state, ChoiceListFlexi ch) throws PrismLangException
	{
		// Create new choice (action index is 0 - not needed)
		ChoiceListFlexi chNew = processUpdatesAndCreateNewChoice(0, ups, state);
		// Build product with existing
		ch.productWith(chNew);
	}
}
