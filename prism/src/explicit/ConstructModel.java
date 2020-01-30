//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford)
//	* Vojtech Forejt <vojtech.forejt@cs.ox.ac.uk> (University of Oxford)
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

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import parser.State;
import parser.Values;
import parser.VarList;
import prism.ModelGenerator;
import prism.ModelType;
import prism.Prism;
import prism.PrismComponent;
import prism.PrismException;
import prism.PrismLog;
import prism.PrismNotSupportedException;
import prism.PrismPrintStreamLog;
import prism.ProgressDisplay;
import prism.UndefinedConstants;
import simulator.ModulesFileModelGenerator;

/**
 * Class to perform explicit-state reachability and model construction.
 * The information about the model to be built is provided via a {@link prism.ModelGenerator} interface.
 * To build a PRISM model, use {@link simulator.ModulesFileModelGenerator}.
 */
public class ConstructModel extends PrismComponent
{
	// Options:

	/** Find deadlocks during model construction? */
	protected boolean findDeadlocks = true;
	/** Automatically fix deadlocks? */
	protected boolean fixDeadlocks = true;
	/** Sort the reachable states before constructing the model? */
	protected boolean sortStates = true;
	/** Build a sparse representation, if possible?
	 *  (e.g. MDPSparse rather than MDPSimple data structure) */
	protected boolean buildSparse = true;
	/** Should actions be attached to distributions (and used to distinguish them)? */
	protected boolean distinguishActions = true;
	/** Should labels be processed and attached to the model? */
	protected boolean attachLabels = true;

	// Details of built model:

	/** Reachable states */
	protected List<State> statesList;

	public ConstructModel(PrismComponent parent) throws PrismException
	{
		super(parent);
	}

	/**
	 * Get the list of states associated with the last model construction performed.  
	 */
	public List<State> getStatesList()
	{
		return statesList;
	}

	/**
	 * Automatically fix deadlocks, if needed?
	 * (by adding self-loops in those states)
	 */
	public void setFixDeadlocks(boolean fixDeadlocks)
	{
		this.fixDeadlocks = fixDeadlocks;
	}

	/**
	 * Sort the reachable states before constructing the model?
	 */
	public void setSortStates(boolean sortStates)
	{
		this.sortStates = sortStates;
	}

	/**
	 * Build a sparse representation, if possible?
	 * (e.g. MDPSparse rather than MDPSimple data structure)
	 */
	public void setBuildSparse(boolean buildSparse)
	{
		this.buildSparse = buildSparse;
	}

	/**
	 * Should actions be attached to distributions (and used to distinguish them)?
	 */
	public void setDistinguishActions(boolean distinguishActions)
	{
		this.distinguishActions = distinguishActions;
	}

	/**
	 * Should labels be processed and attached to the model?
	 */
	public void setAttachLabels(boolean attachLabels)
	{
		this.attachLabels = attachLabels;
	}

	/**
	 * Build the set of reachable states for a model and return it.
	 * @param modelGen The ModelGenerator interface providing the model 
	 */
	public List<State> computeReachableStates(ModelGenerator modelGen) throws PrismException
	{
		constructModel(modelGen, true);
		return getStatesList();
	}

	/**
	 * Construct an explicit-state model and return it.
	 * @param modelGen The ModelGenerator interface providing the model 
	 */
	public Model constructModel(ModelGenerator modelGen) throws PrismException
	{
		return constructModel(modelGen, false);
	}

	/**
	 * Construct an explicit-state model and return it.
	 * If {@code justReach} is true, no model is built and null is returned;
	 * the set of reachable states can be obtained with {@link #getStatesList()}.
	 * @param modelGen The ModelGenerator interface providing the model 
	 * @param justReach If true, just build the reachable state set, not the model
	 *
	 */
	public Model constructModel(ModelGenerator modelGen, boolean justReach) throws PrismException
	{
		// Model info
		ModelType modelType;
		// State storage
		StateStorage<State> states;
		LinkedList<State> explore;
		State state, stateNew;
		// Explicit model storage
		ModelSimple modelSimple = null;
		DTMCSimple dtmc = null;
		CTMCSimple ctmc = null;
		MDPSimple mdp = null;
		CTMDPSimple ctmdp = null;
		STPGExplicit stpg = null;
		CSG csg = null;
		SMG smg = null;
		ModelExplicit model = null;
		Distribution distr = null;
		// Misc
		int i, j, nc, nt, src, dest, player;
		long timer;

		// Get model info
		modelType = modelGen.getModelType();

		// Display a warning if there are unbounded vars
		VarList varList = modelGen.createVarList();
		if (modelGen.containsUnboundedVariables())
			mainLog.printWarning("Model contains one or more unbounded variables: model construction may not terminate");

		// Starting reachability...
		mainLog.print("\nComputing reachable states...");
		mainLog.flush();
		ProgressDisplay progress = new ProgressDisplay(mainLog);
		progress.start();
		timer = System.currentTimeMillis();

		HashMap<Integer, String> playerNames = new HashMap<Integer, String>();
		
		String[] players = new String[modelGen.getNumPlayers()];
		int[] indexes = new int[modelGen.getNumPlayers()];
		
		// Create model storage
		if (!justReach) {
			// Create a (simple, mutable) model of the appropriate type
			switch (modelType) {
			case DTMC:
				modelSimple = dtmc = new DTMCSimple();
				dtmc.setVarList(varList);
				break;
			case CTMC:
				modelSimple = ctmc = new CTMCSimple();
				ctmc.setVarList(varList);
				break;
			case CSG:
				// Add player info
				players = new String[modelGen.getNumPlayers()];
				for (i = 0; i < modelGen.getNumPlayers(); i++) {
					players[i] = modelGen.getPlayer(i).getName();
				}
				Arrays.fill(indexes, -1);
				modelSimple = csg = new CSG(players);
				csg.setActions(((ModulesFileModelGenerator) modelGen).getModulesFile().getSynchs()); // not very pretty
				csg.setVarList(varList);
				break;			
			case MDP:
				modelSimple = mdp = new MDPSimple();
				mdp.setVarList(varList);
				break;
			case CTMDP:
				modelSimple = ctmdp = new CTMDPSimple();
				ctmdp.setVarList(varList);
				break;
			case STPG:
				modelSimple = stpg = new STPGExplicit();
				stpg.setVarList(varList);
				break;
			case SMG:
				modelSimple = smg = new SMG();
				smg.setVarList(varList);
				// Add player info
				for (i = 0; i < modelGen.getNumPlayers(); i++) {
					playerNames.put(i + 1, modelGen.getPlayer(i).getName());
				}
				smg.setPlayerInfo(playerNames);
				break;
			case PTA:
			case LTS:
				throw new PrismNotSupportedException("Model construction not supported for " + modelType + "s");
			}
		}
		
		// Initialise states storage
		states = new IndexedSet<State>(true);
		explore = new LinkedList<State>();
		// Add initial state(s) to 'explore', 'states' and to the model
		for (State initState : modelGen.getInitialStates()) {
			explore.add(initState);
			states.add(initState);
			if (!justReach) {
				modelSimple.addState();
				modelSimple.addInitialState(modelSimple.getNumStates() - 1);
			}
		}

		// Explore...
		src = -1;
		while (!explore.isEmpty()) {
			// Pick next state to explore
			// (they are stored in order found so know index is src+1)
			state = explore.removeFirst();
			src++;
			// Explore all choices/transitions from this state
			modelGen.exploreState(state);
			nc = modelGen.getNumChoices();
			// For games, first determine which player owns the state
			if (modelType.multiplePlayers() && !(modelType == ModelType.CSG)) {
				player = -1;
				for (i = 0; i < nc; i++) {
					int iPlayer = modelGen.getPlayerNumberForChoice(i);
					if (player != -1 && iPlayer != -1 && iPlayer != player) {
						throw new PrismException("Choices for both player " + player + " and " + iPlayer + " in state " + state);
					}
					if (iPlayer != -1) {
						player = iPlayer;
					}
				}
				// Assign deadlock states to player 1
				if (nc == 0) {
					player = 1;
				}
				// No assigned player: only allowed when the state is deterministic (one choice)
				// (in which case, assign the state to player 1)
				if (player == -1) {
					if (nc == 1) {
						player = 1;
					} else {
						Set<String> acts = new HashSet<>();
						for (i = 0; i < nc; i++) {
							int iPlayer = modelGen.getPlayerNumberForChoice(i);
							if (iPlayer == -1) {
								acts.add(modelGen.getChoiceAction(i).toString());
							}
						}
						String errMsg = "There are multiple choices (" + String.join(",", acts) +  ") in state " + state + " not assigned to any player"; 
						throw new PrismException(errMsg);
					}
				}
				// Make sure a valid player owns the state
				if (player < 1 || player > modelGen.getNumPlayers()) {
					throw new PrismException("State " + state + " owned by invalid player (" + player + ")");
				}
				if (modelType == ModelType.STPG) {
					stpg.setPlayer(src, player);
				} else if (modelType == ModelType.SMG) {
					smg.setPlayer(src, player);
				}
			}
			else if (modelType == ModelType.CSG) {
				csg.setPlayer(src, 1);
			}
			// Look at each outgoing choice in turn
			for (i = 0; i < nc; i++) {
				// For nondet models, collect transitions in a Distribution
				if (!justReach && modelType.nondeterministic()) {
					distr = new Distribution();
				}
				// Look at each transition in the choice
				nt = modelGen.getNumTransitions(i);
				for (j = 0; j < nt; j++) {
					stateNew = modelGen.computeTransitionTarget(i, j);
					// Is this a new state?
					if (states.add(stateNew)) {
						// If so, add to the explore list
						explore.add(stateNew);
						// And to model
						if (!justReach) {
							modelSimple.addState();
						}
					}
					// Get index of state in state set
					dest = states.getIndexOfLastAdd();
					// Add transitions to model
					if (!justReach) {
						switch (modelType) {
						case DTMC:
							dtmc.addToProbability(src, dest, modelGen.getTransitionProbability(i, j));
							break;
						case CTMC:
							ctmc.addToProbability(src, dest, modelGen.getTransitionProbability(i, j));
							break;
						case CSG:
							indexes = modelGen.getTransitionIndexes(i);
						case MDP:
						case CTMDP:
						case STPG:
						case SMG:
							distr.add(dest, modelGen.getTransitionProbability(i, j));
							break;
						case PTA:
						case LTS:
							throw new PrismNotSupportedException("Model construction not supported for " + modelType + "s");
						}
					}
				}
				// For nondet models, add collated transition to model 
				if (!justReach) {
					if (modelType == ModelType.MDP) {
						if (distinguishActions) {
							mdp.addActionLabelledChoice(src, distr, modelGen.getChoiceAction(i));
						} else {
							mdp.addChoice(src, distr);
						}
					} else if (modelType == ModelType.CTMDP) {
						if (distinguishActions) {
							ctmdp.addActionLabelledChoice(src, distr, modelGen.getChoiceAction(i));
						} else {
							ctmdp.addChoice(src, distr);
						}
					} else if (modelType == ModelType.STPG) {
						if (distinguishActions) {
							stpg.addActionLabelledChoice(src, distr, indexes);
						} else {
							stpg.addChoice(src, distr);
						}
					}
					else if (modelType == ModelType.CSG) {
						if (distinguishActions) {
							csg.addActionLabelledChoice(src, distr, indexes);
						} 
						else {
							csg.addChoice(src, distr, indexes);	
						}
					} else if (modelType == ModelType.SMG) {
						if (distinguishActions) {
							smg.addActionLabelledChoice(src, distr, modelGen.getTransitionAction(i, 0));
						} else {
							smg.addChoice(src, distr);
						}
					}
				}
			}
			// Print some progress info occasionally
			progress.updateIfReady(src + 1);
		}

		// Finish progress display
		progress.update(src + 1);
		progress.end(" states");

		// Reachability complete
		mainLog.print("Reachable states exploration" + (justReach ? "" : " and model construction"));
		mainLog.println(" done in " + ((System.currentTimeMillis() - timer) / 1000.0) + " secs.");
		//mainLog.println(states);
		
		// Find/fix deadlocks (if required)
		if (!justReach && findDeadlocks) {
			if (modelType != ModelType.CSG) {
				modelSimple.findDeadlocks(fixDeadlocks);
			}
			else {
				modelSimple.findDeadlocks(false);
				// Fixes deadlocks for concurrent games...this may have to be changed though
				if(modelSimple.getNumDeadlockStates() > 0) {
					for(Integer s : modelSimple.getDeadlockStates()) {
						csg.fixDeadlock(s);
					}			
				}
			}
		}
				
		int permut[] = null;

		if (sortStates) {
			// Sort states and convert set to list
			mainLog.println("Sorting reachable states list...");
			permut = states.buildSortingPermutation();
			statesList = states.toPermutedArrayList(permut);
			//mainLog.println("Arrays.toString(permut));
		} else {
			statesList = states.toArrayList();
		}
		states.clear();
		states = null;
		//mainLog.println(permut);
		//mainLog.println(statesList);
		
		// Construct new explicit-state model (with correct state ordering)
		if (!justReach) {
			switch (modelType) {
			case DTMC:
				if (buildSparse) {
					model = sortStates ? new DTMCSparse(dtmc, permut) : new DTMCSparse(dtmc);
				} else {
					model = sortStates ? new DTMCSimple(dtmc, permut) : (DTMCSimple) dtmc;
				}
				break;
			case CTMC:
				model = sortStates ? new CTMCSimple(ctmc, permut) : (CTMCSimple) ctmc;
				break;
			case MDP:
				if (buildSparse) {
					model = sortStates ? new MDPSparse(mdp, true, permut) : new MDPSparse(mdp);
				} else {
					model = sortStates ? new MDPSimple(mdp, permut) : mdp;
				}
				break;
			case CTMDP:
				model = sortStates ? new CTMDPSimple(ctmdp, permut) : mdp;
				break;
			case CSG:
				model = sortStates ? new CSG(csg, permut) : csg;
				break;
			case STPG:
				model = sortStates ? new STPGExplicit(stpg, permut) : stpg;
				break;
			case SMG:
				model = sortStates ? new SMG(smg, permut) : smg;
				break;
			case PTA:
			case LTS:
				throw new PrismNotSupportedException("Model construction not supported for " + modelType + "s");
			}
			model.setStatesList(statesList);
			model.setConstantValues(new Values(modelGen.getConstantValues()));
		}

		// Add idle actions
		if (modelType == ModelType.CSG) 
			csg.addIdleIndexes();
			
		// Discard permutation
		permut = null;

		if (!justReach && attachLabels)
			attachLabels(modelGen, model);
				
		return model;
	}

	private void attachLabels(ModelGenerator modelGen, ModelExplicit model) throws PrismException
	{
		// Get state info
		List<State> statesList = model.getStatesList();
		int numStates = statesList.size();
		// Create storage for labels
		int numLabels = modelGen.getNumLabels();
		// No need to continue unless this ModelGenerator uses labels
		if (numLabels == 0) return;
		BitSet bitsets[] = new BitSet[numLabels];
		for (int j = 0; j < numLabels; j++) {
			bitsets[j] = new BitSet();
		}
		// Construct bitsets for labels
		for (int i = 0; i < numStates; i++) {
			State state = statesList.get(i);
			modelGen.exploreState(state);
			for (int j = 0; j < numLabels; j++) {
				if (modelGen.isLabelTrue(j)) {
					bitsets[j].set(i);
				}
			}
		}
		// Attach labels/bitsets
		for (int j = 0; j < numLabels; j++) {
			model.addLabel(modelGen.getLabelName(j), bitsets[j]);
		}
	}

	/**
	 * Test method.
	 */
	public static void main(String[] args)
	{
		try {
			// Simple example: parse a PRISM file from a file, construct the model and export to a .tra file
			PrismLog mainLog = new PrismPrintStreamLog(System.out);
			Prism prism = new Prism(mainLog);
			parser.ast.ModulesFile modulesFile = prism.parseModelFile(new File(args[0]));
			UndefinedConstants undefinedConstants = new UndefinedConstants(modulesFile, null);
			if (args.length > 2)
				undefinedConstants.defineUsingConstSwitch(args[2]);
			modulesFile.setUndefinedConstants(undefinedConstants.getMFConstantValues());
			ConstructModel constructModel = new ConstructModel(prism);
			constructModel.setSortStates(true);
			simulator.ModulesFileModelGenerator modelGen = new simulator.ModulesFileModelGenerator(modulesFile, constructModel);
			Model model = constructModel.constructModel(modelGen, false);
			model.exportToPrismExplicitTra(args[1]);
		} catch (FileNotFoundException e) {
			System.out.println("Error: " + e.getMessage());
			System.exit(1);
		} catch (PrismException e) {
			System.out.println("Error: " + e.getMessage());
			System.exit(1);
		}
	}
}
