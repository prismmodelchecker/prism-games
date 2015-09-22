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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Vector;

import parser.State;
import parser.Values;
import parser.VarList;
import parser.ast.Expression;
import parser.ast.ModulesFile;
import parser.ast.SystemBrackets;
import parser.ast.SystemDefn;
import parser.ast.SystemFullParallel;
import parser.ast.SystemModule;
import parser.ast.SystemParallel;
import parser.ast.SystemReference;
import prism.ModelType;
import prism.Prism;
import prism.PrismComponent;
import prism.PrismException;
import prism.PrismLangException;
import prism.PrismLog;
import prism.PrismPrintStreamLog;
import prism.PrismNotSupportedException;
import prism.ProgressDisplay;
import prism.UndefinedConstants;
import simulator.SimulatorEngine;

public class ConstructModel extends PrismComponent
{
	// The simulator engine
	protected SimulatorEngine engine;

	// Options:
	// Find deadlocks during model construction?
	protected boolean findDeadlocks = true;
	// Automatically fix deadlocks?
	protected boolean fixDeadlocks = true;
        // Automatically check compatibility?
        protected boolean checkCompatibility = false;

	// Details of built model
	protected List<State> statesList;

	public ConstructModel(PrismComponent parent, SimulatorEngine engine) throws PrismException
	{
		super(parent);
		this.engine = engine;
	}

	public List<State> getStatesList()
	{
		return statesList;
	}

	public void setCheckCompatibility(boolean b)
	{
	        checkCompatibility = b;
	}

	public void setFixDeadlocks(boolean b)
	{
		fixDeadlocks = b;
	}

	/**
	 * Build the set of reachable states for a PRISM model language description and return.
	 * @param modulesFile The PRISM model
	 */
    public List<State> computeReachableStates(ModulesFile modulesFile, boolean[] cancel_computation) throws PrismException
	{
	    constructModel(modulesFile, true, false, cancel_computation);
		return statesList;
	}

	/**
	 * Construct an explicit-state model from a PRISM model language description and return.
	 * @param modulesFile The PRISM model
	 */
    public Model constructModel(ModulesFile modulesFile, boolean[] cancel_computation) throws PrismException
	{
	    return constructModel(modulesFile, false, true, cancel_computation);
	}

	/**
	 * Construct an explicit-state model from a PRISM model language description and return.
	 * If {@code justReach} is true, no model is built and null is returned;
	 * the set of reachable states can be obtained with {@link #getStatesList()}.
	 * @param modulesFile The PRISM model
	 * @param justReach If true, just build the reachable state set, not the model
	 * @param buildSparse Build a sparse version of the model (if possible)?
	 */
    public Model constructModel(ModulesFile modulesFile, boolean justReach, boolean buildSparse, boolean[] cancel_computation) throws PrismException
	{
	    return constructModel(modulesFile, justReach, buildSparse, true, false, cancel_computation);
	}

	/**
	 * Construct an explicit-state model from a PRISM model language description and return.
	 * If {@code justReach} is true, no model is built and null is returned;
	 * the set of reachable states can be obtained with {@link #getStatesList()}.
	 * @param modulesFile The PRISM model
	 * @param justReach If true, just build the reachable state set, not the model
	 * @param buildSparse Build a sparse version of the model (if possible)?
	 * @param distinguishActions True if actions should be attached to distributions (and used to distinguish them)
	 * @param compositional True if called for compositional model building (affects action assignments)
	 *
	 */
    public Model constructModel(ModulesFile modulesFile, boolean justReach, boolean buildSparse, boolean distinguishActions, boolean compositional, boolean[] cancel_computation) throws PrismException
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
		SMG smg = null;
		ModelExplicit model = null;
		Distribution distr = null;
		// Misc
		int i, j, nc, nt, src, dest, player;
		long timer;

		// Get model info
		modelType = modulesFile.getModelType();
		
		// For SMGs with a system definition, we build compositionally
		if (modelType == ModelType.SMG && modulesFile.getSystemDefn() != null) {
		    // default is to not check compatibility
		    return constructSMGModelCompositionally(modulesFile, justReach, buildSparse, distinguishActions, cancel_computation);
		}
		
		// Display a warning if there are unbounded vars
		VarList varList = modulesFile.createVarList();
		if (varList.containsUnboundedVariables())
			mainLog.printWarning("Model contains one or more unbounded variables: model construction may not terminate");

		// Starting reachability...
		mainLog.print("\nComputing reachable states...");
		mainLog.flush();
		ProgressDisplay progress = new ProgressDisplay(mainLog);
		progress.start();
		timer = System.currentTimeMillis();

		// Initialise simulator for this model
		engine.createNewOnTheFlyPath(modulesFile);

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
				break;
			case SMG:
			        if(modelType == ModelType.SMG && modulesFile.getSystemDefn() !=  null && modulesFile.getSystemDefnName() != null) {
				        // get top-level system definition name
				        modelSimple = smg = new SMG(modulesFile.getSystemDefnName());
				} else if(modulesFile.getModuleNames() != null && modulesFile.getModuleNames().length > 0) {
				        modelSimple = smg = new SMG(modulesFile.getModuleName(0));
				} else {
				        modelSimple = smg = new SMG();
				}
				// Add player info
				HashMap<Integer, String> playerNames = new HashMap<Integer, String>();
				for (i = 0; i < modulesFile.getNumPlayers(); i++) {
					playerNames.put(i + 1, modulesFile.getPlayer(i).getName());
				}
				smg.setPlayerInfo(playerNames);
				break;
			case PTA:
				throw new PrismNotSupportedException("Model construction not supported for " + modelType + "s");
			}
		}

		// Initialise states storage
		states = new IndexedSet<State>(true);
		explore = new LinkedList<State>();
		// Add initial state(s) to 'explore'
		// Easy (normal) case: just one initial state
		if (modulesFile.getInitialStates() == null) {
			state = modulesFile.getDefaultInitialState();
			explore.add(state);
		}
		// Otherwise, there may be multiple initial states
		// For now, we handle this is in a very inefficient way
		else {
			Expression init = modulesFile.getInitialStates();
			List<State> allPossStates = varList.getAllStates();
			for (State possState : allPossStates) {
				if (init.evaluateBoolean(modulesFile.getConstantValues(), possState)) {
					explore.add(possState);
				}
			}
		}
		// Copy initial state(s) to 'states' and to the model
		for (State initState : explore) {
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

			// Use simulator to explore all choices/transitions from this state
			engine.initialisePath(state);
			nc = engine.getNumChoices();
			// For games, first determine which player owns the state
			if (modelType.multiplePlayers()) {
				player = -1;
				for (i = 0; i < nc; i++) {
				        int iPlayer = determinePlayerForChoice(modulesFile, modelType, i, compositional);
					if (player != -1 && iPlayer != player) {
						throw new PrismException("Choices for both player " + player + " and " + iPlayer + " in state " + state);
					}
					player = iPlayer;
				}
				// Assign deadlock states to player 1
				if (nc == 0) {
					player = 1;
				}
				if (modelType == ModelType.STPG) {
					stpg.setPlayer(src, player);
				} else if (modelType == ModelType.SMG) {
					smg.setPlayer(src, player);
				}
			}
			// Look at each outgoing choice in turn
			for (i = 0; i < nc; i++) {
				// For nondet models, collect transitions in a Distribution
				if (!justReach && modelType.nondeterministic()) {
					distr = new Distribution();
				}
				// Look at each transition in the choice
				nt = engine.getNumTransitions(i);
				for (j = 0; j < nt; j++) {
					stateNew = engine.computeTransitionTarget(i, j);

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
							dtmc.addToProbability(src, dest, engine.getTransitionProbability(i, j));
							break;
						case CTMC:
							ctmc.addToProbability(src, dest, engine.getTransitionProbability(i, j));
							break;
						case MDP:
						case CTMDP:
						case STPG:
						case SMG:
							distr.add(dest, engine.getTransitionProbability(i, j));
							break;
						case PTA:
							throw new PrismNotSupportedException("Model construction not supported for " + modelType + "s");
						}
					}
				}
				// For nondet models, add collated transition to model 
				if (!justReach) {
					if (modelType == ModelType.MDP) {
						if (distinguishActions) {
							mdp.addActionLabelledChoice(src, distr, engine.getTransitionAction(i, 0));
						} else {
							mdp.addChoice(src, distr);
						}
					} else if (modelType == ModelType.CTMDP) {
						if (distinguishActions) {
							ctmdp.addActionLabelledChoice(src, distr, engine.getTransitionAction(i, 0));
						} else {
							ctmdp.addChoice(src, distr);
						}
					} else if (modelType == ModelType.STPG) {
						if (distinguishActions) {
							stpg.addActionLabelledChoice(src, distr, engine.getTransitionAction(i, 0));
						} else {
							stpg.addChoice(src, distr);
						}
					} else if (modelType == ModelType.SMG) {
						if (distinguishActions) {
							smg.addActionLabelledChoice(src, distr, engine.getTransitionAction(i, 0));
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
			modelSimple.findDeadlocks(fixDeadlocks);
		}

		boolean sort = true;
		int permut[] = null;

		if (sort) {
			// Sort states and convert set to list
			mainLog.println("Sorting reachable states list...");
			permut = states.buildSortingPermutation();
			statesList = states.toPermutedArrayList(permut);
			//mainLog.println(permut);
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
				model = sort ? new DTMCSimple(dtmc, permut) : (DTMCSimple) dtmc;
				break;
			case CTMC:
				model = sort ? new CTMCSimple(ctmc, permut) : (CTMCSimple) ctmc;
				break;
			case MDP:
				if (buildSparse) {
					model = sort ? new MDPSparse(mdp, true, permut) : new MDPSparse(mdp);
				} else {
					model = sort ? new MDPSimple(mdp, permut) : mdp;
				}
				break;
			case CTMDP:
				model = sort ? new CTMDPSimple(ctmdp, permut) : mdp;
				break;
			case STPG:
				model = sort ? new STPGExplicit(stpg, permut) : stpg;
				break;
			case SMG:
				model = sort ? new SMG(smg, permut) : smg;
				break;
			case PTA:
				throw new PrismNotSupportedException("Model construction not supported for " + modelType + "s");
			}
			model.setStatesList(statesList);
			model.setConstantValues(new Values(modulesFile.getConstantValues()));
		}

		// Discard permutation
		permut = null;

		return model;
	}

	/**
	 * For game models, determine the player who owns the {@code i}th choice in
	 * the state currently being explored by the simulator. Returns the index
	 * (starting from 1) of the player.
	 * @param compositional True if called for compositional model building (affects action assignments)
	 *
	 */
    private int determinePlayerForChoice(ModulesFile modulesFile, ModelType modelType, int i, boolean compositional) throws PrismException
	{
		int modAct, player;

		try {
		    modAct = engine.getTransitionModuleOrActionIndex(i, 0);
		} catch (PrismLangException e) {
		    throw new PrismException(String.format("%s. Did you forget to assign an action to a player?"));
		}
		// Synchronous action
		if (modAct > 0) {
		    // first try to get action from player specifications
		    String action_name = modulesFile.getSynch(modAct - 1);
		    if(compositional) {
			// get inputs and outputs of subsystem
			Vector<String> inputs = new Vector<String>();
			Vector<String> outputs = new Vector<String>();
			for(int n = 0; n < modulesFile.getNumModules(); n++) {
			    inputs.addAll(modulesFile.getModule(n).getAllInputActions());
			    outputs.addAll(modulesFile.getModule(n).getAllOutputActions());
			}

			if(action_name == null || "".equals(action_name)) {
			    player = 2; // tau controlled by player 2
			} else if(inputs.contains(action_name)) {
			    player = 2; // inputs controlled by player 2
			} else if(outputs.contains(action_name)) {
			    player = 1; // outputs controlled by player 1
			} else {
			    throw new PrismException("Action \"" + engine.getTransitionModuleOrAction(i, 0) + "\" is not assigned to any player");
			}
		    } else {
			player = modulesFile.getPlayerForAction(action_name);
			if (player == -1) {
			    throw new PrismException("Action \"" + modulesFile.getSynch(modAct - 1) + "\" is not assigned to any player");
			}
			// 0-indexed to 1-indexed
			player++;
		    }
		}
		// Asynchronous action
		else {
			if (compositional) {
			    player = 2; // tau controlled by player 2
			} else {
				player = modulesFile.getPlayerForModule(engine.getTransitionModuleOrAction(i, 0));
				if (player == -1) {
					throw new PrismException("Module \"" + engine.getTransitionModuleOrAction(i, 0) + "\" is not assigned to any player");
				}
				// 0-indexed to 1-indexed
				player++;
			}
		}
		return player;
	}

	/**
	 * Construct an explicit-state model for an SMG with a system definition (done compositionally).
	 * If {@code justReach} is true, no model is built and null is returned;
	 * the set of reachable states can be obtained with {@link #getStatesList()}.
	 * @param modulesFile The PRISM model
	 * @param justReach If true, just build the reachable state set, not the model
	 * @param buildSparse Build a sparse version of the model (if possible)?
	 * @param distinguishActions True if actions should be attached to distributions (and used to distinguish them)
	 */
    public Model constructSMGModelCompositionally(ModulesFile modulesFile, boolean justReach, boolean buildSparse, boolean distinguishActions, boolean[] cancel_computation) throws PrismException
        {
	    return constructSMGModelCompositionally(modulesFile, justReach, buildSparse, distinguishActions, null, null, true, cancel_computation);
	}
        /**
	 * Compositionally constructs the explicit-state model for an SMG with a system definition.
	 * Returns the component games in {@code subsystems} and the corresponding modules files in
	 * {@code subsystemModulesFiles}. If {@code buildFullModel} or {@code checkCompatibility} is
	 * {@code true}, the full model is constructed. Otherwise, the full model is not constructed
	 * and and only the components are returned.
	 **/
        public Model constructSMGModelCompositionally(ModulesFile modulesFile, boolean justReach, boolean buildSparse, boolean distinguishActions,
						      List<SMG> subsystems, List<ModulesFile> subsystemModulesFiles, boolean buildFullModel, boolean[] cancel_computation) throws PrismException
	{
	        // if compatibility check or full model build is requested, need to build full model in any case
	        if((checkCompatibility || buildFullModel) && subsystems==null)
		        subsystems = new ArrayList<SMG>();
	        if((checkCompatibility || buildFullModel) && subsystemModulesFiles==null)
		        subsystemModulesFiles = new ArrayList<ModulesFile>();

		// Extract subsystems from system definition
		SystemDefn sys = modulesFile.getSystemDefn();
		if (sys == null) {
			throw new PrismException("Can only construct the model compositionally if there is a system definition");
		}
		while (sys instanceof SystemBrackets) {
			sys = ((SystemBrackets) sys).getOperand();
		}

		ArrayList<SystemReference> sysRefs = new ArrayList<SystemReference>();
		ModulesFile.extractSubsystemRefs(sys, sysRefs);
		
		// Extract modules in each subsystem ...
		int numComps = sysRefs.size();
		ArrayList<List<String>> moduleNameLists = new ArrayList<List<String>>();
		for (int i = 0; i < numComps; i++) {
		        SystemDefn subsys = modulesFile.getSystemDefnByName(sysRefs.get(i).getName());
			if (subsys == null) {
			        throw new PrismException("Unknown system reference" + sysRefs.get(i));
			}
			ArrayList<String> moduleNames = new ArrayList<String>();
			ModulesFile.extractSubsystemModuleNames(subsys, moduleNames);
			moduleNameLists.add(moduleNames);
		}
		
		// ... to build the subsystems individually.
		setFixDeadlocks(false); // deadlocks not fixed
		SMG m = null;
		for (int i = 0; i < numComps; i++) {
		        // extract and store modules files for the composition later
		        ModulesFile modulesFile_i = (ModulesFile) modulesFile.deepCopy(moduleNameLists.get(i));
			if(subsystemModulesFiles != null) subsystemModulesFiles.add(modulesFile_i);
			// construct the models
			m = (SMG) constructModel(modulesFile_i, justReach, buildSparse, distinguishActions, true, cancel_computation);
			// convert to normal form if necessary
			if(m != null) m.toNormalForm();
			// add model to subsystem list
			if(subsystems != null) subsystems.add(m);
		}
		if(!justReach && (checkCompatibility || buildFullModel) &&
		       sys instanceof SystemFullParallel) {
		        // form composition and check compatibility (if requested)
		        double t0 = (double)System.nanoTime();
			SMG smg_compositional = new SMG(modulesFile, subsystemModulesFiles, subsystems, mainLog, checkCompatibility, cancel_computation);
			
			mainLog.print(String.format("product construction took %f s\n", ((double) (System.nanoTime() - t0)) / 1e9));
			return smg_compositional;
		} else {
		        return m;
		}
	}


	/**
	 * Test method.
	 */
	public static void main(String[] args)
	{
		boolean[] cancel_computation = new boolean[1]; // false by default

		try {
			// Simple example: parse a PRISM file from a file, construct the model and export to a .tra file
			PrismLog mainLog = new PrismPrintStreamLog(System.out);
			Prism prism = new Prism(mainLog, mainLog);
			ModulesFile modulesFile = prism.parseModelFile(new File(args[0]));
			UndefinedConstants undefinedConstants = new UndefinedConstants(modulesFile, null);
			if (args.length > 2)
				undefinedConstants.defineUsingConstSwitch(args[2]);
			modulesFile.setUndefinedConstants(undefinedConstants.getMFConstantValues());
			ConstructModel constructModel = new ConstructModel(prism, prism.getSimulator());
			Model model = constructModel.constructModel(modulesFile, cancel_computation);
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
