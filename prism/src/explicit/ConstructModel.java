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
import java.util.BitSet;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import parser.State;
import parser.Values;
import parser.VarList;
import parser.ast.ModulesFile;
import parser.ast.SystemBrackets;
import parser.ast.SystemDefn;
import parser.ast.SystemFullParallel;
import parser.ast.SystemReference;
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
	/** Build a sparse representation, if possible?
	 *  (e.g. MDPSparse rather than MDPSimple data structure) */
	protected boolean buildSparse = true;
	/** Should actions be attached to distributions (and used to distinguish them)? */
	protected boolean distinguishActions = true;
	/** Should labels be processed and attached to the model? */
	protected boolean attachLabels = true;

	/** Automatically check compatibility? */
	protected boolean checkCompatibility = false;

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
	 * Automatically check compatibility?
	 */
	public void setCheckCompatibility(boolean checkCompatibility)
	{
		this.checkCompatibility = checkCompatibility;
	}

	/**
	 * Build the set of reachable states for a model and return it.
	 * @param modelGen The ModelGenerator interface providing the model 
	 */
	public List<State> computeReachableStates(ModelGenerator modelGen) throws PrismException
	{
		constructModel(modelGen, true, null);
		return getStatesList();
	}

	/**
	 * Construct an explicit-state model and return it.
	 * @param modelGen The ModelGenerator interface providing the model 
	 */
	public Model constructModel(ModelGenerator modelGen) throws PrismException
	{
		return constructModel(modelGen, false, null);
	}

	/**
	 * Construct an explicit-state model and return it.
	 * If {@code justReach} is true, no model is built and null is returned;
	 * the set of reachable states can be obtained with {@link #getStatesList()}.
	 * @param modelGen The ModelGenerator interface providing the model 
	 * @param justReach If true, just build the reachable state set, not the model
	 *
	 */
	public Model constructModel(ModelGenerator modelGen, boolean justReach, boolean[] cancel_computation) throws PrismException
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
		modelType = modelGen.getModelType();

		// For PRISM SMGs with a system definition, we build compositionally
		if (modelGen instanceof ModulesFileModelGenerator && modelType == ModelType.SMG
				&& ((ModulesFileModelGenerator) modelGen).getModulesFile().getSystemDefn() != null) {
			// default is to not check compatibility
			return constructSMGModelCompositionally(((ModulesFileModelGenerator) modelGen).getModulesFile(), justReach, cancel_computation);
		}

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
				if (modelGen instanceof ModulesFileModelGenerator && modelType == ModelType.SMG
						&& ((ModulesFileModelGenerator) modelGen).getModulesFile().getModuleNames() != null
						&& ((ModulesFileModelGenerator) modelGen).getModulesFile().getModuleNames().length > 0) {
					modelSimple = smg = new SMG(((ModulesFileModelGenerator) modelGen).getModulesFile().getModuleName(0));
				} else {
					modelSimple = smg = new SMG();
				}
				// Add player info
				HashMap<Integer, String> playerNames = new HashMap<Integer, String>();
				for (i = 0; i < modelGen.getNumPlayers(); i++) {
					playerNames.put(i + 1, modelGen.getPlayer(i).getName());
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
			if (modelType.multiplePlayers()) {
				player = -1;
				for (i = 0; i < nc; i++) {
					int iPlayer = modelGen.getPlayerNumberForChoice(i);
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
						case MDP:
						case CTMDP:
						case STPG:
						case SMG:
							distr.add(dest, modelGen.getTransitionProbability(i, j));
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
							mdp.addActionLabelledChoice(src, distr, modelGen.getTransitionAction(i, 0));
						} else {
							mdp.addChoice(src, distr);
						}
					} else if (modelType == ModelType.CTMDP) {
						if (distinguishActions) {
							ctmdp.addActionLabelledChoice(src, distr, modelGen.getTransitionAction(i, 0));
						} else {
							ctmdp.addChoice(src, distr);
						}
					} else if (modelType == ModelType.STPG) {
						if (distinguishActions) {
							stpg.addActionLabelledChoice(src, distr, modelGen.getTransitionAction(i, 0));
						} else {
							stpg.addChoice(src, distr);
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
			model.setConstantValues(new Values(modelGen.getConstantValues()));
		}

		// Discard permutation
		permut = null;

		if (attachLabels)
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
	 * Construct an explicit-state model for an SMG with a system definition (done compositionally).
	 * If {@code justReach} is true, no model is built and null is returned;
	 * the set of reachable states can be obtained with {@link #getStatesList()}.
	 * @param modulesFile The PRISM model
	 */
	public Model constructSMGModelCompositionally(ModulesFile modulesFile, boolean justReach, boolean[] cancel_computation) throws PrismException
	{
		return constructSMGModelCompositionally(modulesFile, justReach, null, null, true, cancel_computation);
	}

	/**
	* Compositionally constructs the explicit-state model for an SMG with a system definition.
	* Returns the component games in {@code subsystems} and the corresponding modules files in
	* {@code subsystemModulesFiles}. If {@code buildFullModel} or {@code checkCompatibility} is
	* {@code true}, the full model is constructed. Otherwise, the full model is not constructed
	* and and only the components are returned.
	**/
	public Model constructSMGModelCompositionally(ModulesFile modulesFile, boolean justReach, List<SMG> subsystems, List<ModulesFile> subsystemModulesFiles,
			boolean buildFullModel, boolean[] cancel_computation) throws PrismException
	{
		// if compatibility check or full model build is requested, need to build full model in any case
		if ((checkCompatibility || buildFullModel) && subsystems == null)
			subsystems = new ArrayList<SMG>();
		if ((checkCompatibility || buildFullModel) && subsystemModulesFiles == null)
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
			if (subsystemModulesFiles != null)
				subsystemModulesFiles.add(modulesFile_i);
			// construct the models
			ModulesFileModelGenerator modelGen = new ModulesFileModelGenerator(modulesFile_i, this);
			modelGen.setCompositional(true);
			m = (SMG) constructModel(modelGen, justReach, cancel_computation);
			// convert to normal form if necessary
			if (m != null)
				m.toNormalForm();
			// add model to subsystem list
			if (subsystems != null)
				subsystems.add(m);
		}
		if (!justReach && (checkCompatibility || buildFullModel) && sys instanceof SystemFullParallel) {
			// form composition and check compatibility (if requested)
			double t0 = (double) System.nanoTime();
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
			Prism prism = new Prism(mainLog);
			parser.ast.ModulesFile modulesFile = prism.parseModelFile(new File(args[0]));
			UndefinedConstants undefinedConstants = new UndefinedConstants(modulesFile, null);
			if (args.length > 2)
				undefinedConstants.defineUsingConstSwitch(args[2]);
			modulesFile.setUndefinedConstants(undefinedConstants.getMFConstantValues());
			ConstructModel constructModel = new ConstructModel(prism);
			simulator.ModulesFileModelGenerator modelGen = new simulator.ModulesFileModelGenerator(modulesFile, constructModel);
			Model model = constructModel.constructModel(modelGen, false, cancel_computation);
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
