//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford)
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

import java.io.*;
import java.util.*;

import parser.*;
import parser.ast.*;
import prism.*;
import simulator.*;

public class ConstructModel
{
	// The simulator engine and a log for output
	private SimulatorEngine engine;
	private PrismLog mainLog;
	
	// PRISM settings object (optional)
	private PrismSettings settings;

	// Basic info needed about model
	//	private ModelType modelType;

	// Details of built model
	private List<State> statesList;

	// Partition of modules into players
	private Set<String> player1mods;
	private Set<String> player2mods;
	// Partition of synchronised actions
	private Set<String> player1syncs;
	private Set<String> player2syncs;

	public ConstructModel(SimulatorEngine engine, PrismLog mainLog)
	{
		this.engine = engine;
		this.mainLog = mainLog;
		settings = null;
	}

	public void setSettings(PrismSettings settings)
	{
		this.settings = settings;
	}
	
	/** 
	 *	Methods extracts game parameters: division of modules and synchronised actions from a given string.
	 *
	 *	@param params string containing partitions of synchronised actions and 
	 *	module names for two players. Format for player 1: p1s=[action1,action2,...,actionN]
	 *	p1m=[module1,module2,...,moduleN]. For player 2 it is the same but names are p2s p2m. 
	 */
	private void loadGameParams(String params)
	{
		try {
			Properties props = new Properties();

			for (String p : params.split(" "))
				props.load(new StringReader(p));

			if (props.containsKey("p1m"))
				player1mods.addAll(Arrays.asList(((String) props.get("p1m")).substring(1).split("\\]|,| ")));
			if (props.containsKey("p2m"))
				player2mods.addAll(Arrays.asList(((String) props.get("p2m")).substring(1).split("\\]|,| ")));
			if (props.containsKey("p1s"))
				player1syncs.addAll(Arrays.asList(((String) props.get("p1s")).substring(1).split("\\]|,| ")));
			if (props.containsKey("p2s"))
				player2syncs.addAll(Arrays.asList(((String) props.get("p2s")).substring(1).split("\\]|,| ")));

		} catch (Exception e) {
			// Loading of parameters failed.
			e.printStackTrace();
		}
	}

	public List<State> getStatesList()
	{
		return statesList;
	}

	/**
	 * Build the set of reachable states for a PRISM model language description and return.
	 * @param modulesFile The PRISM model
	 * @param initialState The initial state (for reachability)
	 */
	public List<State> computeReachableStates(ModulesFile modulesFile, Values initialState) throws PrismException
	{
		constructModel(modulesFile, initialState, true, false);
		return statesList;
	}

	/**
	 * Construct an explicit-state model from a PRISM model language description and return.
	 * @param modulesFile The PRISM model
	 * @param initialState The initial state (for reachability)
	 */
	public Model constructModel(ModulesFile modulesFile, Values initialState) throws PrismException
	{
		return constructModel(modulesFile, initialState, false, false);
	}

	/**
	 * Construct an explicit-state model from a PRISM model language description and return.
	 * If {@code justReach} is true, no model is built and null is returned;
	 * the set of reachable states can be obtained with {@link #getStatesList()}.
	 * @param modulesFile The PRISM model
	 * @param initialState The initial state (for reachability)
	 * @param justReach If true, just build the reachable state set, not the model
	 * @param buildSparse Build a sparse version of the model (if possible)?
	 */
	public Model constructModel(ModulesFile modulesFile, Values initialState, boolean justReach, boolean buildSparse)
			throws PrismException
	{
		// Model info
		ModelType modelType;
		// State storage
		IndexedSet<State> states;
		LinkedList<State> explore;
		State state, stateNew;
		List<Integer> stateLabels;
		List<State> statesref = null;
		// Explicit model storage
		ModelSimple modelSimple = null;
		DTMCSimple dtmc = null;
		CTMCSimple ctmc = null;
		MDPSimple mdp = null;
		STPGExplicit stpg = null;
		Model model = null;
		Distribution distr = null;
		// Misc
		int i, j, k, nc, nt, src, dest;
		long timer, timerProgress;
		boolean fixdl = false;
		String actionLabel = null;
		int player;
		int id;

		// Process game info
		this.player1mods = new HashSet<String>();
		this.player2mods = new HashSet<String>();
		this.player1syncs = new HashSet<String>();
		this.player2syncs = new HashSet<String>();
		//loadGameParams("p1m=[scheduler,task_generator,sensor1,sensor3,sensor5,sensor7] p2m=[sensor2,sensor4,sensor6] p1s=[initialise,scheduling,str1,str2,str3,str4,str5,str6,str7,fin1,fin2,fin3,fin4,fin5,fin6,fin7] p2s=[]");
		loadGameParams(settings.getString(PrismSettings.PRISM_GAME_OPTIONS));
		
		// For now, don't use sparse (so can use actions) (TODO: fix)
		buildSparse = false;
		
		// Don't support multiple initial states
		if (modulesFile.getInitialStates() != null) {
			throw new PrismException("Cannot do explicit-state reachability if there are multiple initial states");
		}

		// Starting reachability...
		mainLog.print("\nComputing reachable states...\n");
		mainLog.flush();
		timer = timerProgress = System.currentTimeMillis();

		// Initialise simulator for this model
		modelType = modulesFile.getModelType();

		// TODO HACK!
		//modelType = ModelType.STPG;

		engine.createNewOnTheFlyPath(modulesFile);

		// Create model storage
		if (!justReach) {
			// Create a (simple, mutable) model of the appropriate type
			switch (modelType) {
			case DTMC:
				modelSimple = dtmc = new DTMCSimple();
				break;
			case CTMC:
				modelSimple = ctmc = new CTMCSimple();
				break;
			case MDP:
				modelSimple = mdp = new MDPSimple();
				break;
			case STPG:
				modelSimple = stpg = new STPGExplicit();
				break;

			}
		}

		// Initialise states storage
		states = new IndexedSet<State>(true);
		explore = new LinkedList<State>();
		stateLabels = new ArrayList<Integer>();
		// Add initial state to lists/model
		state = new State(modulesFile.getInitialValues());
		states.add(state);
		explore.add(state);
		if (!justReach) {
			if (modelType == ModelType.STPG) {
				engine.initialisePath(state);
				actionLabel = engine.getTransitionModuleOrAction(0, 0);
				if (actionLabel != null) {
					id = stpg.addState(getPlayer(actionLabel));
					stateLabels.add(1);
				} else
					throw new PrismException(
							"Every state must have an action. Cannot determine the player which owns the state.");
			} else
				modelSimple.addState();

			modelSimple.addInitialState(0);
		}

		// Explore...
		src = -1;
		while (!explore.isEmpty()) {
			// Pick next state to explore
			// (they are stored in order found so know index is src+1)
			state = explore.removeFirst();
			src++;
			// System.out.println(src);

			// Use simulator to explore all choices/transitions from this state
			engine.initialisePath(state);
			// Look at each outgoing choice in turn
			nc = engine.getNumChoices();
			for (i = 0; i < nc; i++) {
				if (!justReach && (modelType == ModelType.MDP || modelType == ModelType.STPG)) {
					distr = new Distribution();
				}
				// Look at each transition in the choice
				nt = engine.getNumTransitions(i);
				for (j = 0; j < nt; j++) {
					stateNew = engine.computeTransitionTarget(i, j);

					// labeling the source state with player
					if (modelType == ModelType.STPG && j == 0) {
						actionLabel = engine.getTransitionModuleOrAction(i, j);
						if (actionLabel != null)
							stateLabels.set(src, getPlayer(actionLabel));
						else
							throw new PrismException(
									"Every state must have an action. Cannot determine the player which owns the state.");
					}

					// Is this a new state?
					if (states.add(stateNew)) {
						// If so, add to the explore list
						explore.add(stateNew);
						// And to model
						if (!justReach) {
							if (modelType == ModelType.STPG)
								stateLabels.add(stpg.addState());
							else
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
							distr.add(dest, engine.getTransitionProbability(i, j));
							break;
						case STPG:
							distr.add(dest, engine.getTransitionProbability(i, j));
							break;
						}
					}
				}
				if (!justReach) {
					if (modelType == ModelType.MDP) {
						k = mdp.addChoice(src, distr);
						mdp.setAction(src, k, engine.getTransitionModuleOrAction(i, 0));
					} else if (modelType == ModelType.STPG) {
						k = stpg.addChoice(src, distr);
						stpg.setAction(src, k, engine.getTransitionModuleOrAction(i, 0));
					}
				}
			}

			// Print some progress info occasionally
			if (System.currentTimeMillis() - timerProgress > 3000) {
				mainLog.print(" " + (src + 1));
				mainLog.flush();
				timerProgress = System.currentTimeMillis();
			}
		}

		// labeling states
		if (modelType == ModelType.STPG)
			for (int st = 0; st < stateLabels.size(); st++) {

				stpg.setPlayer(st, stateLabels.get(st));
			}
		// Finish progress display
		mainLog.println(" " + (src + 1));

		// Reachability complete
		mainLog.print("Reachable states exploration" + (justReach ? "" : " and model construction"));
		mainLog.println(" done in " + ((System.currentTimeMillis() - timer) / 1000.0) + " secs.");
		// mainLog.println(states);

		// Fix deadlocks (if required)
		if (!justReach && fixdl) {
			BitSet deadlocks = modelSimple.findDeadlocks(true);
			if (deadlocks.cardinality() > 0) {
				mainLog.println("Added self-loops in " + deadlocks.cardinality() + " states...");
			}
		}

		// Sort states and convert set to list
		mainLog.println("Sorting reachable states list...");
		int permut[] = states.buildSortingPermutation();
		statesList = states.toPermutedArrayList(permut);
		states.clear();
		states = null;
		// mainLog.println(permut);
		// mainLog.println(statesList);

		// Construct new explicit-state model (with correct state ordering)
		if (!justReach) {
			switch (modelType) {
			case DTMC:
				model = new DTMCSimple(dtmc, permut);
				((ModelSimple) model).statesList = statesList;
				((ModelSimple) model).constantValues = new Values(modulesFile.getConstantValues());
				break;
			case CTMC:
				model = new CTMCSimple(ctmc, permut);
				((ModelSimple) model).statesList = statesList;
				((ModelSimple) model).constantValues = new Values(modulesFile.getConstantValues());
				break;
			case MDP:
				if (buildSparse) {
					model = new MDPSparse(mdp, true, permut);
					((ModelSparse) model).statesList = statesList;
					((ModelSparse) model).constantValues = new Values(modulesFile.getConstantValues());
				} else {
					model = new MDPSimple(mdp, permut);
					((ModelSimple) model).statesList = statesList;
					((ModelSimple) model).constantValues = new Values(modulesFile.getConstantValues());
				}
				break;
			case STPG:
				model = new STPGExplicit(stpg, permut);
				((ModelSimple) model).statesList = statesList;
				((ModelSimple) model).constantValues = new Values(modulesFile.getConstantValues());
				break;
			}
			mainLog.println("Model: " + model);
		}

		// Discard permutation
		permut = null;

		return model;
	}

	/**
	 * Extracts player information from the action label
	 * 
	 * @param actionLabel action label
	 * @return STPGExplicit.{PLAYER_1,PLAYER_2}
	 * @throws PrismException when the action label is not assigned to any
	 *           player
	 */
	private int getPlayer(String actionLabel) throws PrismException
	{
		int player;
		if (actionLabel.startsWith("[")) {
			actionLabel = actionLabel.substring(1, actionLabel.length() - 1);
			if (player1syncs.contains(actionLabel))
				player = STPGExplicit.PLAYER_1;
			else if (player2syncs.contains(actionLabel))
				player = STPGExplicit.PLAYER_2;
			else
				throw new PrismException("Synchronised action label '" + actionLabel
						+ "' is not assigned to any player.");
		} else {
			if (player1mods.contains(actionLabel))
				player = STPGExplicit.PLAYER_1;
			else if (player2mods.contains(actionLabel))
				player = STPGExplicit.PLAYER_2;
			else
				throw new PrismException("Module '" + actionLabel + "' is not assigned to any player.");
		}
		return player;
	}

	/**
	 * Test method.
	 */
	public static void main(String[] args)
	{

		// TODO HACK!
		testSTPG(args[0]);
		if (1 == 1)
			return;

		try {
			// Simple example: parse a PRISM file from a file, construct the model and export to a .tra file
			PrismLog mainLog = new PrismPrintStreamLog(System.out);
			Prism prism = new Prism(mainLog, mainLog);
			ModulesFile modulesFile = prism.parseModelFile(new File(args[0]));
			UndefinedConstants undefinedConstants = new UndefinedConstants(modulesFile, null);
			if (args.length > 2)
				undefinedConstants.defineUsingConstSwitch(args[1]);
			modulesFile.setUndefinedConstants(undefinedConstants.getMFConstantValues());
			ConstructModel constructModel = new ConstructModel(prism.getSimulator(), mainLog);
			Model model = constructModel.constructModel(modulesFile, modulesFile.getInitialValues());
			MDP mdp = (MDP) model;
			mainLog.println(mdp);
			mainLog.println(constructModel.getStatesList());
			MDPModelChecker mc = new MDPModelChecker();
			mc.setLog(prism.getMainLog());
			BitSet target = new BitSet();
			target.set(2);
			ModelCheckerResult res = mc.computeReachProbs(mdp, target, true);
			mainLog.println(res.soln);

		} catch (FileNotFoundException e) {
			System.out.println("Error: " + e.getMessage());
			System.exit(1);
		} catch (PrismException e) {
			System.out.println("Error: " + e.getMessage());
			System.exit(1);
		}
	}

	/** Method to test STPG construction */
	public static void testSTPG(String filename)
	{
		try {
			// Simple example: parse a PRISM file from a file, construct the model
			// and export to a .tra file
			PrismLog mainLog = new PrismPrintStreamLog(System.out);
			Prism prism = new Prism(mainLog, mainLog);
			ModulesFile modulesFile = prism.parseModelFile(new File(filename));
			UndefinedConstants undefinedConstants = new UndefinedConstants(modulesFile, null);

			modulesFile.setUndefinedConstants(undefinedConstants.getMFConstantValues());

			ConstructModel constructModel = new ConstructModel(prism.getSimulator(), mainLog);

			Model model = constructModel.constructModel(modulesFile, modulesFile.getInitialValues());

			STPGExplicit stpg = (STPGExplicit) model;
			mainLog.println(stpg);
			mainLog.println(constructModel.getStatesList());
			STPGModelChecker mc = new STPGModelChecker();

			mc.setLog(prism.getMainLog());
			BitSet target = new BitSet();
			target.set(8);
			stpg.exportToDotFile("stpg.dot", target);
			System.out.println("min min: " + mc.computeReachProbs(stpg, target, true, true).soln[0]);
			System.out.println("max min: " + mc.computeReachProbs(stpg, target, false, true).soln[0]);
			System.out.println("min max: " + mc.computeReachProbs(stpg, target, true, false).soln[0]);
			System.out.println("max max: " + mc.computeReachProbs(stpg, target, false, false).soln[0]);

		} catch (FileNotFoundException e) {
			System.out.println("Error: " + e.getMessage());
			System.exit(1);
		} catch (PrismException e) {
			System.out.println("Error: " + e.getMessage());
			System.exit(1);
		}
	}
}
