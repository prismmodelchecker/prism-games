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
	 */
	public List<State> computeReachableStates(ModulesFile modulesFile) throws PrismException
	{
		constructModel(modulesFile, true, false);
		return statesList;
	}

	/**
	 * Construct an explicit-state model from a PRISM model language description and return.
	 * @param modulesFile The PRISM model
	 */
	public Model constructModel(ModulesFile modulesFile) throws PrismException
	{
		return constructModel(modulesFile, false, false);
	}

	/**
	 * Construct an explicit-state model from a PRISM model language description and return.
	 * If {@code justReach} is true, no model is built and null is returned;
	 * the set of reachable states can be obtained with {@link #getStatesList()}.
	 * @param modulesFile The PRISM model
	 * @param justReach If true, just build the reachable state set, not the model
	 * @param buildSparse Build a sparse version of the model (if possible)?
	 */
	public Model constructModel(ModulesFile modulesFile, boolean justReach, boolean buildSparse) throws PrismException
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
		SMG smg = null;
		Model model = null;
		Distribution distr = null;
		// Misc
		int i, j, k, nc, nt, src, dest, prev, tmp;
		long timer, timerProgress;
		boolean fixdl = false;
		String actionLabel = null;
		int player;
		int id;
		int nPlayers = 0;
		State tempState;

		// SMG vars
		Map<Integer, List<Integer>> playerChoices = new HashMap<Integer, List<Integer>>();
		List<Integer> originalStates = new ArrayList<Integer>();
		List<Integer> choices;
		boolean lazy = false;
		int schedIndx = 0;;

		// Process game info
		this.player1mods = new HashSet<String>();
		this.player2mods = new HashSet<String>();
		this.player1syncs = new HashSet<String>();
		this.player2syncs = new HashSet<String>();

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
		modelType = ModelType.SMG;

		if (modelType == ModelType.SMG) {
			// adding a scheduler variable to the model
			lazy = false;
			nPlayers = 2;
			schedIndx = 0;
			modulesFile.addGlobal(new Declaration("_sched", new DeclarationInt(Expression.Int(0), Expression
					.Int(nPlayers))));
			modulesFile.tidyUp();
		}
		if (modelType == ModelType.STPG) {
			loadGameParams(settings.getString(PrismSettings.PRISM_GAME_OPTIONS));
		}

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
			case SMG:
				modelSimple = smg = new SMG();
				break;
			}
		}

		// Initialise states storage
		states = new IndexedSet<State>(true);
		explore = new LinkedList<State>();
		stateLabels = new ArrayList<Integer>();
		// Add initial state to lists/model
		if (modulesFile.getInitialStates() != null) {
			throw new PrismException("Explicit model construction does not support multiple initial states");
		}
		state = modulesFile.getDefaultInitialState();
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

		// constructing stochastic multiplayer game model
		if (modelType == ModelType.SMG) {
			src = -1;
			originalStates.add(0);
			while (!explore.isEmpty()) {

				// Pick next state to explore
				state = explore.removeFirst();
				src++;

				// Use simulator to explore all choices/transitions from this state
				engine.initialisePath(state);

				// get number of choices
				nc = engine.getNumChoices();

				// Group choices by player
				for (i = 0; i < nc; i++) {
					player = Integer.parseInt(engine.getTransitionModuleOrAction(i, 0).substring(6));
					if (!playerChoices.containsKey(player))
						playerChoices.put(player, new ArrayList<Integer>());
					playerChoices.get(player).add(i);
				}

				// adding states for 'lazy' players
				if (lazy && playerChoices.size() < nPlayers) {
					choices = new ArrayList<Integer>(playerChoices.keySet());
					Collections.sort(choices);

					for (i = 0, j = 0; i < nPlayers; i++) {
						if ((j < choices.size() && choices.get(j) != i + 1) || !(j < choices.size())) {

							k = smg.addState(i + 1);
							tempState = new State(state);
							tempState.setValue(schedIndx, i + 1); // last variable is _sched
							states.add(tempState);
							if (!justReach) {

								// add transition to new state
								distr = new Distribution();
								distr.add(k, 1.0);
								smg.setAction(originalStates.get(src), smg.addChoice(originalStates.get(src), distr),
										"player" + (i + 1));

								// add transition back
								distr = new Distribution();
								distr.add(originalStates.get(src), 1.0);
								smg.setAction(k, smg.addChoice(k, distr), "player" + (i + 1));
							}
						}

					}

				}

				// add a state for each player and transition to it from current state
				for (Integer pl : playerChoices.keySet()) {

					// create a state and add transition to it from the current state
					k = smg.addState(pl);

					// creating a new state
					tempState = new State(state);
					tempState.setValue(schedIndx, pl); // last variable is _sched
					states.add(tempState);

					if (!justReach) {
						distr = new Distribution();
						distr.add(k, 1.0);
						smg.setAction(originalStates.get(src), smg.addChoice(originalStates.get(src), distr), "player"
								+ pl);
					}

					// copy the choices by player pl from the original state to the new one
					for (Integer ch : playerChoices.get(pl)) {

						if (!justReach)
							distr = new Distribution();

						// Look at each transition in the choice
						nt = engine.getNumTransitions(ch);
						for (j = 0; j < nt; j++) {
							stateNew = engine.computeTransitionTarget(ch, j);

							// Is this a new state?
							if (states.add(stateNew)) {
								// If so, add to the explore list
								explore.add(stateNew);
								// And to model
								originalStates.add(modelSimple.addState());
							}

							// Get index of state in state set
							dest = states.getIndexOfLastAdd();

							if (!justReach) {
								// Add transitions to model
								distr.add(dest, engine.getTransitionProbability(ch, j));
							}
						}

						if (!justReach) {
							// adding a choice to the state
							smg.setAction(k, smg.addChoice(k, distr), engine.getTransitionModuleOrAction(ch, 0));
						}
					}
				}

				playerChoices.clear();
			}

		} else {
			// constructing all other models

			// Explore...
			src = -1;
			while (!explore.isEmpty()) {
				// Pick next state to explore
				// (they are stored in order found so know index is src+1)
				state = explore.removeFirst();
				src++;

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
						} else if (modelType == ModelType.SMG) {
							k = smg.addChoice(src, distr);
							smg.setAction(src, k, engine.getTransitionModuleOrAction(i, 0));
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
			case SMG:
				model = new SMG(smg, permut);
				((ModelSimple) model).statesList = statesList;
				((ModelSimple) model).constantValues = new Values(modulesFile.getConstantValues());
				break;
			}
			//			mainLog.println("Model: " + model);
		}

		// Discard permutation
		permut = null;

		return model;
	}

	/**
	 * State placeholder
	 * @author aissim
	 *
	 */
	static class DummyState extends State
	{

		public DummyState(int n)
		{
			super(n);
		}

		public int compareTo(State s)
		{
			return -1;
		}

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
		testSMG(args[0]);
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
			Model model = constructModel.constructModel(modulesFile);
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

			Model model = constructModel.constructModel(modulesFile);

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

	/**
	 * Method to test SMG construction
	 * @param filename
	 */
	public static void testSMG(String filename)
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

			Model model = constructModel.constructModel(modulesFile);


			System.out.println("\nStates list: ");
			mainLog.println(constructModel.getStatesList());

			SMG smg = (SMG) model;

			Set<Integer> player1 = new HashSet<Integer>();
			player1.add(1);
			SMG stpg_nondet = smg.clone().reduceToSTPG(player1, SMG.SCHED_NONDET);
			
			player1 = new HashSet<Integer>();
			player1.add(1);
			SMG stpg_rand = smg.clone().reduceToSTPG(player1, SMG.SCHED_RANDOM);
			
			
			System.out.println("\nSMG: ");
			mainLog.println(smg);
			System.out.println("\nSTPG-NONDET: ");
			mainLog.println(stpg_nondet);
			System.out.println("\nSTPG-RAND: ");
			mainLog.println(stpg_rand);
			
			
			
			smg.exportToDotFile("/auto/users/aissim/Desktop/smg.dot");
			stpg_nondet.exportToDotFile("/auto/users/aissim/Desktop/stpg_nondet.dot");
			stpg_rand.exportToDotFile("/auto/users/aissim/Desktop/stpg_rand.dot");

			
		} catch (FileNotFoundException e) {
			System.out.println("Error: " + e.getMessage());
			System.exit(1);
		} catch (PrismException e) {
			System.out.println("Error: " + e.getMessage());
			System.exit(1);
		}
	}
}
