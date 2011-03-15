// ==============================================================================
//	
// Copyright (c) 2002-
// Authors:
// * Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford)
//	
// ------------------------------------------------------------------------------
//	
// This file is part of PRISM.
//	
// PRISM is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//	
// PRISM is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//	
// You should have received a copy of the GNU General Public License
// along with PRISM; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//	
// ==============================================================================

package explicit;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import parser.State;
import parser.Values;
import parser.ast.ModulesFile;
import prism.ModelType;
import prism.Prism;
import prism.PrismException;
import prism.PrismLog;
import prism.PrismPrintStreamLog;
import prism.UndefinedConstants;
import simulator.SimulatorEngine;

public class ConstructModel
{
	// The simulator engine and a log for output
	private SimulatorEngine engine;
	private PrismLog mainLog;

	// Basic info needed about model
	// private ModelType modelType;

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

		this.player1mods = new HashSet<String>();
		this.player2mods = new HashSet<String>();

		this.player1syncs = new HashSet<String>();
		this.player2syncs = new HashSet<String>();

		initSensors();
	}

	// test method to initialise module and action allocation
	private void initSensors()
	{
		player1mods.add("scheduler");
		player1mods.add("task_generator");
		player1mods.add("sensor1");
		player1mods.add("sensor3");
		player1mods.add("sensor5");
		player1mods.add("sensor7");

		player2mods.add("sensor2");
		player2mods.add("sensor4");
		player2mods.add("sensor6");

		player1syncs.add("initialise");
		player1syncs.add("scheduling");

		player1syncs.add("str1");
		player1syncs.add("str2");
		player1syncs.add("str3");
		player1syncs.add("str4");
		player1syncs.add("str5");
		player1syncs.add("str6");
		player1syncs.add("str7");

		player1syncs.add("fin1");
		player1syncs.add("fin2");
		player1syncs.add("fin3");
		player1syncs.add("fin4");
		player1syncs.add("fin5");
		player1syncs.add("fin6");
		player1syncs.add("fin7");
	}

	public List<State> getStatesList()
	{
		return statesList;
	}

	/**
	 * Build the set of reachable states for a PRISM model language description
	 * and return.
	 * 
	 * @param modulesFile The PRISM model
	 * @param initialState The initial state (for reachability)
	 */
	public List<State> computeReachableStates(ModulesFile modulesFile, Values initialState) throws PrismException
	{
		constructModel(modulesFile, initialState, true, false);
		return statesList;
	}

	/**
	 * Construct an explicit-state model from a PRISM model language description
	 * and return.
	 * 
	 * @param modulesFile The PRISM model
	 * @param initialState The initial state (for reachability)
	 */
	public Model constructModel(ModulesFile modulesFile, Values initialState) throws PrismException
	{
		return constructModel(modulesFile, initialState, false, false);
	}

	/**
	 * Construct an explicit-state model from a PRISM model language description
	 * and return. If {@code justReach} is true, no model is built and null is
	 * returned; the set of reachable states can be obtained with
	 * {@link #getStatesList()}.
	 * 
	 * @param modulesFile The PRISM model
	 * @param initialState The initial state (for reachability)
	 * @param justReach If true, just build the reachable state set, not the
	 *          model
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
		int i, j, nc, nt, src, dest;
		long timer, timerProgress;
		boolean fixdl = false;
		String actionLabel = null;
		int player;
		int id;

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
		modelType = ModelType.STPG;

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
					if (modelType == ModelType.MDP)
						mdp.addChoice(src, distr);
					else if (modelType == ModelType.STPG)
						stpg.addDistribution(src, distr);
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
				model = modelSimple;
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
			// Simple example: parse a PRISM file from a file, construct the model
			// and export to a .tra file
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
			ModelCheckerResult res = mc.probReach(mdp, target, true);
			mainLog.println(res.soln);

		} catch (FileNotFoundException e) {
			System.out.println("Error: " + e.getMessage());
			System.exit(1);
		} catch (PrismException e) {
			System.out.println("Error: " + e.getMessage());
			System.exit(1);
		}
	}

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
			System.out.println("min min: " + mc.probReach(stpg, target, true, true).soln[0]);
			System.out.println("max min: " + mc.probReach(stpg, target, false, true).soln[0]);
			System.out.println("min max: " + mc.probReach(stpg, target, true, false).soln[0]);
			System.out.println("max max: " + mc.probReach(stpg, target, false, false).soln[0]);

		} catch (FileNotFoundException e) {
			System.out.println("Error: " + e.getMessage());
			System.exit(1);
		} catch (PrismException e) {
			System.out.println("Error: " + e.getMessage());
			System.exit(1);
		}
	}
}
