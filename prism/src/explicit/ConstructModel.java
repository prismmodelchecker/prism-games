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

import java.io.*;
import java.util.*;

import parser.*;
import parser.ast.*;
import prism.*;
import simulator.*;

public class ConstructModel {
	// The simulator engine and a log for output
	private SimulatorEngine engine;
	private PrismLog mainLog;

	// PRISM settings object (optional)
	private PrismSettings settings;

	// Options:
	// Find deadlocks during model construction?
	private boolean findDeadlocks = true;
	// Automatically fix deadlocks?
	private boolean fixDeadlocks = true;
	
	// Basic info needed about model
	// private ModelType modelType;

	// Details of built model
	private List<State> statesList;

	public ConstructModel(SimulatorEngine engine, PrismLog mainLog) {
		this.engine = engine;
		this.mainLog = mainLog;
		settings = null;
	}

	public void setSettings(PrismSettings settings) {
		this.settings = settings;
	}

	public List<State> getStatesList() {
		return statesList;
	}

	public void setFixDeadlocks(boolean b)
	{
		fixDeadlocks = b;
	}
	
	/**
	 * Build the set of reachable states for a PRISM model language description
	 * and return.
	 * 
	 * @param modulesFile
	 *            The PRISM model
	 */
	public List<State> computeReachableStates(ModulesFile modulesFile)
			throws PrismException {
		constructModel(modulesFile, true, false);
		return statesList;
	}

	/**
	 * Construct an explicit-state model from a PRISM model language description
	 * and return.
	 * 
	 * @param modulesFile
	 *            The PRISM model
	 */
	public Model constructModel(ModulesFile modulesFile) throws PrismException {
		return constructModel(modulesFile, false, false, true);
	}

	/**
	 * Construct an explicit-state model from a PRISM model language description
	 * and return. If {@code justReach} is true, no model is built and null is
	 * returned; the set of reachable states can be obtained with
	 * {@link #getStatesList()}.
	 * 
	 * @param modulesFile
	 *            The PRISM model
	 * @param justReach
	 *            If true, just build the reachable state set, not the model
	 * @param buildSparse
	 *            Build a sparse version of the model (if possible)?
	 */
	public Model constructModel(ModulesFile modulesFile, boolean justReach,
			boolean buildSparse) throws PrismException {
		return constructModel(modulesFile, justReach, buildSparse, true);
	}

	/**
	 * Construct an explicit-state model from a PRISM model language description
	 * and return. If {@code justReach} is true, no model is built and null is
	 * returned; the set of reachable states can be obtained with
	 * {@link #getStatesList()}.
	 * 
	 * @param modulesFile
	 *            The PRISM model
	 * @param justReach
	 *            If true, just build the reachable state set, not the model
	 * @param buildSparse
	 *            Build a sparse version of the model (if possible)?
	 * @param distinguishActions
	 *            True if actions should be attached to distributions (and used
	 *            to distinguish them)
	 */
	public Model constructModel(ModulesFile modulesFile, boolean justReach,
			boolean buildSparse, boolean distinguishActions)
			throws PrismException {
		// Model info
		ModelType modelType;
		// State storage
		IndexedSet<State> states;
		LinkedList<State> explore;
		State state, stateNew;
		// Explicit model storage
		ModelSimple modelSimple = null;
		DTMCSimple dtmc = null;
		CTMCSimple ctmc = null;
		MDPSimple mdp = null;
		STPGExplicit stpg = null;
		SMG smg = null;
		ModelExplicit model = null;
		Distribution distr = null;
		// Misc
		int i, j, nc, nt, src, dest, player;
		long timer, timerProgress;
		// int id;
		int nPlayers = 0;
		State tempState;

		// SMG vars
		Map<Integer, List<Integer>> playerChoices = new HashMap<Integer, List<Integer>>();
		List<Integer> originalStates = new ArrayList<Integer>();
		List<Integer> choices;
		boolean lazy = false;
		int schedIndx = 0;

		// Don't support multiple initial states
		if (modulesFile.getInitialStates() != null) {
			throw new PrismException(
					"Cannot do explicit-state reachability if there are multiple initial states");
		}

		// Starting reachability...
		mainLog.print("\nComputing reachable states...\n");
		mainLog.flush();
		timer = timerProgress = System.currentTimeMillis();

		// Initialise simulator for this model
		modelType = modulesFile.getModelType();

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
		// Add initial state to lists/model
		if (modulesFile.getInitialStates() != null) {
			throw new PrismException(
					"Explicit model construction does not support multiple initial states");
		}
		state = modulesFile.getDefaultInitialState();
		states.add(state);
		explore.add(state);
		if (!justReach) {
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

			// Use simulator to explore all choices/transitions from this state
			engine.initialisePath(state);
			nc = engine.getNumChoices();
			// For an STPG, first determine which player owns the state
			if (modelType == ModelType.STPG || modelType == ModelType.SMG) {
				player = -1;
				for (i = 0; i < nc; i++) {
					int iPlayer = determinePlayerForChoice(modulesFile,
							modelType, i);
					if (player != -1 && iPlayer != player) {
						throw new PrismException("Choices for both player "
								+ player + " and " + iPlayer + " in state "
								+ state);
					}
					player = iPlayer;
				}
				if (modelType == ModelType.STPG)
					stpg.setPlayer(src, player);
				else
					smg.setPlayer(src, player);
			}
			// Look at each outgoing choice in turn
			for (i = 0; i < nc; i++) {

				if (!justReach
						&& (modelType == ModelType.MDP
								|| modelType == ModelType.STPG || modelType == ModelType.SMG)) {
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
							dtmc.addToProbability(src, dest,
									engine.getTransitionProbability(i, j));
							break;
						case CTMC:
							ctmc.addToProbability(src, dest,
									engine.getTransitionProbability(i, j));
							break;
						case MDP:
						case STPG:
						case SMG:
							distr.add(dest,
									engine.getTransitionProbability(i, j));
							break;
						}
					}
				}
				if (!justReach) {
					if (modelType == ModelType.MDP) {
						if (distinguishActions) {
							mdp.addActionLabelledChoice(src, distr,
									engine.getTransitionAction(i, 0));
						} else {
							mdp.addChoice(src, distr);
						}
					} else if (modelType == ModelType.STPG) {
						// TODO: need addActionLabelledChoice
						int k = stpg.addChoice(src, distr);
						stpg.setAction(src, k, engine.getTransitionAction(i, 0));
					} else if (modelType == ModelType.SMG) {
						int k = smg.addChoice(src, distr);
						smg.setAction(src, k, engine.getTransitionAction(i, 0));
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

		// Finish progress display
		mainLog.println(" " + (src + 1));

		// Reachability complete
		mainLog.print("Reachable states exploration"
				+ (justReach ? "" : " and model construction"));
		mainLog.println(" done in "
				+ ((System.currentTimeMillis() - timer) / 1000.0) + " secs.");
		// mainLog.println(states);

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
			// mainLog.println(permut);
		} else {
			statesList = states.toArrayList();
		}
		states.clear();
		states = null;
		// mainLog.println(permut);
		// mainLog.println(statesList);

		// Construct new explicit-state model (with correct state ordering)
		if (!justReach) {
			switch (modelType) {
			case DTMC:
				model = new DTMCSimple(dtmc, permut);
				break;
			case CTMC:
				model = new CTMCSimple(ctmc, permut);
				break;
			case MDP:
				if (buildSparse) {
					model = new MDPSparse(mdp, true, permut);
				} else {
					model = new MDPSimple(mdp, permut);
				}
				break;
			case STPG:
				model = new STPGExplicit(stpg, permut);
				break;
			case SMG:
				model = new SMG(smg, permut);
				break;
			}
			model.setStatesList(statesList);
			model.setConstantValues(new Values(modulesFile.getConstantValues()));
			// System.out.println(model);
			// mainLog.println("Model: " + model);
		}

		// Discard permutation
		permut = null;

		return model;
	}

	/**
	 * For game models, determine the player who owns the {@code i}th choice in
	 * the state currently being explored by the simulator. Returns the index
	 * (starting from 1) of the player.
	 */
	private int determinePlayerForChoice(ModulesFile modulesFile,
			ModelType modelType, int i) throws PrismException {
		int modAct, player;

		modAct = engine.getTransitionModuleOrActionIndex(i, 0);
		// Synchronous action
		if (modAct > 0) {
			player = modulesFile.getPlayerForAction(modulesFile
					.getSynch(modAct - 1));
			if (player == -1) {
				throw new PrismException("Action \""
						+ modulesFile.getSynch(modAct - 1)
						+ "\" is not assigned to any player");
			}
			// 0-indexed to 1-indexed
			player++;
		}
		// Asynchronous action
		else {
			player = modulesFile.getPlayerForModule(engine
					.getTransitionModuleOrAction(i, 0));
			if (player == -1) {

				// for backwards compatibility trying to parse player from the
				// module name (e.g., playerX)
				try {
					player = Integer.parseInt(engine
							.getTransitionModuleOrAction(i, 0).substring(6));
				} catch (Exception e) {
					throw new PrismException("Module \""
							+ engine.getTransitionModuleOrAction(i, 0)
							+ "\" is not assigned to any player");
				}
			}
			// 0-indexed to 1-indexed
			player++;
		}

		return player;
	}

	/**
	 * State placeholder
	 * 
	 * @author aissim
	 * 
	 */
	static class DummyState extends State {

		public DummyState(int n) {
			super(n);
		}

		public int compareTo(State s) {
			return -1;
		}

	}

	/**
	 * Test method.
	 */
	public static void main(String[] args) {

		// TODO HACK!
		testSMG(args[0]);
		if (1 == 1)
			return;

		try {
			// Simple example: parse a PRISM file from a file, construct the
			// model and export to a .tra file
			PrismLog mainLog = new PrismPrintStreamLog(System.out);
			Prism prism = new Prism(mainLog, mainLog);
			ModulesFile modulesFile = prism.parseModelFile(new File(args[0]));
			UndefinedConstants undefinedConstants = new UndefinedConstants(
					modulesFile, null);
			if (args.length > 2)
				undefinedConstants.defineUsingConstSwitch(args[1]);
			modulesFile.setUndefinedConstants(undefinedConstants
					.getMFConstantValues());
			ConstructModel constructModel = new ConstructModel(
					prism.getSimulator(), mainLog);
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
	public static void testSTPG(String filename) {
		try {
			// Simple example: parse a PRISM file from a file, construct the
			// model
			// and export to a .tra file
			PrismLog mainLog = new PrismPrintStreamLog(System.out);
			Prism prism = new Prism(mainLog, mainLog);
			ModulesFile modulesFile = prism.parseModelFile(new File(filename));
			UndefinedConstants undefinedConstants = new UndefinedConstants(
					modulesFile, null);

			modulesFile.setUndefinedConstants(undefinedConstants
					.getMFConstantValues());

			ConstructModel constructModel = new ConstructModel(
					prism.getSimulator(), mainLog);

			Model model = constructModel.constructModel(modulesFile);

			STPGExplicit stpg = (STPGExplicit) model;
			mainLog.println(stpg);
			mainLog.println(constructModel.getStatesList());
			STPGModelChecker mc = new STPGModelChecker();

			mc.setLog(prism.getMainLog());
			BitSet target = new BitSet();
			target.set(8);
			stpg.exportToDotFile("stpg.dot", target);
			System.out.println("min min: "
					+ mc.computeReachProbs(stpg, target, true, true).soln[0]);
			System.out.println("max min: "
					+ mc.computeReachProbs(stpg, target, false, true).soln[0]);
			System.out.println("min max: "
					+ mc.computeReachProbs(stpg, target, true, false).soln[0]);
			System.out.println("max max: "
					+ mc.computeReachProbs(stpg, target, false, false).soln[0]);

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
	 * 
	 * @param filename
	 */
	public static void testSMG(String filename) {
		try {
			// Simple example: parse a PRISM file from a file, construct the
			// model
			// and export to a .tra file
			PrismLog mainLog = new PrismPrintStreamLog(System.out);
			Prism prism = new Prism(mainLog, mainLog);
			ModulesFile modulesFile = prism.parseModelFile(new File(filename));
			UndefinedConstants undefinedConstants = new UndefinedConstants(
					modulesFile, null);

			modulesFile.setUndefinedConstants(undefinedConstants
					.getMFConstantValues());

			ConstructModel constructModel = new ConstructModel(
					prism.getSimulator(), mainLog);

			Model model = constructModel.constructModel(modulesFile);

			System.out.println("\nStates list: ");
			mainLog.println(constructModel.getStatesList());

			SMG smg = (SMG) model;

			Set<Integer> player1 = new HashSet<Integer>();
			player1.add(1);
			SMG stpg_nondet = smg.clone();

			player1 = new HashSet<Integer>();
			player1.add(1);
			SMG stpg_rand = smg.clone();

			System.out.println("\nSMG: ");
			mainLog.println(smg);
			System.out.println("\nSTPG-NONDET: ");
			mainLog.println(stpg_nondet);
			System.out.println("\nSTPG-RAND: ");
			mainLog.println(stpg_rand);

			// smg.exportToDotFile("smg.dot");
			// stpg_nondet.exportToDotFile("stpg_nondet.dot");
			// stpg_rand.exportToDotFile("stpg_rand.dot");

		} catch (FileNotFoundException e) {
			System.out.println("Error: " + e.getMessage());
			System.exit(1);
		} catch (PrismException e) {
			System.out.println("Error: " + e.getMessage());
			System.exit(1);
		}
	}
}
