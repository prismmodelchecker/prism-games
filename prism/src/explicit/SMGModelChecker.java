//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford)
// * Clemens Wiltsche <clemens.wiltsche@cs.ox.ac.uk> (University of Oxford)
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

import java.math.BigInteger;
import java.util.AbstractMap.SimpleEntry;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.apache.commons.math3.fraction.BigFraction;

import acceptance.AcceptanceReach;
import explicit.rewards.ConstructRewards;
import explicit.rewards.Rewards;
import explicit.rewards.SMGRewards;
import explicit.rewards.SMGRewardsSimple;
import explicit.rewards.StateRewardsConstant;
import parma_polyhedra_library.C_Polyhedron;
import parma_polyhedra_library.Coefficient;
import parma_polyhedra_library.Constraint;
import parma_polyhedra_library.Constraint_System;
import parma_polyhedra_library.Generator;
import parma_polyhedra_library.Generator_System;
import parma_polyhedra_library.Generator_Type;
import parma_polyhedra_library.Linear_Expression;
import parma_polyhedra_library.Linear_Expression_Coefficient;
import parma_polyhedra_library.Linear_Expression_Sum;
import parma_polyhedra_library.Linear_Expression_Times;
import parma_polyhedra_library.Polyhedron;
import parma_polyhedra_library.Relation_Symbol;
import parma_polyhedra_library.Variable;
import parser.ast.Coalition;
import parser.ast.Expression;
import parser.ast.ExpressionConstant;
import parser.ast.ExpressionProb;
import parser.ast.ExpressionQuant;
import parser.ast.ExpressionReward;
import parser.ast.ExpressionTemporal;
import parser.ast.ExpressionVar;
import parser.ast.RelOp;
import parser.ast.RewardStruct;
import prism.ModelType;
import prism.OpRelOpBound;
import prism.PointList;
import prism.PrismComponent;
import prism.PrismException;
import prism.PrismLangException;
import prism.PrismSettings;
import prism.PrismUtils;
import strat.StochasticUpdateStrategy;

/**
 * Explicit-state model checker for multi-player stochastic games (SMGs).
 */
public class SMGModelChecker extends ProbModelChecker
{
	// Multi-Objective Synthesis
	protected long maxCIter = 500;
	protected int maxDIter = 100;
	protected int dIterOffset = 1;
	protected int maxRIter = 500;
	protected int minM = 2;
	protected int maxM = 16;
	protected double increase_factor = 1.01;
	protected long max_accuracy = Integer.MAX_VALUE / 4;
	protected boolean gaussSeidel = true;

	// logging options for Pareto sets and Strategy
	protected boolean logCPareto = false;
	protected boolean logDPareto = false;
	protected boolean logRPareto = false;
	protected boolean logStrategy = false;

	// relative termination criterion
	protected double varepsilon = 0.0001;

	// tracking for issuing warnings in batch
	private List<String> strictToNonstrict = new ArrayList<String>();
	private List<String> unfolded = new ArrayList<String>();

	/**
	 * Create a new SMGModelChecker, inherit basic state from parent (unless null).
	 */
	public SMGModelChecker(PrismComponent parent) throws PrismException
	{
		super(parent);
		if (settings != null) {
			varepsilon = settings.getDouble(PrismSettings.PRISM_PARETO_EPSILON);
			if (varepsilon < 0)
				throw new PrismException(String.format("Must have %s >= 0", settings.getSettingName(PrismSettings.PRISM_PARETO_EPSILON)));
			gaussSeidel = settings.getBoolean(PrismSettings.PRISM_MULTI_GAUSS_SEIDEL);
			maxCIter = settings.getInteger(PrismSettings.PRISM_MULTI_MAX_C_ITER);
			if (maxCIter < 1)
				throw new PrismException(String.format("Must have %s >= 1", settings.getSettingName(PrismSettings.PRISM_MULTI_MAX_C_ITER)));
			maxRIter = settings.getInteger(PrismSettings.PRISM_MULTI_MAX_R_ITER);
			if (maxRIter < 1)
				throw new PrismException(String.format("Must have %s >= 1", settings.getSettingName(PrismSettings.PRISM_MULTI_MAX_R_ITER)));
			maxDIter = settings.getInteger(PrismSettings.PRISM_MULTI_MAX_D_ITER);
			if (maxDIter < 1)
				throw new PrismException(String.format("Must have %s >= 1", settings.getSettingName(PrismSettings.PRISM_MULTI_MAX_D_ITER)));
			dIterOffset = settings.getInteger(PrismSettings.PRISM_MULTI_D_ITER_OFFSET);
			if (dIterOffset < 1)
				throw new PrismException(String.format("Must have %s >= 1", settings.getSettingName(PrismSettings.PRISM_MULTI_D_ITER_OFFSET)));
			minM = settings.getInteger(PrismSettings.PRISM_MULTI_MIN_M);
			maxM = settings.getInteger(PrismSettings.PRISM_MULTI_MAX_M);
			if (maxM < minM || minM < 2)
			        throw new PrismException("Box size parameters invalid");
			logCPareto = settings.getBoolean(PrismSettings.LOG_MULTI_C_PARETO);
			logDPareto = settings.getBoolean(PrismSettings.LOG_MULTI_D_PARETO);
			logRPareto = settings.getBoolean(PrismSettings.LOG_MULTI_R_PARETO);
			logStrategy = settings.getBoolean(PrismSettings.LOG_MULTI_STRATEGY);
			increase_factor = settings.getDouble(PrismSettings.PRISM_MULTI_INCREASE_FACTOR);
			if (increase_factor < 1)
				throw new PrismException(String.format("Must have %s >= 1", settings.getSettingName(PrismSettings.PRISM_MULTI_INCREASE_FACTOR)));
			max_accuracy = Integer.MAX_VALUE / 4;
		}
	}

	@Override
	public StateValues checkExpressionMultiObjective(Model model, List<List<Expression>> cnf, Coalition coalition) throws PrismException
        {
	        // initialise the Parma Polyhedra Library
	        PPLSupport.initPPL();

		// extract simple expression from MQ
		MultiParameters params = initialiseRewards(model, cnf);

		// direct method [QEST'13, MFCS'13, TACAS'15]
		return checkExpressionMultiDirect(model, params, coalition);
	}
	
	// Model checking functions
	
	/**
	 * Gets the states in the maximal irreducible component (MIC) induced by pivot {@code t}.
	 * 
	 * @param smg The game for which to compute the MIC
	 * @param reach The set of states reached from t while remaining in states of {@code reach}.
	 * @param t The pivot.
	 * 
	 * @return Whether a MIC is induced from t. The set of states in the MIC is in {@code reach}.
	 */
	private boolean getMIC(SMG smg, BitSet reach, int t) throws PrismException
	{
		throw new UnsupportedOperationException();
		
		// Disabled for now because support disabling transitions has been removed
		// during refactoring of explicit engine model classes
		
		/*int gameSize = smg.numStates;

		Map<Integer, BitSet> disabledChoices;
		boolean someChoicesDisabled;

		BitSet reach_temp = new BitSet();
		reach_temp.or(reach);

		BitSet target = new BitSet(gameSize);
		target.set(t);

		go_through_states: for (int s = 0; s < gameSize; s++) {
			if (reach.get(s) && s != t) { // don't remove pivot
				disabledChoices = new HashMap<Integer, BitSet>(smg.disabledChoices); // keep for backtracking
				someChoicesDisabled = smg.someChoicesDisabled;

				BitSet disabled = smg.disableState(s, false, true); // close for P2

				if (disabled.get(t)) { // if pivot was removed, pick a different s, need to restore disabled transitions
					smg.disabledChoices = disabledChoices;
					smg.someChoicesDisabled = someChoicesDisabled;
					continue go_through_states;
				}

				reach.andNot(disabled); // remove states that were disabled

				double[] minmaxreach = createSTPGModelChecker().computeReachProbs((STPG) smg, reach, target, true, false, null, null).soln; // for all P1 strategies, for some P2 strategy

				for (int r = 0; r < gameSize; r++) {
					if (PrismUtils.doublesAreEqual(minmaxreach[r], 0.0)) { // pivot t not reached from s ... not an IC
						// try making IC by removing another state
						if (getMIC(smg, reach, t)) {
							return true;
						} else {
							smg.disabledChoices = disabledChoices;
							smg.someChoicesDisabled = someChoicesDisabled;
							continue go_through_states;
						}
					}
				}

				// at this point we have an IC
				return true;
			}
		}

		// if no more states can be removed, and not an IC
		// reset reach
		reach.clear();
		reach.or(reach_temp);
		return false;*/
	}

	/**
	* Checks whether the game is controllable multichain (CM),
	* Note that the current implementation is sound but not complete, in the sense
	* that it may return false even if the game is CM.
	*
	* @param smg The game to be checked.
	*
	* @return True if the game is controllable multichain.
	**/
	private boolean isControllableMultichain(SMG smg) throws PrismException
	{
		int gameSize = smg.numStates;
		int initial = smg.getFirstInitialState();

		// first need all strongly connected components (SCCs)
		// in the game induced by uniformly randomising strategies
		// i.e., the transients should be avoided
		SCCConsumerStore sccStore = new SCCConsumerStore();
		SCCComputer sccComputer = SCCComputer.createSCCComputer(this, smg, sccStore);
		sccComputer.computeSCCs();
		List<BitSet> sccs = sccStore.getSCCs();
		BitSet inSCC = new BitSet(); // states not in some SCC
		int n = sccs.size();
		for (int i = 0; i < n; i++)
			inSCC.or(sccs.get(i));

		// then want a list of P2-closed subtrees of the game,
		// that exclusively lie in the sccs
		Set<BitSet> subtrees = new HashSet<BitSet>();
		for (BitSet scc : sccs) {
			// go through all P1 states in the SCC
			boolean p1state_exists = false;
			BitSet ui = new BitSet(gameSize);
			BitSet uj = new BitSet(gameSize);
			BitSet covered = new BitSet(gameSize);
			for (int i = 0; i < gameSize; i++) {
				if (scc.get(i) && smg.getPlayer(i) == 1) {
					p1state_exists = true;
					// first get all states reachable in one step from i
					ui.clear();
					ui.set(i);
					uj.clear();
					((STPGExplicit) smg).reachpositivestep(ui, false, true, uj);
					for (int j = uj.nextSetBit(0); j != -1; j = uj.nextSetBit(j + 1)) { // go through all states reachable in one step from i
						if (smg.getPlayer(j) != 1 && !covered.get(j)) {
							BitSet subtree = getSubtree((STPG) smg, j, 2); // get P2-closed subtree
							subtree.set(i); // put the P1 state (the initial root) into the subtree
							subtrees.add(subtree); // get P2-closed subtree
							covered.or(subtree);
						}
					}

				}
			}
			if (!p1state_exists) { // if no P1 state was found, add the full SCC
				subtrees.add(scc);
			}
		}

		// keep whether CM, to be able to print all unreachable trees
		boolean result = true;
		HashSet printedSubtrees = new HashSet<BitSet>(subtrees.size());
		boolean print_all_counterexamples = false;
		// once the subtrees are found, can perform reachability on them
		int tmpverbosity = this.getVerbosity();
		boolean genStrat = generateStrategy;
		try {
			generateStrategy = false; // turn off strategy generation for STPGModelChecker
			this.setVerbosity(0); // temporarily turn off logger
			go_through_subtrees: for (BitSet subtree : subtrees) {
				double[] maxminreach = createSTPGModelChecker().computeReachProbs((STPG) smg, null, subtree, false, true, null, null).soln; // for some P1 strategy (max), for all P2 strategies (min)
				// need to chech reachability from each state
				for (int s = inSCC.nextSetBit(0); s != -1; s = inSCC.nextSetBit(s + 1)) { // go through all non-transient states
					if (!PrismUtils.doublesAreEqual(maxminreach[s], 1.0)) { // cannot reach the subtree
						if (!printedSubtrees.contains(subtree))
							mainLog.print(String.format("Unreachable subtree: %s\n", subtree));
						printedSubtrees.add(subtree);
						result = false;
						if (!print_all_counterexamples)
							break go_through_subtrees;
					}
				}
			}
		} finally { // do "cleanup" by switching logger on again, and set strategy generation to previous setting
			generateStrategy = genStrat;
			this.setVerbosity(tmpverbosity);
		}

		// here result is only true if no subtree was unreachable
		return result;
	}

	/**
	 * Determine the states of an STPG which, with min/max probability strictly greater than 0,
	 * are reached from a state in {@code origin}, while remaining in sates of {@code remain}.
	 * @param stpg The STPG
	 * @param remain The states to remain in
	 * @param origin Origin states
	 * @param min1 Min or max probabilities for player 1 (true=lower bound, false=upper bound)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 */
	public BitSet reachpositive(STPG stpg, BitSet remain, BitSet origin, boolean min1, boolean min2)
	{
		int n, iters;
		BitSet u, soln;
		boolean u_done;
		long timer;

		// Start precomputation
		timer = System.currentTimeMillis();
		if (verbosity >= 1)
			mainLog.println("Starting ReachPositive (" + (min1 ? "min" : "max") + (min2 ? "min" : "max") + ")...");

		// Special case: no origin states
		if (origin.cardinality() == 0) {
			soln = new BitSet(stpg.getNumStates());
			return soln; // all zeroes
		}

		// Initialise vectors
		n = stpg.getNumStates();
		u = new BitSet(n);
		soln = new BitSet(n);

		// Fixed point loop
		iters = 0;
		u_done = false;
		// Least fixed point - should start from 0 but we optimise by
		// starting from 'origin', thus bypassing first iteration
		u.or(origin);
		soln.or(origin);
		while (!u_done) {
			iters++;
			// Single step of ReachPositive
			((STPGExplicit) stpg).reachpositivestep(u, min1, min2, soln);
			// ensure to remain in set
			if (remain != null)
				soln.and(remain);
			// Check termination
			u_done = soln.equals(u);
			// u = soln
			u.clear();
			u.or(soln);
		}

		// Finished precomputation
		timer = System.currentTimeMillis() - timer;
		if (verbosity >= 1) {
			mainLog.print("ReachPositive (" + (min1 ? "min" : "max") + (min2 ? "min" : "max") + ")");
			mainLog.println(" took " + iters + " iterations and " + timer / 1000.0 + " seconds.");
		}

		return u;
	}

	/**
	* @param u The subtree so far
	* @param closedPlayer Player for which subtree is closed
	* @param result The subtree after extending
	 **/
	public BitSet getSubtree(STPG stpg, int s, int closedPlayer)
	{
		int n;
		BitSet u, soln;
		boolean u_done;

		// Initialise vectors
		n = stpg.getNumStates();
		u = new BitSet(n);
		soln = new BitSet(n);

		// iterative loop, starting with only s in the subtree
		u_done = false;
		u.set(s);
		soln.set(s);
		while (!u_done) {
			// Single step of ReachPositive
			((STPGExplicit) stpg).subtreeStep(u, closedPlayer, soln);
			// Check termination
			u_done = soln.equals(u);
			// u = soln
			u.clear();
			u.or(soln);
		}

		return u;
	}

	/**
	* Extracts simple expressions from a multi-objective query
	* and puts the structures required for the multi-objective engine
	* in a MultiParameters object that is returned.
	*
	* @param model The SMG model.
	* @param expr The expression to be checked.
	**/
	private MultiParameters initialiseRewards(Model model, List<List<Expression>> expr) throws PrismException
        {
		MultiParameters params = new MultiParameters(); // return value

		params.expr = expr;
		params.rewards = new ArrayList<SMGRewards>(); // the reward structures
		params.reward_names = new ArrayList<String>();
		params.divisors = new ArrayList<SMGRewards>(); // the divisors of ratio rewards (null if not present)
		params.divisor_names = new ArrayList<String>();
		params.shifts = new ArrayList<Double>();
		params.reward_types = new ArrayList<Integer>();
		params.directions = new ArrayList<Integer>(); // direction of inequality: 1 max, -1 min
		params.expressions = new ArrayList<Expression>(); // for displaying the Pareto set
		params.bounds = new ArrayList<Double>(); // the required bounds
		params.nested = new ArrayList<Boolean>();
		params.varepsilon = varepsilon;
		params.increase_factor = increase_factor;
		params.dIterOffset = dIterOffset;


		// reset warnings
		strictToNonstrict.clear();
		unfolded.clear();

		mainLog.print(String.format("expr: %s\n", expr));

		// dimensionality of the problem
		int CONJUNCTS = expr.size();
		int[] DISJUNCTS = new int[CONJUNCTS];
		int i = 0;
		for(List<Expression> expr_i : expr) {
		    DISJUNCTS[i] = expr_i.size();
		    int j = 0;
		    for(Expression expr_ij : expr_i) {
			if(expr_ij instanceof ExpressionQuant) {
			    parseExpressionQuant(model, (ExpressionQuant) expr_ij, params);
			    params.expressions.add(expr_ij);
			} else {
			    throw new PrismException("Only the R operator is supported.");
			}
			j++;
		    }
		    i++;
		}

		// issue new warnings
		if (strictToNonstrict.size() > 0) {
			String warning = "Strict inequalities ignored and turned into nonstrict inequalities:\n";
			for (String w : strictToNonstrict)
				warning += String.format("\t%s\n", w);
			mainLog.printWarning(warning);
		}
		strictToNonstrict.clear();

		if (unfolded.size() > 0) {
			String warning = "Expressions not folded in Pareto set:\n";
			for (String w : unfolded)
				warning += String.format("\t%s\n", w);
			mainLog.printWarning(warning);
		}
		unfolded.clear();

		// engine dimensions of objective
		params.CONJUNCTS = CONJUNCTS;
		params.DISJUNCTS = DISJUNCTS;

		// set iteration counts, unless overridden in specification
		params.maxCIter = maxCIter;
		params.maxDIter = maxDIter;
		params.maxRIter = maxRIter;
		params.rounding = settings.getBoolean(PrismSettings.PRISM_MULTI_ROUNDING);
		params.baseline_accuracy = settings.getInteger(PrismSettings.PRISM_MULTI_BASELINE_ACCURACY);

		// set objective type
		setObjectiveType(params);

		// assume consistent objectives
		if (params.objective_type == MultiParameters.ETCR) { // total reward
			// check if stopping game
		} else if (params.objective_type == MultiParameters.EAR ||
			   params.objective_type == MultiParameters.ERCR) { // expected average or ratio
		    mainLog.printWarning("Ensure that game is controllable multichain. Not tested.");
		    //if(!isControllableMultichain((SMG)model))
		    //mainLog.printWarning("Cannot determine whether game is controllable multichain.");
		} else if (params.objective_type == MultiParameters.PAR ||
			   params.objective_type == MultiParameters.PRCR) { // satisfaction average or ratio
		    if(!isConjunction(params)) // query must be conjunction
			throw new PrismException("Only conjunctions supported for satisfaction objectives");
		    if(computePareto) // if Pareto computation requested, need equivalence with expectations
			mainLog.printWarning("Ensure that game is controllable multichain. Not tested.");
		    //if(!isControllableMultichain((SMG)model))
		    //    mainLog.printWarning("Cannot determine whether game is controllable multichain.");
		} else if (params.objective_type == -1) { // mixed objectives
			throw new PrismException("Cannot mix objective types");
		}

		return params;
	    
	}


	protected StateValues checkExpressionMultiDirect(Model model, MultiParameters params, Coalition coalition) throws PrismException
	{
		return checkExpressionMultiDirect(model, params, coalition, false);
	}

	/**
	 * arguments:
	 *    @param model The stochastic game
	 *    @param params The parameters with the parsed expression
	 *    @param checkBounds If true then the iteration stops after the values in double[] bound are hit,
	 *                    further, if true, then no Pareto set is displayed
	 **/
	protected StateValues checkExpressionMultiDirect(Model model, MultiParameters params, Coalition coalition, boolean checkBounds) throws PrismException
	{
 	        int init = model.getFirstInitialState();

		mainLog.print(String.format("/////////////////   NEW (DIRECT) MODEL CHECKING TASK     /////////////////////\n"));
		mainLog.print(String.format("Property:\n\t%s\n\n", params.expr));
		mainLog.print(String.format("initial state: %d\n", init));
		mainLog.print(String.format("operation: %s\n", computePareto ? "Pareto set computation" : generateStrategy ? "Strategy generation" : "Verification"));
		mainLog.flush();

		// check model type
		if (model.getModelType() != ModelType.SMG)
			throw new PrismException("Only SMGs supported by multi-objective engine.");

		// make sure there is only one initial state
		if (model.getNumInitialStates() != 1)
			throw new PrismException("Multi-objective engine supports only models with a single initial state.");

		// set coalition:
		if (coalition != null)
			((SMG) model).setCoalition(coalition);

		// If no info, player 1 is alone by default
		else {
			Coalition coalitionNew = new Coalition();
			List<String> coalitionNewPlayers = new ArrayList<String>();
			coalitionNewPlayers.add("1");
			coalitionNew.setPlayers(coalitionNewPlayers);
			((SMG) model).setCoalition(coalitionNew);
		}

		if (computePareto) {
			// COMPUTE PARETO SET
			double t0 = (double) System.nanoTime();
			setRewardBrackets(params, model); // find brackets around rewards to start iteration and/or rounding
			Pareto[] pareto = computeRatioMQParetoSet((SMG) model, params);
			mainLog.print(String.format("Pareto set computation took %f s\n", ((double) (System.nanoTime() - t0)) / 1e9));

			// after computation need to turn the negative directions around
			interiorPareto(pareto);
			turnParetoAndBounds(pareto, params);

			// display Pareto set
			pareto_set = pareto[init]; // register Pareto set - for compositional mainly
			parsed_params = params; // register parameters
			if (!checkBounds) {
				PointList.addStoredPointList("M", new PointList(pareto_set, params.expressions, params.bounds));
			}
			mainLog.print("Resulting Pareto set:\n");
			PPLSupport.printReachabilityPolyhedron(pareto, params.rewards.size(), init, mainLog);
			mainLog.flush();
			return StateValues.createFromParetoArray(pareto, model);
		} else {
			// CHECK MQ AND COMPUTE STRATEGY (if requested)
			double t0 = (double) System.nanoTime();
			// convert the ratio objectives in the MQ to average objectives
			MultiParameters new_params = convertRatioMQToMQ((SMG) model, params);

			setRewardBrackets(new_params, model); // find brackets around rewards to start iteration and/or rounding
			Entry<StateValues, StochasticUpdateStrategy> SvS = checkMQ((SMG) model, new_params, generateStrategy);
			mainLog.print(String.format("%s took %f s\n", generateStrategy ? "Synthesis" : "Verification", ((double) (System.nanoTime() - t0)) / 1e9));

			parsed_params = params; // register parameters
			result.setStrategy(SvS.getValue()); // register strategy
			mainLog.print(String.format("strategy: %s\n", result.getStrategy()));

			// optional: simulate strategy
			//testStrategy_QEST((SMG) model, strategy, params);
			mainLog.flush();
			return SvS.getKey(); // return state values
		}
	}

	/**
	* Convert the ratio rewards in the expression to average rewards:
	* every ratio reward E[r/c] >= w is converted to an average reward A[r - wc] >=0.
	* This is used for strategy computation, but does not work for Pareto set computations,
	* as the parameterisation on w is no longer explicit.
	* Note: this operation is destructive.
	*
	* @param smg The stochastic game.
	* @param params The parameters specifying the MQ
	*
	* @return The converted MQ
	**/
	private MultiParameters convertRatioMQToMQ(SMG smg, MultiParameters params) throws PrismException
	{
		int gameSize = smg.getNumStates();
		int n = params.rewards.size();

		MultiParameters new_params = new MultiParameters();
		new_params.shallow_copy(params);

		// for every ratio objective, convert it to an average objective, assuming the bounds are fixed
		for (int i = 0; i < n; i++) {
			if (new_params.reward_types.get(i) == MultiParameters.ERCR
			    || new_params.reward_types.get(i) == MultiParameters.PRCR) {
				// build new reward structure: r - lambda c
				double lambda = new_params.bounds.get(i);
				SMGRewardsSimple rlc = new SMGRewardsSimple(gameSize);
				SMGRewards r = params.rewards.get(i);
				SMGRewards c = params.divisors.get(i);
				for (int s = 0; s < gameSize; s++) {
					rlc.setStateReward(s, r.getStateReward(s) - lambda * c.getStateReward(s));
					for (int t = 0; t < smg.getNumChoices(s); t++)
						rlc.setTransitionReward(s, t, r.getTransitionReward(s, t) - lambda * c.getTransitionReward(s, t));
				}
				new_params.rewards.set(i, rlc);
				new_params.reward_types.set(i, MultiParameters.EAR);
				new_params.divisors.set(i, null); // no longer needed
				new_params.shifts.set(i, getShiftFromReward(smg, rlc)); // recalculate here!
				new_params.bounds.set(i, 0.0);
			}
		}

		setObjectiveType(new_params); // reset objective type

		return new_params;
	}

	private void interiorPareto(Pareto[] pareto)
	{
		if (pareto == null)
			return;

		for (int s = 0; s < pareto.length; s++) {
			pareto[s].interior(varepsilon);
		}
	}

	/**
	* Given a collection of Pareto sets, one for each state,
	* together with the parameters for an MQ, turns around the directions appropriately,
	* i.e. minimised directions are turned, maximisations are left alone.
	* Modifies the given Pareto set and parameters.
	**/
	public static void turnParetoAndBounds(Pareto[] pareto, MultiParameters params)
	{
		int dim = params.rewards.size();

		// if Pareto set empty, don't do anything
		if (pareto == null)
			return;

		// minimised rewards need to be turned around
		List<Variable> negative_dimensions = new ArrayList<Variable>();
		for (int i = 0; i < dim; i++) {
			if (params.directions.get(i) < 0) { // if direction is minimised, turn it around now
				negative_dimensions.add(new Variable(i));
				params.bounds.set(i, -params.bounds.get(i));
			}
		}

		// now turn the directions
		for (int s = 0; s < pareto.length; s++)
			pareto[s].turn(negative_dimensions);
	}

	/**
	* Computes brackets around the rewards in each direction - the biggest and the smallest.
	* Used to initialise the Pareto set computation for total rewards,
	* and in any case if rounding is enabled.
	*
	* the results are directly set in the parameters params
	**/
	private void setRewardBrackets(MultiParameters params, Model model) throws PrismException
	{
		int tmpverbosity = this.getVerbosity();
		boolean genStrat = generateStrategy;

		int gameSize = model.getNumStates();
		int dim = params.rewards.size(); // full dimensionality of the problem

		// set baseline bracket
		double[] baseline_biggest_reward = params.rounding ? new double[dim] : null;
		double[] baseline_smallest_reward = params.rounding ? new double[dim] : null;
		params.baseline_biggest_reward = baseline_biggest_reward;
		params.baseline_smallest_reward = baseline_smallest_reward;
		params.no_union_with_previous = false;

		double[][] MIN = new double[dim][gameSize];
		params.MIN = MIN;

		int i = 0;
		try {
			generateStrategy = false; // turn off strategy generation for STPGModelChecker
			this.setVerbosity(0); // temporarily turn off logger

			for (i = 0; i < dim; i++) {
				SMGRewards reward = params.rewards.get(i);
				switch (params.reward_types.get(i)) {
				case MultiParameters.ETCR: // expected total cumulative reward
					// note: computeReachRewardsValIter doesn't check for infinity
					double[] maxminreward = null;
					double[] minminreward = null;
					// check if rewards are all non-negative
					boolean hasNegative = false;
					search_for_negative: for (int s = 0; s < gameSize; s++) {
						// state rewards
						if (reward.getStateReward(s) < 0.0) {
							hasNegative = true;
							break search_for_negative;
						}
						// transition rewards
						for (int d = 0; d < ((SMG) model).getNumChoices(s); d++) {
							if (reward.getTransitionReward(s, d) < 0.0) {
								hasNegative = true;
								break search_for_negative;
							}
						}
					}
					if (hasNegative) {
						try {
							minminreward = createSTPGModelChecker().computeReachRewardsValIter((STPG) model, reward, new BitSet(gameSize), new BitSet(gameSize), true, true,
									null, null).soln;
						} catch (PrismException e) {
							// value iteration did not converge
							// if computing Pareto set, can still approximate,
							// so start with zero for this dimension,
							// but disable union with previous
							if (computePareto) {
								minminreward = new double[gameSize];
								for (int s = 0; s < gameSize; s++)
									minminreward[s] = 0.0;
								params.no_union_with_previous = true;
								mainLog.printWarning("Could not initialise value iteration, because the reward for objective " + i
										+ " does not converge. Pareto set computation started at 0 and safety not guaranteed.");
							} else {
								throw e;
							}
						}
					} else {
						// all rewards non-negative, can set minmin to zero
						minminreward = new double[gameSize];
						for (int s = 0; s < gameSize; s++)
							minminreward[s] = 0.0;
					}
					if (params.rounding) {
					        // defaults
					        double tmp_biggest = 1;
						double tmp_smallest = 0;

						try {
						        maxminreward = createSTPGModelChecker().computeReachRewardsValIter((STPG) model, reward, new BitSet(gameSize), new BitSet(gameSize), false, true, null,
								null).soln;
							tmp_biggest = Double.NEGATIVE_INFINITY;
							tmp_smallest = Double.POSITIVE_INFINITY;
							for (int s = 0; s < gameSize; s++) {
							    if (tmp_biggest < maxminreward[s])
								    tmp_biggest = maxminreward[s];
							    if (tmp_smallest > minminreward[s])
								    tmp_smallest = minminreward[s];
							}
						} catch (PrismException e) {
						    // value iteration did not converge
						    // retain default values
						    // and issue a warning
						    tmp_biggest = 1;
						    tmp_smallest = 0;
						    mainLog.printWarning("Could not initialise value iteration, because the reward for objective " + i
									 + " does not converge. Rounding not sensitive to maximal rewards in the dimension " + i + ".");
						}

						baseline_biggest_reward[i] = tmp_biggest;
						baseline_smallest_reward[i] = tmp_smallest;
					}
					System.arraycopy(minminreward, 0, MIN[i], 0, gameSize);
					break;
				case MultiParameters.EAR: // expected average reward
				case MultiParameters.PAR: // probability (satisfaction) average reward
				case MultiParameters.ERCR: // expected ratio cumulative reward
				case MultiParameters.PRCR: // probability (satisfaction) ratio reward
					// default value for average- (and, by extension, ratio-) rewards
					if (params.rounding) {
					        baseline_biggest_reward[i] = 100;
						baseline_smallest_reward[i] = 0;
					}
					break;
				default:
					throw new PrismException("Reward type not supported.");
				}
			}
		} catch (PrismException e) {
			// if minmin or maxmin reward diverges, treat set of solutions as empty
			throw new PrismException("Could not initialise value iteration, because the reward for objective " + i + " does not converge. Pareto set empty");
		} finally { // do "cleanup" by switching logger on again, and set strategy generation to previous setting
			generateStrategy = genStrat;
			this.setVerbosity(tmpverbosity);
		}

		// fill in biggest reward
		params.biggest_reward = params.rounding ? new double[params.CONJUNCTS] : null;
		if (params.rounding)
			for (i = 0; i < params.CONJUNCTS; i++)
				params.biggest_reward[i] = Math.max(Math.abs(params.baseline_biggest_reward[i]), Math.abs(params.baseline_smallest_reward[i]));

	}

	/**
	 * Resolves StateRewardsConstant in reward structure.
	 * 
	 * @param model The SMG.
	 * @param rewStruct The reward structure to be resolved.
	 * @return
	 * @throws PrismException
	 */
	private SMGRewardsSimple constructSMGRewards(SMG model, RewardStruct rewStruct) throws PrismException
	{
		if (rewStruct == null)
			return null;

		int gameSize = model.numStates;
		ConstructRewards constructRewards = new ConstructRewards(this);
		SMGRewards smgreward = constructRewards.buildSMGRewardStructure(model, rewStruct, constantValues);
		SMGRewardsSimple reward = null;
		if (smgreward instanceof SMGRewardsSimple) {
			reward = (SMGRewardsSimple) smgreward;
		} else if (smgreward instanceof StateRewardsConstant) {
			// construct reward structure manually
			reward = new SMGRewardsSimple(gameSize);
			for (int s = 0; s < gameSize; s++) {
				reward.setStateReward(s, ((StateRewardsConstant) smgreward).getStateReward(s));
			}
		}
		return reward;
	}

	private void parseExpressionQuant(Model model, ExpressionQuant expr_i, MultiParameters params) throws PrismException
	{
		// First check if this is an R(path) nested in a P>=1 (rather than just an R)
		ExpressionReward exprReward = null;
		OpRelOpBound opInfo = expr_i.getRelopBoundInfo(constantValues);
		boolean nested;
		if (expr_i instanceof ExpressionProb) {
			if (opInfo.getBound() != 1.0) {
				throw new PrismLangException("Expression cannot occur in SMG multi-objective query", expr_i);
			}
			if (!((expr_i.getExpression() instanceof ExpressionReward)) && ((ExpressionReward) expr_i.getExpression()).getModifier().equals("path")) {
				throw new PrismLangException("Expression cannot occur in SMG multi-objective query", expr_i);
			}
			exprReward = (ExpressionReward) expr_i.getExpression();
			opInfo = exprReward.getRelopBoundInfo(constantValues);
			nested = true;

		} else if (expr_i instanceof ExpressionReward) {
			exprReward = (ExpressionReward) expr_i;
			nested = false;
		} else {
			throw new PrismLangException("Expression cannot occur in SMG multi-objective query", expr_i);
		}
		params.nested.add(nested);
		
		int gameSize = ((SMG) model).numStates;

		RelOp relOp = exprReward.getRelOp(); // direction and strictness of operator
		boolean minimize = false;
		if (relOp == RelOp.LT || relOp == RelOp.GT) {
			strictToNonstrict.add(exprReward.toString());
		} else if (relOp != RelOp.GEQ && relOp != RelOp.LEQ && relOp != RelOp.MIN && relOp != RelOp.MAX) {
			throw new PrismException("Only minimization or maximization supported.");
		}
		if (relOp == RelOp.LT || relOp == RelOp.LEQ || relOp == RelOp.MIN) {
			minimize = true;
			params.directions.add(-1);
		} else {
			params.directions.add(1);
		}

		// Get reward structures from expression
		
		
		RewardStruct reward_struct = modulesFile.getRewardStruct(exprReward.getRewardStructIndexByIndexObject(modulesFile.getRewardStructNames(), constantValues));
		SMGRewardsSimple reward = constructSMGRewards((SMG) model, reward_struct);
		params.reward_names.add(reward_struct.getName());
		int divisor_struct_index = exprReward.getRewardStructDivIndexByIndexObject(modulesFile.getRewardStructNames(), constantValues);
		RewardStruct divisor_struct = divisor_struct_index == -1 ? null : modulesFile.getRewardStruct(divisor_struct_index);
		SMGRewardsSimple divisor = constructSMGRewards((SMG) model, divisor_struct);
		params.divisor_names.add(divisor_struct == null ? null : divisor_struct.getName());
		// register reward structures
		params.rewards.add(reward);
		params.divisors.add(divisor);

		// check type of reward
		Expression e = exprReward.getExpression();
		if (e instanceof ExpressionTemporal) {
			switch (((ExpressionTemporal) e).getOperator()) {
			case ExpressionTemporal.R_C: // expected cumulative reward
				if (divisor == null)
				    if(nested)
					throw new PrismException("Cannot nest total reward");
				    else
					params.reward_types.add(MultiParameters.ETCR); // expected total
				else
				    throw new PrismException("Transient ratio rewards not supported");
				break;
			case ExpressionTemporal.R_S: // expected average reward
				if (divisor == null)
				    if(nested)
				        params.reward_types.add(MultiParameters.PAR); // satisfaction mean
				    else
				        params.reward_types.add(MultiParameters.EAR); // expected mean
				else
				    if(nested)
					params.reward_types.add(MultiParameters.PRCR); // satisfaction ratio
				    else
					params.reward_types.add(MultiParameters.ERCR); // ratio expectation
				break;
			default:
				throw new PrismException("Only total, average and ratio rewards supported so far.");
			}
		} else {
			throw new PrismException("Only temporal expressions supported in rewards so far.");
		}

		// for minimization need to revert sign of reward (only the numerator, not the divisor):
		if (minimize) {
			for (int s = 0; s < gameSize; s++) {
				reward.setStateReward(s, -reward.getStateReward(s));
				for (int c = 0; c < ((NondetModel) model).getNumChoices(s); c++) {
					reward.setTransitionReward(s, c, -reward.getTransitionReward(s, c));
				}
			}
		}

		// if it is an long-run average reward, need to evaluate a shift as well
		// note: always non-positive: this is what needs to be added after completion
		if (((ExpressionTemporal) e).getOperator() == ExpressionTemporal.R_S && divisor==null)
			params.shifts.add(getShiftFromReward(model, reward));
		else
			params.shifts.add(0.0);

		// evaluate reward bound
		Expression rb = exprReward.getReward();
		if (!(rb instanceof ExpressionVar || rb instanceof ExpressionConstant) && computePareto) {
			unfolded.add(exprReward.toString());
		}
		double r = 0.0; // reward bound
		if (rb != null) {
			r = rb.evaluateDouble(constantValues);
		} else if (relOp != RelOp.MIN && relOp != RelOp.MAX) {
			throw new PrismException("Reward bound required");
		}
		params.bounds.add(minimize ? -r : r); // add bound to vector
	}

	private double getShiftFromReward(Model model, SMGRewards reward)
	{
		int gameSize = ((SMG) model).getNumStates();

		double shift = 0.0;
		for (int s = 0; s < gameSize; s++) {
			// state
			double rew = reward.getStateReward(s);
			if (rew < 0 && shift > rew)
				shift = rew;
			// transitions
			for (int c = 0; c < ((NondetModel) model).getNumChoices(s); c++) {
				rew = reward.getTransitionReward(s, c);
				if (rew < 0 && shift > rew)
					shift = rew;
			}
		}
		return shift;
	}

	protected StochasticUpdateStrategy constructStrategy(SMG G, List<Double> bounds, Pareto[] X, List<Pareto>[] Y, MultiParameters params,
			boolean energy_objective) throws PrismException
	{
		double[] d_bounds = new double[bounds.size()];
		for (int i = 0; i < bounds.size(); i++)
			d_bounds[i] = bounds.get(i);
		return new StochasticUpdateStrategy(G, d_bounds, X, Y, params.rewards, params.biggest_reward, params.baseline_accuracy, true, params.rounding,
				energy_objective ? varepsilon : 0.0, logStrategy, mainLog);
	}

	@Override
	protected StateValues checkProbPathFormulaCosafeLTL(Model model, Expression expr, boolean qual, MinMax minMax, BitSet statesOfInterest) throws PrismException
	{
		// Temporarily make SMG into an STPG by setting coalition and do computation on STPG
		SMG smg = (SMG) model;
		smg.setCoalition(minMax.getCoalition());
		StateValues probs = createSTPGModelChecker().checkProbPathFormulaCosafeLTL(model, expr, qual, minMax, statesOfInterest);
		smg.setCoalition(null);
		return probs;
	}
	
	/**
	 * Compute rewards for a co-safe LTL reward operator.
	 */
	protected StateValues checkRewardCoSafeLTL(Model model, Rewards modelRewards, Expression expr, MinMax minMax, BitSet statesOfInterest) throws PrismException
	{
		// Build product of SMG and DFA for the LTL formula, convert rewards and do any required exports
		LTLModelChecker mcLtl = new LTLModelChecker(this);
		LTLModelChecker.LTLProduct<SMG> product = mcLtl.constructDFAProductForCosafetyReward(this, (SMG) model, expr, statesOfInterest);
		SMGRewards productRewards = ((SMGRewards) modelRewards).liftFromModel(product);
		doProductExports(product);

		// Find accepting states + compute reachability rewards
		BitSet acc = ((AcceptanceReach)product.getAcceptance()).getGoalStates();
		mainLog.println("\nComputing reachability rewards...");
		SMGModelChecker mcProduct = new SMGModelChecker(this);
		mcProduct.inheritSettings(this);
		ModelCheckerResult res = mcProduct.computeReachRewards(product.getProductModel(), productRewards, acc, STPGModelChecker.R_INFINITY, minMax.isMin1(), minMax.isMin2(), minMax.getCoalition());
		StateValues rewardsProduct = StateValues.createFromDoubleArrayResult(res, product.getProductModel());
		
		// Mapping rewards in the original model
		StateValues rewards = product.projectToOriginalModel(rewardsProduct);
		rewardsProduct.clear();
		
		return rewards;
	}
	
	/**
	* Convert the multi-objective query (MQ) from CNF into a CQ
	* as suggested in MFCS'13, given that we have finite rewards in all dimensions.
	* The MQ may contain ratio objectives, which are translated to
	* ratio objectives in the resulting CQ.
	*
	* @param x The weight vector.
	* @param smg The game.
	* @param params The parameters specifying the current MQ (in CNF).
	*
	* @return The parameters specifying the transformed CQ.
	**/
	private MultiParameters convertRatioMQToRatioCQ(double[][] x, SMG smg, MultiParameters params) throws PrismException
	{
		int gameSize = smg.getNumStates();

		MultiParameters new_params = new MultiParameters(); // the return value
		new_params.shallow_copy(params); // shallow copy of parameters, sufficient here, as relevant things are changed now ...

		// compute a new reward and new bound for each conjunct
		List<SMGRewards> new_rewards = new ArrayList<SMGRewards>();
		List<SMGRewards> new_divisors = new ArrayList<SMGRewards>();
		List<Integer> new_types = new ArrayList<Integer>();
		List<Double> new_shifts = new ArrayList<Double>();
		List<Double> new_bounds = new ArrayList<Double>();
		List<Integer> new_ji = new ArrayList<Integer>();
		// set up parameters of new CQ
		new_params.rewards = new_rewards;
		new_params.divisors = new_divisors;
		new_params.reward_types = new_types;
		new_params.bounds = new_bounds;
		new_params.shifts = new_shifts;
		params.ji = new_ji; // sic. needed as sort of return value

		// if rounding enabled, have to recompute biggest_reward and smallest_reward
		double[] new_baseline_biggest_reward = params.rounding ? new double[x.length] : null;
		double[] new_baseline_smallest_reward = params.rounding ? new double[x.length] : null;
		new_params.baseline_biggest_reward = new_baseline_biggest_reward;
		new_params.baseline_smallest_reward = new_baseline_smallest_reward;
		// the biggest reward then is the maximum bound in either direction
		double[] new_biggest_reward = params.rounding ? new double[x.length] : null;
		new_params.biggest_reward = new_biggest_reward;

		// also recompute MIN, for each state - to initialise sets
		double[][] new_MIN = new double[x.length][gameSize];
		new_params.MIN = new_MIN;

		int count = 0; // index of reward (flat, i.e. in the list of all rewards, not in the CNF form)
		int tmpcount = 0; // used to store index during loops so that we start at the same index
		for (int i = 0; i < params.CONJUNCTS; i++) { // iterate through conjuncts

			// check type
			int type_i = -1;
			count = tmpcount;
			for (int j = 0; j < params.DISJUNCTS[i]; j++) { // iterate through disjuncts
				if (type_i < 0)
					type_i = params.reward_types.get(count);
				else if (type_i != params.reward_types.get(count))
					throw new PrismException("Cannot mix reward types within a disjunction");
				count++;
			}
			new_types.add(type_i);

			// calculate bounds, shifts and rewards
			if (type_i == MultiParameters.ERCR || type_i == MultiParameters.PRCR) { // ratio rewards
				// assert: only arrive here if Pareto set computation is requested
				/**
				 * Ratio rewards: /\_i \/_j E[r_ij/c_ij] >= v_ij
				 * 1. transformed to average rewards: /\_i \/_j A[r_ij - v_ij c_ij] >= 0
				 * 2. weighted by x: if j(i) is such that x_ij(i) != 0, then have (@ is componentwise product)
				 *      /\_i A[x_i r_i - ((1-e_j(i)) @ x_i) c_i - k_i x_ij(i) c_ij(i)] >= 0
				 * 3. denote: rho_i = x_i r_i - ((1-e_j(i)) @ x_i) c_i
				 *            gamma_i = x_ij(i) c_ij(i)
				 *    and rewrite as 
				 *      /\_i A[rho_i - k_i gamma_i] >= 0
				 * 4. compute Pareto set, the backwards transformation is done in transformCQParetoToMQPareto
				 **/

				// I. find ji, s.t. x[i][ji] != 0 
				int ji = -1;
				find_ji: for (int j = 0; j < params.DISJUNCTS[i]; j++) {
					if (Math.abs(x[i][j]) >= PrismUtils.epsilonDouble) {
						ji = j;
						break find_ji;
					}
				}
				new_ji.add(ji);
				if (ji < 0)
					throw new PrismException("Weight vector for conjunct " + i + " is zero");

				// II: compute rho_i and gamma_i
				SMGRewardsSimple numerator = new SMGRewardsSimple(gameSize);
				SMGRewardsSimple divisor = new SMGRewardsSimple(gameSize);
				for (int s = 0; s < gameSize; s++) {
					// compute new reward for state
					double rho_i = 0.0;
					count = tmpcount;
					SMGRewards c_iji = null; // c_ij(i)
					for (int j = 0; j < params.DISJUNCTS[i]; j++) { // iterate through disjuncts
						SMGRewards r_ij = params.rewards.get(count);
						SMGRewards c_ij = params.divisors.get(count);
						// x_i r_i
						rho_i += r_ij.getStateReward(s) * x[i][j];
						// - ((1-e_j(i)) @ x_i) c_i
						if (j != ji)
							rho_i -= c_ij.getStateReward(s) * x[i][j];
						else
							// j == ji
							c_iji = c_ij;
						count++;
					}
					numerator.setStateReward(s, rho_i);
					divisor.setStateReward(s, x[i][ji] * c_iji.getStateReward(s)); // gamma_i = x_ij(i) c_ij(i)

					// compute new reward for transitions
					for (int c = 0; c < smg.getNumChoices(s); c++) {
						rho_i = 0.0;
						count = tmpcount;
						c_iji = null; // c_ij(i)
						for (int j = 0; j < params.DISJUNCTS[i]; j++) { // iterate through disjuncts
							SMGRewards r_ij = params.rewards.get(count);
							SMGRewards c_ij = params.divisors.get(count);
							// x_i r_i
							rho_i += r_ij.getTransitionReward(s, c) * x[i][j];
							// - ((1-e_j(i)) @ x_i) c_i
							if (j != ji)
								rho_i -= c_ij.getTransitionReward(s, c) * x[i][j];
							else
								// j == ji
								c_iji = c_ij;
							count++;
						}
						numerator.setTransitionReward(s, c, rho_i);
						divisor.setTransitionReward(s, c, x[i][ji] * c_iji.getTransitionReward(s, c)); // gamma_i = x_ij(i) c_ij(i)
					}
				}

				// III: now transform to /\_i A[rho_i - k_i gamma_i] >= 0
				new_bounds.add(0.0);
				new_rewards.add(numerator);
				new_divisors.add(divisor);

				new_shifts.add(0.0); // is recomputed later!

				// MIN stays all zeros
			} else {
				/**
				 * Total cumulative or Average rewards: /\_i \/_j R[r_ij] >= v_ij
				 * 1. Weighted by x: 
				 *      /\_i R[x_i r_i] >= x_i v_i
				 * 2. Compute Pareto set, the backwards transformation is done in transformCQParetoToMQPareto
				 **/
				// the bound
				double bound = 0.0;
				count = tmpcount;
				for (int j = 0; j < params.DISJUNCTS[i]; j++) {
					bound += params.bounds.get(count) * x[i][j];
					count++;
				}
				new_bounds.add(bound);
				// the reward
				count = tmpcount;
				SMGRewardsSimple new_reward = new SMGRewardsSimple(gameSize);
				for (int s = 0; s < gameSize; s++) {
					// compute new reward for state
					double reward = 0.0;
					count = tmpcount;
					for (int j = 0; j < params.DISJUNCTS[i]; j++) { // iterate through disjuncts
						reward += params.rewards.get(count).getStateReward(s) * x[i][j];
						count++;
					}
					new_reward.setStateReward(s, reward);

					// compute new reward for transitions
					for (int c = 0; c < smg.getNumChoices(s); c++) {
						reward = 0.0;
						count = tmpcount;
						for (int j = 0; j < params.DISJUNCTS[i]; j++) {
							reward += params.rewards.get(count).getTransitionReward(s, c) * x[i][j];
							count++;
						}
						new_reward.setTransitionReward(s, c, reward);
					}

				}
				new_rewards.add(new_reward);
				// the shift
				if (type_i == MultiParameters.EAR)
					new_shifts.add(getShiftFromReward(smg, new_reward));
				else
					new_shifts.add(0.0);

				// recompute MIN
				for (int s = 0; s < gameSize; s++) {
					count = tmpcount;
					for (int j = 0; j < x[i].length; j++) {
						new_MIN[i][s] += params.MIN[count][s] * x[i][j];
						count++;
					}
				}

				// dummy ji, not used
				new_ji.add(0);
			}

			// if rounding, recompute accuracy parameters
			if (params.rounding) {
				count = tmpcount;
				for (int j = 0; j < x[i].length; j++) {
					new_baseline_biggest_reward[i] += params.baseline_biggest_reward[count] * x[i][j];
					new_baseline_smallest_reward[i] += params.baseline_smallest_reward[count] * x[i][j];
					count++;
				}
				// the biggest reward then is the maximum bound in either direction
				new_biggest_reward[i] = Math.max(Math.abs(new_baseline_biggest_reward[i]), Math.abs(new_baseline_smallest_reward[i]));
			}

			tmpcount = count;
		}

		new_params.DISJUNCTS = new int[new_params.CONJUNCTS];
		Arrays.fill(new_params.DISJUNCTS, 1); // only one "disjunct" per conjunct now

		return new_params;
	}

	private Pareto[] convertCQParetoToMQPareto(double[][] x, SMG smg, Pareto[] Px, MultiParameters params) throws PrismException
	{
		int gameSize = smg.getNumStates();
		int dim = params.rewards.size();
		long denominator = (long) (1.0 / PrismUtils.epsilonDouble);

		Pareto[] result = new Pareto[gameSize]; // stores the result to be returned
		// for average and cumulative rewards:
		//     U_{x \in X} cvx( U_{v \in Px} intersect_{i \le n} {p | p x <= v_i})
		// for ratio rewards:
		//     U_{x \in X} cvx( U_{v \in Px} intersect_{i \le n} {p | p x <= [1, 1, ... v_i, ..., 1] x_i})
		// note: sufficient to take vertices of Px
		// add the hyperplane for each state
		for (int s = 0; s < gameSize; s++) {
			for (Generator vertex : Px[s].get().generators()) {
				if (vertex.type() != Generator_Type.POINT)
					continue; // only operate on points, not rays

				// construct a polyhedron for each vertex, and take the union of such polyhedra
				Polyhedron poly = null;

				// get coefficients vertex_i, to make up the right hand side of the constraint
				Map<Integer, BigInteger> vertexi = PPLSupport.getCoefficients(vertex.linear_expression());
				int count = 0;
				for (int i = 0; i < x.length; i++) { // for each conjunct
					// left hand side of the constraint
					Linear_Expression lnumj = null;
					for (int j = 0; j < x[i].length; j++) {
						BigInteger numij = BigInteger.valueOf((long) (x[i][j] * ((double) denominator)));
						if (lnumj != null)
							lnumj = new Linear_Expression_Sum(lnumj, new Linear_Expression_Times(new Coefficient(numij), new Variable(count)));
						else
							lnumj = new Linear_Expression_Times(new Coefficient(numij), new Variable(count));
						count++;
					}
					Linear_Expression lhs = new Linear_Expression_Times(lnumj, vertex.divisor());

					// right hand side of the constraint
					Linear_Expression rhs = null;
					BigInteger v_i = BigInteger.ZERO;
					if (vertexi.size() != 0 && vertexi.get(i) != null)
						v_i = vertexi.get(i);

					if (params.reward_types.get(i) == MultiParameters.ERCR
					    || params.reward_types.get(i) == MultiParameters.PRCR) { // ratio reward
						BigInteger k_i = BigInteger.ZERO;
						for (int j = 0; j < x[i].length; j++) {
							double xij = x[i][j] * ((double) denominator);
							if (params.ji.get(i) == j) {
								k_i = k_i.add(v_i.multiply(BigInteger.valueOf((long) xij)));
							} else {
								k_i = k_i.add(BigInteger.valueOf((long) xij).multiply(vertex.divisor().getBigInteger()));
							}
						}
						rhs = new Linear_Expression_Coefficient(new Coefficient(k_i));
					} else { // expected average or total cumulative reward
						rhs = new Linear_Expression_Coefficient(new Coefficient(v_i.multiply(BigInteger.valueOf(denominator))));
					}
					// form new constraint: lhs <= rhs
					Constraint con = new Constraint(lhs, Relation_Symbol.LESS_OR_EQUAL, rhs);

					// add constraint to new poly
					if (poly != null) {
						poly.add_constraint(con);
					} else {
						Constraint_System cs = new Constraint_System();
						cs.add(con);
						poly = new C_Polyhedron(cs);
						// add full number of dimensions to poly
						poly.add_space_dimensions_and_embed(dim - poly.space_dimension());
					}
				}
				if (poly == null)
					throw new PrismException("Error in computing polyhedra from CQ evaluation");

				if (result[s] != null)
					result[s].get().upper_bound_assign(poly); // union with previous polys
				else
					result[s] = new Pareto(poly);
			}
		}

		return result;
	}

	/**
	* The set for the weighted CQ is computed as well for player and stochastic states,
	* which can be used to then construct the strategy.
	*
	* @param x The weight vector
	* @param smg The game
	* @param params Parameters for the computation
	* @param Px Reference to return CQ value for states
	* @param stochasticStates Reference to return CQ value for moves (i.e. stochastic states)
	* @param checkBounds Abort computation when bounds are met (bounds are in params)
	* @param energy_objective Is the objective an energy objective?
	* @param cq_params Reference to parameters for the transformed CQ
	*
	* @return Whether CQ Pareto set computation converged
	**/
	private boolean iterateMQParetoSet(double[][] x, SMG smg, MultiParameters params, Pareto[] Px, List<Pareto>[] stochasticStates, boolean checkBounds,
			boolean energy_objective, MultiParameters cq_params) throws PrismException
	{
		// transform MQ (which is assumed to be in CNF) to CQ
		cq_params.shallow_copy(convertRatioMQToRatioCQ(x, smg, params));

		// CQ Pareto set computation - returns Pareto sets in Px, to be used below to get the MQ pareto sets
		boolean converged = computeCQParetoSet(smg, cq_params, Px, stochasticStates, checkBounds, energy_objective);

		// log
		if (logDPareto)
			PPLSupport.printReachabilityPolyhedron(convertCQParetoToMQPareto(x, smg, Px, params), params.rewards.size(), smg.getFirstInitialState(), mainLog);

		// return whether CQ Pareto set computation converged
		return converged;
	}

	/**
	* Computes the Pareto set of an MQ for a specific choice of weight vectors,
	* where the MQ may contain ratio objectives.
	* Shifting is delegated to {@code computeRatioCQParetoSet}.
	*
	* @param x The weight vector
	* @param smg The game
	* @param params Parameters for the computation
	*
	* @return A convex Pareto set for each state for the MQ.
	**/
	private Pareto[] iterateRatioMQParetoSet(double[][] x, SMG smg, MultiParameters params) throws PrismException
	{
		// CQ Pareto set computation - returns Pareto sets in Px, to be used below to get the MQ pareto sets
		Pareto[] Px = computeRatioCQParetoSet(smg, convertRatioMQToRatioCQ(x, smg, params));

		// log
		if (logRPareto)
			PPLSupport.printReachabilityPolyhedron(Px, params.CONJUNCTS, smg.getFirstInitialState(), mainLog);

		// compute Pareto sets for MQ and return
		return convertCQParetoToMQPareto(x, smg, Px, params);
	}

	/**
	* Shifts the rewards sets in {@code params.rewards} by the amount provided
	* in {@code params.shifts}.
	*
	* @param revert_direction Revert the direction of shifting.
	**/
	private void shiftRewards(Model model, MultiParameters params, boolean revert_direction)
	{
		// default is to use the shifts specified in params, i.e. no strategy construction for energy objective
		shiftRewards(model, params, false, revert_direction);
	}

	/**
	* Shifts the rewards sets in {@code params.rewards} by the amount provided
	* in {@code params.shifts} if {@code construct_energy_strategy} is false,
	* otherwise, if true, the shift is made by {@code params.bounds}.
	*
	* @param revert_direction Revert the direction of shifting.
	**/
	private void shiftRewards(Model model, MultiParameters params, boolean revert_direction, boolean construct_energy_strategy)
	{
		int gameSize = ((SMG) model).numStates;

		// apply shift to rewards
		// if not constructing strategy for an energy objective, use the standard shifts
		// if constructing a strategy, use the bounds as shifts
		Iterator<Double> shifts_iterator = params.shifts.iterator();
		int i = 0;
		for (SMGRewards reward_i : params.rewards) {
			SMGRewardsSimple reward = (SMGRewardsSimple) reward_i;
			Double shift = construct_energy_strategy ? params.bounds.get(i) : shifts_iterator.next();
			for (int s = 0; s < gameSize; s++) {
				reward.setStateReward(s, revert_direction ? reward.getStateReward(s) + shift : reward.getStateReward(s) - shift);
				for (int c = 0; c < ((NondetModel) model).getNumChoices(s); c++) {
					reward.setTransitionReward(s, c, revert_direction ? reward.getTransitionReward(s, c) + shift : reward.getTransitionReward(s, c) - shift);
				}
			}
			i++;
		}
	}

	/**
	* Shifts the bounds sets in {@code params.bounds} by the amount provided
	* in {@code params.shifts} if {@code construct_energy_strategy} is false,
	* otherwise, if true, no shift is made.
	*
	* @param revert_direction Revert the direction of shifting.
	**/
	private void shiftBounds(Model model, MultiParameters params, boolean revert_direction, boolean construct_energy_strategy)
	{
		// if constructing a strategy for an energy objective, don't shift at all
		if (construct_energy_strategy)
			return;
		// apply shift to bounds
		Iterator<Double> shifts_iterator = params.shifts.iterator();
		for (int i = 0; i < params.bounds.size(); i++) {
			Double shift_i = shifts_iterator.next();
			Double bound_i = params.bounds.get(i);
			params.bounds.set(i, revert_direction ? bound_i + shift_i : bound_i - shift_i);
		}
	}

	/**
	* Shifts the Pareto sets in {@code Px} by the amount provided
	* in {@code params.shifts}.
	**/
	private void shiftPareto(SMG smg, MultiParameters params, Pareto[] Px) throws PrismException
	{
		int gameSize = smg.getNumStates();
		int n = params.rewards.size();
		double[] extra_rewards = new double[n];
		for (int i = 0; i < n; i++)
			extra_rewards[i] = params.shifts.get(i);
		for (int s = 0; s < gameSize; s++)
			Px[s].replace(0, PPLSupport.add_rewards(Px[s].get(), s, Integer.MIN_VALUE, null, extra_rewards));
	}

	/**
	 * Checks if the query is a conjunction, and if so returns true.
	 **/
	private boolean isConjunction(MultiParameters params)
	{
		for (int i = 0; i < params.CONJUNCTS; i++)
			if (params.DISJUNCTS[i] != 1)
				return false;
		return true;
	}

	/**
	* Computes the Pareto set of an MQ in CNF that potentially contains ratio objectives.
	* The rewards can be negative as well, but shifting is delegated to the subroutines.
	*
	* @param smg The stochastic game.
	* @param params The parameters specifying the MQ.
	*
	* @return The Pareto set, a union of convex sets.
	**/
	public Pareto[] computeRatioMQParetoSet(SMG smg, MultiParameters params) throws PrismException
	{
		int gameSize = smg.getNumStates();

		// sanity check for conjunction valiter count
		if (params.maxCIter < 1)
			throw new PrismException("Iteration count for conjunctions has to be greater or equal to one.");

		if (isConjunction(params)) { // IF CONJUNCTION
			return computeRatioCQParetoSet(smg, params);

		} else { // IF MQ

			// initialise polyhedron lists
			Pareto[] result = new Pareto[gameSize]; // stores the result to be returned - a union of convex sets
			for (int s = 0; s < gameSize; s++)
				result[s] = new Pareto(params.maxDIter);

			// sanity check for disjunction valiter count
			if (params.dIterOffset < 1)
				throw new PrismException("Iteration offset for disjunctions has to be greater or equal to one");
			if (params.maxDIter < 1)
				throw new PrismException("Iteration count for disjunctions has to be greater or equal to one.");

			// selector for hyperplanes / weight vectors
			HyperplaneSelector hyperplaneSelector = new HyperplaneSelector(params.CONJUNCTS, params.DISJUNCTS);

			// ITERATE THROUGH HYPERPLANES
			for (int iter = 1; iter < params.maxDIter + params.dIterOffset; iter++) {
				// get new choice of hyperplane
				double[][] x = hyperplaneSelector.next_hyperplane();
				if (iter < params.dIterOffset)
					continue;
				if (logDPareto)
					mainLog.print(String.format("D-ITER (%d/%d): x=%s\texclude=%s\n", iter, params.maxDIter + params.dIterOffset, Arrays.deepToString(x),
							Arrays.toString(hyperplaneSelector.getExclude())));

				// evalute Pareto set for this choice of x (i.e. the hyperplanes)
				Pareto[] Ux = iterateRatioMQParetoSet(x, smg, params);

				// log
				if (logDPareto)
					PPLSupport.printReachabilityPolyhedron(Ux, params.rewards.size(), smg.getFirstInitialState(), mainLog);

				// put result in list of convex sets per state
				for (int s = 0; s < gameSize; s++)
					result[s].add(Ux[s]);
			}

			return result;
		}
	}

	/**
	 * if the objective consists of only one type of objectives,
	* set objective type of parameters to this type, otherwise to -1;
	**/
	private void setObjectiveType(MultiParameters params) throws PrismException
	{
		int objective_type = -1;
		int count = 0;
		for (int i = 0; i < params.CONJUNCTS; i++) {
			for (int j = 0; j < params.DISJUNCTS[i]; j++) {
				if (objective_type < 0) {
					objective_type = params.reward_types.get(count);
				} else if (objective_type != params.reward_types.get(count)) {
					params.objective_type = -1;
					return;
				}
				count++;
			}
		}
		params.objective_type = objective_type;
	}

	/**
	 * Checks MQ and optionally also constructs a winning strategy if requested using {@code construct_strategy}.
	 *
	 * @param smg The game.
	 * @param params  The parameters specifying the MQ.
	 * @param construct_strategy Whether to construct a strategy.
	 * 
	 * @return State values with bit set of states satisfying the objective,
	 * and a winning strategy for initial state if objective is satisfied
	 **/
	public Entry<StateValues, StochasticUpdateStrategy> checkMQ(SMG smg, MultiParameters params, boolean construct_strategy) throws PrismException
	{
		int gameSize = smg.getNumStates();
		int initialState = smg.getFirstInitialState();

		// return values
		StateValues sv = null;
		StochasticUpdateStrategy strategy = null;

		// for strategy construction, check if all objectives are of the same type
		if (construct_strategy && params.objective_type < 0)
			throw new PrismException(
					"Strategy construction only available if all objectives are either expected cumulative total rewards, or expected cumulative ratio and expected average rewards");
		// strategy construction for average rewards is via energy objectives (cf. Chatterjee et al. FSTTCS'10)
		boolean energy_objective = (params.objective_type == MultiParameters.EAR || params.objective_type == MultiParameters.PAR);

		// sets for player and stochastic states and bounds for the (potentially tranformed) CQ
		Pareto[] Px = new Pareto[gameSize]; // sets at states
		List<Pareto>[] stochasticStates = construct_strategy ? new List[gameSize] : null;

		List<Double> cq_bounds = null; // bounds used to construct the strategy and test the objective
		MultiParameters cq_params = null; // parameters used to construct the strategy

		// only check bounds during iterating if the objective is puerly expected total cumulative rewards
		boolean checkBounds = (params.objective_type == MultiParameters.ETCR);

		// sanity check for conjunction valiter count
		if (params.maxCIter < 1)
			throw new PrismException("Iteration count for conjunctions has to be greater or equal to one.");

		boolean satisfied = false; // stores if objective satisfied for BOX_M iterations

		shiftRewards(smg, params, false, energy_objective); // apply shifts to rewards
		shiftBounds(smg, params, false, energy_objective); // apply shifts to bounds
		try {
			if (isConjunction(params)) { // IF CONJUNCTION
				// OUTER LOOP IN TACAS'15 - INCREASE M
				params.M = minM;
				increaseM: do { // terminate outer loop below if not energy objective in any case
					boolean converged = computeCQParetoSet(smg, params, Px, stochasticStates, checkBounds, energy_objective);

					// ENERGY: increase M if initial state is still empty or not converged
					if (energy_objective && (Px[initialState].get().is_empty() ||
								 !converged)) {
						params.M = params.M * params.M;
						mainLog.print(String.format("Increasing M to %d\n", params.M));
					} else if (energy_objective && converged) {
					        satisfied = true;
						break increaseM; // converged 
					} else if (!energy_objective) {
						break increaseM; // if not an energy objective, terminate outer loop in any case
					}
				} while (params.M <= maxM);
				cq_params = params;
				cq_bounds = energy_objective ? initialCreditVector(Px[initialState]) : cq_params.bounds;
			} else { // IF MQ

				// sanity check for disjunction valiter count
				if (params.dIterOffset < 1)
					throw new PrismException("Iteration offset for disjunctions has to be greater or equal to zero");
				if (params.maxDIter < 1)
					throw new PrismException("Iteration count for disjunctions has to be greater or equal to one.");

				// OUTER LOOP IN TACAS'15 - INCREASE M
				params.M = minM;
				increaseM: do { // terminate outer loop below if not energy objective in any case
					// selector for hyperplanes / weight vectors
					HyperplaneSelector hyperplaneSelector = new HyperplaneSelector(params.CONJUNCTS, params.DISJUNCTS);
					cq_params = new MultiParameters();

					// ITERATE THROUGH HYPERPLANES
					iterate_disj: for (int iter = 1; iter < params.maxDIter + params.dIterOffset; iter++) {
						// get new choice of hyperplane
						double[][] x = hyperplaneSelector.next_hyperplane();
						if (iter < params.dIterOffset)
							continue;
						if (logDPareto)
							mainLog.print(String.format("D-ITER (%d/%d): x=%s\n", iter, params.maxDIter + params.dIterOffset, Arrays.deepToString(x)));

						// evaluate Pareto set for this choice of x (i.e. the hyperplanes)
						boolean converged = iterateMQParetoSet(x, smg, params, Px, stochasticStates, checkBounds, energy_objective, cq_params);
						cq_bounds = energy_objective ? initialCreditVector(Px[initialState]) : cq_params.bounds;
						if (logDPareto)
							PPLSupport.printReachabilityPolyhedron(Px, params.CONJUNCTS, initialState, mainLog);

						if (energy_objective && !converged) {
							continue iterate_disj;
						}

						if (logDPareto) mainLog.print(String.format("checking bounds %s on %d\n", cq_bounds, initialState));
						// check if point satisfies the objective
						if (PPLSupport.checkBound(Px[initialState], cq_bounds, cq_params)) {
							satisfied = true;
							break iterate_disj;
						}
					}

					// ENERGY: increase M if initial state is still empty
					if (energy_objective && !satisfied) {
						params.M = params.M * params.M;
						mainLog.print(String.format("Increasing M to %d\n", params.M));
					} else if (energy_objective && satisfied) {
						break increaseM; // converged 
					} else if (!energy_objective) {
						break increaseM; // if not an energy objective, terminate outer loop in any case
					}

				} while (params.M <= maxM);
			}

			// vector of whether states satisfy the objective
			if (cq_bounds == null || (energy_objective && !satisfied))
				sv = StateValues.createFromBitSet(new BitSet(gameSize), smg);
			else
				sv = StateValues.createFromBitSet(PPLSupport.checkBounds(Px, cq_bounds, cq_params), smg);

			// construct strategy if requested and objective met
			if (construct_strategy && sv.getBitSet().get(initialState)) {
				double t0 = System.nanoTime();
				strategy = constructStrategy(smg, cq_bounds, Px, stochasticStates, cq_params, energy_objective);
				mainLog.print(String.format("Strategy construction took %f s\n", ((double) (System.nanoTime() - t0)) / 1e9));
				strategy.setInfo(params.getParameterString());
			}

		} finally {
			shiftRewards(smg, params, true, energy_objective); // shift rewards back
			shiftBounds(smg, params, true, energy_objective); // shift bounds back
		}

		return new SimpleEntry<StateValues, StochasticUpdateStrategy>(sv, strategy);
	}

	/**
	* Returns an arbitrary initial credit vector from the Pareto set for energy objectives.
	* If set is empty, return null.
	**/
	private List<Double> initialCreditVector(Pareto Pxt) throws PrismException
	{
	        List<Generator> gens = Pxt.get().generators();
	    
		for(int gi = gens.size()-1; gi >= 0; gi--) {
		        Generator g = gens.get(gi);
			if (g.type() == Generator_Type.POINT) {
				List<Double> result = PPLSupport.getGeneratorAsVector(g, (int) Pxt.get().space_dimension());
				for (int i = 0; i < result.size(); i++)
					result.set(i, result.get(i) - varepsilon); // can be lenient since already converged
				return result;
			}
		}
		return null;
	}

	/**
	 * Initialises the Pareto sets for CQ queries
	 * so that at state s, the Pareto set will contain
	 * the downwards closure of {@code MIN[][s]}, which is an n-dimensional vector.
	 **/
	private Pareto[] initialiseCQParetoSet(int gameSize, int n, double[][] MIN) throws PrismException
	{
		Pareto[] Qx = new Pareto[gameSize]; // the return value

		for (int s = 0; s < gameSize; s++) {
			Generator_System gs = new Generator_System();

			double[] x = new double[n];
			for (int i = 0; i < n; i++) {
				x[i] = MIN[i][s];
				// generate ray for downward closure
				Linear_Expression ray = new Linear_Expression_Times(new Coefficient((BigInteger.ONE).negate()), new Variable(i));
				gs.add(Generator.ray(ray));
			}
			gs.add(PPLSupport.generatorFromPoint(x));

			C_Polyhedron cp = new C_Polyhedron(gs); // generate initial polyhedra: X^0_s
			Qx[s] = new Pareto(cp);
		}

		return Qx;
	}

	/**
	 * If all rewards are positive, the Pareto set can be initialised
	 * by taking one generator per dimension, and taking the convex hull.
	 *
	 * Note: currently not used in direct method.
	 **/
	public Pareto[] initialiseCQParetoSetPositive(int gameSize, int n, double[][] maxmin)
	{
		Pareto[] Qx = new Pareto[gameSize];
		for (int s = 0; s < gameSize; s++) {
			Generator_System gs = new Generator_System();
			// step through rewards and add to r
			for (int i = 0; i < n; i++) {
				BigFraction ri = new BigFraction(maxmin[i][s]);
				BigInteger num = ri.getNumerator();
				BigInteger den = ri.getDenominator();
				Linear_Expression r_num = new Linear_Expression_Times(new Coefficient(num), new Variable(i));
				gs.add(Generator.point(r_num, new Coefficient(den)));
				// generate ray for downward closure
				Linear_Expression ray = new Linear_Expression_Times(new Coefficient((BigInteger.ONE).negate()), new Variable(i));
				gs.add(Generator.ray(ray));
			}
			// generate initial polyhedra: X^0_s
			C_Polyhedron cp = new C_Polyhedron(gs);
			Qx[s] = new Pareto(cp);
		}
		return Qx;
	}

	/**
	* Increases the baseline accuracy by the increase factor specified in the Prism options.
	**/
	private long increaseBaselineAccuracy(long baseline_accuracy)
	{
		long tmp_ba;
		if (increase_factor > 1.0 && (long) (((double) baseline_accuracy) * increase_factor) <= baseline_accuracy) {
			tmp_ba = baseline_accuracy + 1;
		} else {
			tmp_ba = (long) (((double) baseline_accuracy) * increase_factor);
		}
		if (tmp_ba < max_accuracy && tmp_ba > 0) { // sanity check to prevent overflow
			return tmp_ba;
		} else {
			return max_accuracy;
		}
	}

	/**
	* Convert CQ with ratio rewards to CQ without ratio rewards:
	* for every ratio reward E[r/c] >= w that you find,
	* A[xr - yc] >= v
	* and then for every vector v >= 0, take w = y/x.
	**/
	private MultiParameters convertRatioCQToCQ(SMG smg, double[][] x, MultiParameters params) throws PrismException
	{
		int gameSize = smg.getNumStates();

		MultiParameters new_params = new MultiParameters();
		new_params.shallow_copy(params);

		List<SMGRewards> new_rewards = new ArrayList<SMGRewards>(params.CONJUNCTS);
		List<SMGRewards> new_divisors = new ArrayList<SMGRewards>(params.CONJUNCTS);
		List<Integer> new_types = new ArrayList<Integer>(params.CONJUNCTS);
		List<Double> new_shifts = new ArrayList<Double>(params.CONJUNCTS);
		List<Double> new_bounds = new ArrayList<Double>(params.CONJUNCTS);
		// set up parameters of new CQ
		new_params.rewards = new_rewards;
		new_params.divisors = new_divisors;
		new_params.reward_types = new_types;
		new_params.bounds = new_bounds;
		new_params.shifts = new_shifts;

		int ratio_count = 0;
		for (int i = 0; i < params.CONJUNCTS; i++) {
			switch (params.reward_types.get(i)) {
			case MultiParameters.ETCR: // expected total cumulative reward
			case MultiParameters.EAR: // expected average reward		    
				// just copy over
				new_rewards.add(params.rewards.get(i));
				new_divisors.add(null); // no divisors needed
				new_types.add(params.reward_types.get(i));
				new_bounds.add(params.bounds.get(i));
				new_shifts.add(params.shifts.get(i));
				break;
			case MultiParameters.ERCR: // ratio rewards
			case MultiParameters.PRCR:
				// build new reward structure rmc: xr - yc
				SMGRewardsSimple xrmyc = new SMGRewardsSimple(gameSize);
				SMGRewards r = params.rewards.get(i);
				SMGRewards c = params.divisors.get(i);
				for (int s = 0; s < gameSize; s++) {
					xrmyc.setStateReward(s, x[ratio_count][0] * r.getStateReward(s) - x[ratio_count][1] * c.getStateReward(s));
					for (int t = 0; t < smg.getNumChoices(s); t++) {
						xrmyc.setTransitionReward(s, t, x[ratio_count][0] * r.getTransitionReward(s, t) - x[ratio_count][1] * c.getTransitionReward(s, t));
					}
				}
				new_rewards.add(xrmyc);
				new_divisors.add(null); // no divisors needed
				new_types.add(MultiParameters.EAR); // convert to average reward
				new_bounds.add(0.0);
				new_shifts.add(getShiftFromReward(smg, xrmyc));
				ratio_count++;
				break;
			default:
				throw new PrismException("Reward type not supported.");
			}
		}

		setObjectiveType(new_params);

		return new_params;
	}

	/**
	* Computes the Pareto set of a CQ which potentially includes ratio objectives.
	* Takes care of reward shifting if they are negative in some dimensions (note that {@code computeCQParetoSet} doesn't shift,
	* but assumes rewards are non-negative.)
	* 
	* @param smg The game.
	* @param params The CQ with potentially some ratio objectives.
	**/
	private Pareto[] computeRatioCQParetoSet(SMG smg, MultiParameters params) throws PrismException
	{
		int gameSize = smg.getNumStates();
		int initial_state = smg.getFirstInitialState();
		Pareto[] result = new Pareto[gameSize];

		// evaluate number of ratio rewards
		int ratio_rewards = 0;
		for (Integer reward_type : params.reward_types)
			if (reward_type == MultiParameters.ERCR || reward_type == MultiParameters.PRCR)
				ratio_rewards++;

		// if no ratio rewards present, just call standard routine
		if (ratio_rewards == 0) {
			// apply shifts to rewards
			shiftRewards(smg, params, false);
			try {
				computeCQParetoSet(smg, params, result, null, false, false);
				shiftPareto(smg, params, result); // apply shift back to Pareto set
				return result;
			} finally {
				shiftRewards(smg, params, true); // shift rewards back
			}
		}

		// if there are ratio rewards, need to build up the Pareto set through another value iteration
		// initialise downward-closed sets
		double min = -100.0; // TODO: an issue for disjunctions?
		double[] mins = new double[ratio_rewards];
		Arrays.fill(mins, min);
		Generator_System[] gs = new Generator_System[gameSize];
		for (int s = 0; s < gameSize; s++) {
			gs[s] = new Generator_System();
			for (int i = 0; i < params.CONJUNCTS; i++)
				gs[s].add(Generator.ray(new Linear_Expression_Times(new Coefficient((BigInteger.ONE).negate()), new Variable(i))));
			gs[s].add(PPLSupport.generatorFromPoint(mins));
		}

		// keep track of which vectors are guaranteed to be good and which are guaranteed to be bad
		Generator_System gl = new Generator_System();
		for (int i = 0; i < ratio_rewards; i++)
			gl.add(Generator.ray(new Linear_Expression_Times(new Coefficient((BigInteger.ONE).negate()), new Variable(i))));
		gl.add(PPLSupport.generatorFromPoint(mins));
		Polyhedron good_lambdas = new C_Polyhedron(gl);
		Set<double[]> bad_lambdas = new HashSet<double[]>();

		// sanity check for ratio valiter count
		if (params.maxRIter < 1)
			throw new PrismException("Iteration count for ratios has to be greater or equal to one.");

		// use the hyperplane selection to get values for lambda
		int[] DISJUNCTS = new int[ratio_rewards];
		Arrays.fill(DISJUNCTS, 2); // two for each ratio reward
		double[][] x = new double[ratio_rewards][2];
		HyperplaneSelector hyperplaneSelectorPos = new HyperplaneSelector(ratio_rewards, DISJUNCTS);
		iterate_through_hyperplanes: for (int iter = 0; iter < params.maxRIter; iter++) {
			double[][] y = hyperplaneSelectorPos.next_hyperplane();
			double start = hyperplaneSelectorPos.getStart();
			// translate and rotate
			for (int i = 0; i < ratio_rewards; i++) {
				x[i][0] = y[i][0] + y[i][1] - (1 - start);
				x[i][1] = y[i][0] - y[i][1];
			}

			// check if this weight selection is well defined in the first place
			int i = 0;
			double[] lambda = new double[ratio_rewards];
			for (Integer reward_type : params.reward_types) {
				if (reward_type == MultiParameters.ERCR || reward_type == MultiParameters.PRCR) {
					if (Math.abs(x[i][0]) < PrismUtils.epsilonDouble) {
						iter--; // nothing computed yet
						continue iterate_through_hyperplanes;
					}
					lambda[i] = x[i][1] / x[i][0];
					i++;
				}
			}

			// test if weight selection already covered and found good enough
			Generator_System tgs = new Generator_System();
			tgs.add(PPLSupport.generatorFromPoint(lambda));
			Polyhedron tp = new C_Polyhedron(tgs);
			if (good_lambdas.contains(tp)) {
				continue iterate_through_hyperplanes;
			}
			// test if weight selection already covered and found too bad
			look_for_bad_lambda: for (double[] bad_lambda : bad_lambdas) {
				for (i = 0; i < ratio_rewards; i++) {
					if (lambda[i] < bad_lambda[i]) // this lambda is less in at least one coordinate than an already bad lambda
						continue look_for_bad_lambda;
				}
				// fall through if lambda is larger in all coordinates than some already bad lambda
				continue iterate_through_hyperplanes;
			}

			// if all ok with the weights, convert the ratio objectives to average rewards
			MultiParameters new_params = convertRatioCQToCQ(smg, x, params);

			// apply shifts to rewards
			shiftRewards(smg, new_params, false);
			try {
				Pareto[] Qx = new Pareto[gameSize];
				computeCQParetoSet(smg, new_params, Qx, null, false, false);
				shiftPareto(smg, new_params, Qx); // apply shift back to Pareto set
				if(logRPareto) {
				        mainLog.print(String.format("QX[%d/%d](%s):", iter, params.maxRIter, Arrays.toString(lambda)));
					PPLSupport.printReachabilityPolyhedron(Qx, params.CONJUNCTS, initial_state, mainLog);
				}
				for (int s = 0; s < gameSize; s++) {
					boolean good_vertex_found = false;
					go_through_vertices: for (Generator vertex : Qx[s].get().generators()) {
						if (vertex.type() != Generator_Type.POINT)
							continue;
						Map<Integer, BigInteger> vertexi = PPLSupport.getCoefficients(vertex.linear_expression());
						int ratio_count = 0;
						i = 0;
						double[] vector = new double[params.CONJUNCTS]; // vector of values to be added to computed Pareto set
						for (Integer reward_type : params.reward_types) {
							double vali = vertexi.get(i) == null ? 0.0 : vertexi.get(i).doubleValue();
							if (reward_type == MultiParameters.ERCR || reward_type == MultiParameters.PRCR) { // ratio reward
								if (vali < 0.0) // do not add point, lambda too large
									continue go_through_vertices;
								else
									// lambda still good so far
									vector[i] = lambda[ratio_count];
								ratio_count++;
							} else { // average or cumulative total reward
								vector[i] = vali / vertex.divisor().getBigInteger().doubleValue();
							}
							i++;
						}
						if (s == initial_state)
							good_vertex_found = true;
						gs[s].add(PPLSupport.generatorFromPoint(vector));
					}
					if (s == initial_state && good_vertex_found) {
						good_lambdas.add_generator(PPLSupport.generatorFromPoint(lambda));
						if(logRPareto)
						        mainLog.print(String.format("R-ITER[%d/%d]: lambda=%s\n", iter + 1, params.maxRIter, Arrays.toString(lambda)));
					} else if (s == initial_state)
						bad_lambdas.add(lambda);
				}
			} finally {
				shiftRewards(smg, new_params, true); // shift rewards back
			}
		}

		// fill in results
		for (int s = 0; s < gameSize; s++)
			result[s] = new Pareto(new C_Polyhedron(gs[s]));

		return result;
	}

	/**
	 * Compute CQ Pareto sets for cumulative total and average rewards.
	 * Result is returned in {@code Px} for states and in {@code stochasticStates} for moves.
	 * The iteration is stopped on hitting {@code params.bounds} if {@code checkBounds} is true.
	 * The assumption is that all rewards are non-negative if {@code energy_objective} is false.
	 *
	 * @param smg The stochastic game
	 * @param params The parameters of the query
	 * @param Px Pareto sets for player states, used as return value
	 * @param stochasticStates Pareto sets for stochastic states, used as return value
	 * @params checkBounds Whether the computation is stopped once the bounds are met
	 * @params energy_objective Whether the objective is evaluated as an energy objective (for synthesis)
	 *
	 * @return Whether the value iteration converged.
	 */
	public boolean computeCQParetoSet(SMG smg, MultiParameters params, Pareto[] Px, List<Pareto>[] stochasticStates, boolean checkBounds,
			boolean energy_objective) throws PrismException
	{
		int gameSize = smg.getNumStates();
		int n = params.rewards.size();
		int init = smg.getFirstInitialState();

		// first check if we can use the single-objective engine to speed things up
		if (n == 1 && !energy_objective && !checkBounds) {
			// Create solution vector(s)
			double[] soln = new double[gameSize];
			double[] soln2 = new double[gameSize];
			double[] tmpsoln;

			BitSet unknown = new BitSet();
			unknown.set(0, gameSize);

			// Start iterations
			for (int iters = 0; iters < params.maxCIter; iters++) {
				// Matrix-vector multiply and min/max ops
				((STPG) smg).mvMultRewMinMax(soln, params.rewards.get(0), false, true, soln2, unknown, false, null, 1.0);

				// Swap vectors for next iter
				tmpsoln = soln;
				soln = soln2;
				soln2 = tmpsoln;
			}

			if (params.objective_type == MultiParameters.EAR) // evaluate average
				for (int s = 0; s < gameSize; s++)
					soln[s] /= ((double) (2 * params.maxCIter));

			// build Pareto set
			double[][] maxmin = new double[1][];
			maxmin[0] = soln;
			Pareto[] temp = initialiseCQParetoSet(gameSize, 1, maxmin);

			// copy to result
			System.arraycopy(temp, 0, Px, 0, temp.length);

			return true; // if not energy objective, will speak of convergence in any case
		}

		// only allow Gauss-Seidel when all dimensions are total cumulative rewards,
		// or if we have an energy objective
		boolean localGaussSeidel = gaussSeidel && ((params.objective_type == MultiParameters.ETCR) || energy_objective);

		// INITIALISATION: compute polyhedra X_s^0
		Pareto[] Qx = initialiseCQParetoSet(gameSize, n, params.MIN);

		// set up arrays for average reward (needed to check bounds and convergence and later rescale the sets)
		int[] step = new int[n];
		Arrays.fill(step, 1); // default is 1
		double[] bounds = new double[n];
		double[] base_bounds = new double[n];
		// apply shift to bounds - note: this is not used for energy objectives
		for (int i = 0; i < n; i++)
			base_bounds[i] = params.bounds.get(i) - params.shifts.get(i);

		// ITERATE FUNCTIONAL APPLICATION: compute X_s^k+1 = F(X_s^k), cf. MFCS'13 / TACAS'15
		boolean converged = false;
		long baseline_accuracy = params.baseline_accuracy;
		iterate_cq: for (int k = 0; k < params.maxCIter; k++) {
			if (logCPareto)
			        mainLog.print(String.format("C-ITER %d/%s, %s", k + 1, params.maxCIter,
						params.rounding ? String.format("acc = %d, ", baseline_accuracy) : ""));
			mainLog.flush();
			// set up factors for average reward
			for (int i = 0; i < n; i++) {
				if (!energy_objective && params.reward_types.get(i) == MultiParameters.EAR) {
					// take step count times two, because every iteration the functional is applied twice!
					step[i] = (k + 1) * 2;
					bounds[i] = base_bounds[i] * ((double) ((k + 1) * 2));
				} else {
					bounds[i] = base_bounds[i];
				}
			}

			// VALUE ITERATION STEP
			Pareto[] temp = smg.pMultiObjective(Qx, params.rewards, localGaussSeidel, baseline_accuracy, params.biggest_reward,
					stochasticStates, params.rounding, !params.no_union_with_previous & !energy_objective, energy_objective, params.M);
			System.arraycopy(temp, 0, Px, 0, temp.length); // copy to result

			if (logCPareto)
			    PPLSupport.printReachabilityPolyhedron(Px, params.CONJUNCTS, init, mainLog);
			    //PPLSupport.printReachabilityPolyhedra(Px, stochasticStates, params.CONJUNCTS, mainLog);

			// test varepsilon-convergence
			if (convergeNorm(Px, Qx, n, step, energy_objective, init)) {
				if (logCPareto)
					mainLog.print("CQ value iteration converged.\n");
				converged = true;
				break iterate_cq; // if converged, break cq iteration
			}

			// test if target met
			if (checkBounds && !energy_objective && PPLSupport.checkBound(Px[init], bounds, params))
				break iterate_cq; // if target met, break cq iteration 

			// increase accuracy
			baseline_accuracy = increaseBaselineAccuracy(baseline_accuracy);

			// keep current as previous Pareto (for convergence check)
			System.arraycopy(Px, 0, Qx, 0, Px.length);
		}

		// MEAN/TOTAL/RATIO: rescale if required by average reward if not energy objective
		if (!energy_objective) {
			double[] alpha = new double[n]; // scaling factor
			Arrays.fill(alpha, 1.0); // default is 1.0
			for (int i = 0; i < n; i++)
				if (params.reward_types.get(i) == MultiParameters.EAR)
					alpha[i] = 1.0 / ((double) step[i]);
			PPLSupport.discountPareto(Px, alpha);
			PPLSupport.discountPareto(stochasticStates, alpha);

		}
		// finished - Pareto sets now in Px and stochasticStates

		// return whether converged
		return converged;
	}

	/**
	 * Tests convergence using epsilon-growth criterion (relative!),
	 * that is, test whether (prev \cap current) + epsilon \supseteq (prev \cup current).
	 * If the sets are monotonically increasing, this reduces to whether prev + epsilon \supseteq current.
	 * Note that stochastic states are ignored here.
	 *
	 * arguments:
	 * @param result Current Pareto sets
	 * @param prev_result Previous Pareto sets
	 * @param n Dimension
	 * @param k Current step number (if positive, used for average reward)
	 * @param energy_objective If dealing with an energy objective
	 * @param init Index of initial state
	 *
	 * @return Whether all (or initial if energy objective) sets have converged.
	 **/
	private boolean convergeNorm(final Pareto[] result, final Pareto[] prev_result, int n, int[] k, boolean energy_objective, int init) throws PrismException
	{
		for (int s = 0; s < result.length; s++) {
			Polyhedron ck1 = new C_Polyhedron((C_Polyhedron) result[s].get()); // deep copy - current
			Polyhedron ck = new C_Polyhedron((C_Polyhedron) prev_result[s].get()); // deep copy - previous

			// add step-discount if required
			for (int i = 0; i < n; i++) {
				if (k[i] > 1) {
					Variable var = new Variable(i);
					Linear_Expression expr = new Linear_Expression_Times(new Coefficient(1), var);
					Coefficient den1 = new Coefficient(k[i]);
					Coefficient den = new Coefficient(k[i] - 1);
					ck1.affine_image(var, expr, den1);
					ck.affine_image(var, expr, den);
				}
			}

			Polyhedron ck_prime;

			// if not monotonically increasing anyway, aply union and intersection
			if (!ck1.contains(ck)) {
				ck_prime = new C_Polyhedron(ck.generators()); // deep copy
				ck.intersection_assign(ck1); // the set that's supposed to be smaller holds the intersection
				ck1.upper_bound_assign(ck_prime); // the set that's supposed to be larger holds the union
			}

			Generator_System ngs = new Generator_System();
			// first set up the reward vector that should be added to each point generator

			BigFraction r = new BigFraction(varepsilon);
			BigInteger num = r.getNumerator();
			BigInteger den = r.getDenominator();

			// prepare vector pointing in direction (varepsilon, varepsilon, ...)
			Linear_Expression le = new Linear_Expression_Times(new Coefficient(num), new Variable(0));
			Coefficient c = new Coefficient(den);
			for (int i = 1; i < n; i++) {
				le = new Linear_Expression_Sum(le, new Linear_Expression_Times(new Coefficient(num), new Variable(i)));
			}

			// now add reward vector to each point generator 
			for (Generator g : ck.generators()) {
				if (g.type() == Generator_Type.POINT) {
					Linear_Expression nle = new Linear_Expression_Sum(le.times(g.divisor()), g.linear_expression().times(c));
					Coefficient nc = new Coefficient(g.divisor().getBigInteger().multiply(c.getBigInteger()));
					ngs.add(Generator.point(nle, nc));
				} else {
					ngs.add(g);
				}
			}
			ck_prime = new C_Polyhedron(ngs);
			// now test containment
			if (ck1.is_empty())
				continue; // converged for this state
			if (ck_prime.is_empty() && !ck1.is_empty())
				return false; // not converged yet
			if (!ck_prime.is_empty() && !ck_prime.contains(ck1))
				return false; // not converged yet
		}
		return true; // only fall through if all polyhedra converge
	}
	
	// Numerical computation functions
	
	/**
	 * Compute next-state probabilities.
	 * i.e. compute the probability of being in a state in {@code target} in the next step.
	 * @param smg The SMG
	 * @param target Target states
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 * @param coalition The coalition of players which define player 1
	 */
	public ModelCheckerResult computeNextProbs(SMG smg, BitSet target, boolean min1, boolean min2, Coalition coalition) throws PrismException
	{
		// Temporarily make SMG into an STPG by setting coalition and do computation on STPG
		smg.setCoalition(coalition);
		ModelCheckerResult res = createSTPGModelChecker().computeNextProbs(smg, target, min1, min2);
		smg.setCoalition(null);
		return res;
	}

	/**
	 * Compute bounded reachability/until probabilities.
	 * i.e. compute the min/max probability of reaching a state in {@code target},
	 * within k steps, and while remaining in states in @{code remain}.
	 * @param smg The SMG
	 * @param remain Remain in these states (optional: null means "all")
	 * @param target Target states
	 * @param k Bound
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 * @param coalition The coalition of players which define player 1
	 */
	public ModelCheckerResult computeBoundedUntilProbs(SMG smg, BitSet remain, BitSet target, int k, boolean min1, boolean min2, Coalition coalition) throws PrismException
	{
		// Temporarily make SMG into an STPG by setting coalition and do computation on STPG
		smg.setCoalition(coalition);
		ModelCheckerResult res = createSTPGModelChecker().computeBoundedUntilProbs(smg, remain, target, k, min1, min2);
		smg.setCoalition(null);
		return res;
	}

	/**
	 * Compute until probabilities.
	 * i.e. compute the min/max probability of reaching a state in {@code target},
	 * while remaining in those in @{code remain}.
	 * @param smg The SMG
	 * @param remain Remain in these states (optional: null means "all")
	 * @param target Target states
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 * @param coalition The coalition of players which define player 1
	 */
	public ModelCheckerResult computeUntilProbs(SMG smg, BitSet remain, BitSet target, boolean min1, boolean min2, Coalition coalition) throws PrismException
	{
		// Temporarily make SMG into an STPG by setting coalition and do computation on STPG
		smg.setCoalition(coalition);
		ModelCheckerResult res = createSTPGModelChecker().computeUntilProbs(smg, remain, target, min1, min2, -1);
		smg.setCoalition(null);
		return res;
	}

	/**
	 * Compute expected reachability rewards, where the runs that don't reach
	 * the final state get infinity. i.e. compute the min/max reward accumulated
	 * to reach a state in {@code target}.
	 * @param smg The SMG
	 * @param rewards The rewards
	 * @param target Target states
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 * @param coalition The coalition of players which define player 1
	 */
	public ModelCheckerResult computeReachRewards(SMG smg, SMGRewards rewards, BitSet target, int unreachingSemantics, boolean min1, boolean min2, Coalition coalition) throws PrismException
	{
		// Temporarily make SMG into an STPG by setting coalition and do computation on STPG
		smg.setCoalition(coalition);
		ModelCheckerResult res = createSTPGModelChecker().computeReachRewards(smg, rewards, target, min1, min2, null, null, unreachingSemantics);
		smg.setCoalition(null);
		return res;
	}
	
	// Utility methods
	
	/**
	 * Create a new STPG model checker with the same settings as this one. 
	 */
	private STPGModelChecker createSTPGModelChecker() throws PrismException
	{
		STPGModelChecker mcSTPG = new STPGModelChecker(this);
		mcSTPG.inheritSettings(this);
		return mcSTPG;
	}
}
