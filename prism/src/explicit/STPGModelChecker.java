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

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import acceptance.AcceptanceReach;
import common.IterableBitSet;
import explicit.rewards.MDPRewardsSimple;
import explicit.rewards.STPGRewards;
import explicit.rewards.STPGRewardsSimple;
import explicit.rewards.StateRewardsConstant;
import parser.ast.Expression;
import prism.AccuracyFactory;
import prism.PrismComponent;
import prism.PrismException;
import prism.PrismFileLog;
import prism.PrismLog;
import prism.PrismNotSupportedException;
import prism.PrismUtils;
import strat.BoundedRewardDeterministicStrategy;
import strat.MemorylessDeterministicStrategy;
import strat.StepBoundedDeterministicStrategy;

/**
 * Explicit-state model checker for two-player stochastic games (STPGs).
 */
public class STPGModelChecker extends ProbModelChecker
{
	/**
	 * Used when calling methods computing reachability rewards. Says that the
	 * runs which don't reach the target get reward infinity.
	 */
	public static final int R_INFINITY = 0;
	/**
	 * Used when calling methods computing reachability rewards. Says that the
	 * runs which don't reach the target get the reward cumulated along the run
	 * (i.e. possibly infinite, possibly finite)
	 */
	public static final int R_CUMULATIVE = 1;
	/**
	 * Used when calling methods computing reachability rewards. Says that the
	 * runs which don't reach the target get reward zero.
	 */
	public static final int R_ZERO = 2;

	/**
	 * Create a new STPGModelChecker, inherit basic state from parent (unless null).
	 */
	public STPGModelChecker(PrismComponent parent) throws PrismException
	{
		super(parent);
	}

	// Model checking functions

	@Override
	protected StateValues checkProbPathFormulaLTL(Model model, Expression expr, boolean qual, MinMax minMax, BitSet statesOfInterest) throws PrismException
	{
		throw new PrismNotSupportedException("Full LTL model checking not yet supported for stochastic games");
	}

	@Override
	protected StateValues checkProbPathFormulaCosafeLTL(Model model, Expression expr, boolean qual, MinMax minMax, BitSet statesOfInterest) throws PrismException
	{
		// Build product of STPG and DFA for the LTL formula, and do any required exports
		LTLModelChecker mcLtl = new LTLModelChecker(this);
		LTLModelChecker.LTLProduct<STPG> product = mcLtl.constructDFAProductForCosafetyProbLTL(this, (STPG) model, expr, statesOfInterest);
		doProductExports(product);
		
		// Find accepting states + compute reachability probabilities
		BitSet acc = ((AcceptanceReach)product.getAcceptance()).getGoalStates();
		mainLog.println("\nComputing reachability probabilities...");
		STPGModelChecker mcProduct = new STPGModelChecker(this);
		mcProduct.inheritSettings(this);
		ModelCheckerResult res = mcProduct.computeReachProbs(product.getProductModel(), acc, minMax.isMin1(), minMax.isMin2());
		StateValues probsProduct = StateValues.createFromDoubleArrayResult(res, product.getProductModel());

		// Output vector over product, if required
		if (getExportProductVector()) {
				mainLog.println("\nExporting product solution vector matrix to file \"" + getExportProductVectorFilename() + "\"...");
				PrismFileLog out = new PrismFileLog(getExportProductVectorFilename());
				probsProduct.print(out, false, false, false, false);
				out.close();
		}

		// Mapping probabilities in the original model
		StateValues probs = product.projectToOriginalModel(probsProduct);
		probsProduct.clear();

		return probs;
	}

	// Numerical computation functions

	/**
	 * Compute next=state probabilities.
	 * i.e. compute the probability of being in a state in {@code target} in the next step.
	 * @param stpg The STPG
	 * @param target Target states
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 */
	public ModelCheckerResult computeNextProbs(STPG stpg, BitSet target, boolean min1, boolean min2) throws PrismException
	{
		ModelCheckerResult res = null;
		int n;
		double soln[], soln2[];
		boolean genAdv = exportAdv || generateStrategy;
		int[] adv = null;
		long timer;

		timer = System.currentTimeMillis();

		// Store num states
		n = stpg.getNumStates();

		// Create/initialise solution vector(s)
		soln = Utils.bitsetToDoubleArray(target, n);
		soln2 = new double[n];

		// Create/initialise adversary storage
		if (genAdv) {
			adv = new int[n];
			for (int i = 0; i < n; i++) {
				adv[i] = -1;
			}
		}

		// Next-step probabilities 
		stpg.mvMultMinMax(soln, min1, min2, soln2, null, false, adv);

		// Store results/strategy
		res = new ModelCheckerResult();
		res.accuracy = AccuracyFactory.boundedNumericalIterations();
		res.soln = soln2;
		res.numIters = 1;
		res.timeTaken = timer / 1000.0;
		if (generateStrategy) {
			res.strat = new MemorylessDeterministicStrategy(adv);
		}

		return res;
	}

	/**
	 * Compute reachability probabilities.
	 * i.e. compute the min/max probability of reaching a state in {@code target}.
	 * @param stpg The STPG
	 * @param target Target states
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 */
	public ModelCheckerResult computeReachProbs(STPG stpg, BitSet target, boolean min1, boolean min2) throws PrismException
	{
		return computeReachProbs(stpg, target, min1, min2, -1);
	}

	/**
	 * Compute reachability probabilities.
	 * i.e. compute the min/max probability of reaching a state in {@code target}.
	 * @param stpg The STPG
	 * @param target Target states
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 */
	public ModelCheckerResult computeReachProbs(STPG stpg, BitSet target, boolean min1, boolean min2, double bound) throws PrismException
	{
		return computeReachProbs(stpg, null, target, min1, min2, null, null, bound);
	}

	/**
	 * Compute until probabilities.
	 * i.e. compute the min/max probability of reaching a state in {@code target},
	 * while remaining in those in {@code remain}.
	 * @param stpg The STPG
	 * @param remain Remain in these states (optional: null means "all")
	 * @param target Target states
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 */
	public ModelCheckerResult computeUntilProbs(STPG stpg, BitSet remain, BitSet target, boolean min1, boolean min2, double bound) throws PrismException
	{
		return computeReachProbs(stpg, remain, target, min1, min2, null, null, bound);
	}

	/**
	 * Compute reachability/until probabilities.
	 * i.e. compute the min/max probability of reaching a state in {@code target},
	 * while remaining in those in {@code remain}.
	 * @param stpg The STPG
	 * @param remain Remain in these states (optional: null means "all")
	 * @param target Target states
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 * @param init Optionally, an initial solution vector (may be overwritten) 
	 * @param known Optionally, a set of states for which the exact answer is known
	 * Note: if 'known' is specified (i.e. is non-null, 'init' must also be given and is used for the exact values.  
	 */
	public ModelCheckerResult computeReachProbs(STPG stpg, BitSet remain, BitSet target, boolean min1, boolean min2, double init[], BitSet known)
			throws PrismException
	{
		return computeReachProbs(stpg, remain, target, min1, min2, init, known, -1);
	}

	/**
	 * Compute reachability/until probabilities.
	 * i.e. compute the min/max probability of reaching a state in {@code target},
	 * while remaining in those in @{code remain}.
	 * @param stpg The STPG
	 * @param remain Remain in these states (optional: null means "all")
	 * @param target Target states
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 * @param init Optionally, an initial solution vector (may be overwritten) 
	 * @param known Optionally, a set of states for which the exact answer is known
	 * Note: if 'known' is specified (i.e. is non-null, 'init' must also be given and is used for the exact values.  
	 */
	public ModelCheckerResult computeReachProbs(STPG stpg, BitSet remain, BitSet target, boolean min1, boolean min2, double init[], BitSet known, double bound)
			throws PrismException
	{
		ModelCheckerResult res = null;
		BitSet no, yes;
		int n, numYes, numNo;
		long timer, timerProb0, timerProb1;
		boolean genAdv;

		// Check for some unsupported combinations
		if (solnMethod == SolnMethod.VALUE_ITERATION && valIterDir == ValIterDir.ABOVE && !(precomp && prob0)) {
			throw new PrismException("Precomputation (Prob0) must be enabled for value iteration from above");
		}

		// Are we generating an optimal adversary?
		genAdv = exportAdv || generateStrategy;

		// Start probabilistic reachability
		timer = System.currentTimeMillis();
		if (verbosity >= 1)
			mainLog.println("\nStarting probabilistic reachability...");

		// Check for deadlocks in non-target state (because breaks e.g. prob1)
		stpg.checkForDeadlocks(target);

		// Store num states
		n = stpg.getNumStates();

		// Optimise by enlarging target set (if more info is available)
		if (init != null && known != null && !known.isEmpty()) {
			BitSet targetNew = (BitSet) target.clone();
			for (int i : new IterableBitSet(known)) {
				if (init[i] == 1.0) {
					targetNew.set(i);
				}
			}
			target = targetNew;
		}

		// Precomputation
		timerProb0 = System.currentTimeMillis();
		if (precomp && prob0) {
			no = prob0(stpg, remain, target, min1, min2);
		} else {
			no = new BitSet();
		}
		timerProb0 = System.currentTimeMillis() - timerProb0;
		timerProb1 = System.currentTimeMillis();
		if (precomp && prob1 && !genAdv) {
			yes = prob1(stpg, remain, target, min1, min2);
		} else {
			yes = (BitSet) target.clone();
		}
		timerProb1 = System.currentTimeMillis() - timerProb1;

		// Print results of precomputation
		numYes = yes.cardinality();
		numNo = no.cardinality();
		if (verbosity >= 1)
			mainLog.println("target=" + target.cardinality() + ", yes=" + numYes + ", no=" + numNo + ", maybe=" + (n - (numYes + numNo)));
		// do value iteration only if the values needed wasn't handled by
		// precomputation
		if (bound < 1.0 || !(precomp && prob1 && !genAdv)) {
			// Compute probabilities
			switch (solnMethod) {
			case VALUE_ITERATION:
				res = computeReachProbsValIter(stpg, no, yes, min1, min2, init, known);
				break;
			case GAUSS_SEIDEL:
				res = computeReachProbsGaussSeidel(stpg, no, yes, min1, min2, init, known);
				break;
			default:
				throw new PrismException("Unknown STPG solution method " + solnMethod);
			}
		} else {
			res = new ModelCheckerResult();
			res.numIters = 0;
			res.soln = Utils.bitsetToDoubleArray(yes, n);
			res.accuracy = AccuracyFactory.doublesFromQualitative();
			mainLog.println("Bound is 1, hence I am skipping the computation of other values than 1.");
		}

		// Finished probabilistic reachability
		timer = System.currentTimeMillis() - timer;
		if (verbosity >= 1)
			mainLog.println("Probabilistic reachability took " + timer / 1000.0 + " seconds.");

		// Update time taken
		res.timeTaken = timer / 1000.0;
		res.timeProb0 = timerProb0 / 1000.0;
		res.timePre = (timerProb0 + timerProb1) / 1000.0;

		return res;
	}
	
	/**
	 * Prob0 precomputation algorithm.
	 * i.e. determine the states of an STPG which, with min/max probability 0,
	 * reach a state in {@code target}, while remaining in those in {@code remain}.
	 * {@code min}=true gives Prob0E, {@code min}=false gives Prob0A. 
	 * @param stpg The STPG
	 * @param remain Remain in these states (optional: null means "all")
	 * @param target Target states
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 */
	public BitSet prob0(STPG stpg, BitSet remain, BitSet target, boolean min1, boolean min2)
	{
		int n, iters;
		BitSet u, soln, unknown;
		boolean u_done;
		long timer;

		// Start precomputation
		timer = System.currentTimeMillis();
		if (verbosity >= 1)
			mainLog.println("Starting Prob0 (" + (min1 ? "min" : "max") + (min2 ? "min" : "max") + ")...");

		// Special case: no target states
		if (target.cardinality() == 0) {
			soln = new BitSet(stpg.getNumStates());
			soln.set(0, stpg.getNumStates());
			return soln;
		}

		// Initialise vectors
		n = stpg.getNumStates();
		u = new BitSet(n);
		soln = new BitSet(n);

		// Determine set of states actually need to perform computation for
		unknown = new BitSet();
		unknown.set(0, n);
		unknown.andNot(target);
		if (remain != null)
			unknown.and(remain);

		// Fixed point loop
		iters = 0;
		u_done = false;
		// Least fixed point - should start from 0 but we optimise by
		// starting from 'target', thus bypassing first iteration
		u.or(target);
		soln.or(target);
		while (!u_done) {
			iters++;
			// Single step of Prob0
			stpg.prob0step(unknown, u, min1, min2, soln);
			// Check termination
			u_done = soln.equals(u);
			// u = soln
			u.clear();
			u.or(soln);
		}

		// Negate
		u.flip(0, n);

		// Finished precomputation
		timer = System.currentTimeMillis() - timer;
		if (verbosity >= 1) {
			mainLog.print("Prob0 (" + (min1 ? "min" : "max") + (min2 ? "min" : "max") + ")");
			mainLog.println(" took " + iters + " iterations and " + timer / 1000.0 + " seconds.");
		}

		return u;
	}

	/**
	 * Prob1 precomputation algorithm.
	 * i.e. determine the states of an STPG which, with min/max probability 1,
	 * reach a state in {@code target}, while remaining in those in {@code remain}.
	 * @param stpg The STPG
	 * @param remain Remain in these states (optional: null means "all")
	 * @param target Target states
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 */
	public BitSet prob1(STPG stpg, BitSet remain, BitSet target, boolean min1, boolean min2)
	{
		int n, iters;
		BitSet u, v, soln, unknown;
		boolean u_done, v_done;
		long timer;

		// Start precomputation
		timer = System.currentTimeMillis();
		if (verbosity >= 1)
			mainLog.println("Starting Prob1 (" + (min1 ? "min" : "max") + (min2 ? "min" : "max") + ")...");

		// Special case: no target states
		if (target.cardinality() == 0) {
			return new BitSet(stpg.getNumStates());
		}

		// Initialise vectors
		n = stpg.getNumStates();
		u = new BitSet(n);
		v = new BitSet(n);
		soln = new BitSet(n);

		// Determine set of states actually need to perform computation for
		unknown = new BitSet();
		unknown.set(0, n);
		unknown.andNot(target);
		if (remain != null)
			unknown.and(remain);

		// Nested fixed point loop
		iters = 0;
		u_done = false;
		// Greatest fixed point
		u.set(0, n);
		while (!u_done) {
			v_done = false;
			// Least fixed point - should start from 0 but we optimise by
			// starting from 'target', thus bypassing first iteration
			v.clear();
			v.or(target);
			soln.clear();
			soln.or(target);
			while (!v_done) {
				iters++;
				// Single step of Prob1
				stpg.prob1step(unknown, u, v, min1, min2, soln);
				// Check termination (inner)
				v_done = soln.equals(v);
				// v = soln
				v.clear();
				v.or(soln);
			}
			// Check termination (outer)
			u_done = v.equals(u);
			// u = v
			u.clear();
			u.or(v);
		}

		// Finished precomputation
		timer = System.currentTimeMillis() - timer;
		if (verbosity >= 1) {
			mainLog.print("Prob1 (" + (min1 ? "min" : "max") + (min2 ? "min" : "max") + ")");
			mainLog.println(" took " + iters + " iterations and " + timer / 1000.0 + " seconds.");
		}

		return u;
	}

	/**
	 * Compute reachability probabilities using value iteration.
	 * @param stpg The STPG
	 * @param no Probability 0 states
	 * @param yes Probability 1 states
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 * @param init Optionally, an initial solution vector (will be overwritten) 
	 * @param known Optionally, a set of states for which the exact answer is known
	 * Note: if 'known' is specified (i.e. is non-null, 'init' must also be given and is used for the exact values.  
	 */
	protected ModelCheckerResult computeReachProbsValIter(STPG stpg, BitSet no, BitSet yes, boolean min1, boolean min2, double init[], BitSet known)
			throws PrismException
	{
		ModelCheckerResult res = null;
		BitSet unknown;
		int i, n, iters;
		double soln[], soln2[], tmpsoln[], initVal;
		int adv[] = null;
		boolean genAdv, done;
		long timer;

		// Are we generating an optimal adversary?
		genAdv = exportAdv || generateStrategy;

		// Start value iteration
		timer = System.currentTimeMillis();
		if (verbosity >= 1)
			mainLog.println("Starting value iteration (" + (min1 ? "min" : "max") + (min2 ? "min" : "max") + ")...");

		// Store num states
		n = stpg.getNumStates();

		// Create solution vector(s)
		soln = new double[n];
		soln2 = (init == null) ? new double[n] : init;

		// Initialise solution vectors. Use (where available) the following in order of preference:
		// (1) exact answer, if already known; (2) 1.0/0.0 if in yes/no; (3) passed in initial value; (4) initVal
		// where initVal is 0.0 or 1.0, depending on whether we converge from below/above. 
		initVal = (valIterDir == ValIterDir.BELOW) ? 0.0 : 1.0;
		if (init != null) {
			if (known != null) {
				for (i = 0; i < n; i++)
					soln[i] = soln2[i] = known.get(i) ? init[i] : yes.get(i) ? 1.0 : no.get(i) ? 0.0 : init[i];
			} else {
				for (i = 0; i < n; i++)
					soln[i] = soln2[i] = yes.get(i) ? 1.0 : no.get(i) ? 0.0 : init[i];
			}
		} else {
			for (i = 0; i < n; i++)
				soln[i] = soln2[i] = yes.get(i) ? 1.0 : no.get(i) ? 0.0 : initVal;
		}

		// Determine set of states actually need to compute values for
		unknown = new BitSet();
		unknown.set(0, n);
		unknown.andNot(yes);
		unknown.andNot(no);
		if (known != null)
			unknown.andNot(known);

		// Create/initialise adversary storage
		if (genAdv) {
			adv = new int[n];
			for (i = 0; i < n; i++) {
				adv[i] = -1;
			}

			int s;
			for (i = 0; i < no.length(); i++) {
				s = no.nextSetBit(i);
				for (int c = 0; c < stpg.getNumChoices(s); c++) {
					if (stpg.allSuccessorsInSet(s, c, no)) {
						adv[i] = c;
						break;
					}
				}
			}
		}

		// Start iterations
		iters = 0;
		done = false;
		while (!done && iters < maxIters) {
			iters++;
			// Matrix-vector multiply and min/max ops
			stpg.mvMultMinMax(soln, min1, min2, soln2, unknown, false, genAdv ? adv : null);
			// Check termination
			done = PrismUtils.doublesAreClose(soln, soln2, termCritParam, termCrit == TermCrit.ABSOLUTE);
			// Swap vectors for next iter
			tmpsoln = soln;
			soln = soln2;
			soln2 = tmpsoln;
		}

		// Finished value iteration
		timer = System.currentTimeMillis() - timer;
		if (verbosity >= 1) {
			mainLog.print("Value iteration (" + (min1 ? "min" : "max") + (min2 ? "min" : "max") + ")");
			mainLog.println(" took " + iters + " iterations and " + timer / 1000.0 + " seconds.");
		}

		// Non-convergence is an error (usually)
		if (!done && errorOnNonConverge) {
			String msg = "Iterative method did not converge within " + iters + " iterations.";
			msg += "\nConsider using a different numerical method or increasing the maximum number of iterations";
			throw new PrismException(msg);
		}

		// Store results/strategy
		res = new ModelCheckerResult();
		res.soln = soln;
		double maxDiff = PrismUtils.measureSupNorm(soln, soln2, termCrit == TermCrit.ABSOLUTE);
		res.accuracy = AccuracyFactory.valueIteration(termCritParam, maxDiff, termCrit == TermCrit.ABSOLUTE);
		res.numIters = iters;
		res.timeTaken = timer / 1000.0;
		if (generateStrategy) {
			res.strat = new MemorylessDeterministicStrategy(adv);
		}

		// Print adversary
		if (genAdv) {
			PrismLog out = new PrismFileLog(exportAdvFilename);
			if (exportAdvFilename.lastIndexOf('.') != -1 && exportAdvFilename.substring(exportAdvFilename.lastIndexOf('.') + 1).equals("dot")) {
				stpg.exportToDotFileWithStrat(out, null, adv);
			} else {
				for (i = 0; i < n; i++) {
					out.println(i + " " + (adv[i] != -1 ? stpg.getAction(i, adv[i]) : "-"));
				}
				out.println();
			}
			out.close();
		}

		return res;
	}

	/**
	 * Compute reachability probabilities using Gauss-Seidel.
	 * @param stpg The STPG
	 * @param no Probability 0 states
	 * @param yes Probability 1 states
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 * @param init Optionally, an initial solution vector (will be overwritten) 
	 * @param known Optionally, a set of states for which the exact answer is known
	 * Note: if 'known' is specified (i.e. is non-null, 'init' must also be given and is used for the exact values.  
	 */
	protected ModelCheckerResult computeReachProbsGaussSeidel(STPG stpg, BitSet no, BitSet yes, boolean min1, boolean min2, double init[], BitSet known)
			throws PrismException
	{
		ModelCheckerResult res;
		BitSet unknown;
		int i, n, iters;
		double soln[], initVal, maxDiff = Double.POSITIVE_INFINITY;
		boolean done;
		long timer;

		// Start value iteration
		timer = System.currentTimeMillis();
		if (verbosity >= 1)
			mainLog.println("Starting value iteration (" + (min1 ? "min" : "max") + (min2 ? "min" : "max") + ")...");

		// Store num states
		n = stpg.getNumStates();

		// Create solution vector
		soln = (init == null) ? new double[n] : init;

		// Initialise solution vector. Use (where available) the following in order of preference:
		// (1) exact answer, if already known; (2) 1.0/0.0 if in yes/no; (3) passed in initial value; (4) initVal
		// where initVal is 0.0 or 1.0, depending on whether we converge from below/above.
		initVal = (valIterDir == ValIterDir.BELOW) ? 0.0 : 1.0;
		if (init != null) {
			if (known != null) {
				for (i = 0; i < n; i++)
					soln[i] = known.get(i) ? init[i] : yes.get(i) ? 1.0 : no.get(i) ? 0.0 : init[i];
			} else {
				for (i = 0; i < n; i++)
					soln[i] = yes.get(i) ? 1.0 : no.get(i) ? 0.0 : init[i];
			}
		} else {
			for (i = 0; i < n; i++)
				soln[i] = yes.get(i) ? 1.0 : no.get(i) ? 0.0 : initVal;
		}

		// Determine set of states actually need to compute values for
		unknown = new BitSet();
		unknown.set(0, n);
		unknown.andNot(yes);
		unknown.andNot(no);
		if (known != null)
			unknown.andNot(known);

		// Start iterations
		iters = 0;
		done = false;
		while (!done && iters < maxIters) {
			iters++;
			// Matrix-vector multiply and min/max ops
			maxDiff = stpg.mvMultGSMinMax(soln, min1, min2, unknown, false, termCrit == TermCrit.ABSOLUTE);
			// Check termination
			done = maxDiff < termCritParam;
		}

		// Finished Gauss-Seidel
		timer = System.currentTimeMillis() - timer;
		if (verbosity >= 1) {
			mainLog.print("Value iteration (" + (min1 ? "min" : "max") + (min2 ? "min" : "max") + ")");
			mainLog.println(" took " + iters + " iterations and " + timer / 1000.0 + " seconds.");
		}

		// Non-convergence is an error (usually)
		if (!done && errorOnNonConverge) {
			String msg = "Iterative method did not converge within " + iters + " iterations.";
			msg += "\nConsider using a different numerical method or increasing the maximum number of iterations";
			throw new PrismException(msg);
		}

		// Return results
		res = new ModelCheckerResult();
		res.soln = soln;
		res.accuracy = AccuracyFactory.valueIteration(termCritParam, maxDiff, termCrit == TermCrit.ABSOLUTE);
		res.numIters = iters;
		res.timeTaken = timer / 1000.0;
		return res;
	}

	/**
	 * Construct strategy information for min/max reachability probabilities.
	 * (More precisely, list of indices of player 1 choices resulting in min/max.)
	 * (Note: indices are guaranteed to be sorted in ascending order.)
	 * @param stpg The STPG
	 * @param state The state to generate strategy info for
	 * @param target The set of target states to reach
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 * @param lastSoln Vector of probabilities from which to recompute in one iteration 
	 */
	public List<Integer> probReachStrategy(STPG stpg, int state, BitSet target, boolean min1, boolean min2, double lastSoln[]) throws PrismException
	{
		double val = stpg.mvMultMinMaxSingle(state, lastSoln, min1, min2);
		return stpg.mvMultMinMaxSingleChoices(state, lastSoln, min1, min2, val);
	}

	/**
	 * Compute bounded reachability probabilities.
	 * i.e. compute the min/max probability of reaching a state in {@code target} within k steps.
	 * @param stpg The STPG
	 * @param target Target states
	 * @param k Bound
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 */
	public ModelCheckerResult computeBoundedReachProbs(STPG stpg, BitSet target, int k, boolean min1, boolean min2) throws PrismException
	{
		return computeBoundedReachProbs(stpg, null, target, k, min1, min2, null, null);
	}

	/**
	 * Compute bounded until probabilities.
	 * i.e. compute the min/max probability of reaching a state in {@code target},
	 * within k steps, and while remaining in states in {@code remain}.
	 * @param stpg The STPG
	 * @param remain Remain in these states (optional: null means "all")
	 * @param target Target states
	 * @param k Bound
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 */
	public ModelCheckerResult computeBoundedUntilProbs(STPG stpg, BitSet remain, BitSet target, int k, boolean min1, boolean min2) throws PrismException
	{
		return computeBoundedReachProbs(stpg, remain, target, k, min1, min2, null, null);
	}

	/**
	 * Compute bounded reachability/until probabilities.
	 * i.e. compute the min/max probability of reaching a state in {@code target},
	 * within k steps, and while remaining in states in {@code remain}.
	 * @param stpg The STPG
	 * @param remain Remain in these states (optional: null means "all")
	 * @param target Target states
	 * @param k Bound
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 * @param init Initial solution vector - pass null for default
	 * @param results Optional array of size k+1 to store (init state) results for each step (null if unused)
	 */
	public ModelCheckerResult computeBoundedReachProbs(STPG stpg, BitSet remain, BitSet target, int k, boolean min1, boolean min2, double init[],
			double results[]) throws PrismException
	{
		// TODO: implement until

		ModelCheckerResult res = null;
		int i, n, iters;
		double soln[], soln2[], tmpsoln[];
		long timer;
		List<List<Integer>> stratChoices = null;
		int[] adv = null;

		// Start bounded probabilistic reachability
		timer = System.currentTimeMillis();
		if (verbosity >= 1)
			mainLog.println("\nStarting bounded probabilistic reachability...");

		// Store num states
		n = stpg.getNumStates();

		// Create solution vector(s)
		soln = new double[n];
		soln2 = (init == null) ? new double[n] : init;

		// create strategy vectors
		if (generateStrategy) {
			stratChoices = new ArrayList<List<Integer>>(n);
			for (i = 0; i < n; i++)
				stratChoices.add(new LinkedList<Integer>());
		}

		// Initialise solution vectors. Use passed in initial vector, if present
		if (init != null) {
			for (i = 0; i < n; i++)
				soln[i] = soln2[i] = target.get(i) ? 1.0 : init[i];
		} else {
			for (i = 0; i < n; i++)
				soln[i] = soln2[i] = target.get(i) ? 1.0 : 0.0;
		}
		// Store intermediate results if required
		// (compute min/max value over initial states for first step)
		if (results != null) {
			results[0] = Utils.minMaxOverArraySubset(soln2, stpg.getInitialStates(), min2);
		}

		// Start iterations
		iters = 0;
		while (iters < k) {
			iters++;

			if (generateStrategy)
				adv = new int[n];

			// Matrix-vector multiply and min/max ops
			stpg.mvMultMinMax(soln, min1, min2, soln2, target, true, generateStrategy ? adv : null);
			// Store intermediate results if required
			// (compute min/max value over initial states for this step)
			if (results != null) {
				results[iters] = Utils.minMaxOverArraySubset(soln2, stpg.getInitialStates(), min2);
			}

			// Store strategy information
			if (generateStrategy) {
				for (int s = 0; s < n; s++) {
					i = stratChoices.get(s).size();
					// if not yet initialised, or choice has changed, storing
					// initial choice
					if (i == 0 || stratChoices.get(s).get(i - 1) != adv[s]) {
						stratChoices.get(s).add(iters);
						stratChoices.get(s).add(adv[s]);
					} else {
						// increase the count
						stratChoices.get(s).set(stratChoices.get(s).size() - 2, stratChoices.get(s).get(stratChoices.get(s).size() - 2) + 1);
					}
				}
			}
			//			System.out.println(Arrays.toString(soln));
			//			System.out.println(Arrays.toString(adv));

			// Swap vectors for next iter
			tmpsoln = soln;
			soln = soln2;
			soln2 = tmpsoln;
		}

		// Print vector (for debugging)
		//mainLog.println(soln);

		// Finished bounded probabilistic reachability
		timer = System.currentTimeMillis() - timer;
		if (verbosity >= 1) {
			mainLog.print("Bounded probabilistic reachability (" + (min1 ? "min" : "max") + (min2 ? "min" : "max") + ")");
			mainLog.println(" took " + iters + " iterations and " + timer / 1000.0 + " seconds.");
		}

		// Creating strategy object
		int[][] choices = null;
		if (generateStrategy) {
			// converting list into array
			choices = new int[n][];
			for (i = 0; i < n; i++) {
				choices[i] = new int[stratChoices.get(i).size()];
				// reversing the list
				for (int j = stratChoices.get(i).size() - 2, x = 0; j >= 0; j -= 2, x += 2) {
					choices[i][x] = stratChoices.get(i).get(j);
					choices[i][x + 1] = stratChoices.get(i).get(j + 1);
				}
			}
		}

		// Store results/strategy
		res = new ModelCheckerResult();
		res.soln = soln;
		res.lastSoln = soln2;
		res.accuracy = AccuracyFactory.boundedNumericalIterations();
		res.numIters = iters;
		res.timeTaken = timer / 1000.0;
		res.timePre = 0.0;
		if (generateStrategy) {
			res.strat = new StepBoundedDeterministicStrategy(stpg, choices, k);
		}
		
		return res;
	}

	/**
	 * Compute expected reachability rewards.
	 * @param stpg The STPG
	 * @param rewards The rewards
	 * @param target Target states
	 * @param min1 Min or max rewards for player 1 (true=min, false=max)
	 * @param min2 Min or max rewards for player 2 (true=min, false=max)
	 */
	public ModelCheckerResult computeReachRewards(STPG stpg, STPGRewards rewards, BitSet target, boolean min1, boolean min2) throws PrismException
	{
		return computeReachRewards(stpg, rewards, target, min1, min2, null, null);
	}

	/**
	 * Compute expected reachability rewards.
	 * i.e. compute the min/max reward accumulated to reach a state in {@code target}.
	 * @param stpg The STPG
	 * @param rewards The rewards
	 * @param target Target states
	 * @param min1 Min or max rewards for player 1 (true=min, false=max)
	 * @param min2 Min or max rewards for player 2 (true=min, false=max)
	 * @param init Optionally, an initial solution vector (may be overwritten) 
	 * @param known Optionally, a set of states for which the exact answer is known
	 * Note: if 'known' is specified (i.e. is non-null, 'init' must also be given and is used for the exact values.  
	 */
	public ModelCheckerResult computeReachRewards(STPG stpg, STPGRewards rewards, BitSet target, boolean min1, boolean min2, double init[], BitSet known)
			throws PrismException
	{
		return computeReachRewards(stpg, rewards, target, min1, min2, init, known, R_INFINITY);
	}

	/**
	 * Compute expected reachability rewards, where the runs that don't reach
	 * the final state get infinity. i.e. compute the min/max reward accumulated
	 * to reach a state in {@code target}.
	 * @param stpg The STPG
	 * @param rewards The rewards
	 * @param target Target states
	 * @param min1 Min or max rewards for player 1 (true=min, false=max)
	 * @param min2 Min or max rewards for player 2 (true=min, false=max)
	 * @param init Optionally, an initial solution vector (may be overwritten)
	 * @param known Optionally, a set of states for which the exact answer is known
	 * Note: if 'known' is specified (i.e. is non-null, 'init' must also be given and is used for the exact values.  
	 * @param unreachingSemantics Determines how to treat runs that don't reach the target.
	 * One of {@link #R_INFINITY}, {@link #R_CUMULATIVE} and {@link #R_ZERO}.
	 */
	public ModelCheckerResult computeReachRewards(STPG stpg, STPGRewards rewards, BitSet target, boolean min1, boolean min2, double init[], BitSet known,
			int unreachingSemantics) throws PrismException
	{
		switch (unreachingSemantics) {
		case R_INFINITY:
			return computeReachRewardsInfinity(stpg, rewards, target, min1, min2, init, known);
		case R_CUMULATIVE:
			return computeReachRewardsCumulative(stpg, rewards, target, min1, min2, init, known);
		case R_ZERO:
			return computeReachRewardsZero(stpg, rewards, target, min1, min2, init, known);
		default:
			throw new PrismException("Unknown semantics for runs unreaching the target in STPGModelChecker: " + unreachingSemantics);
		}
	}
	
	/**
	 * Compute expected reachability rewards using value iteration.
	 * @param stpg The STPG
	 * @param rewards The rewards
	 * @param target Target states
	 * @param inf States for which reward is infinite
	 * @param min1 Min or max rewards for player 1 (true=min, false=max)
	 * @param min2 Min or max rewards for player 2 (true=min, false=max)
	 * @param init Optionally, an initial solution vector (will be overwritten) 
	 * @param known Optionally, a set of states for which the exact answer is known
	 * Note: if 'known' is specified (i.e. is non-null, 'init' must also be given and is used for the exact values.
	 */
	protected ModelCheckerResult computeReachRewardsValIter(STPG stpg, STPGRewards rewards, BitSet target, BitSet inf, boolean min1, boolean min2,
			double init[], BitSet known) throws PrismException
	{
		ModelCheckerResult res;
		BitSet unknown, notInf;
		int i, n, iters;
		double soln[], soln2[], tmpsoln[];
		int adv[] = null;
		boolean genAdv, done;
		long timer;

		// Are we generating an optimal adversary?
		genAdv = exportAdv || generateStrategy;

		// Start value iteration
		timer = System.currentTimeMillis();
		if (verbosity >= 1)
			mainLog.println("Starting value iteration (" + (min1 ? "min" : "max") + (min2 ? "min" : "max") + ")...");

		// Store num states
		n = stpg.getNumStates();

		// Create solution vector(s)
		soln = new double[n];
		soln2 = (init == null) ? new double[n] : init;

		// Initialise solution vectors. Use (where available) the following in order of preference:
		// (1) exact answer, if already known; (2) 0.0/infinity if in target/inf; (3) passed in initial value; (4) 0.0
		if (init != null) {
			if (known != null) {
				for (i = 0; i < n; i++)
					soln[i] = soln2[i] = known.get(i) ? init[i] : target.get(i) ? 0.0 : inf.get(i) ? Double.POSITIVE_INFINITY : init[i];
			} else {
				for (i = 0; i < n; i++)
					soln[i] = soln2[i] = target.get(i) ? 0.0 : inf.get(i) ? Double.POSITIVE_INFINITY : init[i];
			}
		} else {
			for (i = 0; i < n; i++)
				soln[i] = soln2[i] = target.get(i) ? 0.0 : inf.get(i) ? Double.POSITIVE_INFINITY : 0.0;
		}

		// Determine set of states actually need to compute values for
		unknown = new BitSet();
		unknown.set(0, n);
		unknown.andNot(target);
		unknown.andNot(inf);
		if (known != null)
			unknown.andNot(known);

		// constructing not infinity set
		notInf = (BitSet) inf.clone();
		notInf.flip(0, n);

		// Create/initialise adversary storage
		if (genAdv) {
			adv = new int[n];
			for (i = 0; i < n; i++) {
				adv[i] = -1;
			}

			int s;
			for (i = 0; i < inf.length(); i++) {
				s = inf.nextSetBit(i);
				for (int c = 0; c < stpg.getNumChoices(s); c++) {
					// for player 1 check
					if (stpg.getPlayer(s) == 1 && !stpg.allSuccessorsInSet(s, c, notInf)) {
						adv[i] = c;
						break;
					}
				}
			}
		}

		// Start iterations
		iters = 0;
		done = false;
		while (!done && iters < maxIters) {
			
		        //mainLog.println(soln);
			//mainLog.println(rewards);
			//mainLog.println(min1);
			//mainLog.println(min2);
			//mainLog.println(soln2);
			//mainLog.println(unknown);
			//mainLog.println(genAdv);
			
			iters++;
			// Matrix-vector multiply and min/max ops
			stpg.mvMultRewMinMax(soln, rewards, min1, min2, soln2, unknown, false, genAdv ? adv : null, useDiscounting ? discountFactor : 1.0);

			// Check termination
			done = PrismUtils.doublesAreClose(soln, soln2, termCritParam, termCrit == TermCrit.ABSOLUTE);
			// Swap vectors for next iter
			tmpsoln = soln;
			soln = soln2;
			soln2 = tmpsoln;
		}

		// Finished value iteration
		timer = System.currentTimeMillis() - timer;
		if (verbosity >= 1) {
			mainLog.print("Value iteration (" + (min1 ? "min" : "max") + (min2 ? "min" : "max") + ")");
			mainLog.println(" took " + iters + " iterations and " + timer / 1000.0 + " seconds.");
		}

		// Print adversary
		if (genAdv) {
			PrismLog out = new PrismFileLog(exportAdvFilename);
			if (exportAdvFilename.lastIndexOf('.') != -1 && exportAdvFilename.substring(exportAdvFilename.lastIndexOf('.') + 1).equals("dot")) {
				stpg.exportToDotFileWithStrat(out, null, adv);
			} else {
				for (i = 0; i < n; i++) {
					out.println(i + " " + (adv[i] != -1 ? stpg.getAction(i, adv[i]) : "-"));
				}
				out.println();
			}
			out.close();
		}

		// Non-convergence is an error (usually)
		if (!done && errorOnNonConverge) {
			String msg = "Iterative method did not converge within " + iters + " iterations.";
			msg += "\nConsider using a different numerical method or increasing the maximum number of iterations";
			throw new PrismException(msg);
		}

		// Store results/strategy
		res = new ModelCheckerResult();
		res.soln = soln;
		double maxDiff = PrismUtils.measureSupNorm(soln, soln2, termCrit == TermCrit.ABSOLUTE);
		res.accuracy = AccuracyFactory.valueIteration(termCritParam, maxDiff, termCrit == TermCrit.ABSOLUTE);
		res.numIters = iters;
		res.timeTaken = timer / 1000.0;
		if (generateStrategy) {
			res.strat = new MemorylessDeterministicStrategy(adv);
		}

		return res;
	}

	/**
	 * Computes the reachability reward under the semantics where nonreaching
	 * runs get infinity.
	 * 
	 * @param stpg
	 * @param rewards
	 * @param target
	 * @param min1
	 * @param min2
	 * @param init
	 * @param known
	 * @return
	 * @throws PrismException
	 */
	public ModelCheckerResult computeReachRewardsInfinity(STPG stpg, STPGRewards rewards, BitSet target, boolean min1, boolean min2, double init[], BitSet known)
			throws PrismException
	{
		ModelCheckerResult res = null;
		BitSet inf;
		int n, numTarget, numInf;
		long timer, timerProb1, timerApprox;

		// Start expected reachability
		timer = System.currentTimeMillis();
		if (verbosity >= 1)
			mainLog.println("\nStarting expected reachability...");

		// Check for deadlocks in non-target state (because breaks e.g. prob1)
		stpg.checkForDeadlocks(target);

		// Store num states
		n = stpg.getNumStates();

		// Optimise by enlarging target set (if more info is available)
		if (init != null && known != null && !known.isEmpty()) {
			BitSet targetNew = (BitSet) target.clone();
			for (int i : new IterableBitSet(known)) {
				if (init[i] == 1.0) {
					targetNew.set(i);
				}
			}
			target = targetNew;
		}

		timerProb1 = System.currentTimeMillis();
		
		// identify infinite values
		inf = prob1(stpg, null, target, !min1, !min2);
		inf.flip(0, n);

		timerProb1 = System.currentTimeMillis() - timerProb1;

		// Print results of precomputation
		numTarget = target.cardinality();
		numInf = inf.cardinality();
		if (verbosity >= 1)
			mainLog.println("target=" + numTarget + ", inf=" + numInf + ", rest=" + (n - (numTarget + numInf)));

		// Compute rewards with epsilon instead of zero. This is used to get the
		// over-approximation
		// of the real result, which deals with the problem of staying in zero
		// components for free
		// when infinity should be gained.

		// first, get the minimum nonzero reward and maximal reward, will be
		// used as a basis for epsilon
		// also, check if by any chance all rewards are nonzero, then we don't
		// need to precompute
		double minimumReward = Double.POSITIVE_INFINITY;
		double maximumReward = 0.0;
		boolean allNonzero = true;
		double r;
		for (int i = 0; i < n; i++) {
			r = rewards.getStateReward(i);
			if (r > 0.0 && r < minimumReward)
				minimumReward = r;
			if (r > maximumReward)
				maximumReward = r;
			allNonzero = allNonzero && r > 0;

			for (int j = 0; j < stpg.getNumChoices(i); j++) {
				r = rewards.getTransitionReward(i, j);
				if (r > 0.0 && r < minimumReward)
					minimumReward = r;
				if (r > maximumReward)
					maximumReward = r;
				allNonzero = allNonzero && rewards.getTransitionReward(i, j) > 0;

				for (int k = 0; k < stpg.getNumNestedChoices(i, j); k++) {
					r = rewards.getNestedTransitionReward(i, j, k);
					if (r > 0.0 && r < minimumReward)
						minimumReward = r;
					if (r > maximumReward)
						maximumReward = r;
					allNonzero = allNonzero && r > 0;
				}
			}
		}

		if (!allNonzero && !(rewards instanceof StateRewardsConstant)) {
			timerApprox = System.currentTimeMillis();
			// A simple heuristic that gives small epsilon, but still is
			// hopefully safe floating-point-wise
			double epsilon = Math.min(minimumReward, maximumReward * 0.01);
			;

			if (verbosity >= 1) {
				mainLog.println("Computing the upper bound where " + epsilon + " is used instead of 0.0");
			}

			// Modify the rewards
			double origZeroReplacement;
			if (rewards instanceof MDPRewardsSimple) {
				origZeroReplacement = ((MDPRewardsSimple) rewards).getZeroReplacement();
				((MDPRewardsSimple) rewards).setZeroReplacement(epsilon);
			} else {
				throw new PrismException("To compute expected reward I need to modify the reward structure. But I don't know how to modify"
						+ rewards.getClass().getName());
			}

			// Compute the value when rewards are nonzero
			switch (solnMethod) {
			case VALUE_ITERATION:
				res = computeReachRewardsValIter(stpg, rewards, target, inf, min1, min2, init, known);
				break;
			default:
				throw new PrismException("Unknown STPG solution method " + solnMethod);
			}

			// Set the value iteration result to be the initial solution for the
			// next part
			// in which "proper" zero rewards are used
			init = res.soln;

			// Return the rewards to the original state
			if (rewards instanceof MDPRewardsSimple) {
				((MDPRewardsSimple) rewards).setZeroReplacement(origZeroReplacement);
			}

			timerApprox = System.currentTimeMillis() - timerApprox;

			if (verbosity >= 1) {
				mainLog.println("Computed an over-approximation of the solution (in " + timerApprox / 1000
						+ " seconds), this will now be used to get the solution");
			}
		}

		// Compute real rewards
		switch (solnMethod) {
		case VALUE_ITERATION:
			res = computeReachRewardsValIter(stpg, rewards, target, inf, min1, min2, init, known);
			break;
		default:
			throw new PrismException("Unknown STPG solution method " + solnMethod);
		}

		// Finished expected reachability
		timer = System.currentTimeMillis() - timer;
		if (verbosity >= 1)
			mainLog.println("Expected reachability took " + timer / 1000.0 + " seconds.");

		// Update time taken
		res.timeTaken = timer / 1000.0;
		res.timePre = timerProb1 / 1000.0;

		return res;
	}

	/**
	 * Computes the reachability reward under the semantics where nonreaching
	 * runs get their total cumulative reward (i.e. anything between 0 and
	 * infinity).
	 * 
	 * @param stpg
	 * @param rewards
	 * @param target
	 * @param min1
	 * @param min2
	 * @param init
	 * @param known
	 * @return
	 * @throws PrismException
	 */
	public ModelCheckerResult computeReachRewardsCumulative(STPG stpg, STPGRewards rewards, BitSet target, boolean min1, boolean min2, double init[],
			BitSet known) throws PrismException
	{
		ModelCheckerResult res = null;
		BitSet inf;
		int i, n, numTarget, numInf;
		long timer, timerProb1, timerApprox;

		// Start expected reachability
		timer = System.currentTimeMillis();
		if (verbosity >= 1)
			mainLog.println("\nStarting expected reachability...");

		// Check for deadlocks in non-target state (because breaks e.g. prob1)
		stpg.checkForDeadlocks(target);

		// Store num states
		n = stpg.getNumStates();

		// Optimise by enlarging target set (if more info is available)
		if (init != null && known != null) {
			BitSet targetNew = new BitSet(n);
			for (i = 0; i < n; i++) {
				targetNew.set(i, target.get(i) || (known.get(i) && init[i] == 0.0));
			}
			target = targetNew;
		}

		timerProb1 = System.currentTimeMillis();
		// identify infinite values
		BitSet aRew = new BitSet();

		if (!useDiscounting) {

			for (i = 0; i < n; i++) {
				// skipping target states
				if (target.get(i))
					continue;

				// check for state reward
				if (rewards.getStateReward(i) > 0.0)
					aRew.set(i);

				// check for transition rewards
				int nonZeroRewards = 0;
				int inftyRewards = 0;
				double trp;
				for (int k = 0; k < stpg.getNumChoices(i); k++) {
					trp = rewards.getTransitionReward(i, k);
					// ignoring infinite rewards as these transitions will neven be
					// taken
					if (trp > 0.0 && trp != Double.POSITIVE_INFINITY && trp != Double.NEGATIVE_INFINITY) {
						nonZeroRewards++;
						aRew.set(i);
					} else if (trp == Double.POSITIVE_INFINITY || trp == Double.NEGATIVE_INFINITY)
						inftyRewards++;
				}

				if (nonZeroRewards != 0 && nonZeroRewards != stpg.getNumChoices(i) - inftyRewards)
					throw new PrismException("If transition reward is nonzero, all transitions going from the state must be.");
			}
		}
		BitSet b1 = aRew;
		BitSet b2 = new BitSet();

		BitSet all = new BitSet(n);
		all.flip(0, n);
		// BitSet none = new BitSet();

		while (true) {
			b2 = prob1(stpg, null, b1, min1, min2);

			BitSet b3 = new BitSet();
			stpg.prob1step(all, b2, all, min1, min2, b3);
			b3.and(b1);

			// check if the alg is correct
			for (i = 0; i < n; i++) {
				if (b3.get(i) && !b1.get(i)) {
					throw new PrismException("There is some error in the implementation");
				}
			}

			if (b3.equals(b1))
				break;

			BitSet tmp = b3;
			b3 = b1;
			b1 = tmp;
		}

		inf = prob0(stpg, null, b1, min1, min2);
		inf.flip(0, n);

		timerProb1 = System.currentTimeMillis() - timerProb1;

		// Print results of precomputation
		numTarget = target.cardinality();
		numInf = inf.cardinality();
		if (verbosity >= 1)
			mainLog.println("target=" + numTarget + ", inf=" + numInf + ", rest=" + (n - (numTarget + numInf)));

		// Compute real rewards
		switch (solnMethod) {
		case VALUE_ITERATION:
			res = computeReachRewardsValIter(stpg, rewards, target, inf, min1, min2, init, known);
			break;
		default:
			throw new PrismException("Unknown STPG solution method " + solnMethod);
		}

		// Finished expected reachability
		timer = System.currentTimeMillis() - timer;
		if (verbosity >= 1)
			mainLog.println("Expected reachability took " + timer / 1000.0 + " seconds.");

		// Update time taken
		res.timeTaken = timer / 1000.0;
		res.timePre = timerProb1 / 1000.0;

		return res;
	}

	/**
	 * Computes the reachability reward under the semantics where nonreaching
	 * runs get 0.
	 * 
	 * @param stpg
	 * @param rewards
	 * @param target
	 * @param min1
	 * @param min2
	 * @param init
	 * @param known
	 * @param unreachingAsInfinity
	 * @return
	 * @throws PrismException
	 */
	public ModelCheckerResult computeReachRewardsZero(STPG stpg, STPGRewards rewards, BitSet target, boolean min1, boolean min2, double init[], BitSet known)
			throws PrismException
	{
		ModelCheckerResult res = null;
		BitSet inf;
		int i, n, numTarget, numInf;
		long timer, timerProb1, timerApprox;
		List<List<Integer>> stratChoices = null;
		int[] adv = null;
		boolean updateChoice;

		// Start expected reachability
		timer = System.currentTimeMillis();
		if (verbosity >= 1)
			mainLog.println("\nStarting expected reachability...");

		// Check for deadlocks in non-target state (because breaks e.g. prob1)
		stpg.checkForDeadlocks(target);

		// Store num states
		n = stpg.getNumStates();

		// Currently we only allow integer rewards, check if all rewards are
		// (close to) an integer.
		// While traversing rewards, get largest reward, too
		boolean hasNonInt = false;
		double nonInt = 0; // will be used for output if there is non-integer
		// reward
		int maxReward = 0;
		checkrewards: for (int s = 0; s < n; s++) {
			double sr = rewards.getStateReward(s);
			if (sr != Math.floor(sr)) {
				hasNonInt = true;
				nonInt = sr;
				break;
			}

			if (sr > maxReward)
				maxReward = (int) sr;

			for (int c = 0; c < stpg.getNumChoices(s); c++) {
				double tr = rewards.getTransitionReward(s, c);
				if (tr != Math.floor(tr)) {
					hasNonInt = true;
					nonInt = tr;
					break checkrewards;
				}

				if (tr > maxReward)
					maxReward = (int) tr;
			}
		}

		if (verbosity >= 1)
			mainLog.println("Maximal reward is " + maxReward);

		if (hasNonInt)
			throw new PrismException("For 'zero' semantics reachability reward all rewards must be integers." + "There is at least one non-integer reward: "
					+ nonInt);

		// Optimise by enlarging target set (if more info is available)
		if (init != null && known != null) {
			BitSet targetNew = new BitSet(n);
			for (i = 0; i < n; i++) {
				targetNew.set(i, target.get(i) || (known.get(i) && init[i] == 0.0));
			}
			target = targetNew;
		}

		// TODO identify "dead" states, i.e. those from which F can't be reached
		// with >0 prob
		// and those from which the bad player can ensure 0. This is optional,
		// but
		// should bring some speedup.

		BitSet zeroProb = prob0(stpg, null, target, min1, min2);

		BitSet positiveProb = new BitSet();
		for (int k = 0; k < n; k++)
			positiveProb.set(k, !zeroProb.get(k));

		// ...and those from which the bad player can ensure 0.
		// BitSet zeroReward = zeroRewards(stpg, rewards, positiveProb, target,
		// !min1, !min2);

		// Identify states that get infinity.
		// First, remove the rewards which are gained at places from which the
		// target can't be reached
		STPGRewards rewardsRestricted;
		if (rewards instanceof MDPRewardsSimple) {
			// And make sure only the best actions are used
			STPGRewardsSimple rewardsRestrictedSimple = new STPGRewardsSimple((MDPRewardsSimple) rewards);

			for (int s = 0; s < n; s++) {
				for (int c = 0; c < stpg.getNumChoices(s); c++) {
					Iterator<Entry<Integer, Double>> iterator = stpg.getTransitionsIterator(s, c);
					boolean hasPositiveSuccessor = false;
					while (iterator.hasNext()) {
						if (positiveProb.get(iterator.next().getKey())) {
							hasPositiveSuccessor = true;
							break;
						}
					}
					if (!hasPositiveSuccessor)
						rewardsRestrictedSimple.setTransitionReward(s, c, 0);
				}
				if (!positiveProb.get(s))
					rewardsRestrictedSimple.setStateReward(s, 0);
			}
			rewardsRestricted = rewardsRestrictedSimple;
		} else if (rewards instanceof StateRewardsConstant) {
			// And make sure only the best actions are used
			STPGRewardsSimple rewardsRestrictedSimple = new STPGRewardsSimple(n);

			for (int s = 0; s < n; s++) {
				if (positiveProb.get(s))
					rewardsRestrictedSimple.setStateReward(s, rewards.getStateReward(s));
			}
			rewardsRestricted = rewardsRestrictedSimple;
		} else {
			throw new PrismException("To compute expected reward I need to modify the reward structure. But I don't know how to modify"
					+ rewards.getClass().getName());
		}

		BitSet aRew = new BitSet();

		for (i = 0; i < n; i++) {
			// skipping target states
			if (target.get(i))
				continue;

			// check for state reward
			if (rewardsRestricted.getStateReward(i) > 0.0)
				aRew.set(i);

			// check for transition rewards
			int nonZeroRewards = 0;
			int inftyRewards = 0;
			double trp;
			for (int k = 0; k < stpg.getNumChoices(i); k++) {
				trp = rewards.getTransitionReward(i, k);
				// ignoring infinite rewards as these transitions will neven be
				// taken
				if (trp > 0.0 && trp != Double.POSITIVE_INFINITY && trp != Double.NEGATIVE_INFINITY) {
					nonZeroRewards++;
					aRew.set(i);
				} else if (trp == Double.POSITIVE_INFINITY || trp == Double.NEGATIVE_INFINITY)
					inftyRewards++;
			}

			if (nonZeroRewards != 0 && nonZeroRewards != stpg.getNumChoices(i) - inftyRewards)
				throw new PrismException("If transition reward is nonzero, all transitions going from the state must be.");
		}

		BitSet b1 = aRew;
		BitSet b2 = new BitSet();

		BitSet all = new BitSet(n);
		all.flip(0, n);
		// BitSet none = new BitSet();

		while (true) {
			b2 = prob1(stpg, null, b1, min1, min2);

			BitSet b3 = new BitSet();
			stpg.prob1step(all, b2, all, min1, min2, b3);
			b3.and(b1);

			// check if the alg is correct
			for (i = 0; i < n; i++) {
				if (b3.get(i) && !b1.get(i)) {
					throw new PrismException("There is some error in the implementation");
				}
			}

			if (b3.equals(b1))
				break;

			BitSet tmp = b3;
			b3 = b1;
			b1 = tmp;
		}

		// excluding states from which cannot reach the target
		b1.and(positiveProb);

		// computing infty states
		inf = prob0(stpg, null, b1, min1, min2);
		inf.flip(0, n);

		// Get the rich man's strategy and its values
		// Start with computing optimal probabilities to reach the final state
		ModelCheckerResult mcrprob = computeReachProbs(stpg, target, min1, min2);

		// Next, reweigh the rewards and make sure that only optimal actions are
		// taken
		if (rewards instanceof MDPRewardsSimple) {
			// And make sure only the best actions are used
			STPGRewardsSimple rewardsRestrictedSimple = new STPGRewardsSimple((MDPRewardsSimple) rewards);

			for (int s = 0; s < n; s++) {
				for (int c = 0; c < stpg.getNumChoices(s); c++) {
					double prob = 0.0;
					Iterator<Entry<Integer, Double>> it = stpg.getTransitionsIterator(s, c);
					while (it.hasNext()) {
						Entry<Integer, Double> e = it.next();
						prob += e.getValue() * mcrprob.soln[e.getKey()];
					}

					// as a hack, set the transition reward of nonoptimal
					// transitions
					// to something extreme so they are never chosen

					if (stpg.getNumChoices(s) > 1 && prob < mcrprob.soln[s] && ((stpg.getPlayer(s) == 1 && !min1) || (stpg.getPlayer(s) == 2 && !min2))) {
						rewardsRestrictedSimple.setTransitionReward(s, c, Double.NEGATIVE_INFINITY);
					} else if (stpg.getNumChoices(s) > 1 && prob > mcrprob.soln[s] && ((stpg.getPlayer(s) == 1 && min1) || (stpg.getPlayer(s) == 2 && min2))) {
						rewardsRestrictedSimple.setTransitionReward(s, c, Double.POSITIVE_INFINITY);
					} else {
						double newReward = rewards.getTransitionReward(s, c) * mcrprob.soln[s];
						rewardsRestrictedSimple.setTransitionReward(s, c, newReward);
					}
				}
				double newReward = rewards.getStateReward(s) * mcrprob.soln[s];
				rewardsRestrictedSimple.setStateReward(s, newReward);
			}
			rewardsRestricted = rewardsRestrictedSimple;
		} else if (rewards instanceof StateRewardsConstant) {
			STPGRewardsSimple rewardsRestrictedSimple = new STPGRewardsSimple(n);

			for (int s = 0; s < n; s++) {
				for (int c = 0; c < stpg.getNumChoices(s); c++) {
					double prob = 0.0;
					Iterator<Entry<Integer, Double>> it = stpg.getTransitionsIterator(s, c);
					while (it.hasNext()) {
						Entry<Integer, Double> e = it.next();
						prob += e.getValue() * mcrprob.soln[e.getKey()];
					}

					// as a hack, set the transition reward of nonoptimal
					// transitions
					// to something extreme so they are never chosen
					if (stpg.getNumChoices(s) > 1 && prob < mcrprob.soln[s] && ((stpg.getPlayer(s) == 1 && !min1) || (stpg.getPlayer(s) == 2 && !min2))) {
						rewardsRestrictedSimple.setTransitionReward(s, c, Double.NEGATIVE_INFINITY);
					} else if (stpg.getNumChoices(s) > 1 && prob > mcrprob.soln[s] && ((stpg.getPlayer(s) == 1 && min1) || (stpg.getPlayer(s) == 2 && min2))) {
						rewardsRestrictedSimple.setTransitionReward(s, c, Double.POSITIVE_INFINITY);
					} // else the reward remains 0
				}
				double newReward = rewards.getStateReward(s) * mcrprob.soln[s];
				rewardsRestrictedSimple.setStateReward(s, newReward);
			}
			rewardsRestricted = rewardsRestrictedSimple;
		} else {
			throw new PrismException("To compute expected reward I need to modify the reward structure. But I don't know how to modify"
					+ rewards.getClass().getName());
		}
		// Next, compute the value for the rich man's strategy.
		ModelCheckerResult mcrrich = computeReachRewards(stpg, rewardsRestricted, target, min1, min2, init, known, R_CUMULATIVE);
		//System.out.println("maximal rews for rich man's strategy: " + Arrays.toString(mcrrich.soln));

		// TODO generate rich man strategy.

		// compute B from the values for the rich man's strategy
		int lastSwitch = 0;
		for (int s = 0; s < n; s++) {
			// for all choices c, find maximal B such that
			// sum_{s'} prob(c,s')(r(s) + r(c) + B + rewRich(s')) >
			// probRich(s)*B + rewRich(s)
			for (int c = 0; c < stpg.getNumChoices(s); c++) {
				double numerator = mcrrich.soln[s];
				double denominator = -mcrprob.soln[s];
				double tRew = rewards.getTransitionReward(s, c);
				Iterator<Entry<Integer, Double>> it = stpg.getTransitionsIterator(s, c);
				while (it.hasNext()) {
					Entry<Integer, Double> e = it.next();
					int ts = e.getKey();
					double sRew = rewards.getStateReward(s);
					double p = e.getValue();

					numerator -= p * (mcrprob.soln[ts] * (sRew + tRew) + mcrrich.soln[ts]);
					denominator += p * mcrprob.soln[ts];

					int b = (denominator == 0) ? 0 : (int) Math.floor(numerator / denominator);

					if (lastSwitch < b)
						lastSwitch = b;
				}
			}
		}

		if (verbosity >= 1)
			mainLog.println("Last switching point is when the reward cumulated in the past becomes " + lastSwitch);

		// TODO using gcd of rewards could save us iterations in many cases
		int kSize = maxReward + 1;
		double[][] rews = new double[kSize][n];
		int iters = 0;

		double[] tmp;

		// fill in the initial values from the rich man's strategy
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < kSize; k++) {
				if (inf.get(j))
					rews[k][j] = Double.POSITIVE_INFINITY;
				else
					rews[k][j] = mcrrich.soln[j] + (lastSwitch + k) * mcrprob.soln[j];
			}
		}

		// create strategy vectors
		if (generateStrategy) {
			stratChoices = new ArrayList<List<Integer>>(n);
			for (i = 0; i < n; i++)
				stratChoices.add(new LinkedList<Integer>());
		}

		for (int x = lastSwitch; x >= 0; x--) {
			// reward[s,x] =
			// opt_c sum_{s'} p(s ->c s')(prob_F[s,x+r(s)+r(c)]*(r(s)+r(c))
			// +
			// reward[s'][x+r(s)+r(c)])
			// where opt is either max or min, depending on the owner of s
			// probs[s,x] is the probability of reaching F under the choice
			// c
			// chosen by rews
			boolean done = false;

			double difference = 0;
			do {
				difference = 0;

				iters++;

				for (int s = 0; s < n; s++) {
					if (target.get(s)) {
						rews[0][s] = x;
					}
				}

				if (generateStrategy)
					adv = new int[n];
				for (int s = 0; s < n; s++) {
					if (target.get(s)) {
						continue;
					}

					// non-target states
					boolean min = (stpg.getPlayer(s) == 1) ? min1 : min2;
					double stateRew = -1;
					for (int c = 0; c < stpg.getNumChoices(s); c++) {
						double choiceRew = 0;
						double r = rewards.getStateReward(s) + rewards.getTransitionReward(s, c);
						int index = (int) r; // the reward determines in
						// which
						// array we will look

						// System.out.println("s="+s + " c="+c +"index="+index);

						Iterator<Entry<Integer, Double>> it = stpg.getTransitionsIterator(s, c);
						while (it.hasNext()) {
							Entry<Integer, Double> e = it.next();
							int ts = e.getKey();
							double p = e.getValue();
							// choiceRew += p*(probs[index][ts]*r +
							// rews[index][ts]);
							choiceRew += p * (rews[index][ts]);
						}

						updateChoice = false;
						if (stateRew < 0) {
							stateRew = choiceRew;
							if (generateStrategy)
								adv[s] = c;
						} else if (min && stateRew > choiceRew) {
							stateRew = choiceRew;
							if (generateStrategy)
								adv[s] = c;
						} else if (!min && stateRew < choiceRew) {
							stateRew = choiceRew;
							if (generateStrategy)
								adv[s] = c;
						}

					}
					double cDif = Math.abs(rews[0][s] - stateRew);
					if (cDif > difference)
						difference = cDif;

					rews[0][s] = stateRew;

				}

			} while (difference > 10e-6); // TODO some smarter convergence
			// test

			// Store strategy information
			if (generateStrategy) {
				for (int s = 0; s < n; s++) {
					i = stratChoices.get(s).size();
					// if not yet initialised, or choice has changed, storing
					// initial choice
					if (i == 0 || stratChoices.get(s).get(i - 1) != adv[s]) {
						stratChoices.get(s).add(iters);
						stratChoices.get(s).add(adv[s]);
					} else {
						// increase the count
						stratChoices.get(s).set(stratChoices.get(s).size() - 2, stratChoices.get(s).get(stratChoices.get(s).size() - 2) + 1);
					}
				}
			}

			// shift the array
			double[] tmpRews = rews[kSize - 1];
			for (i = kSize - 1; i >= 1; i--) {
				rews[i] = rews[i - 1];
			}
			rews[0] = tmpRews;

		}
		timer = System.currentTimeMillis() - timer;

		// Creating strategy object
		int[][] choices = null;
		if (generateStrategy) {
			// converting list into array
			choices = new int[n][];
			for (i = 0; i < n; i++) {
				choices[i] = new int[stratChoices.get(i).size()];

				// reversing the list
				for (int j = stratChoices.get(i).size() - 2, x = 0; j >= 0; j -= 2, x += 2) {
					choices[i][x] = stratChoices.get(i).get(j);
					choices[i][x + 1] = stratChoices.get(i).get(j + 1);
				}
			}
		}

		// Store results/strategy
		res = new ModelCheckerResult();
		res.soln = (rews.length > 1) ? rews[1] : rews[0];
		res.lastSoln = (rews.length > 2) ? rews[2] : null;
		res.numIters = lastSwitch;
		res.timeTaken = timer / 1000;
		res.numIters = iters;
		if (generateStrategy) {
			res.strat = new BoundedRewardDeterministicStrategy(stpg, choices, lastSwitch, rewards);
		}
		
		return res;
	}

	/**
	 * Zero cummulative reward precomputation algorithm. i.e. determine the
	 * states of an STPG which, with probability 1 get min/max reward equal to
	 * 0.0 before (possibly) reaching a state in {@code target}, while remaining
	 * in those in {@code remain}.
	 * @param stpg The STPG
	 * @param remain Remain in these states (optional: null means "all")
	 * @param target Target states
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 */
	public BitSet zeroRewards(STPG stpg, STPGRewards rewards, BitSet remain, BitSet target, boolean min1, boolean min2)
	{
		int n, iters;
		double[] soln1, soln2;
		BitSet unknown;
		boolean done;
		long timer;

		// Start precomputation
		timer = System.currentTimeMillis();
		if (verbosity >= 1)
			mainLog.println("Starting zeroRewards (" + (min1 ? "min" : "max") + (min2 ? "min" : "max") + ")...");

		// Initialise vectors
		n = stpg.getNumStates();
		soln1 = new double[n];
		soln2 = new double[n];

		// Determine set of states actually need to perform computation for
		unknown = new BitSet();
		unknown.set(0, n);
		if (target != null)
			unknown.andNot(target);
		if (remain != null)
			unknown.and(remain);

		// initialise the solution so that the forbidden states are penalised
		for (int i = 0; i < n; i++) {
			if (remain != null && !remain.get(i) && target != null && !target.get(i))
				soln1[i] = Double.POSITIVE_INFINITY;
		}

		// Nested fixed point loop
		iters = 0;
		done = false;
		while (!done) {
			iters++;
			// at every iter at least one state must go from zero to nonzero,
			// hence we have
			// at most n iterations
			assert iters <= n + 1;

			stpg.mvMultRewMinMax(soln1, rewards, min1, min2, soln2, unknown, false, null);

			// Check termination (outer)
			done = true;

			double[] tmp = soln2;
			soln2 = soln1;
			soln1 = tmp;

			done = true;
			for (int i = 0; i < n; i++) {
				if (soln1[i] > 0.0 && soln2[i] == 0.0) {
					done = false;
					break;
				}
			}
		}

		// Finished precomputation
		timer = System.currentTimeMillis() - timer;
		if (verbosity >= 1) {
			mainLog.print("Zero Rewards (" + (min1 ? "min" : "max") + (min2 ? "min" : "max") + ")");
			mainLog.println(" took " + iters + " iterations and " + timer / 1000.0 + " seconds.");
		}

		BitSet result = new BitSet(n);
		for (int i = 0; i < n; i++) {
			if (soln1[i] == 0.0)
				result.set(i);
		}

		return result;
	}

	/**
	 * Simple test program.
	 */
	public static void main(String args[])
	{
		STPGModelChecker mc;
		STPGAbstrSimple stpg;
		ModelCheckerResult res;
		BitSet target = new BitSet();
		Map<String, BitSet> labels;
		boolean min1 = true, min2 = true;
		try {
			mc = new STPGModelChecker(null);
			stpg = new STPGAbstrSimple();
			stpg.buildFromPrismExplicit(args[0]);
			stpg.addInitialState(0);
			//System.out.println(stpg);
			labels = StateModelChecker.loadLabelsFile(args[1]);
			//System.out.println(labels);
			target = labels.get(args[2]);
			if (target == null)
				throw new PrismException("Unknown label \"" + args[2] + "\"");
			for (int i = 3; i < args.length; i++) {
				if (args[i].equals("-minmin")) {
					min1 = true;
					min2 = true;
				} else if (args[i].equals("-maxmin")) {
					min1 = false;
					min2 = true;
				} else if (args[i].equals("-minmax")) {
					min1 = true;
					min2 = false;
				} else if (args[i].equals("-maxmax")) {
					min1 = false;
					min2 = false;
				}
			}
			// stpg.exportToDotFile("stpg.dot", target);
			// stpg.exportToPrismExplicit("stpg");
			res = mc.computeReachProbs(stpg, target, min1, min2);
			System.out.println(res.soln[0]);
		} catch (PrismException e) {
			System.out.println(e);
		}
	}
}
