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

import java.util.*;
import java.util.Map.Entry;

import parser.ast.Expression;
import parser.ast.ExpressionTemporal;
import prism.*;
import explicit.rewards.MDPRewards;
import explicit.rewards.MDPRewardsSimple;
import explicit.rewards.STPGRewards;
import explicit.rewards.STPGRewardsSimple;

/**
 * Explicit-state model checker for two-player stochastic games (STPGs).
 */
public class STPGModelChecker extends ProbModelChecker
{
	// Model checking functions

	/**
	 * Compute probabilities for the contents of a P operator.
	 */
	protected StateValues checkProbPathFormula(Model model, Expression expr, boolean min1, boolean min2, double bound) throws PrismException
	{
		// Test whether this is a simple path formula (i.e. PCTL)
		// and then pass control to appropriate method. 
		if (expr.isSimplePathFormula()) {
			return checkProbPathFormulaSimple(model, expr, min1, min2, bound);
		} else {
			throw new PrismException("Explicit engine does not yet handle LTL-style path formulas");
		}
	}

	/**
	 * Compute probabilities for a simple, non-LTL path operator.
	 */
	protected StateValues checkProbPathFormulaSimple(Model model, Expression expr, boolean min1, boolean min2, double bound) throws PrismException
	{
		StateValues probs = null;

		// Temporal operators
		if (expr instanceof ExpressionTemporal) {
			ExpressionTemporal exprTemp = (ExpressionTemporal) expr;
			
			// Next
			if(exprTemp.getOperator() == ExpressionTemporal.P_X)
			{
				probs = checkProbNext(model, exprTemp, min1, min2);
			}
			
			// Until
			else if (exprTemp.getOperator() == ExpressionTemporal.P_U) {
				if (exprTemp.hasBounds()) {
					probs = checkProbBoundedUntil(model, exprTemp, min1, min2);
				} else {
					probs = checkProbUntil(model, exprTemp, min1, min2, bound);
				}
			}
			// Anything else - convert to until and recurse
			else {
				probs = checkProbPathFormulaSimple(model, exprTemp.convertToUntilForm(), min1, min2, bound);
			}
		}

		if (probs == null)
			throw new PrismException("Unrecognised path operator in P operator");

		return probs;
	}

	private StateValues checkProbNext(Model model, ExpressionTemporal expr,
			boolean min1, boolean min2) throws PrismException {

		BitSet b;
		int n;
		double soln[], soln2[];
		STPG stpg;

		
		stpg = (STPG) model;
		
		// model check the operand
		b = checkExpression(model, expr.getOperand2()).getBitSet();
		
		// Store num states
		n = model.getNumStates();
		
		// Create solution vector(s)
		soln = new double[n];
		soln2 = new double[n];
		
		for(int i=0; i<n; i++)
			soln[i] = b.get(i)?1.0:0.0;
		
		stpg.mvMultMinMax(soln, min1, min2, soln2, null, false, null);
		
		// Return results
		return StateValues.createFromDoubleArray(soln2, model);
	}

	/**
	 * Computes a probability that a formula (F G target) is satisfied.
	 * <p/>
	 * This method exploits the fact that under any fixed strategies
	 * this probability is equivalent to 
	 * (1 - probability of satisfying formula (G F not target))
	 * 
	 * @param stpg The model.
	 * @param target The set of states satisfying the innermost subformula.
	 * @param expr The expression to verify (including the F G part)
	 * @param min1 Min or max probabilities for player 1 (true=lower bound, false=upper bound)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 * @return
	 * @throws PrismException
	 */
	protected StateValues checkGF(Model model, Expression expr, boolean min1, boolean min2) throws PrismException
	{
		BitSet b;
		ModelCheckerResult res = null;
		
		//check if the formula is of the FG form
		Expression subformula = null;
		if (expr instanceof ExpressionTemporal) {
			ExpressionTemporal exprT = (ExpressionTemporal) expr;
			// And children, if present, must be state (not path) formulas
			if (exprT.getOperator() == ExpressionTemporal.P_G) {
				Expression expr2 = exprT.getOperand2();
				if (expr2 instanceof ExpressionTemporal) {
					ExpressionTemporal expr2T = (ExpressionTemporal) expr2;
					if (expr2T.getOperator() == ExpressionTemporal.P_F) {
						Expression expr3 = expr2T.getOperand2();
						if (!(expr3 instanceof ExpressionTemporal)) {
							subformula = expr3;
						}
					}
				}
			}
		}
		
		if (subformula == null)
			throw new PrismException("The expression passed to checkGF must be of the form (G F psi) where psi is a state formula");

		// model check operand first, then negate it because we want 'not target'
		b = checkExpression(model, subformula).getBitSet();
		b.flip(0,((STPG)model).getNumStates());

		//model check using FG and swapped mins/maxes
		res = computeFG((STPG) model, b, !min1, !min2);
		
		//get (1 - result) and return it
		for(int i = 0; i < res.soln.length; i++)
			res.soln[i] = 1 - res.soln[i];
		
		StateValues probs = StateValues.createFromDoubleArray(res.soln, model);
		return probs;
	}
	
	/**
	 * Computes a probability that a formula (FG target) is satisfied.
	 * <p/>
	 * This method exploits the fact that this probability is
	 * equivalent to the probability of satisfying formula
	 * (F (P=1 [ G target ])).
	 * @param stpg The model.
	 * @param target The set of states satisfying the innermost subformula.
	 * @param expr The expression to verify (including the F G part)
	 * @param min1 Min or max probabilities for player 1 (true=lower bound, false=upper bound)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 * @return
	 * @throws PrismException
	 */
	protected StateValues checkFG(Model model, Expression expr, boolean min1, boolean min2) throws PrismException
	{
		BitSet b;
		ModelCheckerResult res = null;
		
		System.out.println("sss");
		
		//check if the formula is of the FG form
		Expression subformula = null;
		if (expr instanceof ExpressionTemporal) {
			ExpressionTemporal exprT = (ExpressionTemporal) expr;
			// And children, if present, must be state (not path) formulas
			if (exprT.getOperator() == ExpressionTemporal.P_F) {
				Expression expr2 = exprT.getOperand2();
				if (expr2 instanceof ExpressionTemporal) {
					ExpressionTemporal expr2T = (ExpressionTemporal) expr2;
					if (expr2T.getOperator() == ExpressionTemporal.P_G) {
						Expression expr3 = expr2T.getOperand2();
						if (!(expr3 instanceof ExpressionTemporal)) {
							subformula = expr3;
						}
					}
				}
			}
		}
		
		if (subformula == null)
			throw new PrismException("The expression passed to checkFG must be of the form (F G psi) where psi is a state formula");
		
		// model check operand first
		b = checkExpression(model, subformula).getBitSet();
		
		res = computeFG((STPG) model, b, min1, min2);
		
		StateValues probs = StateValues.createFromDoubleArray(res.soln, model);
		return probs;
	}
	
	/**
	 * Compute probabilities for a bounded until operator.
	 * @param bound The bound of the probability, can be specified if the caller
	 * is only interested whether a probability is greater/lower than {@code bound},
	 * this method can use it to avoid unnecessary computations.
	 */
	protected StateValues checkProbBoundedUntil(Model model, ExpressionTemporal expr, boolean min1, boolean min2) throws PrismException
	{
		int time;
		BitSet b1, b2;
		StateValues probs = null;
		ModelCheckerResult res = null;

		// get info from bounded until
		time = expr.getUpperBound().evaluateInt(constantValues);
		if (expr.upperBoundIsStrict())
			time--;
		if (time < 0) {
			String bound = expr.upperBoundIsStrict() ? "<" + (time + 1) : "<=" + time;
			throw new PrismException("Invalid bound " + bound + " in bounded until formula");
		}

		// model check operands first
		b1 = checkExpression(model, expr.getOperand1()).getBitSet();
		b2 = checkExpression(model, expr.getOperand2()).getBitSet();

		// print out some info about num states
		// mainLog.print("\nb1 = " + JDD.GetNumMintermsString(b1,
		// allDDRowVars.n()));
		// mainLog.print(" states, b2 = " + JDD.GetNumMintermsString(b2,
		// allDDRowVars.n()) + " states\n");

		// Compute probabilities

		// a trivial case: "U<=0"
		if (time == 0) {
			// prob is 1 in b2 states, 0 otherwise
			probs = StateValues.createFromBitSetAsDoubles(b2, model);
		} else {
			res = computeBoundedUntilProbs((STPG) model, b1, b2, time, min1, min2);
			probs = StateValues.createFromDoubleArray(res.soln, model);
		}

		return probs;
	}

	/**
	 * Compute probabilities for an (unbounded) until operator.
	 */
	protected StateValues checkProbUntil(Model model, ExpressionTemporal expr, boolean min1, boolean min2, double bound) throws PrismException
	{
		BitSet b1, b2;
		StateValues probs = null;
		ModelCheckerResult res = null;

		// model check operands first
		b1 = checkExpression(model, expr.getOperand1()).getBitSet();
		b2 = checkExpression(model, expr.getOperand2()).getBitSet();

		// print out some info about num states
		// mainLog.print("\nb1 = " + JDD.GetNumMintermsString(b1,
		// allDDRowVars.n()));
		// mainLog.print(" states, b2 = " + JDD.GetNumMintermsString(b2,
		// allDDRowVars.n()) + " states\n");

		res = computeUntilProbs((STPG) model, b1, b2, min1, min2, bound);
		probs = StateValues.createFromDoubleArray(res.soln, model);

		return probs;
	}

	/**
	 * Compute rewards for the contents of an R operator.
	 */
	protected StateValues checkRewardFormula(Model model, STPGRewards modelRewards, Expression expr, boolean min1, boolean min2) throws PrismException
	{
		StateValues rewards = null;
		
		if (expr instanceof ExpressionTemporal) {
			ExpressionTemporal exprTemp = (ExpressionTemporal) expr;
			switch (exprTemp.getOperator()) {
			case ExpressionTemporal.R_F:
				rewards = checkRewardReach(model, modelRewards, exprTemp, min1, min2);
				break;
			default:
				throw new PrismException("Explicit engine does not yet handle the " + exprTemp.getOperatorSymbol() + " operator in the R operator");
			}
		}
		
		if (rewards == null)
			throw new PrismException("Unrecognised operator in R operator");

		return rewards;
	}

	/**
	 * Compute rewards for a reachability reward operator, where the runs not reaching the final state get infinity.
	 */
	protected StateValues checkRewardReach(Model model, STPGRewards modelRewards, ExpressionTemporal expr, boolean min1, boolean min2) throws PrismException
	{
		return checkRewardReach(model, modelRewards, expr, min1, min2, true);
	}

	/**
	 * Compute rewards for a reachability reward operator.
	 * @param unreachingAsInfinity If true, the runs not reaching the final state get infinity, if false, these get cumulative reward.
	 */
	protected StateValues checkRewardReach(Model model, STPGRewards modelRewards, ExpressionTemporal expr, boolean min1, boolean min2, boolean unreachingAsInfinity) throws PrismException
	{
		BitSet b;
		StateValues rewards = null;
		ModelCheckerResult res = null;

		// model check operand first
		b = checkExpression(model, expr.getOperand2()).getBitSet();

		// print out some info about num states
		// mainLog.print("\nb = " + JDD.GetNumMintermsString(b1,
		// allDDRowVars.n()));

		res = computeReachRewards((STPG) model, modelRewards, b, min1, min2, null, null, unreachingAsInfinity);
		rewards = StateValues.createFromDoubleArray(res.soln, model);

		return rewards;
	}
	
	protected StateValues checkRewardReachZero(Model model, STPGRewards modelRewards, ExpressionTemporal expr, boolean min1, boolean min2) throws PrismException
	{
		BitSet b;
		StateValues rewards = null;
		ModelCheckerResult res = null;

		// model check operand first
		b = checkExpression(model, expr.getOperand2()).getBitSet();

		// print out some info about num states
		// mainLog.print("\nb = " + JDD.GetNumMintermsString(b1,
		// allDDRowVars.n()));

		res = computeReachRewardsZero((STPG) model, modelRewards, b, min1, min2, null, null);
		rewards = StateValues.createFromDoubleArray(res.soln, model);

		return rewards;
	}
	
	// Numerical computation functions

	/**
	 * Compute reachability probabilities.
	 * i.e. compute the min/max probability of reaching a state in {@code target}.
	 * @param stpg The STPG
	 * @param target Target states
	 * @param min1 Min or max probabilities for player 1 (true=lower bound, false=upper bound)
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
	 * @param min1 Min or max probabilities for player 1 (true=lower bound, false=upper bound)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 */
	public ModelCheckerResult computeReachProbs(STPG stpg, BitSet target, boolean min1, boolean min2, double bound) throws PrismException
	{
		return computeReachProbs(stpg, null, target, min1, min2, null, null, bound);
	}

	/**
	 * Compute until probabilities.
	 * i.e. compute the min/max probability of reaching a state in {@code target},
	 * while remaining in those in @{code remain}.
	 * @param stpg The STPG
	 * @param remain Remain in these states (optional: null means "all")
	 * @param target Target states
	 * @param min1 Min or max probabilities for player 1 (true=lower bound, false=upper bound)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 */
	public ModelCheckerResult computeUntilProbs(STPG stpg, BitSet remain, BitSet target, boolean min1, boolean min2, double bound) throws PrismException
	{
		return computeReachProbs(stpg, remain, target, min1, min2, null, null, bound);
	}

	/**
	 * Compute reachability/until probabilities.
	 * i.e. compute the min/max probability of reaching a state in {@code target},
	 * while remaining in those in @{code remain}.
	 * @param stpg The STPG
	 * @param remain Remain in these states (optional: null means "all")
	 * @param target Target states
	 * @param min1 Min or max probabilities for player 1 (true=lower bound, false=upper bound)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 * @param init Optionally, an initial solution vector (may be overwritten) 
	 * @param known Optionally, a set of states for which the exact answer is known
	 * Note: if 'known' is specified (i.e. is non-null, 'init' must also be given and is used for the exact values.  
	 */
	public ModelCheckerResult computeReachProbs(STPG stpg, BitSet remain, BitSet target, boolean min1, boolean min2, double init[], BitSet known)
			throws PrismException
	{
		// TODO: clean this up
		return computeReachProbs(stpg, remain, target, min1, min2, init, known, -1);
	}
	
	/**
	 * Compute reachability/until probabilities.
	 * i.e. compute the min/max probability of reaching a state in {@code target},
	 * while remaining in those in @{code remain}.
	 * @param stpg The STPG
	 * @param remain Remain in these states (optional: null means "all")
	 * @param target Target states
	 * @param min1 Min or max probabilities for player 1 (true=lower bound, false=upper bound)
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
		int i, n, numYes, numNo;
		long timer, timerProb0, timerProb1;
		boolean genAdv;

		// Check for some unsupported combinations
		if (solnMethod == SolnMethod.VALUE_ITERATION && valIterDir == ValIterDir.ABOVE && !(precomp && prob0)) {
			throw new PrismException("Precomputation (Prob0) must be enabled for value iteration from above");
		}

		// Are we generating an optimal adversary?
		genAdv = exportAdv;

		// Start probabilistic reachability
		timer = System.currentTimeMillis();
		if (verbosity >= 1)
			mainLog.println("Starting probabilistic reachability...");

		// Check for deadlocks in non-target state (because breaks e.g. prob1)
		stpg.checkForDeadlocks(target);

		// Store num states
		n = stpg.getNumStates();

		// Optimise by enlarging target set (if more info is available)
		if (init != null && known != null) {
			BitSet targetNew = new BitSet(n);
			for (i = 0; i < n; i++) {
				targetNew.set(i, target.get(i) || (known.get(i) && init[i] == 1.0));
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
		//do value iteration only if the values needed wasn't handled by precomputation
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
			res.soln = new double[n];
			for(int k = 0; k < n; k++)
				res.soln[k] = (yes.get(k)) ? 1.0 : 0.0;
			System.out.println("Bound is 1, hence I am skipping the computation of other values than 1.");
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
	 * Computes a probability that a formula (FG target) is satisfied.
	 * <p/>
	 * This method exploits the fact that this probability is
	 * equivalent to the probability of satisfying formula
	 * (F (P=1 [ G target ])).
	 * @param stpg The model.
	 * @param target The set of states satisfying the innermost subformula.
	 * @param min1 Min or max probabilities for player 1 (true=lower bound, false=upper bound)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 * @return
	 * @throws PrismException
	 */
	public ModelCheckerResult computeFG(STPG stpg, BitSet target, boolean min1, boolean min2) throws PrismException
	{
		//compute G=1
		// we have G=1 target
		// iff
		// F=0 not target
		
		int n = stpg.getNumStates();
		target.flip(0, n);
		//System.out.println(target);
		BitSet g1 = prob0(stpg, null, target, !min1, !min2);
		
		//System.out.println(g1);
		//g1.flip(0,n);
		
		//do reachability
		return computeReachProbs(stpg, g1, min1, min2);
	}

	/**
	 * Prob0 precomputation algorithm.
	 * i.e. determine the states of an STPG which, with min/max probability 0,
	 * reach a state in {@code target}, while remaining in those in @{code remain}.
	 * {@code min}=true gives Prob0E, {@code min}=false gives Prob0A. 
	 * @param stpg The STPG
	 * @param remain Remain in these states (optional: null means "all")
	 * @param target Target states
	 * @param min1 Min or max probabilities for player 1 (true=lower bound, false=upper bound)
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
	 * reach a state in {@code target}, while remaining in those in @{code remain}.
	 * @param stpg The STPG
	 * @param remain Remain in these states (optional: null means "all")
	 * @param target Target states
	 * @param min1 Min or max probabilities for player 1 (true=lower bound, false=upper bound)
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
	 * Zero cummulative reward precomputation algorithm.
	 * i.e. determine the states of an STPG which, with probability 1 get min/max reward equal to 0.0 
	 * before (possibly) reaching a state in {@code target}, while remaining in those in {@code remain}.
	 * @param stpg The STPG
	 * @param remain Remain in these states (optional: null means "all")
	 * @param target Target states
	 * @param min1 Min or max probabilities for player 1 (true=lower bound, false=upper bound)
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
		if(target != null)
			unknown.andNot(target);
		if (remain != null)
			unknown.and(remain);
		
		//initialise the solution so that the forbidden states are penalised
		for(int i = 0; i < n; i++) {
			if (remain != null && !remain.get(i) && target != null && !target.get(i))
				soln1[i] = Double.POSITIVE_INFINITY;
		}
		

		// Nested fixed point loop
		iters = 0;
		done = false;
		while (!done) {
			iters++;
			//at every iter at least one state must go from zero to nonzero, hence we have
			//at most n iterations
			assert iters <= n+1;
			
			stpg.mvMultRewMinMax(soln1, rewards, min1, min2, soln2, unknown, false, null);
			
			// Check termination (outer)
			done = true;
			
			double[] tmp = soln2;
			soln2 = soln1;
			soln1 = tmp;
			
			done = true;
			for(int i = 0; i < n; i++) {
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
		for(int i = 0; i < n; i++) {
			if (soln1[i] == 0.0)
				result.set(i);
		}
		
		return result;
	}

	/**
	 * Compute reachability probabilities using value iteration.
	 * @param stpg The STPG
	 * @param no Probability 0 states
	 * @param yes Probability 1 states
	 * @param min1 Min or max probabilities for player 1 (true=lower bound, false=upper bound)
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
		genAdv = exportAdv;

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

		// Print adversary
		if (genAdv) {
			PrismLog out = new PrismFileLog(exportAdvFilename);
			for (i = 0; i < n; i++) {
				out.println(i + " " + (adv[i] != -1 ? stpg.getAction(i, adv[i]) : "-"));
			}
			out.println();
		}

		// Return results
		res = new ModelCheckerResult();
		res.soln = soln;
		res.numIters = iters;
		res.timeTaken = timer / 1000.0;
		return res;
	}

	/**
	 * Compute reachability probabilities using Gauss-Seidel.
	 * @param stpg The STPG
	 * @param no Probability 0 states
	 * @param yes Probability 1 states
	 * @param min1 Min or max probabilities for player 1 (true=lower bound, false=upper bound)
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
		double soln[], initVal, maxDiff;
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

		// Return results
		res = new ModelCheckerResult();
		res.soln = soln;
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
	 * @param min1 Min or max probabilities for player 1 (true=lower bound, false=upper bound)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 */
	public ModelCheckerResult computeBoundedReachProbs(STPG stpg, BitSet target, int k, boolean min1, boolean min2) throws PrismException
	{
		return computeBoundedReachProbs(stpg, null, target, k, min1, min2, null, null);
	}

	/**
	 * Compute bounded until probabilities.
	 * i.e. compute the min/max probability of reaching a state in {@code target},
	 * within k steps, and while remaining in states in @{code remain}.
	 * @param stpg The STPG
	 * @param remain Remain in these states (optional: null means "all")
	 * @param target Target states
	 * @param k Bound
	 * @param min1 Min or max probabilities for player 1 (true=lower bound, false=upper bound)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 */
	public ModelCheckerResult computeBoundedUntilProbs(STPG stpg, BitSet remain, BitSet target, int k, boolean min1, boolean min2) throws PrismException
	{
		return computeBoundedReachProbs(stpg, remain, target, k, min1, min2, null, null);
	}

	/**
	 * Compute bounded reachability/until probabilities.
	 * i.e. compute the min/max probability of reaching a state in {@code target},
	 * within k steps, and while remaining in states in @{code remain}.
	 * @param stpg The STPG
	 * @param remain Remain in these states (optional: null means "all")
	 * @param target Target states
	 * @param k Bound
	 * @param min1 Min or max probabilities for player 1 (true=lower bound, false=upper bound)
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

		// Start bounded probabilistic reachability
		timer = System.currentTimeMillis();
		if (verbosity >= 1)
			mainLog.println("Starting bounded probabilistic reachability...");

		// Store num states
		n = stpg.getNumStates();

		// Create solution vector(s)
		soln = new double[n];
		soln2 = (init == null) ? new double[n] : init;

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
			// Matrix-vector multiply and min/max ops
			stpg.mvMultMinMax(soln, min1, min2, soln2, target, true, null);
			// Store intermediate results if required
			// (compute min/max value over initial states for this step)
			if (results != null) {
				results[iters] = Utils.minMaxOverArraySubset(soln2, stpg.getInitialStates(), min2);
			}
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

		// Return results
		res = new ModelCheckerResult();
		res.soln = soln;
		res.lastSoln = soln2;
		res.numIters = iters;
		res.timeTaken = timer / 1000.0;
		res.timePre = 0.0;
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
	 * Compute expected reachability rewards, where the runs that don't reach the final state get infinity.
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
	public ModelCheckerResult computeReachRewards(STPG stpg, STPGRewards rewards, BitSet target, boolean min1, boolean min2, double init[], BitSet known) throws PrismException
	{
		return computeReachRewards(stpg, rewards, target, min1, min2, init, known, true);
	}
	
	/**
	 * Compute expected reachability rewards, where the runs that don't reach the final state get infinity.
	 * i.e. compute the min/max reward accumulated to reach a state in {@code target}.
	 * @param stpg The STPG
	 * @param rewards The rewards
	 * @param target Target states
	 * @param min1 Min or max rewards for player 1 (true=min, false=max)
	 * @param min2 Min or max rewards for player 2 (true=min, false=max)
	 * @param init Optionally, an initial solution vector (may be overwritten) 
	 * @param known Optionally, a set of states for which the exact answer is known
	 * @param unreachingAsInfinity If true, the runs not reaching the final state get infinity, if false, these get cumulative reward.
	 * Note: if 'known' is specified (i.e. is non-null, 'init' must also be given and is used for the exact values.  
	 */
	public ModelCheckerResult computeReachRewards(STPG stpg, STPGRewards rewards, BitSet target, boolean min1, boolean min2, double init[], BitSet known, boolean unreachingAsInfinity) throws PrismException
	{
		
		ModelCheckerResult res = null;
		BitSet inf;
		int i, n, numTarget, numInf;
		long timer, timerProb1, timerApprox;

		// Start expected reachability
		timer = System.currentTimeMillis();
		if (verbosity >= 1)
			mainLog.println("Starting expected reachability...");

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
		//identify infinite values
		if(unreachingAsInfinity) {
			inf = prob1(stpg, null, target, !min1, !min2);
			inf.flip(0, n);
		} else {
			//check that there are no transition rewards
			//TODO fix the precomputation alg for infinity so that it
			//allows transition rewards.
			//At the same time find which states have nonzero reward
			BitSet zeroRew = new BitSet();
			for (int s = 0; s < n; s++) {
				for (int t = 0; t < stpg.getNumChoices(s); t++) {
					if (rewards.getTransitionReward(s, t) > 0) {
						System.out.println (s + "," + t + "," + rewards.getTransitionReward(s, t));
						throw new PrismException("The Fc operator cannot currently work with transition rewards.");
					}
				}
				if (rewards.getStateReward(s) == 0)
					zeroRew.set(s);
			}
			
			inf = new BitSet();
			//TODO the following uses numeric computation, should be changed
			//to something that is purely discrete.
			ModelCheckerResult rm = computeFG(stpg, zeroRew, min1, min2);
			for (int s = 0; s < n; s++)
				if (rm.soln[s] < 1)
					inf.set(s);
		}
		
		timerProb1 = System.currentTimeMillis() - timerProb1;

		// Print results of precomputation
		numTarget = target.cardinality();
		numInf = inf.cardinality();
		if (verbosity >= 1)
			mainLog.println("target=" + numTarget + ", inf=" + numInf + ", rest=" + (n - (numTarget + numInf)));
		
		//Compute rewards with epsilon instead of zero. This is used to get the over-approximation
		//of the real result, which deals with the problem of staying in zero components for free
		//when infinity should be gained.
		if (unreachingAsInfinity) {
			//first, get the minimum nonzero reward and maximal reward, will be used as a basis for epsilon
			//also, check if by any chance all rewards are nonzero, then we don't need to precompute
			double minimumReward = Double.POSITIVE_INFINITY;
			double maximumReward = 0.0;
			boolean allNonzero = true;
			double r;
			for (i = 0; i < n; i++) {
				r = rewards.getStateReward(i);
				if (r > 0.0 && r < minimumReward)
					minimumReward = r;
				if (r > maximumReward)
					maximumReward = r;
				allNonzero = allNonzero && r > 0;
				
				for (int j = 0; j < stpg.getNumChoices(i); j++) {
					r = rewards.getTransitionReward(i,j);
					if (r > 0.0  && r < minimumReward)
						minimumReward = r;
					if (r > maximumReward)
						maximumReward = r;
					allNonzero = allNonzero && rewards.getTransitionReward(i,j) > 0;
					
					for (int k= 0; k < stpg.getNumNestedChoices(i, j); k++) {
						r = rewards.getNestedTransitionReward(i, j, k);
						if (r > 0.0 && r < minimumReward)
							minimumReward = r;
						if (r > maximumReward)
							maximumReward = r;
						allNonzero = allNonzero && r > 0;
					}
				}
			}
			
			if (!allNonzero) {
				timerApprox = System.currentTimeMillis();
				//A simple heuristic that gives small epsilon, but still is hopefully safe floating-point-wise
				double epsilon = Math.min(minimumReward, maximumReward * 0.01);;
				
				if (verbosity >= 1) {
					mainLog.println("Computing the upper bound where " + epsilon + " is used instead of 0.0");
				}
				
				//Modify the rewards
				double origZeroReplacement;
				if (rewards instanceof MDPRewardsSimple) {
					origZeroReplacement = ((MDPRewardsSimple) rewards).getZeroReplacement();
					((MDPRewardsSimple) rewards).setZeroReplacement(epsilon);
				} else {
					throw new PrismException("To compute expected reward I need to modify the reward structure. But I don't know how to modify" + rewards.getClass().getName());
				}
				
				//Compute the value when rewards are nonzero
				switch (solnMethod) {
				case VALUE_ITERATION:
					res = computeReachRewardsValIter(stpg, rewards, target, inf, min1, min2, init, known);
					break;
				default:
					throw new PrismException("Unknown STPG solution method " + solnMethod);
				}
				
				//Set the value iteration result to be the initial solution for the next part
				//in which "proper" zero rewards are used
				init = res.soln;
				
				//Return the rewards to the original state
				if (rewards instanceof MDPRewardsSimple) {
					((MDPRewardsSimple)rewards).setZeroReplacement(origZeroReplacement);
				}
				
				timerApprox = System.currentTimeMillis() - timerApprox;
				
				if (verbosity >= 1) {
					mainLog.println("Computed an over-approximation of the solution (in " + timerApprox / 1000 + " seconds), this will now be used to get the solution");
				}
			}
		}
		
		// Compute real rewards
		switch (solnMethod) {
		case VALUE_ITERATION:
			res = computeReachRewardsValIter(stpg, rewards, target, inf, min1, min2, init , known);
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
	protected ModelCheckerResult computeReachRewardsValIter(STPG stpg, STPGRewards rewards, BitSet target, BitSet inf, boolean min1, boolean min2, double init[], BitSet known)
			throws PrismException
	{
		ModelCheckerResult res;
		BitSet unknown;
		int i, n, iters;
		double soln[], soln2[], tmpsoln[];
		int adv[] = null;
		boolean genAdv, done;
		long timer;
		
		
		// Are we generating an optimal adversary?
		genAdv = exportAdv;

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

		// Create/initialise adversary storage
		if (genAdv) {
			adv = new int[n];
			for (i = 0; i < n; i++) {
				adv[i] = -1;
			}
		}
		
		// Start iterations
		iters = 0;
		done = false;
		while (!done && iters < maxIters) {
			//mainLog.println(soln);
			iters++;
			// Matrix-vector multiply and min/max ops
			stpg.mvMultRewMinMax(soln, rewards, min1, min2, soln2, unknown, false, genAdv ? adv : null);
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
			/*Iterator<Entry<Integer, Double>> it;
			PrismLog out = new PrismFileLog(settings.getString(PrismSettings.PRISM_EXPORT_ADV_FILENAME));
			//out.print("Adv:");
			out.println(n + " ?");
			for (i = 0; i < n; i++) {
				if (adv[i] == -1)
					continue;
				//out.print(" " + i + ":" + stpg.getStatesList().get(i) + ":");
				//out.println(adv[i] != -1 ? stpg.getAction(i, adv[i]) : "-");
				
				it = stpg.getTransitionsIterator(i, adv[i]);
				if (it == null)
					continue;
				while (it.hasNext()) {
					Entry<Integer, Double> next = it.next();
					out.print(i + " 0 " + next.getKey() + " " + next.getValue() + " ");
					out.println(adv[i] != -1 ? stpg.getAction(i, adv[i]) : "-");
				}
			}
			out.println();*/
		}

		// Return results
		res = new ModelCheckerResult();
		res.soln = soln;
		res.numIters = iters;
		res.timeTaken = timer / 1000.0;
		return res;
	}
	
	/**
	 * Computes the reachability reward under the semantics where nonreaching runs get 0.
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
	public ModelCheckerResult computeReachRewardsZero(STPG stpg, STPGRewards rewards, BitSet target, boolean min1, boolean min2, double init[], BitSet known) throws PrismException
	{
		ModelCheckerResult res = null;
		BitSet inf;
		int i, n, numTarget, numInf;
		long timer, timerProb1, timerApprox;
		
		// Start expected reachability
		timer = System.currentTimeMillis();
		if (verbosity >= 1)
			mainLog.println("Starting expected reachability...");

		// Check for deadlocks in non-target state (because breaks e.g. prob1)
		stpg.checkForDeadlocks(target);

		// Store num states
		n = stpg.getNumStates();
				
		//Currently we only allow integer rewards, check if all rewards are (close to) an integer.
		//While traversing rewards, get largest reward, too
		boolean hasNonInt = false;
		double nonInt = 0; //will be used for output if there is non-integer reward
		int maxReward = 0;
		checkrewards: for (int s = 0; s < n; s++) {
			double sr = rewards.getStateReward(s);
			//System.out.println(sr);
			if (sr != Math.floor(sr)) {
				hasNonInt = true;
				nonInt = sr;
				break;
			}
			
			if (sr > maxReward)
				maxReward = (int) sr;
			
			for (int c = 0; c < stpg.getNumChoices(s); c++) { 
				double tr = rewards.getTransitionReward(s, c);
				//System.out.println(tr);
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
			throw new PrismException("For 'zero' semantics reachability reward all rewards must be integers." +
					"There is at least one non-integer reward: " + nonInt);

		// Optimise by enlarging target set (if more info is available)
		if (init != null && known != null) {
			BitSet targetNew = new BitSet(n);
			for (i = 0; i < n; i++) {
				targetNew.set(i, target.get(i) || (known.get(i) && init[i] == 0.0));
			}
			target = targetNew;
		}
		
		//TODO identify "dead" states, i.e. those from which F can't be reached with >0 prob
		//and those from which the bad player can ensure 0.	This is optional, but
		//should bring some speedup.

		
		BitSet zeroProb = prob0(stpg, null, target, min1, min2);
		
		BitSet positiveProb = new BitSet();
		for(int k = 0; k < n; k++)
			positiveProb.set(k, !zeroProb.get(k));
		
		
		//...and those from which the bad player can ensure 0.	
		//BitSet zeroReward = zeroRewards(stpg, rewards, positiveProb, target, !min1, !min2);
		
		//Identify states that get infinity.
		ModelCheckerResult mcri = computeReachRewards(stpg, rewards, target, min1, min2, new double[n], zeroProb, false);
		BitSet infinity = new BitSet();
		
		for (int k = 0; k < n; k++)
			infinity.set(k, mcri.soln[k] == Double.POSITIVE_INFINITY);
		
		//Get the rich man's strategy and its values
		//Start with computing optimal probabilities to reach the final state
		ModelCheckerResult mcrprob = computeReachProbs(stpg, target, min1, min2);
		
		//System.out.println("maximal probs:" + Arrays.toString(mcrprob.soln));
		
		//Next, reweigh the rewards and make sure that only optimal actions are taken 
		STPGRewards rewardsRestricted;
		if (rewards instanceof MDPRewardsSimple) {
			//And make sure only the best actions are used
			STPGRewardsSimple rewardsRestrictedSimple = new STPGRewardsSimple((MDPRewardsSimple) rewards);
			
			for (int s = 0; s < n; s++) {
				for (int c = 0; c < stpg.getNumChoices(s); c++)	{
					double prob = 0.0;
					Iterator<Entry<Integer, Double>> it = stpg.getTransitionsIterator(s, c);
					while(it.hasNext()) {
						Entry<Integer, Double> e = it.next();
						prob += e.getValue() * mcrprob.soln[e.getKey()];
					}
					
					//as a hack, set the transition reward of nonoptimal transitions
					//to something extreme so they are never chosen
					if (prob < mcrprob.soln[s] && ((stpg.getPlayer(s) == 1 && !min1) || (stpg.getPlayer(s) == 2 && !min2))) {
						rewardsRestrictedSimple.setTransitionReward(s, c, Double.NEGATIVE_INFINITY);
						//System.out.println("choice " + s + " " + c + " changing " + rewards.getTransitionReward(s, c) + " to -inf");
					} else if (prob > mcrprob.soln[s] && ((stpg.getPlayer(s) == 1 && min1) || (stpg.getPlayer(s) == 2 && min2))) {
						rewardsRestrictedSimple.setTransitionReward(s, c, Double.POSITIVE_INFINITY);
						//System.out.println("choice " + s + " " + c + " changing " + rewards.getTransitionReward(s, c) + " to +inf");
					} else {
						double newReward = rewards.getTransitionReward(s, c) * mcrprob.soln[s];
						rewardsRestrictedSimple.setTransitionReward(s, c, newReward);
						//System.out.println("choice " + s + " " + c + " changing " + rewards.getTransitionReward(s, c) + " to " + newReward);
					}
					//System.out.println(" " + prob + " " + mcrprob.soln[s] + " " +stpg.getPlayer(s) + " " + min1);
				}
				double newReward = rewards.getStateReward(s) * mcrprob.soln[s];
				rewardsRestrictedSimple.setStateReward(s, newReward);
			}
			rewardsRestricted = rewardsRestrictedSimple;
		} else {
			throw new PrismException("To compute expected reward I need to modify the reward structure. But I don't know how to modify" + rewards.getClass().getName());
		}		
		
		//Next, compute the value for the rich man's strategy.
		ModelCheckerResult mcrrich = computeReachRewards(stpg, rewardsRestricted, target, min1, min2, init, known, false);
		//System.out.println("maximal rews for rich man's strategy: " + Arrays.toString(mcrrich.soln));
		
		//compute B from the values for the rich man's strategy
		int lastSwitch = 0;
		for (int s = 0; s < n; s++) {
			//for all choices c, find maximal B such that
			//sum_{s'} prob(c,s')(r(s) + r(c) + B + rewRich(s'))  > probRich(s)*B + rewRich(s)
			for (int c = 0; c < stpg.getNumChoices(s); c++) {
				double numerator = mcrrich.soln[s];
				double denominator = - mcrprob.soln[s];
				double tRew = rewards.getTransitionReward(s, c);
				Iterator<Entry<Integer, Double>> it = stpg.getTransitionsIterator(s, c);
				while(it.hasNext()) {
					Entry<Integer, Double> e = it.next();
					int ts = e.getKey();
					double sRew = rewards.getStateReward(s);
					double p = e.getValue();
					
					numerator -= p*(mcrprob.soln[ts]*(sRew + tRew) + mcrrich.soln[ts]);
					denominator += p*mcrprob.soln[ts];
					
					int b = (int) Math.floor(numerator/denominator);
					//System.out.println(s + " " + c + " " + b);
					if (lastSwitch < b)
						lastSwitch = b;
				}
			}
		}
		
		if (verbosity >= 1)
			mainLog.println("Last switching point is when the reward cumulated in the past becomes " + lastSwitch);
		
		//TODO using gcd of rewards could save us iterations in many cases
		double[][] rews = new double[maxReward + 1][n];
		double[] tmp;
		
		//fill in the initial values from the rich man's strategy
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < maxReward + 1; k++) {
				rews[k][j] = mcrrich.soln[j] + (lastSwitch + k)*mcrprob.soln[j];
			}
		}
		
		//System.out.println(Arrays.toString(rews[0]));
		//System.out.println(Arrays.deepToString(rews));
		//System.out.println();
		
		
		int iters = 0;
		for(int x = lastSwitch; x >=0; x--) {		
			//reward[s,x] =
			// opt_c sum_{s'} p(s ->c s')(prob_F[s,x+r(s)+r(c)]*(r(s)+r(c)) + reward[s'][x+r(s)+r(c)])
			//where opt is either max or min, depending on the owner of s
			//probs[s,x] is the probability of reaching F under the choice c chosen by rews
			boolean done = false;
			
			double difference;
			do {
				difference = 0;
			
				iters++;
				for (int s = 0; s < n; s++) {
					if (target.get(s))
					{
						rews[0][s] = x;
					}
				}
				
				for (int s = 0; s < n; s++) {
					if (target.get(s))
					{
						continue;
					}
					
					//non-target states
					boolean min = (stpg.getPlayer(s) == 1) ? min1 : min2;
					double stateRew = 0;
					for (int c = 0;  c < stpg.getNumChoices(s); c++) { 
						double choiceRew = 0;
						double r = rewards.getStateReward(s) + rewards.getTransitionReward(s, c);
						int index = (int) r; //the reward determines in which array we will look
						
						Iterator<Entry<Integer, Double>> it = stpg.getTransitionsIterator(s, c);
						while(it.hasNext()) {
							Entry<Integer, Double> e = it.next();
							int ts = e.getKey();
							double p = e.getValue();
							//choiceRew += p*(probs[index][ts]*r + rews[index][ts]);
							choiceRew += p*(rews[index][ts]);
						}
						
						//System.out.println(s + " " + c + " " + choiceRew);
						
						if (min && stateRew > choiceRew) {
							stateRew = choiceRew;
						}
						else if (!min && stateRew < choiceRew) {
							stateRew = choiceRew;
						}
						
					}
					double cDif = Math.abs(rews[0][s] - stateRew);
					if (cDif > difference)
						difference = cDif;
					rews[0][s] = stateRew;
				}
				
				//System.out.println(x + ": " + Arrays.toString(rews[0]));
				//System.out.println(Arrays.deepToString(rews));
				//System.out.println();
			} while (difference > 10e-7); //TODO some smarter convergence test
			
			//shift the array 
			double[] tmpRews = rews[maxReward];
			for (i = maxReward; i >= 1; i--) {
				rews[i] = rews[i-1];
			}
			rews[0] = tmpRews;
		}
		
		timer = System.currentTimeMillis() - timer;
		
		res = new ModelCheckerResult();
		res.soln = rews[1];
		res.lastSoln = (rews.length > 2) ? rews[2] : null;
		res.numIters = lastSwitch;
		res.timeTaken = timer / 1000;
		res.numIters = iters;
		return res;
	}
	
	/**
	 * Simple test program.
	 */
	public static void main(String args[])
	{
		STPGModelChecker mc;
		STPGAbstrSimple stpg;
		ModelCheckerResult res;
		BitSet target;
		Map<String, BitSet> labels;
		boolean min1 = true, min2 = true;
		try {
			mc = new STPGModelChecker();
			stpg = new STPGAbstrSimple();
			stpg.buildFromPrismExplicit(args[0]);
			//System.out.println(stpg);
			labels = mc.loadLabelsFile(args[1]);
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
			//stpg.exportToDotFile("stpg.dot", target);
			//stpg.exportToPrismExplicit("stpg");
			res = mc.computeReachProbs(stpg, target, min1, min2);
			System.out.println(res.soln[0]);
		} catch (PrismException e) {
			System.out.println(e);
		}
	}
}
