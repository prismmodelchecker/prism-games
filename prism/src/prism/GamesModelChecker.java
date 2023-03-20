//==============================================================================
//	
//	Copyright (c) 2016-
//	Authors:
//	* Gabriel Santos <gabriel.santos@cs.ox.uk> (University of Oxford)
//	* Dave Parker <david.parker@cs.ox.ac.uk> (University of Oxford)
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

package prism;

import java.util.List;

import explicit.MinMax;
import jdd.JDD;
import jdd.JDDNode;
import jdd.JDDVars;
import parser.ast.Coalition;
import parser.ast.Expression;
import parser.ast.ExpressionFunc;
import parser.ast.ExpressionProb;
import parser.ast.ExpressionReward;
import parser.ast.ExpressionStrategy;
import parser.ast.ExpressionTemporal;
import parser.ast.ExpressionUnaryOp;
import parser.ast.PropertiesFile;
import parser.type.TypeBool;
import parser.type.TypePathBool;
import parser.type.TypePathDouble;

public class GamesModelChecker extends NonProbModelChecker
{
	// Model (SMG)
	protected GamesModel model;

	// Extra (SMG) model info
	protected JDDNode nondetMask;
	protected JDDVars allDDNondetVars;
	protected JDDVars allDDPlayerVars;

	/** BDD over player DD vars encoding players in coalition (i.e. P1) */
	protected JDDNode ddc1;
	/** BDD over player DD vars encoding players not in coalition (i.e. P2) */
	protected JDDNode ddc2;

	// Options (in addition to those inherited from StateModelChecker):

	// Use 0,1 precomputation algorithms?
	// if 'precomp' is false, this disables all (non-essential) use of prob0/prob1
	// if 'precomp' is true, the values of prob0/prob1 determine what is used
	// (currently prob0/prob are not under user control)
	protected boolean precomp;
	protected boolean prob0;
	protected boolean prob1;
	protected int maxIters;

	// Constructor

	public GamesModelChecker(Prism prism, Model m, PropertiesFile pf) throws PrismException
	{
		// Initialise
		super(prism, m, pf);
		if (!(m instanceof GamesModel)) {
			throw new PrismException("Wrong model type passed to GamesModelChecker.");
		}
		model = (GamesModel) m;
		nondetMask = model.getNondetMask();
		allDDNondetVars = model.getAllDDNondetVars();
		allDDPlayerVars = model.getAllDDPlayerVars();
		
		// Inherit some options from parent Prism object and store locally.
		precomp = prism.getPrecomp();
		prob0 = prism.getProb0();
		prob1 = prism.getProb1();
		maxIters = prism.getMaxIters();
	}

	// Model checking functions

	@Override
	public StateValues checkExpression(Expression expr, JDDNode statesOfInterest) throws PrismException
	{
		StateValues res;

		// <<>> or [[]] operator
		if (expr instanceof ExpressionStrategy) {
			res = checkExpressionStrategy((ExpressionStrategy) expr, statesOfInterest);
		}
		// P operator
		else if (expr instanceof ExpressionProb) {
			res = checkExpressionProb((ExpressionProb) expr, statesOfInterest);
		}
		// R operator
		else if (expr instanceof ExpressionReward) {
			res = checkExpressionReward((ExpressionReward) expr, statesOfInterest);
		}
		// Multi-objective
		else if (expr instanceof ExpressionFunc) {
			// Detect "multi" function
			if (((ExpressionFunc) expr).getName().equals("multi")) {
				throw new PrismNotSupportedException("Not currently supported by the MTBDD engine");
			}
			// For any other function, check as normal
			else {
				res = super.checkExpression(expr, statesOfInterest);
			}
		}
		// Otherwise, use the superclass
		else {
			res = super.checkExpression(expr, statesOfInterest);
		}

		// Filter out non-reachable states from solution
		// (only necessary for symbolically stored vectors)
		if (res instanceof StateValuesMTBDD)
			res.filter(reach);

		return res;
	}

	/**
	 * Model check a <<>> or [[]] operator expression.
	 * The result will have valid results at least for the states of interest (use model.getReach().copy() for all reachable states)
	 * <br>[ REFS: <i>result</i>, DEREFS: statesOfInterest ]
	 */
	protected StateValues checkExpressionStrategy(ExpressionStrategy expr, JDDNode statesOfInterest) throws PrismException
	{
		// Will we be quantifying universally or existentially over strategies/adversaries?
		boolean forAll = !expr.isThereExists();

		// Extract coalition info
		Coalition coalition = expr.getCoalition();
		// Deal with the coalition operator here and then remove it
		if (coalition != null && !model.getModelType().multiplePlayers()) {
			if (coalition.isEmpty()) {
				// An empty coalition negates the quantification ("*" has no effect)
				forAll = !forAll;
			}
			coalition = null;
		}

		// Process operand(s)
		List<Expression> exprs = expr.getOperands();
		// Pass onto relevant method:
		// Single P operator
		if (exprs.size() == 1 && exprs.get(0) instanceof ExpressionProb) {
			return checkExpressionProb((ExpressionProb) exprs.get(0), forAll, statesOfInterest, coalition);
		}
		// Single R operator
		else if (exprs.size() == 1 && exprs.get(0) instanceof ExpressionReward) {
			return checkExpressionReward((ExpressionReward) exprs.get(0), forAll, statesOfInterest, coalition);
		}
		// Equilibria
		/*else if (exprs.size() == 1 && exprs.get(0) instanceof ExpressionMultiNash) {
			return checkExpressionNash((ExpressionMultiNash) exprs.get(0), expr.getCoalitions());
		}*/
		// Anything else is treated as multi-objective 
		else {
			throw new PrismNotSupportedException("Not currently supported by the MTBDD engine");
		}
	}

	/**
	 * Model check a P operator expression.
	 * The result will have valid results at least for the states of interest (use model.getReach().copy() for all reachable states)
	 * <br>[ REFS: <i>result</i>, DEREFS: statesOfInterest ]
	 */
	protected StateValues checkExpressionProb(ExpressionProb expr, JDDNode statesOfInterest) throws PrismException
	{
		// Use the default semantics for a standalone P operator
		// (i.e. quantification over all strategies)
		Coalition coalitionAll = new Coalition();
		coalitionAll.setAllPlayers();
		return checkExpressionProb(expr, true, statesOfInterest, coalitionAll);
	}

	/**
	 * Model check a P operator expression.
	 * The result will have valid results at least for the states of interest (use model.getReach().copy() for all reachable states)
	 * <br>[ REFS: <i>result</i>, DEREFS: statesOfInterest ]
	 *
	 * @param expr The P operator expression
	 * @param forAll Are we checking "for all strategies" (true) or "there exists a strategy" (false)? [irrelevant for numerical (=?) queries]
	 */
	protected StateValues checkExpressionProb(ExpressionProb expr, boolean forAll, JDDNode statesOfInterest, Coalition coalition) throws PrismException
	{
		// Get info from P operator
		OpRelOpBound opInfo = expr.getRelopBoundInfo(constantValues);
		MinMax minMax = opInfo.getMinMax(model.getModelType(), forAll, coalition);

		// Check for trivial (i.e. stupid) cases
		if (opInfo.isTriviallyTrue()) {
			mainLog.printWarning("Checking for probability " + opInfo.relOpBoundString() + " - formula trivially satisfies all states");
			JDD.Ref(reach);
			JDD.Deref(statesOfInterest);
			return new StateValuesMTBDD(reach, model);
		} else if (opInfo.isTriviallyFalse()) {
			mainLog.printWarning("Checking for probability " + opInfo.relOpBoundString() + " - formula trivially satisfies no states");
			JDD.Deref(statesOfInterest);
			return new StateValuesMTBDD(JDD.Constant(0), model);
		}

		// Compute probabilities
		StateValues probs = checkProbPathFormula(expr.getExpression(), minMax, statesOfInterest);

		// Print out probabilities
		if (verbose) {
			mainLog.print("\n" + (minMax.isMin() ? "Minimum" : "Maximum") + " probabilities (non-zero only) for all states:\n");
			probs.print(mainLog);
		}

		// For =? properties, just return values
		if (opInfo.isNumeric()) {
			return probs;
		}
		// Otherwise, compare against bound to get set of satisfying states
		else {
			JDDNode sol = probs.getBDDFromInterval(opInfo.getRelOp(), opInfo.getBound());
			// Remove unreachable states from solution
			JDD.Ref(reach);
			sol = JDD.And(sol, reach);
			// Free vector
			probs.clear();
			return new StateValuesMTBDD(sol, model);
		}
	}

	/**
	 * Model check an R operator expression.
	 * The result will have valid results at least for the states of interest (use model.getReach().copy() for all reachable states)
	 * <br>[ REFS: <i>result</i>, DEREFS: statesOfInterest ]
	 */
	protected StateValues checkExpressionReward(ExpressionReward expr, JDDNode statesOfInterest) throws PrismException
	{
		// Use the default semantics for a standalone R operator
		// (i.e. quantification over all strategies)
		Coalition coalitionAll = new Coalition();
		coalitionAll.setAllPlayers();
		return checkExpressionReward(expr, true, statesOfInterest, coalitionAll);
	}

	/**
	 * Model check an R operator expression.
	 * The result will have valid results at least for the states of interest (use model.getReach().copy() for all reachable states)
	 * <br>[ REFS: <i>result</i>, DEREFS: statesOfInterest ]
	 *
	 * @param expr The R operator expression
	 * @param forAll Are we checking "for all strategies" (true) or "there exists a strategy" (false)? [irrelevant for numerical (=?) queries]
	 */
	protected StateValues checkExpressionReward(ExpressionReward expr, boolean forAll, JDDNode statesOfInterest, Coalition coalition) throws PrismException
	{
		// Get info from R operator
		if (expr.getRewardStructIndexDiv() != null)
			throw new PrismException("Ratio rewards not supported with the selected engine and module type.");
		OpRelOpBound opInfo = expr.getRelopBoundInfo(constantValues);
		MinMax minMax = opInfo.getMinMax(model.getModelType(), forAll, coalition);

		// Get rewards
		Object rs = expr.getRewardStructIndex();
		JDDNode stateRewards = getStateRewardsByIndexObject(rs, model, constantValues);
		JDDNode transRewards = getTransitionRewardsByIndexObject(rs, model, constantValues);

		// Compute rewards
		StateValues rewards = null;
		Expression expr2 = expr.getExpression();

		// Currently ignoring states of interest
		JDD.Deref(statesOfInterest);
		
		if (expr2.getType() instanceof TypePathDouble) {
			ExpressionTemporal exprTemp = (ExpressionTemporal) expr2;
			switch (exprTemp.getOperator()) {
			case ExpressionTemporal.R_C:
				if (!exprTemp.hasBounds()) {
					throw new PrismNotSupportedException("Not currently supported by the MTBDD engine");
				} else {
					int bound = exprTemp.getUpperBound().evaluateInt(constantValues);
					if (bound < 0) {
						throw new PrismException("Invalid bound " + bound + " in cumulative rewards property");
					}
					rewards = computeCumulativeRewards(stateRewards, transRewards, bound, minMax.isMin1(), minMax.isMin2(), minMax.getCoalition());
				}
				break;
			case ExpressionTemporal.R_Fc:
				throw new PrismNotSupportedException("Not currently supported by the MTBDD engine");
			case ExpressionTemporal.R_I:
				int time = exprTemp.getUpperBound().evaluateInt(constantValues);
				if (time < 0) {
					throw new PrismException("Invalid bound " + time + " in instantaneous rewards property");
				}
				rewards = computeInstantaneousRewards(stateRewards, time, minMax.isMin1(), minMax.isMin2(), minMax.getCoalition());
				break;
			default:
				throw new PrismNotSupportedException("Not currently supported by the MTBDD engine");
			}
		} else if (expr2.getType() instanceof TypePathBool || expr2.getType() instanceof TypeBool) {
			if (Expression.isReach(expr2)) {
				rewards = computeReachRewardsInfinity((ExpressionTemporal) expr2, stateRewards, transRewards, minMax.isMin1(), minMax.isMin2(),
						minMax.getCoalition());
			} else if (Expression.isCoSafeLTLSyntactic(expr, true)) {
				throw new PrismNotSupportedException("Cosafe rewards not currently supported by the MTBDD engine");
			} else {
				throw new PrismException("R operator contains a path formula that is not syntactically co-safe: " + expr2);
			}
		} else {
			throw new PrismNotSupportedException("Not currently supported by the MTBDD engine");
		}

		// Print out rewards
		if (verbose) {
			mainLog.print("\n" + (minMax.isMin() ? "Minimum" : "Maximum") + " rewards (non-zero only) for all states:\n");
			rewards.print(mainLog);
		}

		// For =? properties, just return values
		if (opInfo.isNumeric()) {
			return rewards;
		}
		// Otherwise, compare against bound to get set of satisfying states
		else {
			JDDNode sol = rewards.getBDDFromInterval(opInfo.getRelOp(), opInfo.getBound());
			// Remove unreachable states from solution
			JDD.Ref(reach);
			sol = JDD.And(sol, reach);
			// Free vector
			rewards.clear();
			return new StateValuesMTBDD(sol, model);
		}
	}

	/*protected StateValues checkExpressionNash(ExpressionMultiNash expr, List<Coalition> coalitions) throws PrismException
	{
		StateValues result;
		JDDNode values;

		List<ExpressionQuant> formulae;

		JDDNode[] remain = new JDDNode[2];
		JDDNode[] target = new JDDNode[2];
		JDDNode[] stateRewards = new JDDNode[2];
		JDDNode[] transRewards = new JDDNode[2];
		Expression expr1;
		IntegerBound bound;
		int[] bounds = new int[coalitions.size()];

		Object r;
		boolean type = true;

		// Assume only two coalitions and that 2nd is complement of 1st
		setCoalition(coalitions.get(0));

		formulae = expr.getOperands();

		for (int e = 0; e < formulae.size() - 1; e++) {
			type = type && formulae.get(e).getClass() == formulae.get(e + 1).getClass();
		}
		if (!type)
			throw new PrismException("Mixing P and R operators is not yet supported");

		for (int p = 0; p < formulae.size(); p++) {
			expr1 = formulae.get(p);
			if (expr1 instanceof ExpressionMultiNashProb) {
				expr1 = ((ExpressionMultiNashProb) (expr1)).getExpression();
				expr1 = Expression.convertSimplePathFormulaToCanonicalForm(expr1);
				if (expr1 instanceof ExpressionTemporal) {
					ExpressionTemporal exprTemp = (ExpressionTemporal) expr1;
					switch (exprTemp.getOperator()) {
					case ExpressionTemporal.P_U:
						if (exprTemp.hasBounds()) {
							throw new PrismNotSupportedException(
									"The probabilisitic operator " + exprTemp.getOperatorSymbol() + " is not yet supported for equilibria-based properties");

						}
						if (!((ExpressionTemporal) expr1).getOperand1().equals(ExpressionLiteral.True()))
							remain[p] = checkExpressionDD(((ExpressionTemporal) expr1).getOperand1(), model.getReach().copy());
						break;
					default:
						throw new PrismNotSupportedException(
								"The probabilisitic operator " + exprTemp.getOperatorSymbol() + " is not yet supported for equilibria-based properties");
					}
				}
			} else if (expr1 instanceof ExpressionMultiNashReward) {
				r = ((ExpressionReward) expr1).getRewardStructIndex();
				stateRewards[p] = getStateRewardsByIndexObject(r, model, constantValues);
				transRewards[p] = getTransitionRewardsByIndexObject(r, model, constantValues);

				mainLog.println("stateRewards " + p);
				//JDD.PrintMinterms(mainLog, stateRewards[p]);
				mainLog.println("transRewards " + p);
				//JDD.PrintMinterms(mainLog, transRewards[p]);

				expr1 = ((ExpressionMultiNashReward) (expr1)).getExpression();

				ExpressionTemporal exprTemp = (ExpressionTemporal) expr1;
				switch (exprTemp.getOperator()) {
				case ExpressionTemporal.P_F:
					target[p] = checkExpressionDD(((ExpressionTemporal) expr1).getOperand2(), model.getReach().copy());
					mainLog.println("target " + p);
					JDD.PrintMinterms(mainLog, target[p]);
					break;
				case ExpressionTemporal.R_I:
					bounds[p] = exprTemp.getUpperBound().evaluateInt(constantValues);
					break;
				case ExpressionTemporal.R_C:
					if (exprTemp.hasBounds()) {
						bound = IntegerBound.fromExpressionTemporal((ExpressionTemporal) expr1, constantValues, true);
						if (bound.hasUpperBound())
							bounds[p] = bound.getHighestInteger();
						else
							throw new PrismException("Only upper bounds are allowed");

					} else {
						throw new PrismNotSupportedException(
								"Total rewards " + exprTemp.getOperatorSymbol() + " is not yet supported for equilibria-based properties");
					}
					break;
				default:
					throw new PrismNotSupportedException(
							"The reward operator " + exprTemp.getOperatorSymbol() + " is not yet supported for equilibria-based properties");
				}
			}
		}

		values = computeReachEquilibria(target, stateRewards, transRewards, 3, false);
		setCoalition(null);
		result = new StateValuesMTBDD(values, model);

		return result;
	}*/

	/**
	 * Compute probabilities for the contents of a P operator.
	 * The result will have valid results at least for the states of interest (use model.getReach().copy() for all reachable states)
	 * <br>[ REFS: <i>result</i>, DEREFS: statesOfInterest ]
	 */
	protected StateValues checkProbPathFormula(Expression expr, MinMax minMax, JDDNode statesOfInterest) throws PrismException
	{
		// Currently ignoring states of interest
		JDD.Deref(statesOfInterest);
		
		// No support for reward-bounded path formulas (i.e. of the form R(path)~r)
		if (Expression.containsRewardBoundedPathFormula(expr)) {
			throw new PrismException("Reward-bounded path formulas not supported");
		}

		// Test whether this is a simple path formula
		boolean useSimplePathAlgo = expr.isSimplePathFormula();

		if (useSimplePathAlgo) {
			return checkProbPathFormulaSimple(expr, minMax);
		} else {
			throw new PrismNotSupportedException("Not currently supported by the MTBDD engine");
		}
	}

	/**
	 * Compute probabilities for a simple, non-LTL path operator.
	 */
	protected StateValues checkProbPathFormulaSimple(Expression expr, MinMax minMax) throws PrismException
	{
		boolean negated = false;
		StateValues probs = null;

		expr = Expression.convertSimplePathFormulaToCanonicalForm(expr);

		// Negation
		if (expr instanceof ExpressionUnaryOp && ((ExpressionUnaryOp) expr).getOperator() == ExpressionUnaryOp.NOT) {
			negated = true;
			minMax = minMax.negate();
			expr = ((ExpressionUnaryOp) expr).getOperand();
		}

		if (expr instanceof ExpressionTemporal) {
			ExpressionTemporal exprTemp = (ExpressionTemporal) expr;
			// Next
			if (exprTemp.getOperator() == ExpressionTemporal.P_X) {
				probs = checkProbNext(exprTemp, minMax);
			}
			// Until
			else if (exprTemp.getOperator() == ExpressionTemporal.P_U) {
				if (exprTemp.hasBounds()) {
					probs = checkProbBoundedUntil(exprTemp, minMax);

				} else {
					probs = checkProbUntil(exprTemp, minMax);
				}
			}
		}

		if (probs == null)
			throw new PrismException("Unrecognised path operator in P operator");

		if (negated) {
			// Subtract from 1 for negation
			probs.subtractFromOne();
		}

		return probs;
	}

	/**
	 * Compute probabilities for a next operator.
	 */
	protected StateValues checkProbNext(ExpressionTemporal expr, MinMax minMax) throws PrismException
	{
		JDDNode b;
		StateValues probs = null;

		b = checkExpressionDD(expr.getOperand2(), model.getReach().copy());
		probs = computeNextProbs(b, minMax.isMin1(), minMax.isMin2(), minMax.getCoalition());

		return probs;
	}

	/**
	 * Compute probabilities for a bounded until operator.
	 */
	protected StateValues checkProbBoundedUntil(ExpressionTemporal expr, MinMax minMax) throws PrismException
	{
		JDDNode b1, b2;
		StateValues probs = null;
		Integer lowerBound;
		IntegerBound bounds;

		// Get and check bounds information
		bounds = IntegerBound.fromExpressionTemporal(expr, constantValues, true);

		// Model check operands first
		b1 = checkExpressionDD(expr.getOperand1(), model.getReach().copy());

		try {
			b2 = checkExpressionDD(expr.getOperand2(), model.getReach().copy());
		} catch (PrismException e) {
			JDD.Deref(b1);
			throw e;
		}

		if (bounds.hasLowerBound()) {
			lowerBound = bounds.getLowestInteger();
		} else {
			lowerBound = 0;
		}

		Integer bound = null; // unbounded
		if (bounds.hasUpperBound()) {
			bound = bounds.getHighestInteger() - lowerBound;
		}

		// Compute probabilities for Until<=windowSize
		if (bound == null) {

			// Unbounded
			try {
				probs = computeUntilProbs(b1, b2, minMax.isMin1(), minMax.isMin2(), minMax.getCoalition());
			} catch (PrismException e) {
				JDD.Deref(b1);
				JDD.Deref(b2);
				throw e;
			}

		} else if (bound == 0) {

			// The trivial case: windowSize = 0
			// Prob is 1 in b2 states, 0 otherwise
			JDD.Ref(b2);
			probs = new StateValuesMTBDD(b2, model);
		} else {
			try {
				probs = computeBoundedUntilProbs(b1, b2, bound, minMax.isMin1(), minMax.isMin2(), minMax.getCoalition());
			} catch (PrismException e) {
				JDD.Deref(b1);
				JDD.Deref(b2);
				throw e;
			}
		}

		// Derefs
		JDD.Deref(b1);
		JDD.Deref(b2);

		return probs;
	}

	/**
	 * Compute probabilities for an (unbounded) until operator.
	 */
	protected StateValues checkProbUntil(ExpressionTemporal expr, MinMax minMax) throws PrismException
	{
		JDDNode b1, b2;
		StateValues result = null;

		// Model check operands first
		b1 = checkExpressionDD(expr.getOperand1(), model.getReach().copy());
		try {
			b2 = checkExpressionDD(expr.getOperand2(), model.getReach().copy());
		} catch (PrismException e) {
			JDD.Deref(b1);
			throw e;
		}
		try {
			result = computeUntilProbs(b1, b2, minMax.isMin1(), minMax.isMin2(), minMax.getCoalition());
		} catch (PrismException e) {
			JDD.Deref(b1);
			JDD.Deref(b2);
			throw e;
		}

		// Derefs
		JDD.Deref(b1);
		JDD.Deref(b2);

		return result;
	}

	// -----------------------------------------------------------------------------------
	// probability computation methods
	// -----------------------------------------------------------------------------------

	// precomputation algorithms

	/**
	 * <br>[ REFS: <i>result</i>; DEREFS: target, remain ]	 
	 * * @param remain
	 * @param target
	 * @return
	 */
	protected JDDNode prob1(JDDNode remain, JDDNode target, boolean min1, boolean min2)
	{
		JDDNode u, v, tmp, tmp1, tmp2;
		long timer;
		int iters = 0;
		boolean u_done = false;

		// Start precomputation
		mainLog.println("Starting Prob1...");
		timer = System.currentTimeMillis();

		// Greatest fixed point
		u = reach.copy();
		while (!u_done) {
			boolean v_done = false;
			// Least fixed point
			v = JDD.Constant(0);
			while (!v_done) {
				iters++;
				// Find all-to-u, some-to-v
				tmp1 = JDD.SwapVariables(u.copy(), allDDRowVars, allDDColVars);
				tmp1 = JDD.ForAll(JDD.Implies(trans01.copy(), tmp1), allDDColVars);
				tmp2 = JDD.SwapVariables(v.copy(), allDDRowVars, allDDColVars);
				tmp2 = JDD.ThereExists(JDD.And(tmp2, trans01.copy()), allDDColVars);
				tmp = JDD.And(tmp1, tmp2);
				// Exists/forall
				tmp = computeExistsForall(tmp, !min1, !min2);
				// Add/remove target/non-remain
				tmp = JDD.And(remain.copy(), tmp);
				tmp = JDD.Or(target.copy(), tmp);
				// Check termination (inner)
				v_done = (tmp.equals(v));
				JDD.Deref(v);
				v = tmp;
			}
			// Check termination (outer)
			u_done = (v.equals(u));
			JDD.Deref(u);
			u = v;
		}
		// Derefs
		JDD.Deref(remain);
		JDD.Deref(target);

		timer = System.currentTimeMillis() - timer;
		mainLog.println("Prob1 took " + iters + " iterations and " + timer / 1000.0 + " seconds.");

		return u;
	}

	/**
	 * <br>[ REFS: <i>result</i>; DEREFS: target ]
	 * @param b
	 * @return
	 */
	protected JDDNode prob0(JDDNode target, boolean min1, boolean min2)
	{
		JDDNode u, tmp;
		boolean u_done;
		long timer;
		int iters = 0;

		// Start precomputation
		mainLog.println("Starting Prob0...");
		timer = System.currentTimeMillis();

		if (target.equals(JDD.ZERO)) {
			return reach.copy();
		}

		u_done = false;
		u = target.copy();
		while (!u_done) {
			iters++;
			// Find some-to-u
			tmp = JDD.SwapVariables(u.copy(), allDDRowVars, allDDColVars);
			tmp = JDD.ThereExists(JDD.And(tmp, trans01.copy()), allDDColVars);
			// Exists/forall
			tmp = computeExistsForall(tmp, !min1, !min2);
			tmp = JDD.Or(tmp, target.copy());
			u_done = (u.equals(tmp));
			JDD.Deref(u);
			u = tmp;
		}
		JDD.Ref(reach);
		JDD.Ref(target);
		u = JDD.And(JDD.Not(u), JDD.Not(target), reach);
		// Derefs
		JDD.Deref(target);

		timer = System.currentTimeMillis() - timer;
		mainLog.println("Prob0 took " + iters + " iterations and " + timer / 1000.0 + " seconds.");

		return u;
	}

	// compute probabilities for next

	/**
	 * Compute next-state probabilities.
	 * i.e. compute the probability of being in a state in {@code target} in the next step.
	 * @param target Target states
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 * @param coalition The coalition of players which define player 1
	 */
	protected StateValues computeNextProbs(JDDNode target, boolean min1, boolean min2, Coalition coalition) throws PrismException
	{
		StateValues probs = null;

		mainLog.println("Computing next probabilities ...");
		long timer = System.currentTimeMillis();

		setCoalition(coalition);

		// Matrix multiply: trans * target
		target = JDD.SwapVariables(target, allDDRowVars, allDDColVars);
		JDDNode sol = JDD.MatrixMultiply(trans.copy(), target, allDDColVars, JDD.CMU);
		// Max (P1) and min (P2)
		sol = computeMaxMin(sol, min1, min2);

		timer = System.currentTimeMillis() - timer;
		mainLog.println("Computation took " + timer / 1000.0 + " seconds.");

		probs = new StateValuesMTBDD(sol, model);

		setCoalition(null);

		return probs;
	}

	/**
	 * Compute probabilities for a bounded until operator.
	 * @param b1 Remain in these states
	 * @param b2 Target states
	 * @param k Step bound
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 * @param coalition The coalition of players which define player 1
	 */
	protected StateValues computeBoundedUntilProbs(JDDNode b1, JDDNode b2, int k, boolean min1, boolean min2, Coalition coalition) throws PrismException
	{
		JDDNode yes, no, maybe;
		JDDNode sol;
		StateValues result = null;

		setCoalition(coalition);

		// Compute trivial yes/no/maybe states
		if (b2.equals(JDD.ZERO)) {
			yes = JDD.Constant(0);
			JDD.Ref(reach);
			no = reach;
			maybe = JDD.Constant(0);
		} else if (b1.equals(JDD.ZERO)) {
			JDD.Ref(b2);
			yes = b2;
			JDD.Ref(reach);
			JDD.Ref(b2);
			no = JDD.And(reach, JDD.Not(b2));
			maybe = JDD.Constant(0);
		} else {
			// Yes
			JDD.Ref(b2);
			yes = b2;
			// No
			JDD.Ref(reach);
			JDD.Ref(b1);
			JDD.Ref(b2);
			no = JDD.And(reach, JDD.Not(JDD.Or(b1, b2)));
			// Maybe
			JDD.Ref(reach);
			JDD.Ref(yes);
			JDD.Ref(no);
			maybe = JDD.And(reach, JDD.Not(JDD.Or(yes, no)));
		}

		// Print out yes/no/maybe
		mainLog.print("\ntarget=" + JDD.GetNumMintermsString(b2, allDDRowVars.n()));
		mainLog.print(", yes=" + JDD.GetNumMintermsString(yes, allDDRowVars.n()));
		mainLog.print(", no=" + JDD.GetNumMintermsString(no, allDDRowVars.n()));
		mainLog.print(", maybe=" + JDD.GetNumMintermsString(maybe, allDDRowVars.n()) + "\n");

		// If maybe is empty, we have the probabilities already
		if (maybe.equals(JDD.ZERO)) {
			JDD.Ref(yes);
			result = new StateValuesMTBDD(yes, model);
		}
		// Otherwise explicitly compute the remaining probabilities
		else {
			// Compute probabilities
			mainLog.println("\nComputing probabilities...");
			mainLog.println("Engine: " + Prism.getEngineString(engine));
			sol = computeReachProbsValIter(yes, maybe, k, min1, min2);
			result = new StateValuesMTBDD(sol, model);
		}

		// Derefs
		JDD.Deref(no);
		JDD.Deref(yes);
		JDD.Deref(maybe);
		setCoalition(null);

		return result;
	}

	/**
	 * Compute probabilities for an (unbounded) until operator (b1 U b2).
	 * @param b1 Remain in these states
	 * @param b2 Target states
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 * @param coalition The coalition of players which define player 1
	 */
	protected StateValues computeUntilProbs(JDDNode b1, JDDNode b2, boolean min1, boolean min2, Coalition coalition) throws PrismException
	{
		JDDNode yes, no, maybe;
		JDDNode sol;
		StateValues result = null;

		setCoalition(coalition);

		// Compute yes/no/maybe states
		if (b2.equals(JDD.ZERO)) {
			yes = JDD.Constant(0);
			JDD.Ref(reach);
			no = reach;
			maybe = JDD.Constant(0);
		} else if (b1.equals(JDD.ZERO)) {
			JDD.Ref(b2);
			yes = b2;
			JDD.Ref(reach);
			JDD.Ref(b2);
			no = JDD.And(reach, JDD.Not(b2));
			maybe = JDD.Constant(0);
		} else {
			// Precomputation
			long timerPrecomp = System.currentTimeMillis();
			if (precomp && prob0) {
				no = prob0(b2.copy(), min1, min2);
			} else {
				no = JDD.And(reach.copy(), JDD.Not(JDD.Or(b1.copy(), b2.copy())));
			}
			if (precomp && prob1) {
				yes = prob1(b1.copy(), b2.copy(), min1, min2);
			} else {
				yes = b2.copy();
			}
			maybe = JDD.And(reach.copy(), JDD.Not(JDD.Or(yes.copy(), no.copy())));
			timerPrecomp = System.currentTimeMillis() - timerPrecomp;
			mainLog.println("Precomputation took " + timerPrecomp / 1000.0 + " seconds.");
		}

		// Print out yes/no/maybe
		mainLog.print("\ntarget=" + JDD.GetNumMintermsString(b2, allDDRowVars.n()));
		mainLog.print(", yes=" + JDD.GetNumMintermsString(yes, allDDRowVars.n()));
		mainLog.print(", no=" + JDD.GetNumMintermsString(no, allDDRowVars.n()));
		mainLog.print(", maybe=" + JDD.GetNumMintermsString(maybe, allDDRowVars.n()) + "\n");

		// If maybe is empty, we have the answer already...
		if (maybe.equals(JDD.ZERO)) {
			JDD.Ref(yes);
			result = new StateValuesMTBDD(yes, model);
		}
		// Otherwise explicitly compute the remaining probabilities
		else {
			// Compute probabilities
			mainLog.println("\nComputing probabilities...");
			mainLog.println("Engine: " + Prism.getEngineString(engine));
			try {
				sol = computeReachProbsValIter(yes, maybe, -1, min1, min2);
				result = new StateValuesMTBDD(sol, model);
			} catch (PrismException e) {
				JDD.Deref(no);
				JDD.Deref(yes);
				JDD.Deref(maybe);
				setCoalition(null);
				throw e;
			}
		}

		// Derefs
		JDD.Deref(no);
		JDD.Deref(yes);
		JDD.Deref(maybe);
		setCoalition(null);

		return result;
	}

	/**
	 * Compute probabilities for an (unbounded) until operator using value iteration.
	 */
	protected JDDNode computeReachProbsValIter(JDDNode yes, JDDNode maybe, int k, boolean min1, boolean min2) throws PrismException
	{
		JDDNode sol, tmp, tr;
		long timer;
		int iters;
		boolean done = false;

		// Start value iteration
		mainLog.println("Starting value iteration...");
		timer = System.currentTimeMillis();

		// Prepare DDs for transition function
		tr = JDD.Times(trans.copy(), maybe.copy());

		// Value iteration loop
		sol = yes.copy();
		iters = 0;
		while (!done) {
			tmp = sol.copy();
			// Matrix multiply: trans * target
			sol = JDD.MatrixMultiply(tr.copy(), JDD.PermuteVariables(sol, allDDRowVars, allDDColVars), allDDColVars, JDD.CMU);
			// Max (P1) and min (P2)
			sol = computeMaxMin(sol, min1, min2);
			// Yes states
			sol = JDD.Max(sol, yes.copy());
			// Check convergence
			done = JDD.EqualSupNorm(sol, tmp, termCritParam);
			JDD.Deref(tmp);
			iters++;
			if (iters == maxIters) {
				break;
			} else if (iters == k) {
				done = true;
			}
		}

		if (done == true) {
			timer = System.currentTimeMillis() - timer;
			mainLog.println("Value iteration converged after " + iters + " iterations and " + timer / 1000.0 + " seconds.");
		}

		// Strategy synthesis
		/*if (done == true && genStrat) {
			tmp = JDD.MatrixMultiply(tr.copy(), JDD.PermuteVariables(sol.copy(), allDDRowVars, allDDColVars), allDDColVars, JDD.CMU);
			tmp = JDD.Apply(JDD.EQUALS, tmp, sol.copy());
			JDDVars newAllNondetVars = new JDDVars();
			for (int v = 0; v < allDDNondetVars.n(); v++) {
				newAllNondetVars.addVar(model.getModelVariables().allocateVariable("newact" + v));
			}
			JDDNode strat = JDD.SwapVariables(tmp, allDDNondetVars, newAllNondetVars);
			newAllNondetVars.derefAll();
		}*/

		if (!done) {
			JDD.Deref(tr);
			JDD.Deref(sol);
			throw new PrismException("Could not converge after " + iters + " iterations.\nConsider increasing the maximum number of iterations");
		}

		JDD.Deref(tr);

		return sol;
	}

	// compute cumulative rewards

	protected StateValues computeCumulativeRewards(JDDNode sr, JDDNode trr, int k, boolean min1, boolean min2, Coalition coalition) throws PrismException
	{
		mainLog.println("\nComputing bounded cumulative reachability...");
		mainLog.println("Engine: " + Prism.getEngineString(engine));
		setCoalition(coalition);
		JDDNode res = computeReachRewardsValIter(null, null, null, sr, trr, k, min1, min2);
		StateValues result = new StateValuesMTBDD(res, model);
		setCoalition(null);
		return result;
	}

	/**
	 * 
	 * @param strw
	 * @param coalitions
	 * @param k
	 * @return
	 */
	protected StateValues computeInstantaneousRewards(JDDNode strw, int k, boolean min1, boolean min2, Coalition coalition) throws PrismException
	{
		StateValues result = null;
		JDDNode sol;
		long timer;
		int i;

		setCoalition(coalition);

		mainLog.println("\nComputing instantaneous rewards...");
		mainLog.println("Engine: " + Prism.getEngineString(engine));
		timer = System.currentTimeMillis();

		sol = JDD.SwapVariables(strw.copy(), allDDRowVars, allDDColVars);

		JDDNode nondetmaskrew = JDD.ITE(nondetMask.copy(), JDD.PlusInfinity(), JDD.Constant(0));

		for (i = 0; i < k; i++) {

			// Matrix multiply: trans * target
			sol = JDD.MatrixMultiply(trans.copy(), sol, allDDColVars, JDD.CMU);
			// Max (P1) and min (P2)
			sol = computeMaxMin(sol, nondetmaskrew, min1, min2);
			sol = JDD.SwapVariables(sol, allDDRowVars, allDDColVars);

		}
		timer = System.currentTimeMillis() - timer;
		mainLog.println("Instantatenous rewards took " + timer / 1000.0 + " seconds.");

		sol = JDD.SwapVariables(sol, allDDColVars, allDDRowVars);
		result = new StateValuesMTBDD(sol, model);
		setCoalition(null);
		JDD.Deref(nondetmaskrew);

		return result;
	}

	/**
	 * 
	 * @param expr
	 * @param sr
	 * @param trr
	 * @param coalitions
	 * @return
	 * @throws PrismException
	 */
	protected StateValues computeReachRewardsInfinity(ExpressionTemporal expr, JDDNode strw, JDDNode trrw, boolean min1, boolean min2, Coalition coalition)
			throws PrismException
	{
		StateValues result = null;

		JDDNode target, inf, init, res;

		mainLog.println("\nStarting expected reachability...");
		mainLog.println("Engine: " + Prism.getEngineString(engine));

		target = checkExpressionDD(expr.getOperand2(), model.getReach().copy());

		setCoalition(coalition);

		// If <<C>>Rmax=?[F t], builds coalition to compute <<N\C>>P=1[F t] that is, the set of states from which C can be forced to reach t with prob 1. 
		// The complement of this set is then <<C>>P<1[F t] as to maximise rewards C will try not reach the target and accumulate rewards indefinitely. 

		// If <<C>>Rmin=?[F t], builds coalition to compute <<C>>P=1[F t] that is, the set of states from which C can reach the target with prob 1 and thus
		// minimise rewards. The complement of this set is the set of states from which C can potentially be forced to accumulate rewards indefinitely. 

		long timerPrecomp = System.currentTimeMillis();
		JDD.Ref(target);
		inf = prob1(JDD.Constant(1), target, !min1, !min2);
		JDD.Ref(reach);
		inf = JDD.And(JDD.Not(inf), reach);
		timerPrecomp = System.currentTimeMillis() - timerPrecomp;
		mainLog.println("Precomputation took " + timerPrecomp / 1000.0 + " seconds.");

		mainLog.println("target=" + JDD.GetNumMintermsString(target, allDDRowVars.n()) + ", inf=" + JDD.GetNumMintermsString(inf, allDDRowVars.n()));

		double minRewards = Math.min(JDD.FindMinPositive(strw), JDD.FindMinPositive(trrw));
		double maxRewards = Math.max(JDD.FindMax(strw), JDD.FindMax(trrw));
		// Would be good to add an allNonZero check at some point, as precomputation would not be required. //
		double epsilon = Math.min(minRewards, maxRewards * 0.01);

		JDD.Ref(strw);
		JDDNode sttRewardsMod = JDD.Equals(strw, 0.0);
		JDD.Ref(trrw);
		JDDNode actRewardsMod = JDD.Equals(trrw, 0.0);

		JDD.Ref(trans01);
		actRewardsMod = JDD.And(actRewardsMod, trans01);
		actRewardsMod = JDD.Times(actRewardsMod, JDD.Constant(epsilon));
		JDD.Ref(trrw);
		actRewardsMod = JDD.Plus(trrw, actRewardsMod);

		sttRewardsMod = JDD.Times(sttRewardsMod, JDD.Constant(epsilon));
		JDD.Ref(strw);
		JDD.Ref(reach);
		JDD.Ref(target);
		sttRewardsMod = JDD.Times(JDD.Plus(sttRewardsMod, strw), reach, JDD.Not(target), JDD.Not(inf.copy())); // This probably can and should be simplified

		init = null;
		try {
			mainLog.println("Computing the upper bound where " + epsilon + " is used instead of 0.0.");
			init = computeReachRewardsValIter(null, inf, target, sttRewardsMod, actRewardsMod, -1, min1, min2);
			res = computeReachRewardsValIter(init, inf, target, strw, trrw, -1, min1, min2);
		} catch (PrismException e) {
			setCoalition(null);
			if (init != null) {
				JDD.Deref(init);
			}
			JDD.Deref(inf);
			JDD.Deref(target);
			JDD.Deref(sttRewardsMod);
			JDD.Deref(actRewardsMod);

			throw e;
		}

		// Derefs
		setCoalition(null);
		JDD.Deref(init);
		JDD.Deref(inf);
		JDD.Deref(target);
		JDD.Deref(sttRewardsMod);
		JDD.Deref(actRewardsMod);

		result = new StateValuesMTBDD(res, model);
		return result;
	}

	/**
	 * Performs value iteration with rewards for a given number of steps or until convergence.
	 * <br>[ REFS: <i>result</i>; DEREFS:  ]
	 * @param init
	 * @param inf
	 * @param target
	 * @param strw
	 * @param trrw
	 * @param epsilon
	 * @return
	 */
	protected JDDNode computeReachRewardsValIter(JDDNode init, JDDNode inf, JDDNode target, JDDNode strw, JDDNode trrw, int k, boolean min1, boolean min2)
			throws PrismException
	{
		JDDNode nondetmaskrew, sol, tmp, tr, sr, trr;
		long timer;
		int iters;
		boolean done = false;

		// Start value iteration
		mainLog.println("Starting value iteration...");
		timer = System.currentTimeMillis();

		tr = trans.copy();
		sr = strw.copy();

		if (inf != null && target != null) {
			JDD.Ref(inf);
			JDD.Ref(target);
			JDD.Ref(inf);
			JDD.Ref(target);
			tr = JDD.Times(tr, JDD.Not(inf), JDD.Not(target)); // Remember to add information about infinite/target states to solution vector
			sr = JDD.Times(sr, JDD.Not(inf), JDD.Not(target)); // Remember to add information about infinite/target states to solution vector
		} else if (inf != null) {
			JDD.Ref(inf);
			JDD.Ref(inf);
			tr = JDD.Times(tr, JDD.Not(inf));
			sr = JDD.Times(sr, JDD.Not(inf));
		} else if (target != null) {
			JDD.Ref(target);
			JDD.Ref(target);
			tr = JDD.Times(tr, JDD.Not(target));
			sr = JDD.Times(sr, JDD.Not(target));
		}

		if (init == null) {
			if (inf == null) {
				sol = JDD.Constant(0);
			} else {
				sol = JDD.ITE(inf.copy(), JDD.PlusInfinity(), JDD.Constant(0));
			}
		} else {
			JDD.Ref(init);
			sol = JDD.SwapVariables(init, allDDRowVars, allDDColVars);
		}

		trr = JDD.Times(tr.copy(), trrw.copy());
		trr = JDD.SumAbstract(trr, allDDColVars);

		iters = 0;

		nondetmaskrew = JDD.ITE(nondetMask.copy(), JDD.PlusInfinity(), JDD.Constant(0));

		tmp = null;

		if (k == 0) {
			done = true;
		}
		while (!done) {

			tmp = (iters != 0) ? sol.copy() : JDD.PlusInfinity();

			// Matrix multiply: trans * target; add transition rewards
			sol = JDD.MatrixMultiply(tr.copy(), sol, allDDColVars, JDD.CMU);
			sol = JDD.Apply(JDD.PLUS, sol, trr.copy());
			// Max (P1) and min (P2)
			sol = computeMaxMin(sol, nondetmaskrew, min1, min2);
			// Add state rewards
			sol = JDD.Plus(sol, sr.copy());

			if (inf != null) {
				sol = JDD.ITE(inf.copy(), JDD.PlusInfinity(), sol);
			}
			sol = JDD.SwapVariables(sol, allDDRowVars, allDDColVars);
			done = JDD.EqualSupNorm(sol, tmp, termCritParam);

			JDD.Deref(tmp);
			iters++;

			if (iters == maxIters) {
				break;
			} else if (iters == k) {
				done = true;
			}
		}

		if (done == true) {
			timer = System.currentTimeMillis() - timer;
			mainLog.println("Value iteration converged after " + iters + " iterations and " + timer / 1000.0 + " seconds.");
		}
		
		// Strategy synthesis
		/*if (done == true && genStrat) {
			tmp = JDD.MatrixMultiply(tr.copy(), JDD.PermuteVariables(sol.copy(), allDDRowVars, allDDColVars), allDDColVars, JDD.CMU);
			tmp = JDD.Apply(JDD.PLUS, tmp, trrw.copy());
			tmp = JDD.Apply(JDD.EQUALS, tmp, sol.copy());
			JDDVars newAllNondetVars = new JDDVars();
			for (int v = 0; v < allDDNondetVars.n(); v++) {
				newAllNondetVars.addVar(model.getModelVariables().allocateVariable("newact" + v));
			}
			JDDNode strat = JDD.SwapVariables(tmp, allDDNondetVars, newAllNondetVars);
			newAllNondetVars.derefAll();
		}*/

		if (!done) {
			JDD.Deref(tr);
			JDD.Deref(sol);
			JDD.Deref(nondetmaskrew);
			JDD.Deref(sr);
			JDD.Deref(trr);
			throw new PrismException("Could not converge after " + iters + " iterations.\nConsider increasing the maximum number of iterations");
		}

		JDD.Deref(tr);
		JDD.Deref(nondetmaskrew);
		JDD.Deref(sr);
		JDD.Deref(trr);

		sol = JDD.SwapVariables(sol, allDDColVars, allDDRowVars);
		if (inf != null) {
			sol = JDD.Max(sol, JDD.Times(inf.copy(), JDD.PlusInfinity()));
		}

		return sol;
	}

	/*protected JDDNode computeReachEquilibria(JDDNode[] targets, JDDNode srs[], JDDNode trrs[], int k, boolean min)
	{
		JDDNode[] sols = new JDDNode[2];
		JDDNode[] tmps = new JDDNode[2];
		JDDNode[] vals = new JDDNode[2];
		JDDNode arg1, arg2, argS, sum, sol, tmp, tr;
		int iters;
		boolean done = false;

		try {
			NondetModel mdp = model.toMDP();
			NondetModelChecker mdp_mc = new NondetModelChecker(prism, model, propertiesFile);
			vals[0] = mdp_mc.computeReachRewards(mdp.trans, JDD.Constant(0), mdp.trans01, mdp.stateRewards[0], mdp.transRewards[0], targets[0], min)
					.convertToStateValuesMTBDD().getJDDNode();
			vals[1] = mdp_mc.computeReachRewards(mdp.trans, JDD.Constant(0), mdp.trans01, mdp.stateRewards[1], mdp.transRewards[1], targets[1], min)
					.convertToStateValuesMTBDD().getJDDNode();
//			vals[0] = mdp_mc.computeCumulRewards(mdp.trans, mdp.stateRewards[0], 
//															mdp.transRewards[0], k, min).convertToStateValuesMTBDD().getJDDNode();
//			vals[1] = mdp_mc.computeCumulRewards(mdp.trans, mdp.stateRewards[1], 
//															mdp.transRewards[1], k, min).convertToStateValuesMTBDD().getJDDNode();
		} catch (PrismException e) {
			e.printStackTrace();
		}

		tr = trans.copy();
		tr = JDD.Times(tr, JDD.Not(targets[0].copy()));
		tr = JDD.Times(tr, JDD.Not(targets[1].copy()));
		sols[0] = JDD.Times(vals[0], JDD.Not(targets[0].copy()));
		sols[1] = JDD.Times(vals[1], JDD.Not(targets[1].copy()));

		iters = 0;
		sol = JDD.Constant(0);

		System.out.println("-- allDDColVars");
		System.out.println(allDDColVars);
		System.out.println("-- allDDRowVars");
		System.out.println(allDDRowVars);

		System.out.println("ddc1");
		JDD.PrintMinterms(mainLog, ddc1.copy());
		System.out.println("ddc2");
		JDD.PrintMinterms(mainLog, ddc2.copy());

		System.out.println("trrs[0]");
		JDD.PrintMinterms(mainLog, trrs[0].copy());
		System.out.println("trrs[1]");
		JDD.PrintMinterms(mainLog, trrs[1].copy());

		System.out.println("sols[0]");
		JDD.PrintMinterms(mainLog, sols[0].copy());

		System.out.println("sols[1]");
		JDD.PrintMinterms(mainLog, sols[1].copy());

		System.out.println("\n ## Value Iteration ##");

		while (!done) {
			tmp = (iters != 0) ? sol.copy() : JDD.Constant(0);

			System.out.println("\n## " + iters);

			sols[0] = JDD.SwapVariables(sols[0], allDDRowVars, allDDColVars);
			sols[1] = JDD.SwapVariables(sols[1], allDDRowVars, allDDColVars);

			// Multiplication
			sols[0] = JDD.MatrixMultiply(tr.copy(), sols[0], allDDColVars, JDD.CMU);
			sols[1] = JDD.MatrixMultiply(tr.copy(), sols[1], allDDColVars, JDD.CMU);

			// Retrieves the actions
			sols[0] = JDD.Times(trans01.copy(), sols[0]);
			sols[1] = JDD.Times(trans01.copy(), sols[1]);

			System.out.println("sols[0]");
			JDD.PrintMinterms(mainLog, sols[0].copy());

			System.out.println("sols[1]");
			JDD.PrintMinterms(mainLog, sols[1].copy());

			// Sum individual rewards
			sols[0] = JDD.Plus(sols[0].copy(), trrs[0].copy());
			sols[1] = JDD.Plus(sols[1].copy(), trrs[1].copy());

			// Sums solution vectors
			sol = JDD.Plus(sols[0], sols[1]);

			System.out.println("sol");
			JDD.PrintMinterms(mainLog, sol.copy());

			System.out.println("sols[0]");
			JDD.PrintMinterms(mainLog, sols[0].copy());

			System.out.println("sols[1]");
			JDD.PrintMinterms(mainLog, sols[1].copy());

			// Find the actions that maximise the sum 
			argS = JDD.MaxAbstract(sol.copy(), allDDNondetVars);
			System.out.println("argS");
			JDD.PrintMinterms(mainLog, argS.copy());
			argS = JDD.And(JDD.Apply(JDD.EQUALS, argS.copy(), sol.copy()), JDD.Apply(JDD.NOTEQUALS, argS, JDD.Constant(0)));
			System.out.println("argS");
			JDD.PrintMinterms(mainLog, argS.copy());

			System.out.println("sols[0]");
			JDD.PrintMinterms(mainLog, sols[0].copy());

			// Finds the actions that maximise for player 1 in states controlled by him
			arg1 = JDD.MaxAbstract(sols[0].copy(), allDDNondetVars);
			System.out.println("arg1");
			JDD.PrintMinterms(mainLog, arg1.copy());
			arg1 = JDD.And(JDD.Apply(JDD.EQUALS, arg1.copy(), sols[0]), JDD.Apply(JDD.NOTEQUALS, arg1, JDD.Constant(0)));
			System.out.println("arg1");
			JDD.PrintMinterms(mainLog, arg1.copy());

			System.out.println("sols[1]");
			JDD.PrintMinterms(mainLog, sols[1].copy());

			// Finds the actions that maximise for player 2 in states controlled by him
			arg2 = JDD.MaxAbstract(sols[1].copy(), allDDNondetVars);
			System.out.println("arg2");
			JDD.PrintMinterms(mainLog, arg2.copy());
			arg2 = JDD.And(JDD.Apply(JDD.EQUALS, arg2.copy(), sols[1]), JDD.Apply(JDD.NOTEQUALS, arg2, JDD.Constant(0)));
			System.out.println("arg2");
			JDD.PrintMinterms(mainLog, arg2.copy());

			// Gets actions that maximise both
			arg1 = JDD.ITE(JDD.And(arg1, argS.copy()), JDD.And(arg1, argS.copy()), arg1);
			arg2 = JDD.ITE(JDD.And(arg2, argS.copy()), JDD.And(arg2, argS.copy()), arg2);

			System.out.println("arg1");
			JDD.PrintMinterms(mainLog, arg1.copy());
			System.out.println("arg2");
			JDD.PrintMinterms(mainLog, arg2.copy());

			// Get unique actions
			tmps[0] = JDD.MaxAbstract(JDD.Times(sols[0].copy(), arg1, ddc1), allDDNondetVars);
			System.out.println("tmps[0]");
			JDD.PrintMinterms(mainLog, tmps[0].copy());
			sols[0] = JDD.Plus(tmps[0], JDD.MaxAbstract(JDD.Times(sols[0].copy(), arg1, ddc2), allDDNondetVars));
			sols[0] = JDD.ThereExists(sols[0], allDDColVars);
			System.out.println("sols[0]");
			JDD.PrintMinterms(mainLog, sols[0].copy());

			tmps[1] = JDD.MaxAbstract(JDD.Times(sols[1].copy(), arg2, ddc1), allDDNondetVars);
			System.out.println("tmps[1]");
			JDD.PrintMinterms(mainLog, tmps[1].copy());
			sols[1] = JDD.Plus(tmps[1], JDD.MaxAbstract(JDD.Times(sols[1].copy(), arg2, ddc2), allDDNondetVars));
			sols[1] = JDD.ThereExists(sols[1], allDDColVars);
			System.out.println("sols[1]");
			JDD.PrintMinterms(mainLog, sols[1].copy());

			JDD.Deref(argS);

			done = JDD.EqualSupNorm(sol, tmp, termCritParam);

			iters++;
			if (iters == maxIters) {
				break;
			} else if (iters == k) {
				done = true;
			}
		}

		sol = JDD.Plus(sols[0], sols[1]);

		System.out.println("##\nsol");
		JDD.PrintMinterms(mainLog, sol.copy());

		return sol;
	}*/

	// Auxiliary functions
	
	/**
	 * Set a coalition of players for the current SMG
	 * (which effectively makes it an STPG with player 1 representing the coalition and 2 the rest).
	 * Create the required DDs (ddc1, ddc2).
	 * If {@code coalition} is null, the DDs are instead derefed.
	 * @param coalition Coalition info object
	 */
	public void setCoalition(Coalition coalition) throws PrismException
	{
		// Clear info if coalition is null
		if (coalition == null) {
			if (ddc1 != null) {
				JDD.Deref(ddc1);
			}
			if (ddc2 != null) {
				JDD.Deref(ddc2);
			}
			return;
		}

		// Build DDs to represent players in each coalition
		model.getPlayerInfo().setCoalition(coalition);
		ddc1 = JDD.Constant(0);
		ddc2 = JDD.Constant(0);
		for (int p = 0; p < model.getNumPlayers(); p++) {
			if (model.getPlayerInfo().getPlayer(p) == 0) {
				ddc1 = JDD.Or(ddc1, model.getDdPlayerCube(p).copy());
			}
			if (model.getPlayerInfo().getPlayer(p) == 1) {
				ddc2 = JDD.Or(ddc2, model.getDdPlayerCube(p).copy());
			}
		}
	}

	/**
	 * Perform the min/max over actions for a state-action vector of values
	 * (max for player 1 and min for player 2).
	 * <br>[ REFS: <i>result</i>; DEREFS: vals ]
	 */
	protected JDDNode computeMaxMin(JDDNode vals, boolean min1, boolean min2)
	{
		return computeMaxMin(vals, null, min1, min2);
	}

	/**
	 * Perform the min/max over actions for a state-action vector of values
	 * (max for player 1 and min for player 2).
	 * Optionally, provide a nondetmask (if null, normal one is used)
	 * <br>[ REFS: <i>result</i>; DEREFS: vals ]
	 */
	protected JDDNode computeMaxMin(JDDNode vals, JDDNode mask, boolean min1, boolean min2)
	{
		// Min/max for player 1
		JDDNode res1 = JDD.Times(vals.copy(), ddc1.copy());
		if (min1) {
			res1 = JDD.Apply(JDD.MAX, res1, mask == null ? nondetMask.copy() : mask.copy());
			res1 = JDD.MinAbstract(res1, allDDNondetVars);
		} else {
			res1 = JDD.MaxAbstract(res1, allDDNondetVars);
		}
		// Min/max for player 2
		JDDNode res2 = JDD.Times(vals.copy(), ddc2.copy());
		if (min2) {
			res2 = JDD.Apply(JDD.MAX, res2, mask == null ? nondetMask.copy() : mask.copy());
			res2 = JDD.MinAbstract(res2, allDDNondetVars);
		} else {
			res2 = JDD.MaxAbstract(res2, allDDNondetVars);
		}
		// Combine/return
		JDD.Deref(vals);
		return JDD.Apply(JDD.PLUS, res1, res2);
	}

	/**
	 * Perform the exists/forall over actions for a state-action set
	 * (if exists1 is true, exists for player 1 and forall for player 2; otherwise reverse)
	 * <br>[ REFS: <i>result</i>; DEREFS: set ]
	 */
	protected JDDNode computeExistsForall(JDDNode set, boolean exists1, boolean exists2)
	{
		// Exists/forall for player 1
		JDDNode res1 = JDD.And(set.copy(), ddc1.copy());
		if (exists1) {
			res1 = JDD.ThereExists(res1, allDDNondetVars);
		} else {
			res1 = JDD.ForAll(JDD.Or(res1, nondetMask.copy()), allDDNondetVars);
		}
		// Exists/forall for player 2
		JDDNode res2 = JDD.And(set.copy(), ddc2.copy());
		if (exists2) {
			res2 = JDD.ThereExists(res2, allDDNondetVars);
		} else {
			res2 = JDD.ForAll(JDD.Or(res2, nondetMask.copy()), allDDNondetVars);
		}
		// Combine/return
		JDD.Deref(set);
		return JDD.Or(res1, res2);
	}
}
