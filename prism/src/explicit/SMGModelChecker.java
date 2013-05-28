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

import java.util.BitSet;

import parser.ast.Expression;
import parser.ast.ExpressionPATL;
import parser.ast.ExpressionTemporal;
import prism.PrismException;
import strat.ExactValueStrategy;
import strat.Strategy;
import explicit.rewards.SMGRewards;

/**
 * Explicit-state model checker for multi-player stochastic games (SMGs).
 */
public class SMGModelChecker extends STPGModelChecker
{
	/**
	 * Compute probabilities for the contents of a P operator.
	 */
	protected StateValues checkProbPathFormula(Model model, ExpressionPATL exprPATL, boolean min) throws PrismException
	{

		//@clemens : don't change this - this works out the player coalition
		// setting coalition parameter
		((SMG) model).setCoalition(exprPATL.getCoalition());

		Expression expr = exprPATL.getExpressionProb().getExpression();

		// Test whether this is a simple path formula (i.e. PCTL)
		// and then pass control to appropriate method.
		if (expr.isSimplePathFormula()) {
			double p = -1;
			Expression pb = exprPATL.getExpressionProb().getProb();
			if (pb != null) {
				p = pb.evaluateDouble(constantValues);
			}
			return super.checkProbPathFormulaSimple(model, expr, min, !min, p);
		}

		/*
		 * TODO implement FG and GF formulae //Test if this is FG if (expr
		 * instanceof ExpressionTemporal) { ExpressionTemporal exprT =
		 * (ExpressionTemporal) expr; if (exprT.getOperator() ==
		 * ExpressionTemporal.P_F) { Expression expr2 = exprT.getOperand2(); if
		 * (expr2 instanceof ExpressionTemporal) { ExpressionTemporal expr2T =
		 * (ExpressionTemporal) expr2; if (expr2T.getOperator() ==
		 * ExpressionTemporal.P_G) { Expression expr3 = expr2T.getOperand2(); if
		 * (!(expr3 instanceof ExpressionTemporal)) { return
		 * super.checkFG(model, expr, min, !min); } } } } }
		 * 
		 * //Test whether this is GF if (expr instanceof ExpressionTemporal) {
		 * ExpressionTemporal exprT = (ExpressionTemporal) expr; if
		 * (exprT.getOperator() == ExpressionTemporal.P_G) { Expression expr2 =
		 * exprT.getOperand2(); if (expr2 instanceof ExpressionTemporal) {
		 * ExpressionTemporal expr2T = (ExpressionTemporal) expr2; if
		 * (expr2T.getOperator() == ExpressionTemporal.P_F) { Expression expr3 =
		 * expr2T.getOperand2(); if (!(expr3 instanceof ExpressionTemporal)) {
		 * return super.checkGF(model, expr, min, !min); } } } } }
		 */

		// in other case
		throw new PrismException("Explicit engine does not yet handle LTL-style path formulas except for GF and FG");
	}

	/**
	 * Compute rewards for the contents of an R operator.
	 */
	protected StateValues checkRewardFormula(Model model, SMGRewards modelRewards, ExpressionPATL exprPATL, boolean min) throws PrismException
	{
		// setting coalition parameter
		((SMG) model).setCoalition(exprPATL.getCoalition());

		StateValues rewards = null;
		Expression expr = exprPATL.getExpressionRew().getExpression();

		if (expr instanceof ExpressionTemporal) {
			ExpressionTemporal exprTemp = (ExpressionTemporal) expr;
			switch (exprTemp.getOperator()) {
			case ExpressionTemporal.R_F:
				rewards = checkRewardReach(model, modelRewards, exprTemp, min, !min, STPGModelChecker.R_INFINITY);
				break;
			case ExpressionTemporal.R_Fc:
				rewards = checkRewardReach(model, modelRewards, exprTemp, min, !min, STPGModelChecker.R_CUMULATIVE);
				break;
			case ExpressionTemporal.R_F0:
				rewards = checkRewardReach(model, modelRewards, exprTemp, min, !min, STPGModelChecker.R_ZERO);
				break;
			default:
				throw new PrismException("Explicit engine does not yet handle the " + exprTemp.getOperatorSymbol() + " operator in the R operator");
			}
		}

		if (rewards == null)
			throw new PrismException("Unrecognised operator in R operator");

		return rewards;
	}

	protected StateValues checkExactProbabilityFormula(Model model, ExpressionPATL expr, double p) throws PrismException
	{
		if (expr.getExpressionProb().getExpression() instanceof ExpressionTemporal
				&& ((ExpressionTemporal) expr.getExpressionProb().getExpression()).hasBounds()) {
			throw new PrismException("The exact probability queries are not supported for step-bounded properties");
		}

		((SMG) model).setCoalition(expr.getCoalition());
		// 1) check whether the game is stopping, if not - terminate
		// 1.1) find states which have self loops only
		BitSet terminal = findTerminalStates(model);
		// 1.2) check whether the minmin prob to reach those states is
		// 1, if not - terminate, if yes continue to 2)
		double[] res = ((SMGModelChecker) this).computeUntilProbs((STPG) model, null, terminal, true, true, 1.0).soln;

		// System.out.println("Terminal states: " + terminal);
		// System.out.println(Arrays.toString(res));
		for (int i = 0; i < res.length; i++)
			if (res[i] < 1.0 - 1e-6)
				throw new PrismException("The game is not stopping. The exact probability queries only work for stopping games");

		// 2) computing minmax and maxmin values for all states
		double[] minmax = null, maxmin = null; // see the do loop below

		// 3) removing states from the game which have minmax>maxmin
		// model.
		int n = model.getNumStates();
		boolean repeat;
		BitSet removed = new BitSet(n), removedNew = new BitSet(n);
		STPG stpg = ((STPG) model);

		Strategy minStrat = null;
		Strategy maxStrat = null;

		do {
			// computing minmax and maxmin
			minmax = this.checkProbPathFormula(model, expr, true).getDoubleArray();
			if (generateStrategy)
				minStrat = strategy;
			maxmin = this.checkProbPathFormula(model, expr, false).getDoubleArray();
			if (generateStrategy)
				maxStrat = strategy;

			repeat = false;
			// checking which states are marked for removal
			for (int i = 0; i < n; i++)
				if (!removed.get(i) && minmax[i] > maxmin[i]) {
					removed.set(i);
					removedNew.set(i);
					repeat = true;
				}

			// disabling choices that have transitions to those states
			removedNew.flip(0, n);
			for (int i = 0; i < n; i++)
				for (int j = 0; j < model.getNumChoices(i); j++)
					if (!stpg.allSuccessorsInSet(i, j, removedNew))
						stpg.disableChoice(i, j);
			removedNew.clear();
			// 4) repeat 2-3 while the set of states from 3 is empty
		} while (repeat);

		// 5) if bound is null, return the interval, otherwise check
		// whether bound is in the interval.
		BitSet ret = new BitSet(n);
		for (int i = 0; i < n; i++)
			ret.set(i, !removed.get(i) && minmax[i] <= p && maxmin[i] >= p);

		// enabling choices that have been disabled for model checking
		stpg.enableAllChoices();

		if (generateStrategy) {
			strategy = new ExactValueStrategy(minStrat, minmax, maxStrat, maxmin, p, (STPG) model);
		}

		return StateValues.createFromBitSet(ret, model);
	}

	protected StateValues checkExactRewardFormula(Model model, SMGRewards modelRewards, ExpressionPATL expr, double p) throws PrismException
	{
		((SMG) model).setCoalition(expr.getCoalition());
		// check if the reward is Fc
		ExpressionTemporal exprTemp = null;
		if (expr.getExpressionRew().getExpression() instanceof ExpressionTemporal) {
			exprTemp = (ExpressionTemporal) expr.getExpressionRew().getExpression();
			switch (exprTemp.getOperator()) {
			case ExpressionTemporal.R_Fc:
				break;
			case ExpressionTemporal.R_F:
				throw new PrismException("Only cumulative reward type is supported for exact values.");
			case ExpressionTemporal.R_F0:
				throw new PrismException("Only cumulative reward type is supported for exact values.");
			default:
				throw new PrismException("Only cumulative reward type is supported for exact values.");
			}
		} else {
			throw new PrismException("Only temporal expression are supported at the moment");
		}

		// 1) check whether the game is stopping, if not - terminate
		// 1.1) find states which have self loops only
		BitSet terminal = findTerminalStates(model);
		// 1.2) check whether the minmin prob to reach those states is
		// 1, if not - terminate, if yes continue to 2)
		double[] res = ((SMGModelChecker) this).computeUntilProbs((STPG) model, null, terminal, true, true, 1.0).soln;

		// System.out.println("Terminal states: " + terminal);
		// System.out.println(Arrays.toString(res));
		for (int i = 0; i < res.length; i++)
			if (res[i] < 1.0 - 1e-6)
				throw new PrismException("The game is not stopping. The exact probability queries only work for stopping games");

		// 2) computing minmax and maxmin values for all states
		double[] minmax = null, maxmin = null; // see the do loop below 

		// 3) removing states from the game which have minmax>maxmin
		// model.
		int n = model.getNumStates();
		boolean repeat;
		BitSet removed = new BitSet(n), removedNew = new BitSet(n);
		STPG stpg = ((STPG) model);

		Strategy minStrat = null;
		Strategy maxStrat = null;

		do {
			// computing minmax and maxmin
			minmax = this.checkRewardReach(model, modelRewards, exprTemp, true, false, STPGModelChecker.R_CUMULATIVE).valuesD;
			if (generateStrategy)
				minStrat = strategy;
			maxmin = this.checkRewardReach(model, modelRewards, exprTemp, false, true, STPGModelChecker.R_CUMULATIVE).valuesD;
			if (generateStrategy)
				maxStrat = strategy;

			repeat = false;
			// checking which states are marked for removal
			for (int i = 0; i < n; i++)
				if (!removed.get(i) && minmax[i] > maxmin[i]) {
					removed.set(i);
					removedNew.set(i);
					repeat = true;
				}

			// disabling choices that have transitions to those states
			removedNew.flip(0, n);
			for (int i = 0; i < n; i++)
				for (int j = 0; j < model.getNumChoices(i); j++)
					if (!stpg.allSuccessorsInSet(i, j, removedNew))
						stpg.disableChoice(i, j);
			removedNew.clear();
			// 4) repeat 2-3 while the set of states from 3 is empty
		} while (repeat);

		// 5) if bound is null, return the interval, otherwise check
		// whether bound is in the interval.
		BitSet ret = new BitSet(n);
		for (int i = 0; i < n; i++)
			ret.set(i, !removed.get(i) && minmax[i] <= p && maxmin[i] >= p);

		// enabling choices that have been disabled for model checking
		stpg.enableAllChoices();

		if (generateStrategy) {
			strategy = new ExactValueStrategy(minStrat, minmax, maxStrat, maxmin, p, (STPG) model);
		}

		return StateValues.createFromBitSet(ret, model);
	}
}
