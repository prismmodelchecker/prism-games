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

import parser.ast.*;
import prism.*;

/**
 * Super class for explicit-state probabilistic model checkers
 */
public class ProbModelChecker extends StateModelChecker
{
	// Model checking functions

	/**
	 * Model check an expression and return the values for all states.
	 */
	public Object checkExpression(Model model, Expression expr) throws PrismException
	{
		Object res;

		// P operator
		if (expr instanceof ExpressionProb) {
			res = checkExpressionProb(model, (ExpressionProb) expr);
		}
		// R operator
		else if (expr instanceof ExpressionReward) {
			res = checkExpressionReward(model, (ExpressionReward) expr);
		}
		// Otherwise, use the superclass
		else {
			res = super.checkExpression(model, expr);
		}

		return res;
	}

	/**
	 * Model check a P operator expression and return the values for all states.
	 */
	protected StateValues checkExpressionProb(Model model, ExpressionProb expr) throws PrismException
	{
		// Relational operator
		String relOp;
		// For nondeterministic models, are we finding min (true) or max (false) probs
		boolean min1 = false;
		boolean min2 = false;
		ModelType modelType = model.getModelType();

		StateValues probs = null;

		// Get info from prob operator
		relOp = expr.getRelOp();

		// Check for unhandled cases
		if (expr.getProb() != null)
			throw new PrismException("Bounded P operators not yet supported");

		// For nondeterministic models, determine whether min or max probabilities needed
		if (modelType.nondeterministic()) {
			if (modelType != ModelType.STPG) {
				if (relOp.equals(">") || relOp.equals(">=") || relOp.equals("min=")) {
					// min
					min1 = true;
				} else if (relOp.equals("<") || relOp.equals("<=") || relOp.equals("max=")) {
					// max
					min1 = false;
				} else {
					throw new PrismException("Can't use \"P=?\" for nondeterministic models; use \"Pmin=?\" or \"Pmax=?\"");
				}
			} else {
				if (relOp.equals("minmin=")) {
					min1 = true;
					min2 = true;
				} else if (relOp.equals("minmax=")) {
					min1 = true;
					min2 = false;
				} else if (relOp.equals("maxmin=")) {
					min1 = false;
					min2 = true;
				} else if (relOp.equals("maxmax=")) {
					min1 = false;
					min2 = false;
				} else {
					throw new PrismException("Use e.g. \"Pminmax=?\" for stochastic games");
				}
			}
		}

		// Compute probabilities
		switch (modelType) {
		case CTMC:
			probs = ((CTMCModelChecker) this).checkProbPathFormula(model, expr.getExpression());
			break;
		case CTMDP:
			probs = ((CTMDPModelChecker) this).checkProbPathFormula(model, expr.getExpression(), min1);
			break;
		case DTMC:
			probs = ((DTMCModelChecker) this).checkProbPathFormula(model, expr.getExpression());
			break;
		case MDP:
			probs = ((MDPModelChecker) this).checkProbPathFormula(model, expr.getExpression(), min1);
			break;
		case STPG:
			probs = ((STPGModelChecker) this).checkProbPathFormula(model, expr.getExpression(), min1, min2);
			break;
		default:
			throw new PrismException("Cannot model check " + expr + " for a " + modelType);
		}

		// Print out probabilities
		if (getVerbosity() > 5) {
			mainLog.print("\nProbabilities (non-zero only) for all states:\n");
			mainLog.print(probs);
		}

		// For =? properties, just return values
		return probs;
	}
	
	/**
	 * Model check an R operator expression and return the values for all states.
	 */
	protected StateValues checkExpressionReward(Model model, ExpressionReward expr) throws PrismException
	{
		String relOp; // Relational operator
		boolean min1 = false;
		boolean min2 = false;
		ModelType modelType = model.getModelType();

		StateValues rews = null;

		// Get info from reward operator
		relOp = expr.getRelOp();

		// Check for unhandled cases
		// TODO

		// For nondeterministic models, determine whether min or max probabilities needed
		if (modelType.nondeterministic()) {
			if (modelType != ModelType.STPG) {
				if (relOp.equals(">") || relOp.equals(">=") || relOp.equals("min=")) {
					// min
					min1 = true;
				} else if (relOp.equals("<") || relOp.equals("<=") || relOp.equals("max=")) {
					// max
					min1 = false;
				} else {
					throw new PrismException("Can't use \"P=?\" for nondeterministic models; use \"Pmin=?\" or \"Pmax=?\"");
				}
			} else {
				if (relOp.equals("minmin=")) {
					min1 = true;
					min2 = true;
				} else if (relOp.equals("minmax=")) {
					min1 = true;
					min2 = false;
				} else if (relOp.equals("maxmin=")) {
					min1 = false;
					min2 = true;
				} else if (relOp.equals("maxmax=")) {
					min1 = false;
					min2 = false;
				} else {
					throw new PrismException("Use e.g. \"Pminmax=?\" for stochastic games");
				}
			}
		}

		// Compute rewards
		switch (modelType) {
		case CTMC:
			rews = ((CTMCModelChecker) this).checkRewardFormula(model, expr.getExpression());
			break;
		case DTMC:
			rews = ((DTMCModelChecker) this).checkRewardFormula(model, expr.getExpression());
			break;
		case MDP:
			rews = ((MDPModelChecker) this).checkRewardFormula(model, expr.getExpression(), min1);
			break;
		case STPG:
			rews = ((STPGModelChecker) this).checkRewardFormula(model, expr.getExpression(), min1, min2);
			break;
		default:
			throw new PrismException("Cannot model check " + expr + " for a " + modelType);
		}

		// Print out probabilities
		if (getVerbosity() > 5) {
			mainLog.print("\nProbabilities (non-zero only) for all states:\n");
			mainLog.print(rews);
		}

		// For =? properties, just return values
		return rews;
	}
}
