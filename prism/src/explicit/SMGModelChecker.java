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
import java.util.List;

import parser.ast.Expression;
import parser.ast.ExpressionProb;
import parser.ast.ExpressionReward;
import parser.ast.ExpressionTemporal;
import prism.PrismComponent;
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
	 * Create a new SMGModelChecker, inherit basic state from parent (unless null).
	 */
	public SMGModelChecker(PrismComponent parent) throws PrismException
	{
		super(parent);
	}
	
	/**
	 * Model check a P operator expression with an exact probability bound (P=p[...])
	 */
	protected StateValues checkExactProbabilityFormula(SMG smg, ExpressionProb expr, List<String> coalition, double p) throws PrismException
	{
		// Just support untimed , non-LTL path formulas
		if (!expr.isSimplePathFormula())
			throw new PrismException("Exact probability queries cannot contain LTL formulas");
		if (Expression.containsTemporalTimeBounds(expr)) {
			throw new PrismException("Exact probability queries cannot contain time-bounded operators");
		}

		// 1) check whether the game is stopping, if not - terminate
		// 1.1) find states which have self loops only
		BitSet terminal = findTerminalStates(smg);
		// 1.2) check whether the minmin prob to reach those states is
		// 1, if not - terminate, if yes continue to 2)
		double[] res = computeUntilProbs(smg, null, terminal, true, true, 1.0).soln;

		// System.out.println("Terminal states: " + terminal);
		// System.out.println(Arrays.toString(res));
		for (int i = 0; i < res.length; i++)
			if (res[i] < 1.0 - 1e-6)
				throw new PrismException("The game is not stopping. The exact probability queries only work for stopping games");

		// 2) computing minmax and maxmin values for all states
		double[] minmax = null, maxmin = null; // see the do loop below

		// 3) removing states from the game which have minmax>maxmin model.
		int n = smg.getNumStates();
		boolean repeat;
		BitSet removed = new BitSet(n), removedNew = new BitSet(n);
		STPG stpg = ((STPG) smg);

		Strategy minStrat = null;
		Strategy maxStrat = null;

		do {
			// computing minmax and maxmin
			MinMax minMax1 = MinMax.minMin(true, false);
			minMax1.setCoalition(coalition);
			minMax1.setBound(p);
			minmax = checkProbPathFormula(smg, expr.getExpression(), minMax1).getDoubleArray();
			if (generateStrategy)
				minStrat = strategy;
			MinMax minMax2 = MinMax.minMin(false, true);
			minMax2.setCoalition(coalition);
			minMax2.setBound(p);
			maxmin = checkProbPathFormula(smg, expr.getExpression(), minMax2).getDoubleArray();
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
				for (int j = 0; j < smg.getNumChoices(i); j++)
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
			strategy = new ExactValueStrategy(minStrat, minmax, maxStrat, maxmin, p, (STPG) smg);
		}

		return StateValues.createFromBitSet(ret, smg);
	}

	protected StateValues checkExactRewardFormula(SMG smg, SMGRewards modelRewards, ExpressionReward exprRew, List<String> coalition, double p) throws PrismException
	{
		// check if the reward is Fc
		ExpressionTemporal exprTemp = null;
		if (exprRew.getExpression() instanceof ExpressionTemporal) {
			exprTemp = (ExpressionTemporal) exprRew.getExpression();
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
		BitSet terminal = findTerminalStates(smg);
		// 1.2) check whether the minmin prob to reach those states is
		// 1, if not - terminate, if yes continue to 2)
		double[] res = computeUntilProbs(smg, null, terminal, true, true, 1.0).soln;

		// System.out.println("Terminal states: " + terminal);
		// System.out.println(Arrays.toString(res));
		for (int i = 0; i < res.length; i++)
			if (res[i] < 1.0 - 1e-6)
				throw new PrismException("The game is not stopping. The exact probability queries only work for stopping games");

		// 2) computing minmax and maxmin values for all states
		double[] minmax = null, maxmin = null; // see the do loop below 

		// 3) removing states from the game which have minmax>maxmin
		// model.
		int n = smg.getNumStates();
		boolean repeat;
		BitSet removed = new BitSet(n), removedNew = new BitSet(n);

		Strategy minStrat = null;
		Strategy maxStrat = null;

		do {
			// computing minmax and maxmin
			MinMax minMax1 = MinMax.minMin(true, false);
			minMax1.setCoalition(coalition);
			minMax1.setBound(p);
			minmax = checkRewardReach(smg, modelRewards, exprTemp, minMax1).getDoubleArray();
			if (generateStrategy)
				minStrat = strategy;
			MinMax minMax2 = MinMax.minMin(false, true);
			minMax2.setCoalition(coalition);
			minMax2.setBound(p);
			maxmin = checkRewardReach(smg, modelRewards, exprTemp, minMax2).getDoubleArray();
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
				for (int j = 0; j < smg.getNumChoices(i); j++)
					if (!smg.allSuccessorsInSet(i, j, removedNew))
						smg.disableChoice(i, j);
			removedNew.clear();
			// 4) repeat 2-3 while the set of states from 3 is empty
		} while (repeat);

		// 5) if bound is null, return the interval, otherwise check
		// whether bound is in the interval.
		BitSet ret = new BitSet(n);
		for (int i = 0; i < n; i++)
			ret.set(i, !removed.get(i) && minmax[i] <= p && maxmin[i] >= p);

		// enabling choices that have been disabled for model checking
		smg.enableAllChoices();

		if (generateStrategy) {
			strategy = new ExactValueStrategy(minStrat, minmax, maxStrat, maxmin, p, smg);
		}

		return StateValues.createFromBitSet(ret, smg);
	}
}
