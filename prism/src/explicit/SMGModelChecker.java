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

import parser.VarList;
import parser.ast.Coalition;
import parser.ast.Declaration;
import parser.ast.DeclarationIntUnbounded;
import parser.ast.Expression;
import parser.ast.ExpressionProb;
import parser.ast.ExpressionReward;
import parser.ast.ExpressionTemporal;
import prism.Prism;
import prism.PrismComponent;
import prism.PrismException;
import prism.PrismFileLog;
import strat.ExactValueStrategy;
import strat.Strategy;
import acceptance.AcceptanceReach;
import acceptance.AcceptanceType;
import explicit.rewards.Rewards;
import explicit.rewards.SMGRewards;

/**
 * Explicit-state model checker for multi-player stochastic games (SMGs).
 */
public class SMGModelChecker extends ProbModelChecker
{
	/**
	 * Create a new SMGModelChecker, inherit basic state from parent (unless null).
	 */
	public SMGModelChecker(PrismComponent parent) throws PrismException
	{
		super(parent);
	}
	
	// Model checking functions
	
	/**
	 * Compute rewards for a co-safe LTL reward operator.
	 */
	protected StateValues checkRewardCoSafeLTL(Model model, Rewards modelRewards, Expression expr, MinMax minMax, BitSet statesOfInterest) throws PrismException
	{
		LTLModelChecker mcLtl;
		SMGRewards productRewards;
		StateValues rewardsProduct, rewards;
		SMGModelChecker mcProduct;
		LTLModelChecker.LTLProduct<SMG> product;
		ModelCheckerResult res;

		// For LTL model checking routines
		mcLtl = new LTLModelChecker(this);

		// Build product of SMG and automaton
		AcceptanceType[] allowedAcceptance = {
				AcceptanceType.RABIN,
				AcceptanceType.REACH
		};
		product = mcLtl.constructProductSMG(this, (SMG)model, expr, statesOfInterest, allowedAcceptance);
		
		// Adapt reward info to product model
		productRewards = ((SMGRewards) modelRewards).liftFromModel(product);
		
		// Output product, if required
		if (getExportProductTrans()) {
				mainLog.println("\nExporting product transition matrix to file \"" + getExportProductTransFilename() + "\"...");
				product.getProductModel().exportToPrismExplicitTra(getExportProductTransFilename());
		}
		if (getExportProductStates()) {
			mainLog.println("\nExporting product state space to file \"" + getExportProductStatesFilename() + "\"...");
			PrismFileLog out = new PrismFileLog(getExportProductStatesFilename());
			VarList newVarList = (VarList) modulesFile.createVarList().clone();
			String daVar = "_da";
			while (newVarList.getIndex(daVar) != -1) {
				daVar = "_" + daVar;
			}
			newVarList.addVar(0, new Declaration(daVar, new DeclarationIntUnbounded()), 1, null);
			product.getProductModel().exportStates(Prism.EXPORT_PLAIN, newVarList, out);
			out.close();
		}
		
		// Find accepting states + compute reachability rewards
		BitSet acc;
		if (product.getAcceptance() instanceof AcceptanceReach) {
			// For a DFA, just collect the accept states
			mainLog.println("\nSkipping end component detection since DRA is a DFA...");
			acc = ((AcceptanceReach)product.getAcceptance()).getGoalStates();
		} else {
			// Usually, we have to detect end components in the product
			mainLog.println("\nFinding accepting end components...");
			acc = mcLtl.findAcceptingECStates(product.getProductModel(), product.getAcceptance());
		}
		mainLog.println("\nComputing reachability rewards...");
		mcProduct = new SMGModelChecker(this);
		mcProduct.inheritSettings(this);
		res = mcProduct.computeReachRewards(product.getProductModel(), productRewards, acc, STPGModelChecker.R_CUMULATIVE, minMax.isMin1(), minMax.isMin2(), minMax.getCoalition());
		rewardsProduct = StateValues.createFromDoubleArray(res.soln, product.getProductModel());
		
		// Mapping rewards in the original model
		rewards = product.projectToOriginalModel(rewardsProduct);
		rewardsProduct.clear();
		
		return rewards;
	}
	
	/**
	 * Model check a P operator expression with an exact probability bound (P=p[...])
	 */
	protected StateValues checkExactProbabilityFormula(SMG smg, ExpressionProb expr, Coalition coalition, double p, BitSet statesOfInterest) throws PrismException
	{
		// No support for reward-bounded path formulas (i.e. of the form R(path)~r)
		if (Expression.containsRewardBoundedPathFormula(expr)) {
			throw new PrismException("Reward-bounded path formulas not supported");
		}
		
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
		BitSet prob1 = createSTPGModelChecker().prob1(smg, null, terminal, true, true);
		if (prob1.cardinality() != smg.getNumStates()) {
			throw new PrismException("The game is not stopping. Exact probability queries only work for stopping games");
		}

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
			minmax = checkProbPathFormula(smg, expr.getExpression(), minMax1, statesOfInterest).getDoubleArray();
			if (generateStrategy)
				minStrat = strategy;
			MinMax minMax2 = MinMax.minMin(false, true);
			minMax2.setCoalition(coalition);
			minMax2.setBound(p);
			maxmin = checkProbPathFormula(smg, expr.getExpression(), minMax2, statesOfInterest).getDoubleArray();
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

	protected StateValues checkExactRewardFormula(SMG smg, SMGRewards modelRewards, ExpressionReward exprRew, Coalition coalition, double p) throws PrismException
	{
		// check if the reward is Fc
		ExpressionTemporal exprTemp = null;
		if (exprRew.getExpression() instanceof ExpressionTemporal) {
			exprTemp = (ExpressionTemporal) exprRew.getExpression();
			switch (exprTemp.getOperator()) {
			case ExpressionTemporal.R_Fc:
				break;
			case ExpressionTemporal.P_F:
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
		BitSet prob1 = createSTPGModelChecker().prob1(smg, null, terminal, true, true);
		if (prob1.cardinality() != smg.getNumStates()) {
			throw new PrismException("The game is not stopping. Exact reward queries only work for stopping games");
		}

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
			minmax = checkRewardReach(smg, modelRewards, exprTemp, minMax1, null).getDoubleArray();
			if (generateStrategy)
				minStrat = strategy;
			MinMax minMax2 = MinMax.minMin(false, true);
			minMax2.setCoalition(coalition);
			maxmin = checkRewardReach(smg, modelRewards, exprTemp, minMax2, null).getDoubleArray();
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
