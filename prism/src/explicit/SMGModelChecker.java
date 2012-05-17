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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

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

				// the case with arbitrary objectives
				// Expression exp = ((ExpressionTemporal) expr).getOperand1();
				// final List<Expression> exps = new ArrayList<Expression>(5);
				// ((ExpressionUnaryOp) exp).getOperand().accept(
				// new ASTTraverse() {
				// @Override
				// public Object visit(ExpressionBinaryOp e)
				// throws PrismLangException {
				// visitPre(e);
				// // System.out.println(e);
				// // if(!exps.contains(e.getOperand1()))
				//
				// // if(!exps.contains(e.getOperand2()))
				// if (e.getOperand2() instanceof ExpressionUnaryOp
				// && !e.getOperand2().toString()
				// .startsWith("((")) {
				// // System.out.println("Adding " +
				// // e.getOperand2());
				// exps.add(e.getOperand2());
				//
				// }
				// if (e.getOperand1() instanceof ExpressionUnaryOp
				// && !e.getOperand1().toString()
				// .startsWith("((")) {
				// // System.out.println("Adding " +
				// // e.getOperand1());
				// exps.add(e.getOperand1());
				// }
				// e.getOperand1().accept(this);
				// e.getOperand2().accept(this);
				// visitPost(e);
				// return null;
				// }
				// });
				//
				// Collections.reverse(exps);
				//
				// // just the case with 2 objectives
				// // BitSet t1 = checkExpression(model,
				// // ((ExpressionTemporal) expr).getOperand1()).getBitSet();
				// // BitSet t2 = checkExpression(model,
				// // ((ExpressionTemporal) expr).getOperand2()).getBitSet();
				//
				// List<BitSet> targets = new ArrayList<BitSet>(2);
				// for (Expression e : exps) {
				// targets.add(checkExpression(model, e).getBitSet());
				// }
				//
				// List<Set<ReachTuple>> result =
				// this.computeReachabilityTuples(
				// min, !min, (STPG) model, targets);
				//
				// System.out.println("Printing results..");
				// for (int i = 0; i < result.size(); i++) {
				// System.out.println(i + " " + result.get(i));
				// }
				// System.out.println("------------------");

				// this.computeIntervalSet(min, !min, (STPG) model, t1, t2);

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
	protected StateValues checkRewardFormula(Model model, SMGRewards modelRewards, ExpressionPATL exprPATL, boolean min)
			throws PrismException
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
				throw new PrismException("Explicit engine does not yet handle the " + exprTemp.getOperatorSymbol()
						+ " operator in the R operator");
			}
		}

		if (rewards == null)
			throw new PrismException("Unrecognised operator in R operator");

		return rewards;
	}

	protected StateValues checkExactProbabilityFormula(Model model, ExpressionPATL expr, double p)
			throws PrismException
	{
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
				throw new PrismException(
						"The game is not stopping. The exact probability queries only work for stopping games");

		// 2) computing minmax and maxmin values for all states
		double[] minmax = null, maxmin = null; // see the do loop below

		// 3) removing states from the game which have minmax>maxmin
		//model.
		int n = model.getNumStates();
		boolean repeat;
		BitSet removed = new BitSet(n), removedNew = new BitSet(n);
		STPG stpg = ((STPG) model);

		Strategy minStrat = null;
		Strategy maxStrat = null;

		this.storeStrategyValues = generateStrategy;
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
			strategy = new ExactValueStrategy(minStrat, maxStrat, p, (STPG) model);
		}

		return StateValues.createFromBitSet(ret, model);
	}

	/**
	 * Computes convex hull of all achievable reachability tuples for targets
	 * 
	 * @param min1
	 *            false if player 1 maximising, true otherwise
	 * @param min2
	 *            false if player 2 maximising, true otherwise
	 * @param stpg
	 *            stochastic two player game model
	 * @param targets
	 *            set of reachability targets (only disjoint is supprted atm)
	 * @return
	 * @throws PrismException
	 */
	public List<Set<ReachTuple>> computeReachabilityTuples(boolean min1, boolean min2, STPG stpg, List<BitSet> targets)
			throws PrismException
	{

		int gameSize = stpg.getStatesList().size();

		// checking for disjointness
		BitSet bs = new BitSet(gameSize);
		bs.set(0, gameSize - 1, true);
		for (BitSet tar : targets)
			bs.and(tar);
		if (bs.cardinality() != 0)
			throw new PrismException("Target sets have to be disjoint!");

		// intialise results storage
		List<Set<ReachTuple>> result = new ArrayList<Set<ReachTuple>>(gameSize);
		for (int i = 0; i < gameSize; i++)
			result.add(new HashSet<ReachTuple>());

		// compute initial tuples for states
		// System.out.println("Computing initial values for states..");
		double[] maxmin;
		ReachTuple tuple;
		BitSet ones = new BitSet(gameSize);
		ones.set(0, gameSize - 1, true);

		// 0 - from initial intervals
		// 1 - from 1
		int initValues = 1;

		for (int t = 0; t < targets.size(); t++) {

			maxmin = computeUntilProbs(stpg, ones, targets.get(t), false, true, 0).soln;
			// System.out.println(Arrays.toString(maxmin));
			for (int s = gameSize - targets.size(); s < gameSize; s++) {
				// if(maxmin[s] == 0) continue;
				tuple = new ReachTuple(targets.size());
				tuple.setValue(t, maxmin[s]);
				result.get(s).add(tuple);
			}
			for (int s = 0; s < gameSize - targets.size(); s++) {
				// if(maxmin[s] == 0) continue;
				tuple = new ReachTuple(targets.size());
				switch (initValues) {
				case 0:
					tuple.setValue(t, maxmin[s]);
					break;
				case 1:
					tuple.setValue(t, 1);
					break;
				// case 2 : tuple.setValue(t, 0); break;
				}
				result.get(s).add(tuple);
			}
			// TODO consider states which are in multiple targets!
		}

		// System.out.println(result.toString());
		// starting value iterationinit
		int maxIter = 20;
		BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
		Set<ReachTuple> sample;
		while (maxIter-- > 0) {
			System.out.println(result.toString());
			System.out.println("--------------");
			sample = result.get(0);
			System.out.println(sample.size());
			System.out.println(sample);
			System.out.println(result.get(1));
			System.out.println("--------------");

			result = ((SMG) stpg).mvMultiObjective(min1, min2, result);
			//			 try {
			//			 in.readLine();
			//			 } catch (IOException e) {
			//			 // TODO Auto-generated catch block
			//			 e.printStackTrace();
			//			 }

			// perform upward perturbation when half way
			//			if (maxIter == 10)
			//				for (Set<ReachTuple> ts : result)
			//					for (ReachTuple t : ts)
			//						if (initValues == 0)
			//							t.perturbateUp();
			//						else
			//							t.perturbateDown();

		}

		return result;
	}

	public void computeIntervalSet(boolean min1, boolean min2, STPG stpg, BitSet target, BitSet zero)
			throws PrismException
	{
		System.out.println(stpg);
		System.out.println(target);
		System.out.println(zero);

		int it;

		int n = stpg.getNumStates();
		List<List<Interval>> intervals = new ArrayList<List<Interval>>(n);
		List<List<Interval>> intervals_ = new ArrayList<List<Interval>>(n);

		// initialising intervals
		// compute initial tuples for states
		System.out.println("Computing initial values for states..");
		double[] minmin;
		double[] minmax;
		double[] maxmin;
		double[] maxmax;
		BitSet ones = new BitSet(n);
		ones.set(0, n - 1, true);

		minmin = computeUntilProbs(stpg, ones, target, true, true, 0).soln;
		minmax = computeUntilProbs(stpg, ones, target, true, false, 0).soln;
		maxmin = computeUntilProbs(stpg, ones, target, false, true, 0).soln;
		maxmax = computeUntilProbs(stpg, ones, target, false, false, 0).soln;

		List<Interval> init;
		for (int s = 0; s < n; s++) {
			init = new LinkedList<Interval>();
			init.add(new Interval(minmin[s], minmax[s]));
			init.add(new Interval(maxmin[s], maxmax[s]));
			intervals.add(init);
		}
		// TODO consider states which are in multiple targets!

		BufferedReader in = new BufferedReader(new InputStreamReader(System.in));

		// set number of iterations
		it = 15;
		System.out.println(intervals);
		while (true) {
			intervals_ = ((SMG) stpg).mvMultIntervals(min1, min2, intervals);

			try {
				in.readLine();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

			// do some checks...
			System.out.println(intervals_);

			intervals = intervals_;
		}
	}
}
