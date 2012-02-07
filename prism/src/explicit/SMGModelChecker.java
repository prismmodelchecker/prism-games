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
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import explicit.rewards.SMGRewards;
import explicit.rewards.STPGRewards;

import parser.ast.Expression;
import parser.ast.ExpressionPATL;
import parser.ast.ExpressionProb;
import parser.ast.ExpressionTemporal;
import parser.type.TypeBool;
import prism.PrismException;

/**
 * Explicit-state model checker for multi-player stochastic games (SMGs).
 */
public class SMGModelChecker extends STPGModelChecker {
	/**
	 * Compute probabilities for the contents of a P operator.
	 */
	protected StateValues checkProbPathFormula(Model model,
			ExpressionPATL exprPATL, boolean min) throws PrismException {
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


//				BitSet target = checkExpression(model,
//						((ExpressionTemporal) expr).getOperand1()).getBitSet();
//				BitSet zero = checkExpression(model,
//						((ExpressionTemporal) expr).getOperand2()).getBitSet();
//				this.computeIntervalSet(min, !min, (STPG) model, target, zero);

		
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
		throw new PrismException(
				"Explicit engine does not yet handle LTL-style path formulas except for GF and FG");

	}

	/**
	 * Compute rewards for the contents of an R operator.
	 */
	protected StateValues checkRewardFormula(Model model,
			SMGRewards modelRewards, ExpressionPATL exprPATL, boolean min)
			throws PrismException {
		// setting coalition parameter
		((SMG) model).setCoalition(exprPATL.getCoalition());

		StateValues rewards = null;
		Expression expr = exprPATL.getExpressionRew().getExpression();

		if (expr instanceof ExpressionTemporal) {
			ExpressionTemporal exprTemp = (ExpressionTemporal) expr;
			switch (exprTemp.getOperator()) {
			case ExpressionTemporal.R_F:
				rewards = checkRewardReach(model, modelRewards, exprTemp, min,
						!min, STPGModelChecker.R_INFINITY);
				break;
			case ExpressionTemporal.R_Fc:
				rewards = checkRewardReach(model, modelRewards, exprTemp, min,
						!min, STPGModelChecker.R_CUMULATIVE);
				break;
			case ExpressionTemporal.R_F0:
				rewards = checkRewardReach(model, modelRewards, exprTemp, min,
						!min, STPGModelChecker.R_ZERO);
				break;
			default:
				throw new PrismException(
						"Explicit engine does not yet handle the "
								+ exprTemp.getOperatorSymbol()
								+ " operator in the R operator");
			}
		}

		if (rewards == null)
			throw new PrismException("Unrecognised operator in R operator");

		return rewards;
	}

	public void computeIntervalSet(boolean min1, boolean min2, STPG stpg,
			BitSet target, BitSet zero) throws PrismException {
		System.out.println(stpg);
		System.out.println(target);
		System.out.println(zero);

		int it;

		int n = stpg.getNumStates();
		List<List<Interval>> intervals = new ArrayList<List<Interval>>(n);
		List<List<Interval>> intervals_ = new ArrayList<List<Interval>>(n);

		// initialising intervals
		List<Interval> init;
		boolean overapprox = false;
		for (int i = 0; i < n; i++) {
			init = new LinkedList<Interval>();
			if (target.get(i))
				init.add(new Interval(1.0, 1.0));
			else if (zero.get(i))
				init.add(new Interval(0.0, 0.0));
			else
			{
				if(overapprox){
					init.add(new Interval(1.0, 1.0));
					init.add(new Interval(0.0, 0.0));
				}
				else
				{
					init.add(new Interval(0.0, 1.0));
				}
				//				if(i==0){
//					init.add(new Interval(0.2, 0.2));
//					init.add(new Interval(0.8, 0.8));
//				}
//				else if(i==1){
//					init.add(new Interval(0.4, 0.4));
//					init.add(new Interval(0.6, 0.6));
//				}
//				else
//				{
//					init.add(new Interval(0.0, 1.0));
//					//init.add(new Interval(0.0, 0.0));
//				}
				
			}
			intervals.add(init);
		}

		BufferedReader in = new BufferedReader(new InputStreamReader(System.in));

		// set number of iterations
		it = 40;
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
