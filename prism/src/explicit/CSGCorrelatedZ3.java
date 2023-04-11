//==============================================================================
//	
//	Copyright (c) 2020-
//	Authors:
//	* Gabriel Santos <gabriel.santos@cs.ox.ac.uk> (University of Oxford)
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
import java.util.HashMap;

import com.microsoft.z3.AlgebraicNum;
import com.microsoft.z3.ArithExpr;
import com.microsoft.z3.BoolExpr;
import com.microsoft.z3.Context;
import com.microsoft.z3.Expr;
import com.microsoft.z3.IntExpr;
import com.microsoft.z3.Model;
import com.microsoft.z3.Optimize;
import com.microsoft.z3.RatNum;
import com.microsoft.z3.RealExpr;
import com.microsoft.z3.Status;
import com.microsoft.z3.Version;

import prism.PrismException;

/**
 * Z3-based implementation for correlated equilibria computation
 * @author Gabriel Santos
 *
 */
public class CSGCorrelatedZ3 implements CSGCorrelated {

	private RealExpr[] vars;
	private RealExpr[] payoff_vars;
	private RealExpr obj_var;
	private ArithExpr[] payoffs;
	
	private HashMap<String, String> cfg;
	private Context ctx;
	private Optimize s;
	
	private IntExpr zero;
	private IntExpr one;
	
	private String name;
	private int n_coalitions;
	
	/**
	 * Creates a new CSGCorrelatedZ3 (without initialisation)
	 */
	public CSGCorrelatedZ3() {
		cfg = new HashMap<String, String>();
		cfg.put("model", "true");
		ctx = new Context(cfg);
		s = ctx.mkOptimize();
		name = Version.getFullVersion();
	}
	
	/**
	 * Creates a new CSGCorrelatedZ3
	 * @param n_entries Number of entries in the utility table
	 * @param n_coalitions Number of coalitions
	 */
	public CSGCorrelatedZ3(int n_entries, int n_coalitions) {
		cfg = new HashMap<String, String>();
		cfg.put("model", "true");
		cfg.put("auto_config", "true");
		ctx = new Context(cfg);
		s = ctx.mkOptimize();
		name = Version.getFullVersion();
		zero = ctx.mkInt(0);
		one = ctx.mkInt(1);
		vars = new RealExpr[n_entries];
		payoffs = new ArithExpr[n_coalitions];
		this.n_coalitions = n_coalitions;
		for (int i = 0; i < n_entries; i++) {
			vars[i] = ctx.mkRealConst("v" + i);
			s.Add(ctx.mkLe(vars[i], one));
			s.Add(ctx.mkGe(vars[i], zero));
		}
		payoff_vars = new RealExpr[n_coalitions];
		obj_var = ctx.mkRealConst("ob");
	}
	
	public EquilibriumResult computeEquilibrium(HashMap<BitSet, ArrayList<Double>> utilities, 
												ArrayList<ArrayList<HashMap<BitSet, Double>>> ce_constraints,
												ArrayList<ArrayList<Integer>> strategies,
												HashMap<BitSet, Integer> ce_var_map, int type) {
		EquilibriumResult result = new EquilibriumResult();
		ArrayList<Double> payoffs_result = new ArrayList<Double>();
		ArrayList<Distribution<Double>> strategy_result = new ArrayList<>();
		Distribution<Double> d = new Distribution<>();
		ArithExpr expr;
		BitSet is, js;
		double u;
		int c, q, r;
		is = new BitSet();
		js = new BitSet();
		
		s.Push();
		expr = ctx.mkInt(0);
		
		// Builds a payoff expression for each player, and one for the sum of all payoffs
		for (int i = 0; i < n_coalitions; i++) {
			payoffs[i] = zero;
		}
		for (BitSet e : utilities.keySet()) {
			u = 0.0;
			for (c = 0; c < n_coalitions; c++) {
				u += utilities.get(e).get(c);
				payoffs[c] = ctx.mkAdd(payoffs[c], ctx.mkMul(vars[ce_var_map.get(e)], ctx.mkReal(String.valueOf(utilities.get(e).get(c)))));
			}
			expr = ctx.mkAdd(expr, ctx.mkMul(vars[ce_var_map.get(e)], ctx.mkReal(String.valueOf(u))));
		}
		for (c = 0; c < n_coalitions; c++) {
			payoff_vars[c] = ctx.mkRealConst("p " + c);
			s.Add(ctx.mkEq(payoff_vars[c], payoffs[c]));
		}
	
		// In case of the fair variant (minimise the difference between the largest and smallest payoffs)
		// we first need to add auxiliary constraints
		// The default case is social-welfare, in which we just maximise the sum
		switch (type) {
			case CSGModelCheckerEquilibria.FAIR : {
				RealExpr ph;
				RealExpr pl;
				ph = ctx.mkRealConst("ph"); // Highest payoff
				pl = ctx.mkRealConst("pl"); // Lowest payoff
				BoolExpr bh;
				BoolExpr bl;
				// Additional constraints for fair
				for (int i = 0; i < n_coalitions; i++) {
					bh = ctx.mkTrue();
					bl = ctx.mkTrue();
					for (int j = 0; j < n_coalitions; j++) {
						if (i != j) {
							bh = ctx.mkAnd(bh, ctx.mkGe(payoffs[i], payoffs[j]));
							bl = ctx.mkAnd(bl, ctx.mkLe(payoffs[i], payoffs[j]));
						}
					}
					s.Add(ctx.mkImplies(bh, ctx.mkEq(ph, payoffs[i])));
					s.Add(ctx.mkImplies(bl, ctx.mkEq(pl, payoffs[i])));
				}
				s.MkMinimize(ctx.mkSub(ph, pl));
				break;
			}
			default : {
				// Primary objective is maximising the sum
				s.Add(ctx.mkEq(obj_var, expr));
				s.MkMaximize(expr);
			}
		}
	
		// Lower priority objectives of maximising the payoff of each player in decreasing order
		for (c = 0; c < n_coalitions; c++) {
			s.MkMaximize(payoffs[(c)]);
		}
		
		// Constraints for correlated equilibria
		for (c = 0; c < n_coalitions; c++) {
			for (q = 0; q < strategies.get(c).size(); q++) {
				expr = zero;
				is.clear();
				is.set(strategies.get(c).get(q));
				for (r = 0; r < strategies.get(c).size(); r++) {
					js.clear();
					for (BitSet e : ce_constraints.get(c).get(q).keySet()) {
						is.or(e);
						js.or(e);
						if (q != r) {
							js.set(strategies.get(c).get(r));
							expr = ctx.mkAdd(expr, ctx.mkMul(vars[ce_var_map.get(is)], 
															 ctx.mkSub(ctx.mkReal(String.valueOf(utilities.get(is).get(c))), 
																	   ctx.mkReal(String.valueOf(utilities.get(js).get(c))))));
						}
						is.andNot(e);
						js.andNot(e);
					}
					if (q != r) {
						s.Add(ctx.mkGe(expr, zero));
					}
				}
			}
		}
		
		// Sets unused variables to zero
		for (c = utilities.size(); c < vars.length; c++) {
			s.Add(ctx.mkEq(vars[c], zero));
		}
		
		expr = ctx.mkInt(0);
		for (c = 0; c < vars.length; c++) {
			expr = ctx.mkAdd(expr, vars[c]);
		}
		s.Add(ctx.mkEq(expr, one));
				
		/*
		for (BoolExpr e : s.getAssertions()) {
			System.out.println(e);
		} 
		*/
		
		// If and optimal solution is found, set the values and strategies
		if (s.Check() == Status.SATISFIABLE) {
			for (c = 0; c < vars.length; c++) {
				d.add(c, getDoubleValue(s.getModel(), vars[c]));
			}
			for (c = 0; c < payoff_vars.length; c++) {
				payoffs_result.add(getDoubleValue(s.getModel(), payoff_vars[c]));
			}
			result.setStatus(CSGModelCheckerEquilibria.CSGResultStatus.SAT);
			result.setPayoffVector(payoffs_result);
			strategy_result.add(d);
			result.setStrategy(strategy_result);
		}
		else {
			result.setStatus(CSGModelCheckerEquilibria.CSGResultStatus.UNSAT);
		}
		
		s.Pop();
		return result;
	}
	
	/**
	 * Return a double value for a given expression (usually variable), converting from BigInt fractions
	 * @param model The SMT model
	 * @param expr The SMT expression
	 * @return
	 */
	public double getDoubleValue(Model model, Expr expr) {
		RatNum v1;
		AlgebraicNum v2;
		if(model.getConstInterp(expr) instanceof RatNum) {
			v1 = (RatNum) model.getConstInterp(expr);
			return (Double) (v1.getBigIntNumerator().doubleValue() / v1.getBigIntDenominator().doubleValue());
		}
		else if (model.getConstInterp(expr) instanceof AlgebraicNum) {
			v2 = (AlgebraicNum) model.getConstInterp(expr);
			v1 = v2.toUpper(12);
			return (Double) (v1.getBigIntNumerator().doubleValue() / v1.getBigIntDenominator().doubleValue());
		}
		else
			return Double.NaN;
	}

	public String getSolverName() {
		return name;
	}

	@Override
	public void clear() {
		
	}
	
	@Override
	public void printModel() {
		// TODO Auto-generated method stub
		
	}
	
}
