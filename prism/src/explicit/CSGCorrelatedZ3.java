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
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.Map.Entry;
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
import prism.Pair;
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
	
	public HashSet<ArrayList<Integer>> buildAllSupports(ArrayList<ArrayList<Integer>> strategies, BitSet players) {
		HashSet<ArrayList<Integer>> supports = new HashSet<ArrayList<Integer>>();
		ArrayList<Integer> support;
		for(int p = players.nextSetBit(0); p >= 0; p = players.nextSetBit(p + 1)) {
			for (int s : strategies.get(p)) {
				support = new ArrayList<Integer>(); 
				for (int i = 0; i < n_coalitions; i++) 
					support.add(-1);
				support.set(p, s);
				buildAllSupportsAux(supports, strategies, players, support, p);
			}	
		}
		return supports;
	}
	
	public void buildAllSupportsAux(HashSet<ArrayList<Integer>> supports, ArrayList<ArrayList<Integer>> strategies, BitSet players, ArrayList<Integer> support, int p) {
		int cardinality = 0;
		for (int s : support) {
			if (s != -1)
				cardinality = cardinality + 1;
		} 
		if (players.nextSetBit(p + 1) >= 0) { 
			for (int s : strategies.get(p + 1)) {
				ArrayList<Integer> curr = new ArrayList<Integer>(support);
				curr.set(p + 1, s);
				buildAllSupportsAux(supports, strategies, players, curr, players.nextSetBit(p + 1));
			}
		}
		else if (cardinality == players.cardinality())  {
			supports.add(new ArrayList<Integer>(support));
		}
	}
	
	public void buildSubGames(Set<BitSet> games, BitSet sp, int p) {
		BitSet prod = new BitSet();
		prod.set(p);
		games.add((BitSet) prod.clone());
		for(int c = sp.nextSetBit(0); c >= 0; c = sp.nextSetBit(c + 1)) {
			BitSet newprod = new BitSet();
			newprod.or(prod);
			newprod.set(c);
			games.add(newprod);
		}
	}
	
	public BitSet buildBitSetStrategy(ArrayList<Integer> support) {
		BitSet result = new BitSet();
		for (int s : support) {
			if (s != -1)
				result.set(s);
		}
		return result;
	}
	
	public ArrayList<Integer> buildStrategyVector(ArrayList<BitSet> strategy_sets, BitSet support) {
		ArrayList<Integer> result = new ArrayList<Integer>(); 
		BitSet tmp = new BitSet();
		for (int c = 0; c < n_coalitions; c++) {
			tmp.clear();
			tmp.or(strategy_sets.get(c));
			tmp.and(support);
			if (!tmp.isEmpty())
				result.add(tmp.nextSetBit(0));
			else
				result.add(-1);
		}
		return result;
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
			// System.out.println("sat\n");			
			// System.out.println(s.getModel());
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
			// System.out.println(strategy_result);
		}
		else {
			result.setStatus(CSGModelCheckerEquilibria.CSGResultStatus.UNSAT);
		}
		
		s.Pop();
		return result;
	}
	
	public EquilibriumResult __computeEquilibrium(HashMap<BitSet, ArrayList<Double>> utilities, 
												ArrayList<ArrayList<HashMap<BitSet, Double>>> ce_constraints,
												ArrayList<ArrayList<Integer>> strategies,
												HashMap<BitSet, Integer> ce_var_map, int type) {
		EquilibriumResult result = new EquilibriumResult();
		ArrayList<Double> payoffs_result = new ArrayList<Double>();
		ArrayList<Distribution<Double>> strategy_result = new ArrayList<Distribution<Double>>();
		Distribution d = new Distribution();
		ArithExpr expr;
		BitSet is, js, ts;
		double u;
		int i, c, n, q, r;
		is = new BitSet();
		js = new BitSet();
		
		s.Push();
		
		/// Print-outs
		System.out.println("\n# utilities");
		System.out.println(utilities); 
		System.out.println("\n# strategies");
		System.out.println(strategies); 
		System.out.println("\n# ce_var_map");
		System.out.println(ce_var_map);
		System.out.println("\n# ce_constraints");
		System.out.println(ce_constraints);
		
		// Utilities with vectors
		HashMap<ArrayList<Integer>, ArrayList<Double>> vectorised_utilities = 
															new HashMap<ArrayList<Integer>, ArrayList<Double>>();
		ArrayList<BitSet> strategy_sets = new ArrayList<BitSet>();
		ArrayList<Integer> strategy_vector = null;
		for (c = 0; c < n_coalitions; c++) {
			ts = new BitSet();
			for (q = 0; q < strategies.get(c).size(); q++) {
				ts.set( strategies.get(c).get(q));
			}
			strategy_sets.add(ts);
		}
		for (Entry<BitSet, ArrayList<Double>> entry : utilities.entrySet()) {
			strategy_vector = new ArrayList<Integer>();
			for (c = 0; c < n_coalitions; c++) { 
				ts = new BitSet();
				ts.or(strategy_sets.get(c));
				ts.and(entry.getKey());
				strategy_vector.add(ts.nextSetBit(0));
			}
			vectorised_utilities.put(strategy_vector, entry.getValue());
		}
		
		// Constraints with vectors
		ArrayList<ArrayList<HashMap<ArrayList<Integer>, Double>>> vectorised_constraints = 
															new ArrayList<ArrayList<HashMap<ArrayList<Integer>, Double>>>();
		for (c = 0; c < n_coalitions; c++) {
			vectorised_constraints.add(new ArrayList<HashMap<ArrayList<Integer>, Double>>());
			for (q = 0; q < strategies.get(c).size(); q++) {
				vectorised_constraints.get(c).add(new HashMap<ArrayList<Integer>, Double>());
				for (Entry<BitSet, Double> e : ce_constraints.get(c).get(q).entrySet()) {
					strategy_vector = buildStrategyVector(strategy_sets, e.getKey());
					strategy_vector.set(c, strategies.get(c).get(q));
					vectorised_constraints.get(c).get(q).put(strategy_vector, e.getValue());
				}
			}
		}
		
		// Builds subsets S \in N
		BitSet players = new BitSet();
		HashSet<BitSet> ss = new HashSet<BitSet>();
		for (c = 0; c < n_coalitions; c++) {
			players.set(c);
		}
		
		for (c = 0; c < n_coalitions; c++) {
			ts = (BitSet) players.clone();
			ts.clear(c);
			buildSubGames(ss, ts, c);
		}	
		System.out.println("\n# Ss");
		ss.add(new BitSet());
		System.out.println(ss);
		
		// Builds all possible cs in Cs
		HashMap<BitSet, ArrayList<ArrayList<Integer>>> cs = new HashMap<BitSet, ArrayList<ArrayList<Integer>>>();
		
		for (BitSet se : ss) {
			cs.put(se, new ArrayList<ArrayList<Integer>>());
			for (ArrayList<Integer> e : buildAllSupports(strategies, se)) {
				cs.get(se).add(e);
			} 
		}
		cs.get(new BitSet()).add(new ArrayList<Integer>());
		System.out.println(cs);
		
		// Builds a payoff expression for each player, and one for the sum of all payoffs
		for (i = 0; i < n_coalitions; i++) {
			payoffs[i] = zero;
		}
		
		// Builds the set of eta variables
		HashMap<Pair<BitSet, ArrayList<Integer>>, RealExpr> eta = new HashMap<Pair<BitSet, ArrayList<Integer>>, RealExpr>();
		HashMap<String, Pair<BitSet, ArrayList<Integer>>> eta_names = new HashMap<String, Pair<BitSet, ArrayList<Integer>>>();
		RealExpr eta_var;
		n = 0;
		expr = zero;
		for (Entry<BitSet, ArrayList<ArrayList<Integer>>> ce : cs.entrySet()) {
			// System.out.println("# " + ce.getValue());
			for (Entry<BitSet, Integer> joint_action : ce_var_map.entrySet()) {
				for (ArrayList<Integer> e : ce.getValue()) {
					// System.out.println(e);
					eta_var = ctx.mkRealConst("n" + n);
					eta.put(new Pair<BitSet, ArrayList<Integer>>(joint_action.getKey(), e), eta_var);	
					eta_names.put("n" + n, new Pair<BitSet, ArrayList<Integer>>(joint_action.getKey(), e));
					s.Add(ctx.mkLe(eta_var, one));
					s.Add(ctx.mkGe(eta_var, zero));
					expr = ctx.mkAdd(expr, eta_var);
					n++;
					for (c = 0; c < n_coalitions; c++) {
						payoffs[c] = ctx.mkAdd(payoffs[c], ctx.mkMul(eta_var, ctx.mkReal(String.valueOf(utilities.get(joint_action.getKey()).get(c)))));
					}
				}
			}
		}
		s.Add(ctx.mkEq(expr, one));
		
		expr = zero;
		for (c = 0; c < n_coalitions; c++) {
			payoff_vars[c] = ctx.mkRealConst("p " + c);
			expr = ctx.mkAdd(expr, payoffs[c]);
			s.Add(ctx.mkEq(payoff_vars[c], payoffs[c]));
		}
		s.MkMinimize(expr);
		
		System.out.println("\n# eta");
		System.out.println(n + " variables");
		System.out.println(eta);
		System.out.println();
		
		// Adds constraints for epsilon-correlated equilibria
		Pair<BitSet, ArrayList<Integer>> eta_index;
		ArrayList<Integer> is_v;
		ArrayList<Integer> js_v;
		ArrayList<Integer> c1_es;
		ArrayList<Integer> c2_es;
				
		for (c = 0; c < n_coalitions; c++) {
			System.out.println("-- player " + c);
			for (q = 0; q < strategies.get(c).size(); q++) { // c_i
				expr = zero;
				is.clear();
				is.set(strategies.get(c).get(q));				
				for (r = 0; r < strategies.get(c).size(); r++) {// e_i
					js.clear();
					for (BitSet e : ce_constraints.get(c).get(q).keySet()) {
						is.or(e);
						js.or(e);
						
						is_v = buildStrategyVector(strategy_sets, is);
						c1_es = buildStrategyVector(strategy_sets, is);
						c2_es = buildStrategyVector(strategy_sets, is);
						if (q != r) {
							js.set(strategies.get(c).get(r));
							
							js_v =  buildStrategyVector(strategy_sets, js);
							System.out.println("is " + is.toString() + ", js " + js.toString());
							System.out.println("is_v " + is_v.toString() + ", js_v " + js_v.toString());
							for (Entry<BitSet, ArrayList<ArrayList<Integer>>> ce : cs.entrySet()) {
								if (!ce.getKey().get(c)) {
									System.out.println(ce);
									for (ArrayList<Integer> es : ce.getValue()) {
										eta_index = new Pair<BitSet, ArrayList<Integer>>(is, es);
										eta_var = eta.get(eta_index);
										System.out.println(eta_var);
										// "is" is c, the joint action
										// replacement
										
										System.out.println("es " + es.toString());
										for (i = 0; i < n_coalitions; i++) {
											if (!es.isEmpty() && es.get(i) != -1) {
												c1_es.set(i, es.get(i));
												c2_es.set(i, es.get(i));
											}
											if (i == c) 
												c2_es.set(i, strategies.get(c).get(r));
										}
										
										System.out.println("c1_es " + c1_es.toString() + ", c2_es " + c2_es.toString());										
										expr = ctx.mkAdd(expr, ctx.mkMul(eta_var, 
															ctx.mkSub(ctx.mkReal(String.valueOf(vectorised_utilities.get(c1_es).get(c))),
																	  ctx.mkReal(String.valueOf(vectorised_utilities.get(c2_es).get(c))))));
									}
								}
							}
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
		
		// Constraints for 2.4
		BitSet sui;
		for (Entry<BitSet, ArrayList<Double>> o : utilities.entrySet()) {
			for (c = 0; c < n_coalitions; c++) {
				for (Entry<BitSet, ArrayList<ArrayList<Integer>>> ce : cs.entrySet()) {
					if (!ce.getKey().get(c)) {
						sui = new BitSet();
						sui.or(ce.getKey());
						sui.set(c);
						for (ArrayList<Integer> esui : cs.get(sui)) {
							for (ArrayList<Integer> es : ce.getValue()) {
								s.Add(ctx.mkImplies(
									ctx.mkGt(eta.get(new Pair<BitSet, ArrayList<Integer>>(o.getKey(), es)), zero), 
									ctx.mkGt(eta.get(new Pair<BitSet, ArrayList<Integer>>(o.getKey(), esui)), zero)));
							}	
						}
					}
				}
			}
		}
		
		// Constraints for 2.3
		RealExpr epsilon = ctx.mkRealConst("e");
		expr = zero;
		for (Entry<BitSet, ArrayList<Double>> o : utilities.entrySet()) {
			for (c = 0; c < n_coalitions; c++) {
				for (Entry<BitSet, ArrayList<ArrayList<Integer>>> ce : cs.entrySet()) {
					if (!ce.getKey().get(c)) {
						
						for (q = 0; q < strategies.get(c).size(); q++) { // c_i
							
						}
						
					}
				}
			}
		}
		
		// s.Add(ctx.mkEq(eta.get(eta_names.get("n15")), zero));
		// s.Add(ctx.mkEq(eta.get(eta_names.get("n22")), one));
		
		if (s.Check() == Status.SATISFIABLE) {
			System.out.println("sat\n");			
			System.out.println(s.getModel());
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
		
		// System.out.println(eta_names.get("n15"));
		// System.out.println(eta_names.get("n22"));
		
		System.out.println();
		
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
