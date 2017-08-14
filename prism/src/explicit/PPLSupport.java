// ==============================================================================
//	
// Copyright (c) 2013-
// Authors:
// * Clemens Wiltsche <clemens.wiltsche@stx.ox.ac.uk> (University of Oxford)
//	
// ------------------------------------------------------------------------------
//	
// This file is part of PRISM.
//	
// PRISM is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//	
// PRISM is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//	
// You should have received a copy of the GNU General Public License
// along with PRISM; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//	
// ==============================================================================

package explicit;

import java.math.BigInteger;
import java.util.AbstractMap.SimpleEntry;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.apache.commons.math3.fraction.BigFraction;

import parma_polyhedra_library.C_Polyhedron;
import parma_polyhedra_library.Coefficient;
import parma_polyhedra_library.Constraint;
import parma_polyhedra_library.Constraint_System;
import parma_polyhedra_library.Generator;
import parma_polyhedra_library.Generator_System;
import parma_polyhedra_library.Generator_Type;
import parma_polyhedra_library.Linear_Expression;
import parma_polyhedra_library.Linear_Expression_Coefficient;
import parma_polyhedra_library.Linear_Expression_Difference;
import parma_polyhedra_library.Linear_Expression_Sum;
import parma_polyhedra_library.Linear_Expression_Times;
import parma_polyhedra_library.Linear_Expression_Unary_Minus;
import parma_polyhedra_library.Linear_Expression_Variable;
import parma_polyhedra_library.Parma_Polyhedra_Library;
import parma_polyhedra_library.Polyhedron;
import parma_polyhedra_library.Relation_Symbol;
import parma_polyhedra_library.Variable;
import prism.PrismException;
import prism.PrismLog;
import prism.PrismUtils;
import explicit.rewards.SMGRewards;

/**
 * Class containing several methods to augment the Parma Polyhedra Library
 * for example by obtaining the coefficients from a linear expression in a list,
 * printing Polyhedra, checking point containment (i.e. whether the bounds are satisfied for an objective),
 * and scaling (discountPareto) and translating sets (add_rewards)
 **/
public class PPLSupport
{
	// Load and initialise the Parma Polyhedra Library (PPL) - used for polyhedra operations
	public static void initPPL() throws PrismException
	{
		try {
			System.loadLibrary("ppl_java");
			Parma_Polyhedra_Library.initialize_library();
		} catch (UnsatisfiedLinkError e1) {
			System.err.println("\nError loading Parma Polyhedra Library:");
			System.err.println(e1);
			throw new PrismException("Parma Polyhedra Library could not be loaded/initialised");
		} catch (Exception e2) {
			System.err.println("\nError loading Parma Polyhedra Library:");
			System.err.println(e2);
			throw new PrismException("Parma Polyhedra Library could not be loaded/initialised");
		}
	}

	public static List<Double> getGeneratorAsVector(Generator g, int n) throws PrismException
	{
		BigInteger div = g.divisor().getBigInteger();
		Map<Variable, BigInteger> map = new HashMap<Variable, BigInteger>();
		getCoefficientsFromLinearExpression(g.linear_expression(), false, BigInteger.ONE, map);
		List<Double> result = new ArrayList<Double>(n);
		for (int i = 0; i < n; i++)
			// initialise with zeros
			result.add(0.0);
		for (Variable v : map.keySet()) {
			if (v != null) {
				BigFraction val = new BigFraction(map.get(v), div);
				result.set((int) v.id(), val.doubleValue());
			}
		}
		return result;
	}

	public static Generator generatorFromPoint(double[] x) throws PrismException
	{
		int n = x.length;
		if (n == 0)
			throw new PrismException("Generator cannot have zero dimensions.");

		Linear_Expression r_num = null;
		BigInteger r_den = BigInteger.ONE;

		for (int i = 0; i < n; i++) {
			BigFraction ri = new BigFraction(x[i]);
			BigInteger num = ri.getNumerator();
			BigInteger den = ri.getDenominator();

			// add r_num/r_den + num/den:
			if (r_num == null) {
				r_num = new Linear_Expression_Times(new Coefficient(num), new Variable(i));
				r_den = den;
			} else {
				Linear_Expression r_num_to_add = new Linear_Expression_Times(new Coefficient(num.multiply(r_den)), new Variable(i));
				if (den.compareTo(BigInteger.ONE) == 0) {
					// (r_num + num*r_den)/r_den
					r_num = new Linear_Expression_Sum(r_num, r_num_to_add);
				} else {
					// (r_num*den + num*r_den)/(r_den*den)
					r_num = new Linear_Expression_Sum(r_num.times(new Coefficient(den)), r_num_to_add);
					r_den = r_den.multiply(den);
				}
			}
		}

		return Generator.point(r_num, new Coefficient(r_den));
	}

	public static Map<Integer, BigInteger> getCoefficients(Linear_Expression le)
	{
		Map<Variable, BigInteger> map = new HashMap<Variable, BigInteger>();
		getCoefficientsFromLinearExpression(le, false, BigInteger.ONE, map);
		Map<Integer, BigInteger> result = new HashMap<Integer, BigInteger>();
		for (Variable v : map.keySet()) {
			if (v != null) {
				result.put((int) v.id(), map.get(v));
			}
		}
		return result;
	}

	public static void getCoefficientsFromLinearExpression(Linear_Expression le, boolean minus, BigInteger coefficient, Map<Variable, BigInteger> result)
	{
		if (le instanceof Linear_Expression_Coefficient) {
			if (result.containsKey(null)) {
				result.put(null,
						minus ? result.get(null).multiply(((Linear_Expression_Coefficient) le).argument().getBigInteger().multiply(coefficient).negate())
								: result.get(null).multiply(((Linear_Expression_Coefficient) le).argument().getBigInteger().multiply(coefficient)));
			} else {
				result.put(null, minus ? ((Linear_Expression_Coefficient) le).argument().getBigInteger().multiply(coefficient).negate()
						: ((Linear_Expression_Coefficient) le).argument().getBigInteger().multiply(coefficient));
			}
		} else if (le instanceof Linear_Expression_Difference) {
			getCoefficientsFromLinearExpression(((Linear_Expression_Difference) le).left_hand_side(), minus, coefficient, result);
			getCoefficientsFromLinearExpression(((Linear_Expression_Difference) le).right_hand_side(), !minus, coefficient, result);
		} else if (le instanceof Linear_Expression_Sum) {
			getCoefficientsFromLinearExpression(((Linear_Expression_Sum) le).left_hand_side(), minus, coefficient, result);
			getCoefficientsFromLinearExpression(((Linear_Expression_Sum) le).right_hand_side(), minus, coefficient, result);
		} else if (le instanceof Linear_Expression_Times) {
			getCoefficientsFromLinearExpression(((Linear_Expression_Times) le).linear_expression(), minus,
					coefficient.multiply(((Linear_Expression_Times) le).coefficient().getBigInteger()), result);
		} else if (le instanceof Linear_Expression_Unary_Minus) {
			getCoefficientsFromLinearExpression(((Linear_Expression_Unary_Minus) le).argument(), !minus, coefficient, result);
		} else if (le instanceof Linear_Expression_Variable) {
			if (result.containsKey(((Linear_Expression_Variable) le).argument())) {
				result.put(((Linear_Expression_Variable) le).argument(), minus ? result.get(((Linear_Expression_Variable) le).argument()).multiply(coefficient)
						.negate() : result.get(((Linear_Expression_Variable) le).argument()).multiply(coefficient));
			} else {
				result.put(((Linear_Expression_Variable) le).argument(), minus ? coefficient.negate() : coefficient);
			}
		}

	}

	public static Entry<double[], double[][]> getConstraintSystem(Constraint_System cs, int n) throws PrismException
	{
		// get the constraint system for P in the form for A x <= b
		int m = cs.size();
		double[][] A = new double[m][n]; // note: initially zero
		double[] b = new double[m]; // note: initially zero

		// note: throw exception on equality constraint
		//       and ignore strict contraints
		int i = 0; // count number of constraint
		for (Constraint c : cs) {
			// get left and right hand sides of constraints A x >/=/< b
			Map<Variable, BigInteger> lhs = new HashMap<Variable, BigInteger>();
			PPLSupport.getCoefficientsFromLinearExpression(c.left_hand_side(), false, BigInteger.ONE, lhs);
			Map<Variable, BigInteger> rhs = new HashMap<Variable, BigInteger>();
			PPLSupport.getCoefficientsFromLinearExpression(c.right_hand_side(), false, BigInteger.ONE, rhs);

			switch (c.kind()) {
			case EQUAL:
				// lhs = rhs
				throw new PrismException("Equality constraints not supported");
			case LESS_THAN:
			case LESS_OR_EQUAL:
				// lhs <(=) rhs
				for (Entry<Variable, BigInteger> iv : lhs.entrySet())
					if (iv.getKey() != null)
						A[i][(int) iv.getKey().id()] += iv.getValue().doubleValue();
					else
						b[i] -= iv.getValue().doubleValue();
				for (Entry<Variable, BigInteger> iv : rhs.entrySet())
					if (iv.getKey() != null)
						A[i][(int) iv.getKey().id()] -= iv.getValue().doubleValue();
					else
						b[i] += iv.getValue().doubleValue();
				break;
			case GREATER_THAN:
			case GREATER_OR_EQUAL:
				// lhs (=)> rhs, same as -lhs <(=) -rhs
				for (Entry<Variable, BigInteger> iv : lhs.entrySet())
					if (iv.getKey() != null)
						A[i][(int) iv.getKey().id()] -= iv.getValue().doubleValue();
					else
						b[i] += iv.getValue().doubleValue();
				for (Entry<Variable, BigInteger> iv : rhs.entrySet())
					if (iv.getKey() != null)
						A[i][(int) iv.getKey().id()] += iv.getValue().doubleValue();
					else
						b[i] -= iv.getValue().doubleValue();
				break;
			default:
				throw new PrismException("Unexpected relation symbol");
			}
			i++;
		}

		return new SimpleEntry(b, A);
	}

	public static boolean checkBound(Pareto Ps, double[] bounds, MultiParameters params)
	{
		for (Polyhedron Xs : Ps.getSets())
			if (checkBound(Xs, bounds, params))
				return true; // union
		// fall through if none of the sets contains the bounds
		return false;
	}

	public static boolean checkBound(Pareto Ps, List<Double> bounds, MultiParameters params)
	{
		for (Polyhedron Xs : Ps.getSets())
			if (checkBound(Xs, bounds, params))
				return true; // union
		// fall through if none of the sets contains the bounds
		return false;
	}

	public static boolean checkBound(Polyhedron Xs, double[] bounds, MultiParameters params)
	{
		if (Xs.is_empty())
			return false; // an empty polyhedron never satisfies any bound

		Linear_Expression le = null;
		BigInteger d = null;

		for (int i = 0; i < bounds.length; i++) {
			BigFraction v_i = new BigFraction(bounds[i]);
			BigInteger num = v_i.getNumerator();
			BigInteger den = v_i.getDenominator();
			if (le == null) {
				le = new Linear_Expression_Times(new Coefficient(num), new Variable(i));
				d = den;
			} else {
				le = new Linear_Expression_Sum(le.times(new Coefficient(den)), new Linear_Expression_Times(new Coefficient(num.multiply(d)), new Variable(i)));
				d = den.multiply(d);
			}
		}
		Generator v = Generator.point(le, new Coefficient(d));
		Generator_System gs = new Generator_System();
		gs.add(v);
		return Xs.contains(new C_Polyhedron(gs));

	}

	public static boolean checkBound(Polyhedron Xs, List<Double> bounds, MultiParameters params)
	{
		if (Xs.is_empty())
			return false; // an empty polyhedron never satisfies any bound

		Linear_Expression le = null;
		BigInteger d = null;

		for (int i = 0; i < bounds.size(); i++) {
			BigFraction v_i = new BigFraction(bounds.get(i));
			BigInteger num = v_i.getNumerator();
			BigInteger den = v_i.getDenominator();
			if (le == null) {
				le = new Linear_Expression_Times(new Coefficient(num), new Variable(i));
				d = den;
			} else {
				le = new Linear_Expression_Sum(le.times(new Coefficient(den)), new Linear_Expression_Times(new Coefficient(num.multiply(d)), new Variable(i)));
				d = den.multiply(d);
			}
		}
		Generator v = Generator.point(le, new Coefficient(d));
		Generator_System gs = new Generator_System();
		gs.add(v);
		return Xs.contains(new C_Polyhedron(gs));

	}

	public static BitSet checkBounds(Pareto[] P, List<Double> bounds, MultiParameters params)
	{
		BitSet result = new BitSet(P.length);
		for (int s = 0; s < P.length; s++)
			result.set(s, checkBound(P[s], bounds, params));
		return result;
	}

	// additionally have a map from concrete states to abstract states f
	public static BitSet checkBounds(Pareto[] P, List<Double> bounds, Map<Integer, Integer> f, MultiParameters params)
	{
		BitSet result = new BitSet(f.size());
		for (int s = 0; s < f.size(); s++)
			result.set(s, checkBound(P[f.get(s)], bounds, params));
		return result;
	}

	// print a single polyhedron with max-generators
	public static void printReachabilityPolyhedron(Pareto[] polyhedra, int dim, int s, PrismLog mainLog) throws PrismException
	{
		if (polyhedra == null || polyhedra.length == 0)
			return;

		int max_points = 0;
		for (int t = 0; t < polyhedra.length; t++) {
			Pareto pt = polyhedra[t];
			if (pt != null && pt.size() > 0) {
				inner: for (int u = 0; u < pt.size(); u++) {
					if (pt.get(u).is_empty())
						continue inner;
					int points_t = polyhedra[t].get(u).minimized_generators().size();
					if (points_t > max_points) {
						max_points = points_t;
					}
				}
			}
		}
		Pareto ps = polyhedra[s];
		if (mainLog != null) {
			mainLog.print(String.format("maxcorners=%d. state %d:\n", max_points, s));
			if (ps == null) {
				mainLog.print("[]\n");
			} else {
				if (ps.size() > 1)
					mainLog.print("[\n");
				for (int u = 0; u < ps.size(); u++)
					mainLog.print(String.format("%s\n", reachabilityPolyhedronToString(ps.get(u), dim)));
				if (ps.size() > 1)
					mainLog.print("]\n");
			}
			mainLog.flush();
		} else {
			System.out.printf("maxcorners=%d. state %d:\n", max_points, s);
			if (ps == null) {
				System.out.printf("[]\n");
			} else {
				if (ps.size() > 1)
					System.out.printf("[\n");
				for (int u = 0; u < ps.size(); u++)
					System.out.printf("%s\n", reachabilityPolyhedronToString(ps.get(u), dim));
				if (ps.size() > 1)
					System.out.printf("]\n");
			}
		}
	}

	// print all polyhedra, including stochastic states, only works for CQs
	public static void printReachabilityPolyhedra(Pareto[] polyhedra, List<Pareto>[] stochasticStates, int dim, PrismLog mainLog) throws PrismException
	{
		for (int s = 0; s < polyhedra.length; s++) {
			mainLog.print(String.format("state=%d: ", s));
			mainLog.print(reachabilityPolyhedronToString(polyhedra[s].get(), dim) + "\n");
			if (stochasticStates != null && stochasticStates[s] != null) {
				for (int t = 0; t < stochasticStates[s].size(); t++) {
					mainLog.print(String.format("    -> %d:" + reachabilityPolyhedronToString(stochasticStates[s].get(t).get(), dim) + "\n", t));
				}
			}
		}
	}

	// print a single polyhedron
	public static String reachabilityPolyhedronToString(Polyhedron polyhedron, int dim) throws PrismException
	{
		if (polyhedron == null)
			return "[]";

		String result = String.format("[");
		for (Generator g : polyhedron.minimized_generators()) {
			if (g.type() == Generator_Type.RAY) {
				result += "r:";
			} else if (g.type() == Generator_Type.LINE) {
				result += "l:";
			}
			result += "[";
			BigInteger den = null;

			Map<Variable, BigInteger> num = new HashMap<Variable, BigInteger>();
			PPLSupport.getCoefficientsFromLinearExpression(g.linear_expression(), false, BigInteger.ONE, num);

			// first evaluate the denominator
			if (g.type() == Generator_Type.RAY | g.type() == Generator_Type.LINE) {
				// normalise rays and lines to one
				den = BigInteger.ONE;
				for (Variable j : num.keySet()) {
					if (j != null & den.compareTo(num.get(j).abs()) < 0) {
						den = num.get(j).abs();
					}
				}
			} else {
				den = g.divisor().getBigInteger();
			}

			boolean init = true;
			for (int i = 0; i < dim; i++) {
				if (!init) {
					result += ", ";
				}
				init = false;
				boolean foundvalue = false;
				for (Variable j : num.keySet()) {
					if (j != null && i == j.id()) {
						BigFraction val = new BigFraction(num.get(j), den);
						result += String.format("%.4f", val.doubleValue());
						foundvalue = true;
						break;
					}
				}
				if (!foundvalue) {
					result += String.format("%.4f", 0.0);
				}
			}
			result += "]";
		}
		result += "]";

		return result;
	}

	public static void printMatlab(Map<Integer, Polyhedron> polyhedra, long dim, int iter)
	{

		int max_points = 0;
		for (int s = 0; s < polyhedra.size(); s++) {
			int points_s = polyhedra.get(s).minimized_generators().size();
			if (points_s > max_points) {
				max_points = points_s;
			}

		}
		System.out.printf("maxpoints{%d} = %d;\n", iter + 1, max_points - dim); // -dim because of rays

		for (int s = 0; s < polyhedra.size(); s++) {
			//System.out.printf("points{%d, %d} = %d;\n", iter+1, s+1, polyhedra.get(s).minimized_generators().size());
			System.out.printf("m{%d, %d} = [", iter + 1, s + 1); // indices must be greater than zero
			boolean init1 = true;
			for (Generator g : polyhedra.get(s).generators()) {
				// ignore rays and lines (for now)
				if (g.type() == Generator_Type.RAY | g.type() == Generator_Type.LINE) {
					continue;
				}
				if (!init1) {
					System.out.printf(" ; ");
				}
				init1 = false;
				BigInteger den = g.divisor().getBigInteger();
				Map<Variable, BigInteger> num = new HashMap<Variable, BigInteger>();
				PPLSupport.getCoefficientsFromLinearExpression(g.linear_expression(), false, BigInteger.ONE, num);
				boolean init2 = true;
				for (int i = 0; i < dim; i++) {
					if (!init2) {
						System.out.printf(", ");
					}
					init2 = false;
					boolean foundvalue = false;
					for (Variable j : num.keySet()) {
						if (j != null && i == j.id()) {
							BigFraction val = new BigFraction(num.get(j), den);
							System.out.printf("%.4f", val.doubleValue());
							foundvalue = true;
							break;
						}
					}
					if (!foundvalue) {
						System.out.printf("%.4f", 0.0);
					}
				}
			}
			System.out.printf("];\n");
		}
	}

	/**
	 * apply discount factor to the Pareto sets
	 * modifies the sets (hence void return value)
	 * scales by beta
	 */
	public static void discountPareto(Pareto[] Xk, double[] beta) throws PrismException
	{
		for (int i = 0; i < beta.length; i++) {
			if (Math.abs(beta[i] - 1.0) >= PrismUtils.epsilonDouble) {
				Variable var = new Variable(i);
				BigFraction val = new BigFraction(beta[i]);
				Linear_Expression expr = new Linear_Expression_Times(new Coefficient(val.getNumerator()), var);
				Coefficient den = new Coefficient(val.getDenominator());
				for (int s = 0; s < Xk.length; s++) {
					Xk[s].get().affine_image(var, expr, den);
				}
			}
		}
	}

	public static void discountPareto(List<Pareto>[] Yk, double[] beta) throws PrismException
	{
		if (Yk == null)
			return;
		for (int i = 0; i < beta.length; i++) {
			if (Math.abs(beta[i] - 1.0) >= PrismUtils.epsilonDouble) {
				Variable var = new Variable(i);
				BigFraction val = new BigFraction(beta[i]);
				Linear_Expression expr = new Linear_Expression_Times(new Coefficient(val.getNumerator()), var);
				Coefficient den = new Coefficient(val.getDenominator());
				for (int s = 0; s < Yk.length; s++) {
					for (Pareto p : Yk[s]) {
						p.get().affine_image(var, expr, den);
					}
				}
			}
		}
	}

	/**
	 * translates polyhedron X by the rewards given in rewards and extra_rewards
	 * does not modify the given polyhedron X
	 *
	 * arguments:
	 * X ... the Polyhedron
	 * s ... state (which state reward to pick if d < 0)
	 * d ... move from s (which transition reward to pick if d >= 0)
	 * rewards ... reward structure of appropriate dimension
	 * extra_rewards ... extra rewards to add at the respective dimensions
	 */
	public static Polyhedron add_rewards(Polyhedron X, int s, int d, List<SMGRewards> rewards, double[] extra_rewards)
	{
		int n = rewards != null ? rewards.size() : extra_rewards.length;

		// check if reward is zero and if so, ignore
		boolean zero_reward = true;
		Iterator<SMGRewards> rewards_iterator = rewards != null ? rewards.iterator() : null;
		for (int i = 0; i < n; i++) {
			SMGRewards reward = rewards_iterator != null ? rewards_iterator.next() : null;
			if (d < 0) {
				double rew = (reward == null ? 0.0 : reward.getStateReward(s)) + (extra_rewards == null ? 0.0 : extra_rewards[i]);
				if (!PrismUtils.doublesAreEqual(rew, 0.0)) {
					zero_reward = false;
					break;
				}
			} else {
				double rew = (reward == null ? 0.0 : reward.getTransitionReward(s, d)) + (extra_rewards == null ? 0.0 : extra_rewards[i]);
				if (!PrismUtils.doublesAreEqual(rew, 0.0)) {
					zero_reward = false;
					break;
				}
			}
		}
		if (zero_reward) {
			return X;
		}

		Generator_System ngs = new Generator_System();
		// first set up the reward vector that should be added to each point generator
		Linear_Expression le = null;
		Coefficient c = new Coefficient(BigInteger.ONE);
		rewards_iterator = rewards != null ? rewards.iterator() : null;
		for (int i = 0; i < n; i++) {
			SMGRewards reward = rewards_iterator != null ? rewards_iterator.next() : null;
			BigFraction r;
			if (d < 0) {
				double rew = (reward == null ? 0.0 : reward.getStateReward(s)) + (extra_rewards == null ? 0.0 : extra_rewards[i]);
				r = new BigFraction(rew);
			} else {
				double rew = (reward == null ? 0.0 : reward.getTransitionReward(s, d)) + (extra_rewards == null ? 0.0 : extra_rewards[i]);
				r = new BigFraction(rew);
			}
			BigInteger num = r.getNumerator();
			BigInteger den = r.getDenominator();
			if (le == null) {
				le = new Linear_Expression_Times(new Coefficient(num), new Variable(i));
				c = new Coefficient(den);
			} else {
				le = new Linear_Expression_Sum(le.times(new Coefficient(den)), new Linear_Expression_Times(new Coefficient(num.multiply(c.getBigInteger())),
						new Variable(i)));
				c = new Coefficient(den.multiply(c.getBigInteger()));
			}
		}
		// now add reward vector to each point generator
		for (Generator g : X.generators()) {
			if (g.type() == Generator_Type.POINT) {
				Linear_Expression nle = new Linear_Expression_Sum(le.times(g.divisor()), g.linear_expression().times(c));
				Coefficient nc = new Coefficient(g.divisor().getBigInteger().multiply(c.getBigInteger()));
				ngs.add(Generator.point(nle, nc));
			} else {
				ngs.add(g);
			}
		}

		return new C_Polyhedron(ngs);
	}

	/**
	 * cut away everything except the negative orthant
	 **/
	public static void cutPositive(Polyhedron p)
	{
		long n = p.space_dimension();
		for (int i = 0; i < n; i++) {
			Linear_Expression lhs = new Linear_Expression_Times(new Coefficient(BigInteger.ONE), new Variable(i));
			Linear_Expression rhs = new Linear_Expression_Coefficient(new Coefficient(BigInteger.ZERO));
			Constraint c = new Constraint(lhs, Relation_Symbol.LESS_OR_EQUAL, rhs);
			p.add_constraint(c);
		}
	}

	/**
	 * cut away everything except the box in negative orthant of maximum extent -M in all dimensions
	 **/
	public static void cutBox(Polyhedron p, long M) throws PrismException
	{
		long n = p.space_dimension();
		for (int i = 0; i < n; i++) {
			// cut positive values
			Linear_Expression lhs_0 = new Linear_Expression_Times(new Coefficient(BigInteger.ONE), new Variable(i));
			Linear_Expression rhs_0 = new Linear_Expression_Coefficient(new Coefficient(BigInteger.ZERO));
			Constraint c_0 = new Constraint(lhs_0, Relation_Symbol.LESS_OR_EQUAL, rhs_0);
			p.add_constraint(c_0);

			// cut values less than -M
			Linear_Expression lhs_M = new Linear_Expression_Times(new Coefficient(BigInteger.ONE), new Variable(i));
			Linear_Expression rhs_M = new Linear_Expression_Coefficient(new Coefficient(BigInteger.valueOf(-M)));
			Constraint c_M = new Constraint(lhs_M, Relation_Symbol.GREATER_OR_EQUAL, rhs_M);
			p.add_constraint(c_M);
		}

		// down-close again
		if (!p.is_empty()) {
			for (int i = 0; i < n; i++) {
				Linear_Expression ray = new Linear_Expression_Times(new Coefficient((BigInteger.ONE).negate()), new Variable(i));
				p.add_generator(Generator.ray(ray));
			}
		}
	}

}
