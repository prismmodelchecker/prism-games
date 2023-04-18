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

import com.microsoft.z3.Model;
import com.microsoft.z3.*;
import explicit.CSGModelCheckerEquilibria.CSGResultStatus;
import prism.Pair;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;

public class CSGSupportEnumerationZ3 implements CSGSupportEnumeration
{
	private HashMap<String, String> cfg;
	private Context ctx;
	private Solver s;

	private HashMap<Integer, HashMap<Integer, ArithExpr>> assertions;
	private ArrayList<ArrayList<ArithExpr>> strategies;
	private ArrayList<ArrayList<Integer>> indexes;
	private ArrayList<ArithExpr> payoffs;

	private IntNum zero;
	private IntNum one;

	private BitSet players;
	private BoolExpr support;

	private int numCoalitions;
	private Params params;

	private ArrayList<Distribution<Double>> strat;

	public CSGSupportEnumerationZ3()
	{

	}

	public CSGSupportEnumerationZ3(int[] maxActions, int numCoalitions)
	{
		cfg = new HashMap<String, String>();
		cfg.put("model", "true");
		ctx = new Context(cfg);
		s = ctx.mkSolver();

		params = ctx.mkParams();
		params.add("timeout", 20);
		s.setParameters(params);

		assertions = new HashMap<Integer, HashMap<Integer, ArithExpr>>();
		strategies = new ArrayList<ArrayList<ArithExpr>>();
		payoffs = new ArrayList<ArithExpr>();
		players = new BitSet();
		strat = new ArrayList<>();

		players.set(0, numCoalitions);
		zero = ctx.mkInt(0);
		one = ctx.mkInt(1);
		for (int c = 0; c < numCoalitions; c++) {
			strategies.add(c, new ArrayList<ArithExpr>());
			payoffs.add(c, ctx.mkRealConst("v_" + c));
			for (int a = 0; a < maxActions[c]; a++) {
				strategies.get(c).add(a, ctx.mkRealConst("p_" + c + "_" + a));
				s.add(ctx.mkGe(strategies.get(c).get(a), zero));
				s.add(ctx.mkLe(strategies.get(c).get(a), one));
			}
		}
	}

	public void computeConstraints(BitSet supp)
	{
		for (int c = 0; c < indexes.size(); c++) {
			for (int a = 0; a < indexes.get(c).size(); a++) {
				if (supp.get(indexes.get(c).get(a))) {
					s.add(ctx.mkEq(assertions.get(c).get(a), payoffs.get(c)));
				} else {
					s.add(ctx.mkLe(assertions.get(c).get(a), payoffs.get(c)));
				}
			}
		}
	}

	public void computeSupport(BitSet supp, HashMap<Integer, int[]> map)
	{
		support = ctx.mkTrue();
		for (int c = 0; c < indexes.size(); c++) {
			for (int a = 0; a < indexes.get(c).size(); a++) {
				if (supp.get(indexes.get(c).get(a))) {
					support = ctx.mkAnd(support, ctx.mkGt(strategies.get(c).get(a), zero));
				} else {
					support = ctx.mkAnd(support, ctx.mkEq(strategies.get(c).get(a), zero));
				}
			}
		}
		s.add(support);
	}

	public EquilibriumResult computeEquilibria(BitSet supp, HashMap<Integer, int[]> map)
	{
		Model model;
		EquilibriumResult result = new EquilibriumResult();
		ArithExpr cnstr1;
		Distribution<Double> d;
		ArrayList<Double> eq = new ArrayList<>();
		double v;
		int a, c;

		s.push();
		computeConstraints(supp);
		computeSupport(supp, map);
		// Sum variables in the supports to one (for each player)
		for (c = 0; c < numCoalitions; c++) {
			cnstr1 = zero;
			for (a = 0; a < indexes.get(c).size(); a++) {
				cnstr1 = ctx.mkAdd(cnstr1, strategies.get(c).get(a));
			}
			s.add(ctx.mkEq(cnstr1, one));
		}
		/*
		System.out.println("-- to prove");
		for(int i = 0; i < s.getNumAssertions(); i++) {
			System.out.println(s.getAssertions()[i]);
		}
		*/
		//System.out.println("\n%% " + supp);
		Status status = s.check();
		strat.clear();
		if (status == Status.SATISFIABLE) {
			model = s.getModel();
			for (c = 0; c < numCoalitions; c++) {
				d = new Distribution<>();
				eq.add(c, getDoubleValue(model, payoffs.get(c)));
				//System.out.println("p" + p + ": " + getDoubleValue(model, payoffs.get(p)));
				for (a = 0; a < strategies.get(c).size(); a++) {
					v = getDoubleValue(model, strategies.get(c).get(a));
					//System.out.println("p" + p + " -> " + a + " " + getDoubleValue(model, strategies.get(p).get(a)));
					if (v > 0)
						d.add(a, v);
				}
				strat.add(d);
			}
		}
		s.pop();
		result.setStatus((status == Status.SATISFIABLE) ? CSGResultStatus.SAT : (status == Status.UNKNOWN) ? CSGResultStatus.UNKNOWN : CSGResultStatus.UNSAT);
		result.setPayoffVector(new ArrayList<>(eq));
		result.setStrategy(new ArrayList<>(strat));
		return result;
	}

	public void translateAssertions(HashMap<Integer, HashMap<Integer, ArrayList<Pair<BitSet, Double>>>> assertionsIdx, HashMap<Integer, int[]> map)
	{
		assertions.clear();
		ArithExpr sum, prod;
		int a, c;
		for (c = 0; c < assertionsIdx.keySet().size(); c++) {
			assertions.put(c, new HashMap<Integer, ArithExpr>());
			for (a = 0; a < assertionsIdx.get(c).keySet().size(); a++) {
				sum = zero;
				for (Pair<BitSet, Double> product : assertionsIdx.get(c).get(a)) {
					prod = getRealValue(product.second);
					for (int id = product.first.nextSetBit(0); id >= 0; id = product.first.nextSetBit(id + 1)) {
						prod = ctx.mkMul(prod, strategies.get(map.get(id)[0]).get(map.get(id)[1]));
					}
					sum = ctx.mkAdd(sum, prod);
				}
				assertions.get(c).put(a, sum);
			}
		}
	}

	public double getDoubleValue(Model model, Expr expr)
	{
		RatNum v1;
		AlgebraicNum v2;
		if (model.getConstInterp(expr) instanceof RatNum) {
			v1 = (RatNum) model.getConstInterp(expr);
			return v1.getBigIntNumerator().doubleValue() / v1.getBigIntDenominator().doubleValue();
		} else if (model.getConstInterp(expr) instanceof AlgebraicNum) {
			v2 = (AlgebraicNum) model.getConstInterp(expr);
			v1 = v2.toUpper(9);
			return v1.getBigIntNumerator().doubleValue() / v1.getBigIntDenominator().doubleValue();
		} else
			return Double.NaN;
	}

	public ArithExpr getRealValue(double v)
	{
		ArithExpr result = null;
		try {
			result = ctx.mkReal(String.valueOf(v));
		} catch (Exception e) {
			System.out.println("Error converting value " + v + " to a real.");
			System.out.println("assertions: " + assertions);
			e.printStackTrace();
		}
		return result;
	}

	public ArrayList<Distribution<Double>> getStrat()
	{
		return new ArrayList<>(this.strat);
	}

	public void setNumPlayers(int n)
	{
		this.numCoalitions = n;
	}

	public void setIndexes(ArrayList<ArrayList<Integer>> a)
	{
		this.indexes = a;
	}

	@Override
	public void setGradient
			(HashMap<Integer, HashMap<Integer, ArrayList<Pair<BitSet, Double>>>> gradient)
	{
		// TODO Auto-generated method stub

	}

	@Override
	public void setAssertions
			(HashMap<Integer, HashMap<Integer, ArrayList<Pair<BitSet, Double>>>> assertions)
	{
		// TODO Auto-generated method stub

	}

	@Override
	public void setMap(HashMap<Integer, int[]> map)
	{
		// TODO Auto-generated method stub

	}

	@Override
	public void init()
	{
		// TODO Auto-generated method stub
	}
}
