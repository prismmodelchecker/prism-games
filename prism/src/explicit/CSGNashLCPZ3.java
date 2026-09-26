//==============================================================================
//
//	Copyright (c) 2026-
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

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.ArrayList;
import java.util.HashMap;

import com.microsoft.z3.ArithExpr;
import com.microsoft.z3.Context;
import com.microsoft.z3.Expr;
import com.microsoft.z3.Model;
import com.microsoft.z3.Optimize;
import com.microsoft.z3.Params;
import com.microsoft.z3.RatNum;
import com.microsoft.z3.RealExpr;
import com.microsoft.z3.Status;
import com.microsoft.z3.Version;

import prism.PrismException;

/**
 * Optimal Nash equilibrium of a bimatrix game with Z3's optimiser: one lexicographic optimisation call
 * over the LCP encoding (see {@link CSGNashLCP}).
 */
public class CSGNashLCPZ3 implements CSGNashLCP
{
	private final Context ctx;
	private final String name;
	private double[] payoffs;
	private ArrayList<Distribution<Double>> strategies;

	public CSGNashLCPZ3() throws PrismException
	{
		try {
			HashMap<String, String> cfg = new HashMap<>();
			cfg.put("model", "true");
			ctx = new Context(cfg);
			name = Version.getFullVersion();
		} catch (UnsatisfiedLinkError e) {
			throw new PrismException("Could not initialise Z3: " + e.getMessage());
		}
	}

	@Override
	public String getSolverName()
	{
		return name;
	}

	@Override
	public void compute(int nrows, int ncols, double[][] a, double[][] b, int crit) throws PrismException
	{
		if (opt == null) {
			opt = ctx.mkOptimize();
			Params params = ctx.mkParams();
			params.add("priority", "lex");
			opt.setParameters(params);
		}
		Optimize o = opt;
		o.Push();
		try {
			solve(o, nrows, ncols, a, b, crit);
		} finally {
			o.Pop();
		}
	}

	/** One optimiser reused for all games (a scope is pushed per game) */
	private Optimize opt = null;

	private void solve(Optimize o, int nrows, int ncols, double[][] a, double[][] b, int crit) throws PrismException
	{

		RealExpr[] x = new RealExpr[nrows];
		RealExpr[] y = new RealExpr[ncols];
		for (int i = 0; i < nrows; i++)
			x[i] = ctx.mkRealConst("x" + i);
		for (int j = 0; j < ncols; j++)
			y[j] = ctx.mkRealConst("y" + j);
		RealExpr u = ctx.mkRealConst("u");
		RealExpr v = ctx.mkRealConst("v");
		ArithExpr zero = ctx.mkReal(0);
		ArithExpr one = ctx.mkReal(1);

		for (RealExpr t : x)
			o.Add(ctx.mkGe(t, zero));
		for (RealExpr t : y)
			o.Add(ctx.mkGe(t, zero));
		o.Add(ctx.mkEq(ctx.mkAdd(x), one));
		o.Add(ctx.mkEq(ctx.mkAdd(y), one));
		for (int i = 0; i < nrows; i++) {
			ArithExpr[] terms = new ArithExpr[ncols];
			for (int j = 0; j < ncols; j++)
				terms[j] = ctx.mkMul(ctx.mkReal(CSGNashLCP.decimal(a[i][j])), y[j]);
			ArithExpr ay = ctx.mkAdd(terms);
			o.Add(ctx.mkGe(u, ay));
			o.Add(ctx.mkOr(ctx.mkEq(x[i], zero), ctx.mkEq(u, ay)));
		}
		for (int j = 0; j < ncols; j++) {
			ArithExpr[] terms = new ArithExpr[nrows];
			for (int i = 0; i < nrows; i++)
				terms[i] = ctx.mkMul(ctx.mkReal(CSGNashLCP.decimal(b[i][j])), x[i]);
			ArithExpr xb = ctx.mkAdd(terms);
			o.Add(ctx.mkGe(v, xb));
			o.Add(ctx.mkOr(ctx.mkEq(y[j], zero), ctx.mkEq(v, xb)));
		}

		if (crit == CSGModelCheckerEquilibria.FAIR) {
			RealExpr g = ctx.mkRealConst("g");
			o.Add(ctx.mkGe(g, ctx.mkSub(u, v)));
			o.Add(ctx.mkGe(g, ctx.mkSub(v, u)));
			o.MkMinimize(g);
		}
		o.MkMaximize(ctx.mkAdd(u, v));
		o.MkMaximize(u);

		if (o.Check() != Status.SATISFIABLE)
			throw new PrismException("Z3 could not find an optimal Nash equilibrium: " + o.getReasonUnknown());
		Model m = o.getModel();
		payoffs = new double[] { value(m, u), value(m, v) };
		strategies = new ArrayList<>();
		strategies.add(distribution(m, x));
		strategies.add(distribution(m, y));
	}

	private Distribution<Double> distribution(Model m, RealExpr[] vars)
	{
		Distribution<Double> d = new Distribution<>();
		for (int i = 0; i < vars.length; i++) {
			double p = value(m, vars[i]);
			if (p > 0.0)
				d.add(i, p);
		}
		return d;
	}

	static double value(Model m, Expr e)
	{
		RatNum r = (RatNum) m.eval(e, true);
		return new BigDecimal(r.getBigIntNumerator()).divide(new BigDecimal(r.getBigIntDenominator()), MathContext.DECIMAL64).doubleValue();
	}

	@Override
	public double[] getPayoffs()
	{
		return payoffs;
	}

	@Override
	public ArrayList<Distribution<Double>> getStrategies()
	{
		return strategies;
	}
}
