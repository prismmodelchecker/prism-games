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
import java.util.ArrayList;

import com.sri.yices.BigRational;
import com.sri.yices.Config;
import com.sri.yices.Context;
import com.sri.yices.Model;
import com.sri.yices.Status;
import com.sri.yices.Terms;
import com.sri.yices.Yices;
import com.sri.yices.YicesException;

import prism.PrismException;

/**
 * Optimal Nash equilibrium of a bimatrix game with Yices (which has no optimiser), over the LCP encoding
 * (see {@link CSGNashLCP}):
 * <ol>
 * <li> Degeneracy check (one query): are there two equilibria with the same support pair but different values (u, v)?
 * <li> If not (UNSAT), the set of equilibrium values is finite, so a bound search is exact and terminates: for each
 *      objective in priority order, find an equilibrium, then repeatedly ask for a strictly better value until UNSAT,
 *      keeping the objective at its optimum for the next ones.
 * <li> Otherwise (SAT, or the check did not finish in time), the bound search could creep along a continuum of
 *      equilibria, so the game is solved with Z3's optimiser ({@link CSGNashLCPZ3}) instead.
 * </ol>
 */
public class CSGNashLCPYices implements CSGNashLCP
{
	/** Time limit (seconds) for the degeneracy check; if exceeded, Z3's optimiser is used */
	public static final int DEGENERACY_CHECK_TIMEOUT = 10;

	private final String name;
	private CSGNashLCPZ3 fallback = null;
	private double[] payoffs;
	private ArrayList<Distribution<Double>> strategies;

	/** Number of games solved by the bound search / by the Z3 fallback (degenerate or check timed out) */
	private long nBound = 0, nFallback = 0;

	public CSGNashLCPYices() throws PrismException
	{
		try {
			name = "Yices " + Yices.version();
		} catch (UnsatisfiedLinkError e) {
			throw new PrismException("Could not initialise Yices: " + e.getMessage());
		}
	}

	@Override
	public String getSolverName()
	{
		return name;
	}

	public long getNumBoundSearch()
	{
		return nBound;
	}

	public long getNumFallback()
	{
		return nFallback;
	}

	/** The LCP encoding for one copy of the variables */
	private static class Encoding
	{
		final int[] x, y;
		final int u, v;
		final int[] formulas;

		Encoding(int nrows, int ncols, int[][] a, int[][] b)
		{
			int real = Yices.realType();
			x = new int[nrows];
			y = new int[ncols];
			for (int i = 0; i < nrows; i++)
				x[i] = Terms.newUninterpretedTerm(real);
			for (int j = 0; j < ncols; j++)
				y[j] = Terms.newUninterpretedTerm(real);
			u = Terms.newUninterpretedTerm(real);
			v = Terms.newUninterpretedTerm(real);
			int one = Terms.intConst(1);
			ArrayList<Integer> f = new ArrayList<>();
			for (int t : x)
				f.add(Terms.arithGeq0(t));
			for (int t : y)
				f.add(Terms.arithGeq0(t));
			f.add(Terms.arithEq(Terms.add(x), one));
			f.add(Terms.arithEq(Terms.add(y), one));
			for (int i = 0; i < nrows; i++) {
				int[] terms = new int[ncols];
				for (int j = 0; j < ncols; j++)
					terms[j] = Terms.mul(a[i][j], y[j]);
				int ay = Terms.add(terms);
				f.add(Terms.arithGeq(u, ay));
				f.add(Terms.or(Terms.arithEq0(x[i]), Terms.arithEq(u, ay)));
			}
			for (int j = 0; j < ncols; j++) {
				int[] terms = new int[nrows];
				for (int i = 0; i < nrows; i++)
					terms[i] = Terms.mul(b[i][j], x[i]);
				int xb = Terms.add(terms);
				f.add(Terms.arithGeq(v, xb));
				f.add(Terms.or(Terms.arithEq0(y[j]), Terms.arithEq(v, xb)));
			}
			formulas = f.stream().mapToInt(Integer::intValue).toArray();
		}
	}

	@Override
	public void compute(int nrows, int ncols, double[][] a, double[][] b, int crit) throws PrismException
	{
		boolean exact;
		try {
			int[][] ta = constants(a), tb = constants(b);
			exact = !degenerate(nrows, ncols, ta, tb);
			if (exact)
				boundSearch(nrows, ncols, ta, tb, crit);
		} catch (YicesException e) {
			throw new PrismException("Yices error: " + e.getMessage());
		} finally {
			// frees the terms created for this game (all contexts using them are closed)
			Yices.yicesGarbageCollect();
		}
		if (exact) {
			nBound++;
		} else {
			nFallback++;
			if (fallback == null)
				fallback = new CSGNashLCPZ3();
			fallback.compute(nrows, ncols, a, b, crit);
			payoffs = fallback.getPayoffs();
			strategies = fallback.getStrategies();
		}
	}

	private static int[][] constants(double[][] m)
	{
		int[][] t = new int[m.length][];
		for (int i = 0; i < m.length; i++) {
			t[i] = new int[m[i].length];
			for (int j = 0; j < m[i].length; j++)
				t[i][j] = Terms.rationalConst(new BigRational(new BigDecimal(CSGNashLCP.decimal(m[i][j]))));
		}
		return t;
	}

	private static Context newContext()
	{
		Config cfg = new Config("QF_LRA");
		cfg.set("mode", "push-pop");
		Context ctx = new Context(cfg);
		cfg.close();
		return ctx;
	}

	/**
	 * Two copies of the encoding with the same support pair (x_i = 0 iff x'_i = 0, same for y) and different values.
	 * @return false if UNSAT (values are constant on each support pair), true if SAT or unknown
	 */
	private boolean degenerate(int nrows, int ncols, int[][] a, int[][] b)
	{
		Encoding p = new Encoding(nrows, ncols, a, b);
		Encoding q = new Encoding(nrows, ncols, a, b);
		try (Context ctx = newContext()) {
			ctx.assertFormulas(p.formulas);
			ctx.assertFormulas(q.formulas);
			for (int i = 0; i < nrows; i++)
				ctx.assertFormula(Terms.iff(Terms.arithEq0(p.x[i]), Terms.arithEq0(q.x[i])));
			for (int j = 0; j < ncols; j++)
				ctx.assertFormula(Terms.iff(Terms.arithEq0(p.y[j]), Terms.arithEq0(q.y[j])));
			ctx.assertFormula(Terms.or(Terms.arithNeq(p.u, q.u), Terms.arithNeq(p.v, q.v)));
			return ctx.check(DEGENERACY_CHECK_TIMEOUT) != Status.UNSAT;
		}
	}

	/** Exact when the game is not degenerate (the set of equilibrium values is finite). */
	private void boundSearch(int nrows, int ncols, int[][] a, int[][] b, int crit) throws PrismException
	{
		Encoding e = new Encoding(nrows, ncols, a, b);
		// objectives in priority order: SF (the gap |u - v|, minimised), then u + v and u (maximised).
		// Each objective is a function of the equilibrium values only, so every strict improvement moves to
		// different values, of which there are finitely many (the game is not degenerate).
		int diff = Terms.sub(e.u, e.v);
		int sum = Terms.add(e.u, e.v);
		boolean fair = crit == CSGModelCheckerEquilibria.FAIR;

		try (Context ctx = newContext()) {
			ctx.assertFormulas(e.formulas);
			if (ctx.check() != Status.SAT)
				throw new PrismException("Yices could not find a Nash equilibrium");
			BigRational[] sol = read(ctx, e);
			for (int k = fair ? 0 : 1; k < 3; k++) {
				while (true) {
					int n = nrows + ncols;
					int keep, better;
					if (k == 0) {
						// gap: |u - v| <= c, then < c
						int c = gap(sol[n], sol[n + 1]);
						keep = Terms.and(Terms.arithLeq(diff, c), Terms.arithGeq(diff, Terms.neg(c)));
						better = Terms.and(Terms.arithLt(diff, c), Terms.arithGt(diff, Terms.neg(c)));
					} else {
						int f = k == 1 ? sum : e.u;
						int c = k == 1 ? Terms.add(Terms.rationalConst(sol[n]), Terms.rationalConst(sol[n + 1])) : Terms.rationalConst(sol[n]);
						keep = Terms.arithGeq(f, c);
						better = Terms.arithGt(f, c);
					}
					// the optimum is at least as good as the incumbent (kept for the lower priorities)
					ctx.assertFormula(keep);
					ctx.push();
					ctx.assertFormula(better);
					Status st = ctx.check();
					if (st == Status.UNSAT) {
						ctx.pop();
						break;
					}
					if (st != Status.SAT) {
						ctx.pop();
						throw new PrismException("Yices returned " + st + " in the bound search");
					}
					sol = read(ctx, e);
					ctx.pop();
				}
			}
			payoffs = new double[] { sol[nrows + ncols].doubleValue(), sol[nrows + ncols + 1].doubleValue() };
			strategies = new ArrayList<>();
			Distribution<Double> dx = new Distribution<>(), dy = new Distribution<>();
			for (int i = 0; i < nrows; i++)
				if (sol[i].getNumerator().signum() > 0)
					dx.add(i, sol[i].doubleValue());
			for (int j = 0; j < ncols; j++)
				if (sol[nrows + j].getNumerator().signum() > 0)
					dy.add(j, sol[nrows + j].doubleValue());
			strategies.add(dx);
			strategies.add(dy);
		}
	}

	/** Constant |u - v| (exact) */
	private static int gap(BigRational u, BigRational v)
	{
		int d = Terms.sub(Terms.rationalConst(u), Terms.rationalConst(v));
		return u.getNumerator().multiply(v.getDenominator()).compareTo(v.getNumerator().multiply(u.getDenominator())) >= 0 ? d : Terms.neg(d);
	}

	/** Current model values of x, y, u, v */
	private static BigRational[] read(Context ctx, Encoding e)
	{
		int n = e.x.length + e.y.length;
		BigRational[] r = new BigRational[n + 2];
		try (Model m = ctx.getModel()) {
			for (int i = 0; i < e.x.length; i++)
				r[i] = m.bigRationalValue(e.x[i]);
			for (int j = 0; j < e.y.length; j++)
				r[e.x.length + j] = m.bigRationalValue(e.y[j]);
			r[n] = m.bigRationalValue(e.u);
			r[n + 1] = m.bigRationalValue(e.v);
		}
		return r;
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
