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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map.Entry;

import prism.PrismException;
import soplex.SoPlex;

/**
 * Correlated equilibria (social welfare or social fairness) of a normal-form game as a sequence of LPs in SoPlex,
 * an alternative to {@link CSGCorrelatedZ3} with the same objectives and priorities:
 * <ul>
 * <li> SW: maximise the sum of payoffs, then the payoff of coalition 0, 1, ..., n-1;
 * <li> SF: minimise the gap t - s with t >= p_c >= s for all c (epigraph form, an LP), then the payoff of coalition 0, 1, ..., n-1.
 * </ul>
 * SoPlex has a single objective, so the priorities are solved in sequence: after each solve the objective just optimised
 * is kept within a small tolerance of its optimum by an extra row, the next objective is set and the LP is re-solved
 * from the current basis (warm start).
 */
public class CSGCorrelatedSoPlex implements CSGCorrelated
{
	/** Relative tolerance when fixing an objective at its optimum before optimising the next one */
	public static final double FIX_TOLERANCE = 1e-9;
	/** Primal/dual feasibility tolerances passed to SoPlex (its default is 1e-6) */
	public static final double LP_TOLERANCE = 1e-9;

	private final SoPlex lp;
	private final int nEntries;
	private final int nCoalitions;
	private String name;
	/** Coalitions in the objectives (null: all) */
	private BitSet active = null;

	@Override
	public void setActiveCoalitions(BitSet active)
	{
		this.active = active;
	}

	private boolean isActive(int c)
	{
		return active == null || active.get(c);
	}

	/**
	 * Creates a new CSGCorrelatedSoPlex
	 * @param n_entries Maximum number of entries in the utility table (joint actions)
	 * @param n_coalitions Number of coalitions
	 */
	public CSGCorrelatedSoPlex(int n_entries, int n_coalitions) throws PrismException
	{
		try {
			lp = new SoPlex();
		} catch (UnsatisfiedLinkError | IllegalStateException e) {
			throw new PrismException("Could not load SoPlex (native library soplexj): " + e.getMessage());
		}
		lp.setRealParam(SoPlex.FEASTOL, LP_TOLERANCE);
		lp.setRealParam(SoPlex.OPTTOL, LP_TOLERANCE);
		nEntries = n_entries;
		nCoalitions = n_coalitions;
		name = "SoPlex";
	}

	/** Set the SoPlex scaling method (one of the SoPlex.SCALER_* constants). */
	public void setScaler(int scaler)
	{
		lp.setScaler(scaler);
	}

	@Override
	public EquilibriumResult computeEquilibrium(HashMap<BitSet, ArrayList<Double>> utilities, ArrayList<ArrayList<HashMap<BitSet, Double>>> ce_constraints,
			ArrayList<ArrayList<Integer>> strategies, HashMap<BitSet, Integer> ce_var_map, int type)
	{
		EquilibriumResult result = new EquilibriumResult();
		boolean fair = type == CSGModelCheckerEquilibria.FAIR;
		int n = nCoalitions;

		// Columns: one probability per joint action (index given by ce_var_map); for SF also t (highest) and s (lowest payoff)
		int nvars = nEntries;
		for (BitSet e : utilities.keySet())
			nvars = Math.max(nvars, ce_var_map.get(e) + 1);
		int ncols = fair ? nvars + 2 : nvars;
		int tcol = nvars, scol = nvars + 1;

		// Utility of each coalition for each probability variable
		double[][] util = new double[n][nvars];
		boolean[] used = new boolean[nvars];
		for (Entry<BitSet, ArrayList<Double>> e : utilities.entrySet()) {
			int k = ce_var_map.get(e.getKey());
			used[k] = true;
			for (int c = 0; c < n; c++)
				util[c][k] = e.getValue().get(c);
		}

		lp.clear();
		double[] obj = new double[ncols];
		double[] lb = new double[ncols];
		double[] ub = new double[ncols];
		for (int k = 0; k < nvars; k++)
			ub[k] = used[k] ? 1.0 : 0.0;
		if (fair) {
			lb[tcol] = lb[scol] = -SoPlex.INFINITY;
			ub[tcol] = ub[scol] = SoPlex.INFINITY;
		}
		lp.addCols(obj, lb, ub);

		// sum of probabilities = 1
		int[] allIdx = new int[nvars];
		double[] ones = new double[nvars];
		for (int k = 0; k < nvars; k++) {
			allIdx[k] = k;
			ones[k] = 1.0;
		}
		lp.addRow(allIdx, ones, 1.0, 1.0);

		// Correlated equilibrium constraints: for each coalition c, recommended action q and deviation r != q,
		//   sum over the others' joint actions e of P(q, e) * (u_c(q, e) - u_c(r, e)) >= 0
		BitSet is = new BitSet(), js = new BitSet();
		for (int c = 0; c < n; c++) {
			ArrayList<Integer> acts = strategies.get(c);
			for (int q = 0; q < acts.size(); q++) {
				HashMap<BitSet, Double> others = ce_constraints.get(c).get(q);
				int m = others.size();
				int[] idx = new int[m];
				double[] base = new double[m];
				BitSet[] profiles = new BitSet[m];
				int j = 0;
				for (BitSet e : others.keySet()) {
					is.clear();
					is.or(e);
					is.set(acts.get(q));
					idx[j] = ce_var_map.get(is);
					base[j] = utilities.get(is).get(c);
					profiles[j++] = e;
				}
				for (int r = 0; r < acts.size(); r++) {
					if (r == q)
						continue;
					double[] vals = new double[m];
					for (j = 0; j < m; j++) {
						js.clear();
						js.or(profiles[j]);
						js.set(acts.get(r));
						vals[j] = base[j] - utilities.get(js).get(c);
					}
					lp.addRow(idx, vals, 0.0, SoPlex.INFINITY);
				}
			}
		}

		// Epigraph rows for SF: t - p_c >= 0 and p_c - s >= 0
		if (fair) {
			int[] idx = new int[nvars + 1];
			double[] vals = new double[nvars + 1];
			System.arraycopy(allIdx, 0, idx, 0, nvars);
			for (int c = 0; c < n; c++) {
				if (!isActive(c))
					continue; // not part of the objectives (cooperates)
				for (int k = 0; k < nvars; k++)
					vals[k] = -util[c][k];
				idx[nvars] = tcol;
				vals[nvars] = 1.0;
				lp.addRow(idx, vals.clone(), 0.0, SoPlex.INFINITY);
				for (int k = 0; k < nvars; k++)
					vals[k] = util[c][k];
				idx[nvars] = scol;
				vals[nvars] = -1.0;
				lp.addRow(idx, vals.clone(), 0.0, SoPlex.INFINITY);
			}
		}

		// Objectives in priority order: primary (sum or gap, -1), then the payoff of each active coalition
		List<Integer> objectives = new ArrayList<>();
		objectives.add(-1);
		for (int c = 0; c < n; c++)
			if (isActive(c))
				objectives.add(c);
		int nobj = objectives.size();
		for (int o = 0; o < nobj; o++) {
			int oc = objectives.get(o);
			boolean maximise = !(fair && oc == -1);
			Arrays.fill(obj, 0.0);
			if (oc == -1 && fair) {
				obj[tcol] = 1.0;
				obj[scol] = -1.0;
			} else if (oc == -1) {
				for (int c = 0; c < n; c++)
					if (isActive(c))
						for (int k = 0; k < nvars; k++)
							obj[k] += util[c][k];
			} else {
				System.arraycopy(util[oc], 0, obj, 0, nvars);
			}
			lp.setMaximise(maximise);
			lp.changeObj(obj);
			if (lp.optimize() != SoPlex.OPTIMAL) {
				result.setStatus(CSGModelCheckerEquilibria.CSGResultStatus.UNSAT);
				return result;
			}
			if (o < nobj - 1) {
				// keep this objective (within tolerance) at its optimum for the lower priorities
				double opt = lp.getObjValue();
				double tol = FIX_TOLERANCE * Math.max(1.0, Math.abs(opt));
				int[] idx = new int[ncols];
				int nz = 0;
				for (int k = 0; k < ncols; k++)
					if (obj[k] != 0.0)
						idx[nz++] = k;
				int[] fidx = Arrays.copyOf(idx, nz);
				double[] fvals = new double[nz];
				for (int k = 0; k < nz; k++)
					fvals[k] = obj[fidx[k]];
				if (maximise)
					lp.addRow(fidx, fvals, opt - tol, SoPlex.INFINITY);
				else
					lp.addRow(fidx, fvals, -SoPlex.INFINITY, opt + tol);
			}
		}

		// Joint distribution and payoffs
		double[] x = lp.getPrimal();
		Distribution<Double> d = new Distribution<>();
		ArrayList<Double> payoffs = new ArrayList<>();
		// remove numerical noise: clip to [0,1] and renormalise
		double sum = 0.0;
		for (int k = 0; k < nvars; k++) {
			x[k] = Math.min(1.0, Math.max(0.0, x[k]));
			sum += x[k];
		}
		for (int k = 0; k < nvars; k++) {
			x[k] /= sum;
			d.add(k, x[k]);
		}
		for (int c = 0; c < n; c++) {
			double p = 0.0;
			for (int k = 0; k < nvars; k++)
				p += util[c][k] * x[k];
			payoffs.add(p);
		}
		ArrayList<Distribution<Double>> strategy = new ArrayList<>();
		strategy.add(d);
		result.setStatus(CSGModelCheckerEquilibria.CSGResultStatus.SAT);
		result.setPayoffVector(payoffs);
		result.setStrategy(strategy);
		return result;
	}

	@Override
	public String getSolverName()
	{
		return name;
	}

	@Override
	public void clear()
	{
		lp.clear();
	}

	@Override
	public void printModel()
	{
	}
}
