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
import java.util.List;

import prism.PrismException;

/**
 * Optimal Nash equilibrium of a bimatrix game with SCIP (floating point): the same encoding as a MILP
 * (complementarity via SOS1 constraints), solved by {@link StageGameSolverScip} in lexicographic stages with the
 * objectives of {@link CSGNashLCP} (for SF: gap, then u + v, then u). Earlier optima are kept within a relative
 * tolerance of 1e-9, so results agree with the exact SMT solvers up to about that tolerance.
 */
public class CSGNashLCPScip implements CSGNashLCP
{
	private final StageGameSolverScip scip;
	private double[] payoffs;
	private ArrayList<Distribution<Double>> strategies;

	public CSGNashLCPScip() throws PrismException
	{
		scip = new StageGameSolverScip();
		scip.setFairnessWelfareTieBreak(true);
	}

	@Override
	public String getSolverName()
	{
		return scip.getSolverName();
	}

	@Override
	public void compute(int nrows, int ncols, double[][] a, double[][] b, int crit) throws PrismException
	{
		StageGame<Double> game = new StageGame<>(new int[] { nrows, ncols });
		for (int i = 0; i < nrows; i++) {
			for (int j = 0; j < ncols; j++) {
				int joint = game.jointIndex(new int[] { i, j });
				game.setPayoff(0, joint, a[i][j]);
				game.setPayoff(1, joint, b[i][j]);
			}
		}
		StageGameSolver.Criterion criterion = crit == CSGModelCheckerEquilibria.FAIR ? StageGameSolver.Criterion.SOCIAL_FAIRNESS
				: StageGameSolver.Criterion.SOCIAL_WELFARE;
		StageGameResult<Double> res = scip.solve(game, StageGameSolver.Concept.NASH, criterion, false);
		if (res.getStatus() != StageGameResult.Status.OPTIMAL)
			throw new PrismException("SCIP could not find an optimal Nash equilibrium (" + res.getStatus() + (res.getMessage() != null ? ": " + res.getMessage() : "") + ")");
		payoffs = new double[] { res.getPayoffs().get(0), res.getPayoffs().get(1) };
		strategies = new ArrayList<>();
		for (int p = 0; p < 2; p++) {
			List<Double> st = res.getStrategies().get(p);
			Distribution<Double> d = new Distribution<>();
			for (int k = 0; k < st.size(); k++)
				if (st.get(k) > 0.0)
					d.add(k, st.get(k));
			strategies.add(d);
		}
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
