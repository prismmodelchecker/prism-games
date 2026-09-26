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

import prism.PrismException;

/**
 * Optimal Nash equilibrium of a bimatrix game (A, B) via the LCP encoding with the values u, v as variables:
 * <pre>
 *   x, y >= 0,  sum x = sum y = 1,
 *   u >= (A y)_i  and  (x_i = 0  or  u = (A y)_i)   for each row i,
 *   v >= (x^T B)_j  and  (y_j = 0  or  v = (x^T B)_j)   for each column j.
 * </pre>
 * Every model is a Nash equilibrium with values (u, v), and the feasible set is closed, so the optimum over all
 * equilibria (including continua of equilibria on one support) is attained. Objectives, in priority order:
 * <ul>
 * <li> social welfare: maximise u + v, then u;
 * <li> social fairness: minimise |u - v|, then maximise u + v, then u.
 * </ul>
 * Payoffs are maximised: for costs (min), the caller passes the negated game.
 */
public interface CSGNashLCP
{
	public String getSolverName();

	/**
	 * Computes the optimal equilibrium of the game with payoff matrices a (row player) and b (column player).
	 * @param crit CSGModelCheckerEquilibria.SWEQ or CSGModelCheckerEquilibria.FAIR
	 */
	public void compute(int nrows, int ncols, double[][] a, double[][] b, int crit) throws PrismException;

	/** Values (u, v) of the last equilibrium computed */
	public double[] getPayoffs();

	/** Strategies (row player, column player) of the last equilibrium computed, over row/column indices */
	public ArrayList<Distribution<Double>> getStrategies();

	/**
	 * Exact decimal string for a matrix entry: the shortest decimal that rounds to the double (as Double.toString),
	 * without exponent notation, so that solvers read e.g. 0.1 as 1/10 and 1.0E-5 correctly.
	 */
	public static String decimal(double d)
	{
		return new BigDecimal(Double.toString(d)).toPlainString();
	}
}
