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

import java.util.List;

/**
 * Result of solving a {@link StageGame}: a single selected equilibrium (or game value).
 *
 * @param <V> value type
 */
public class StageGameResult<V>
{
	public enum Status
	{
		/** Optimal with respect to the criterion and all tie-breaking objectives */
		OPTIMAL,
		/** Time limit reached: the result is a valid equilibrium but not proven optimal (see {@link #getGap()}) */
		TIME_LIMIT,
		/** No solution (should not happen for Nash/correlated equilibria, which always exist) */
		INFEASIBLE,
		/** Solver failure */
		FAILED
	}

	private final Status status;
	private final List<V> payoffs;
	private final List<List<V>> strategies;
	private final List<V> jointDistribution;
	private final double gap;
	private final String message;

	public StageGameResult(Status status, List<V> payoffs, List<List<V>> strategies, List<V> jointDistribution, double gap, String message)
	{
		this.status = status;
		this.payoffs = payoffs;
		this.strategies = strategies;
		this.jointDistribution = jointDistribution;
		this.gap = gap;
		this.message = message;
	}

	public static <V> StageGameResult<V> failed(Status status, String message)
	{
		return new StageGameResult<>(status, null, null, null, Double.NaN, message);
	}

	public Status getStatus()
	{
		return status;
	}

	public boolean hasSolution()
	{
		return payoffs != null;
	}

	/** Payoff of each coalition under the selected equilibrium (in the original, un-negated scale). */
	public List<V> getPayoffs()
	{
		return payoffs;
	}

	/** Mixed strategy of each coalition (Nash equilibria and zero-sum games); null for correlated equilibria. */
	public List<List<V>> getStrategies()
	{
		return strategies;
	}

	/** Distribution over joint actions (correlated equilibria); null otherwise. Indexed as in {@link StageGame#jointIndex(int[])}. */
	public List<V> getJointDistribution()
	{
		return jointDistribution;
	}

	/** Relative optimality gap reported by the solver for the stage at which it stopped (0 if optimal). */
	public double getGap()
	{
		return gap;
	}

	public String getMessage()
	{
		return message;
	}

	@Override
	public String toString()
	{
		return status + " payoffs=" + payoffs + (strategies != null ? " strategies=" + strategies : "")
				+ (jointDistribution != null ? " joint=" + jointDistribution : "") + (message != null ? " (" + message + ")" : "");
	}
}
