//==============================================================================
//
//	Copyright (c) 2026-
//	Authors:
//	* Dave Parker <david.parker@cs.ox.ac.uk> (University of Oxford)
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

import common.IntSet;
import explicit.rewards.STPGRewards;
import prism.PrismException;
import prism.PrismNotSupportedException;

/**
 * Value iteration objects for turn-based stochastic games (STPGs, and hence SMGs),
 * for use with the {@link IterationMethod} framework,
 * e.g. {@link IterationMethod#doValueIteration}.
 * <p>
 * This is kept separate from the {@link IterationMethod} class hierarchy
 * (which is shared with PRISM) to simplify merging.
 * The supported iteration methods are {@link IterationMethodPower}
 * and (forwards) {@link IterationMethodGS}.
 */
public class IterationMethodGames
{
	/**
	 * Create a value iteration object for matrix-vector multiplication followed by min/max,
	 * i.e., for reachability probabilities.
	 * @param method The iteration method (power or Gauss-Seidel)
	 * @param stpg The STPG
	 * @param min1 Min or max for player 1 (true=min, false=max)
	 * @param min2 Min or max for player 2 (true=min, false=max)
	 * @param strat Storage for (memoryless) strategy choice indices (ignored if null)
	 */
	public static IterationMethod.IterationValIter forMvMultMinMax(IterationMethod method, STPG<Double> stpg, boolean min1, boolean min2, int[] strat) throws PrismException
	{
		if (method instanceof IterationMethodPower) {
			return method.new TwoVectorIteration(stpg, null)
			{
				@Override
				public void doIterate(IntSet states)
				{
					stpg.mvMultMinMax(soln, min1, min2, soln2, states.iterator(), strat);
				}
			};
		} else if (method instanceof IterationMethodGS) {
			return method.new SingleVectorIterationValIter(stpg)
			{
				@Override
				public boolean iterateAndCheckConvergence(IntSet states)
				{
					error = stpg.mvMultGSMinMax(soln, min1, min2, states.iterator(), method.absolute, strat);
					return error < method.termCritParam;
				}
			};
		}
		throw new PrismNotSupportedException("Iteration method \"" + method.getDescriptionShort() + "\" is not supported for games");
	}

	/**
	 * Create a value iteration object for (discounted) matrix-vector multiplication and sum of rewards
	 * followed by min/max, i.e., for reachability/total rewards.
	 * @param method The iteration method (power or Gauss-Seidel)
	 * @param stpg The STPG
	 * @param rewards The rewards
	 * @param min1 Min or max for player 1 (true=min, false=max)
	 * @param min2 Min or max for player 2 (true=min, false=max)
	 * @param strat Storage for (memoryless) strategy choice indices (ignored if null)
	 * @param disc Discount factor (1.0 = no discounting)
	 */
	public static IterationMethod.IterationValIter forMvMultRewMinMax(IterationMethod method, STPG<Double> stpg, STPGRewards<Double> rewards, boolean min1, boolean min2, int[] strat, double disc) throws PrismException
	{
		if (method instanceof IterationMethodPower) {
			return method.new TwoVectorIteration(stpg, null)
			{
				@Override
				public void doIterate(IntSet states)
				{
					stpg.mvMultRewMinMax(soln, rewards, min1, min2, soln2, states.iterator(), strat, disc);
				}
			};
		} else if (method instanceof IterationMethodGS) {
			return method.new SingleVectorIterationValIter(stpg)
			{
				@Override
				public boolean iterateAndCheckConvergence(IntSet states)
				{
					error = stpg.mvMultRewGSMinMax(soln, rewards, min1, min2, states.iterator(), method.absolute, strat, disc);
					return error < method.termCritParam;
				}
			};
		}
		throw new PrismNotSupportedException("Iteration method \"" + method.getDescriptionShort() + "\" is not supported for games");
	}
}
