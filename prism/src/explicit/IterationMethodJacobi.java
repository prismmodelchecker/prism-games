//==============================================================================
//	
//	Copyright (c) 2016-
//	Authors:
//	* Joachim Klein <klein@tcs.inf.tu-dresden.de> (TU Dresden)
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
import explicit.rewards.MCRewards;
import explicit.rewards.MDPRewards;
import prism.PrismException;
import prism.PrismNotSupportedException;

/**
 * IterationMethod that encapsulates the functionality of the Jacobi method.
 */
class IterationMethodJacobi extends IterationMethod {

	/**
	 * Constructor.
	 * @param absolute For convergence check, perform absolute comparison?
	 * @param termCritParam For convergence check, the epsilon value to use
	 */
	public IterationMethodJacobi(boolean absolute, double epsilon)
	{
		super(absolute, epsilon);
	}

	@Override
	public IterationValIter forMvMult(DTMC<Double> dtmc)
	{
		return new TwoVectorIteration(dtmc, null) {
			@Override
			public void doIterate(IntSet states)
			{
				dtmc.mvMultJac(soln, soln2, states.iterator());
			}
		};
	}

	@Override
	public IterationIntervalIter forMvMultInterval(DTMC<Double> dtmc, boolean fromBelow, boolean enforceMonotonicity, boolean checkMonotonicity)
	{
		IterationPostProcessor post = (soln, soln2, states) -> {
			twoVectorPostProcessing(soln, soln2, states, fromBelow, enforceMonotonicity, checkMonotonicity);
		};

		return new TwoVectorIteration(dtmc, post) {
			@Override
			public void doIterate(IntSet states)
			{
				dtmc.mvMultJac(soln, soln2, states.iterator());
			}
		};
	}

	@Override
	public IterationValIter forMvMultRew(DTMC<Double> dtmc, MCRewards<Double> rew)
	{
		return new TwoVectorIteration(dtmc, null) {
			@Override
			public void doIterate(IntSet states)
			{
				dtmc.mvMultRewJac(soln, rew, soln2, states.iterator());
			}
		};
	}

	@Override
	public IterationIntervalIter forMvMultRewInterval(DTMC<Double> dtmc, MCRewards<Double> rew, boolean fromBelow, boolean enforceMonotonicity, boolean checkMonotonicity)
	{
		IterationPostProcessor post = (soln, soln2, states) -> {
			twoVectorPostProcessing(soln, soln2, states, fromBelow, enforceMonotonicity, checkMonotonicity);
		};

		return new TwoVectorIteration(dtmc, post) {
			@Override
			public void doIterate(IntSet states)
			{
				dtmc.mvMultRewJac(soln, rew, soln2, states.iterator());
			}
		};
	}

	@Override
	public IterationValIter forMvMultMinMax(MDP<Double> mdp, boolean min, int[] strat) throws PrismException
	{
		throw new PrismNotSupportedException("Jacobi not supported for MDPs");
	}

	@Override
	public IterationIntervalIter forMvMultMinMaxInterval(MDP<Double> mdp, boolean min, int[] strat, boolean fromBelow, boolean enforceMonotonicity,
			boolean checkMonotonicity) throws PrismException
	{
		throw new PrismNotSupportedException("Jacobi not supported for MDPs");
	}

	@Override
	public IterationValIter forMvMultRewMinMax(MDP<Double> mdp, MDPRewards<Double> rewards, boolean min, int[] strat) throws PrismException
	{
		throw new PrismNotSupportedException("Jacobi not supported for MDPs");
	}

	@Override
	public IterationIntervalIter forMvMultRewMinMaxInterval(MDP<Double> mdp, MDPRewards<Double> rewards, boolean min, int[] strat, boolean fromBelow,
			boolean enforceMonotonicity, boolean checkMonotonicity) throws PrismException
	{
		throw new PrismNotSupportedException("Jacobi not supported for MDPs");
	}

	@Override
	public IterationValIter forMvMultMinMaxUnc(UDTMC<Double> udtmc, MinMax minMax) throws PrismException
	{
		throw new PrismNotSupportedException("Jacobi not supported for UDTMCs");
	}

	@Override
	public IterationValIter forMvMultRewMinMaxUnc(UDTMC<Double> udtmc, MCRewards<Double> mcRewards, MinMax minMax) throws PrismException
	{
		throw new PrismNotSupportedException("Jacobi not supported for UDTMCs");
	}

	@Override
	public IterationValIter forMvMultMinMaxUnc(UMDP<Double> imdp, MinMax minMax, int[] strat) throws PrismException
	{
		throw new PrismNotSupportedException("Jacobi not supported for IMDPs");
	}

	@Override
	public IterationValIter forMvMultRewMinMaxUnc(UMDP<Double> imdp, MDPRewards<Double> mdpRewards, MinMax minMax, int[] strat) throws PrismException
	{
		throw new PrismNotSupportedException("Jacobi not supported for IMDPs");
	}

	@Override
	public String getDescriptionShort()
	{
		return "Jacobi";
	}
}
