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

/**
 * IterationMethod that encapsulates the functionality of the Power method.
 */
public class IterationMethodPower extends IterationMethod {

	/**
	 * Constructor.
	 * @param absolute For convergence check, perform absolute comparison?
	 * @param termCritParam For convergence check, the epsilon value to use
	 */
	public IterationMethodPower(boolean absolute, double epsilon)
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
				dtmc.mvMult(soln, soln2, states.iterator());
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
				dtmc.mvMult(soln, soln2, states.iterator());
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
				dtmc.mvMultRew(soln, rew, soln2, states.iterator());
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
				dtmc.mvMultRew(soln, rew, soln2, states.iterator());
			}
		};
	}

	@Override
	public IterationValIter forMvMultMinMax(MDP<Double> mdp, boolean min, int[] strat)
	{
		return new TwoVectorIteration(mdp, null) {
			@Override
			public void doIterate(IntSet states)
			{
				mdp.mvMultMinMax(soln, min, soln2, states.iterator(), strat);
			}
		};
	}

	@Override
	public IterationIntervalIter forMvMultMinMaxInterval(MDP<Double> mdp, boolean min, int[] strat, boolean fromBelow, boolean enforceMonotonicity,
			boolean checkMonotonicity) throws PrismException
	{
		IterationPostProcessor post = (soln, soln2, states) -> {
			twoVectorPostProcessing(soln, soln2, states, fromBelow, enforceMonotonicity, checkMonotonicity);
		};

		return new TwoVectorIteration(mdp, post) {
			@Override
			public void doIterate(IntSet states)
			{
				mdp.mvMultMinMax(soln, min, soln2, states.iterator(), strat);
			}
		};
	}

	@Override
	public String getDescriptionShort()
	{
		return "Power method";
	}

	@Override
	public IterationValIter forMvMultRewMinMax(MDP<Double> mdp, MDPRewards<Double> rewards, boolean min, int[] strat) throws PrismException
	{
		return new TwoVectorIteration(mdp, null) {
			@Override
			public void doIterate(IntSet states)
			{
				mdp.mvMultRewMinMax(soln, rewards, min, soln2, states.iterator(), strat);
			}
		};
	}

	@Override
	public IterationIntervalIter forMvMultRewMinMaxInterval(MDP<Double> mdp, MDPRewards<Double> rewards, boolean min, int[] strat, boolean fromBelow,
			boolean enforceMonotonicity, boolean checkMonotonicity) throws PrismException
	{
		IterationPostProcessor post = (soln, soln2, states) -> {
			twoVectorPostProcessing(soln, soln2, states, fromBelow, enforceMonotonicity, checkMonotonicity);
		};

		return new TwoVectorIteration(mdp, post) {
			@Override
			public void doIterate(IntSet states)
			{
				mdp.mvMultRewMinMax(soln, rewards, min, soln2, states.iterator(), strat);
			}
		};
	}

	@Override
	public IterationValIter forMvMultMinMaxUnc(UDTMC<Double> udtmc, MinMax minMax)
	{
		return new TwoVectorIteration(udtmc, null) {
			@Override
			public void doIterate(IntSet states)
			{
				udtmc.mvMultUnc(soln, minMax, soln2, states.iterator());
			}
		};
	}

	@Override
	public IterationValIter forMvMultRewMinMaxUnc(UDTMC<Double> udtmc, MCRewards<Double> mcRewards, MinMax minMax)
	{
		return new TwoVectorIteration(udtmc, null) {
			@Override
			public void doIterate(IntSet states)
			{
				udtmc.mvMultRewUnc(soln, mcRewards, minMax, soln2, states.iterator());
			}
		};
	}

	@Override
	public IterationValIter forMvMultMinMaxUnc(UMDP<Double> imdp, MinMax minMax, int[] strat)
	{
		return new TwoVectorIteration(imdp, null) {
			@Override
			public void doIterate(IntSet states)
			{
				imdp.mvMultUnc(soln, minMax, soln2, states.iterator(), strat);
			}
		};
	}
	
	@Override
	public IterationValIter forMvMultRewMinMaxUnc(UMDP<Double> imdp, MDPRewards<Double> mdpRewards, MinMax minMax, int[] strat)
	{
		return new TwoVectorIteration(imdp, null) {
			@Override
			public void doIterate(IntSet states)
			{
				imdp.mvMultRewUnc(soln, mdpRewards, minMax, soln2, states.iterator(), strat);
			}
		};
	}
}
