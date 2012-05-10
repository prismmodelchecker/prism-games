package strat;

import explicit.rewards.STPGRewards;

/**
 *
 * @author aissim
 * @version
 */
public class BoundedRewardDeterministicStrategy extends StepBoundedDeterministicStrategy
{

	private STPGRewards rewards;

	public BoundedRewardDeterministicStrategy(int[][] choices, int bound, STPGRewards rewards)
	{
		super(choices, bound);
		this.rewards = rewards;
	}

	@Override
	public void init(int state) throws InvalidStrategyStateException
	{
		memory = bound - (int) rewards.getStateReward(state);
	}

	@Override
	public void updateMemory(int action, int state) throws InvalidStrategyStateException
	{
		memory -= rewards.getStateReward(state);
		if (memory < 0)
			memory = 0;
	}

}
