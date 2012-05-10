package strat;


/**
 *
 * @author aissim
 * @version
 */
public class BoundedRewardDeterministicStrategy extends StepBoundedDeterministicStrategy
{

	public BoundedRewardDeterministicStrategy(int[][] choices, int bound)
	{
		super(choices, bound);
	}

	/**
	 *
	 * @param action
	 * @param state
	 * @throws InvalidStrategyStateException
	 */
	@Override
	public void updateMemory(int action, int state) throws InvalidStrategyStateException
	{
		// TODO Auto-generated method stub

	}

}
