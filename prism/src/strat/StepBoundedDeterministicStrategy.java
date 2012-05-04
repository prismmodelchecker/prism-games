package strat;

import explicit.Distribution;

public class StepBoundedDeterministicStrategy extends FiniteMemoryStrategy {

	// memory: the number of steps currently made
	private int memory;

	// bound: maximum number of steps to be made
	private int bound;

	// choices of a strategy
	private int[][] choices;
	private int chSize;

	/**
	 * Initialises the strategy
	 * 
	 * @param choices
	 *            strategy choice table: for every state contains an array of
	 *            integers structured as follows {bk, ck, bk-1, ck-1,..., b0 c0}
	 *            where ck represents the choice to be made by the strategy when
	 *            B-k steps have elapsed.
	 */
	public StepBoundedDeterministicStrategy(int[][] choices, int bound) {
		this.choices = choices;

		if (bound < 0)
			throw new IllegalArgumentException("The bound should be positive.");

		this.bound = bound;

		// computing the size of the choice function and validating the format
		chSize = 0;
		int prev = -1;
		for (int i = 0; i < choices.length; i++)
			for (int j = 0; j < choices[i].length; j++) {
				chSize++;

				// performing validation
				if (choices[i][j] < 0)
					throw new IllegalArgumentException(
							"The format of choices is invalid: array cannot contain negative numbers.");

				// adjusting the choices to be at most the bound
				if (j % 2 == 0 && choices[i][j] > bound) {
					choices[i][j] = bound;
					prev = bound;
				} else if (j == 0 && choices[i][j] < bound) {
					throw new IllegalArgumentException(
							"The format of choices is invalid: the first pivot has to be >= than the bound.");
				}

				// checking if ordering is correct
				if (j % 2 == 0)
					if (choices[i][j] > prev)
						throw new IllegalArgumentException(
								"The format of choices is invalid: pivots have to be in decreasing order.");
			}
	}

	@Override
	public void init(int state) throws InvalidStrategyStateException {
		memory = bound;
	}

	@Override
	public void updateMemory(int action, int state)
			throws InvalidStrategyStateException {
		memory--;
	}

	@Override
	public Distribution getNextMove(int state)
			throws InvalidStrategyStateException {

		if (state > choices.length)
			throw new InvalidStrategyStateException(
					"The strategy undefined for state " + state + ".");

		// determining the action
		int[] actions = choices[state];
		int c = 0;
		for (int i = 0; i < actions.length; i += 2)
			if (actions[i] >= memory)
				c = actions[i + 1];
			else
				break;

		Distribution dist = new Distribution();
		dist.add(c, 1);

		return dist;
	}

	@Override
	public void reset() {
		memory = bound;
	}

	@Override
	public int getMemorySize() {
		return memory;
	}

	@Override
	public Object getCurrentMemoryElement() {
		return memory;
	}

	@Override
	public void setMemory(Object memory) throws InvalidStrategyStateException {
		if (memory instanceof Integer) {
			this.memory = (Integer) memory;
		} else {
			throw new InvalidStrategyStateException(
					"Memory has to integer for this strategy.");
		}
	}

	@Override
	public String getStateDescription() {
		String desc = "";
		desc += "Finite memory deterministic strategy\n";
		desc += "Size of memory: " + bound + "\n";
		desc += "Size of next move function: " + chSize + " \n";
		return desc;
	}

	public static void main(String[] args) throws InvalidStrategyStateException {
		int[][] choices = { { 30, 1, 28, 2 }, { 25, 1, 24, 2 } };
		int bound = 25;
		
		StepBoundedDeterministicStrategy strat = new StepBoundedDeterministicStrategy(
				choices, bound);
		strat.init(0);

		for (int i = 0; i < 25; i++) {
			System.out.println("i = " + i);
			System.out.println(strat.getNextMove(0) + ", "
					+ strat.getNextMove(1));
			strat.updateMemory(0, 0);
		}
	}

}
