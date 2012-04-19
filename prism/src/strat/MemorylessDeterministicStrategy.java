package strat;

import explicit.Distribution;

/**
 * Implementation of the memoryless deterministic strategy
 * 
 * @author aistis
 *
 */
public class MemorylessDeterministicStrategy implements Strategy {

	private Distribution[] choices;

	public MemorylessDeterministicStrategy(int[] choices) {
		this.choices = new Distribution[choices.length];
		Distribution dist;
		for (int i = 0; i < choices.length; i++) {
			dist = new Distribution();
			dist.add(choices[i], 1);
			this.choices[i] = dist;
		}
	}

	@Override
	public void init(int state) throws InvalidStrategyStateException {
		// do nothing
	}

	@Override
	public void updateMemory(int action, int state)
			throws InvalidStrategyStateException {
		// do nothing
	}

	@Override
	public Distribution getNextMove(int state)
			throws InvalidStrategyStateException {
		
		if(choices == null || state >= choices.length || state < 0)
			throw new InvalidStrategyStateException("Strategy not defined for state " + state +".");
		
		return choices[state];
	}

	@Override
	public void reset() {
		// do nothing
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}
}
