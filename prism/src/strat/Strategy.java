package strat;

import explicit.Distribution;

/**
 * 
 * Generic interface to use the strategies based on stochastic strategy automata [?].
 * 
 * The strategy should be used interactively. The typical usage pattern would be:
 * 1) In the initial state of the model call init method.
 * 2) Call nextMove method with the initial state as parameter.
 * 3) Choose the successor according to the distribution.
 * 4) Once entered the successor state, call updateMemory method with the action taken
 * in the previous state and the successor state.
 * 5) Repeat steps 2-4 forever.
 * 
 * Strategy could also used to obtain the choice after a (strategy-compatible) history:
 * 1) Call reset method to reset the strategy.
 * 2) Follow the previous instructions 1-5 except in step 3, choose the successor
 * based on the given history. 
 * 
 * @author aistis
 *
 */
public interface Strategy {

	/**
	 * Initialises memory based on a state
	 * 
	 * @param state initial state
	 * @throws InvalidStrategyStateException if the initial distribution
	 * function is undefined for the given state
	 */
	public void init(int state) throws InvalidStrategyStateException;
	
	/**
	 * Updates memory
	 * 
	 * @param action action taken in the previous states
	 * @param state the current state
	 * @throws InvalidStrategyStateException if memory update function is not defined
	 * for the given action, state and the current strategy's memory state.
	 */
	public void updateMemory(int action, int state)  throws InvalidStrategyStateException;

	/**
	 * Next move function
	 * 
	 * @param state current state
	 * @return the distribution on actions prescribed by the strategy in a state.
	 * @throws InvalidStrategyStateException if next move function is undefined
	 * for the given state in current strategy's memory state.
	 */
	public Distribution getNextMove(int state) throws InvalidStrategyStateException;
	
	/**
	 * Resets the strategy to uninitialised state
	 */
	public void reset();
}
