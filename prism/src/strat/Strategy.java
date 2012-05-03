package strat;

import explicit.Distribution;
import explicit.Model;

/**
 * 
 * Generic interface to use the strategies based on stochastic strategy automata
 * [?].
 * 
 * The strategy should be used interactively. The typical usage pattern would
 * be: 1) In the initial state of the model call init method. 2) Call nextMove
 * method with the initial state as parameter. 3) Choose the successor according
 * to the distribution. 4) Once entered the successor state, call updateMemory
 * method with the action taken in the previous state and the successor state.
 * 5) Repeat steps 2-4 forever.
 * 
 * Strategy could also used to obtain the choice after a (strategy-compatible)
 * history: 1) Call reset method to reset the strategy. 2) Follow the previous
 * instructions 1-5 except in step 3, choose the successor based on the given
 * history.
 * 
 * @author aistis
 * 
 */
public interface Strategy
{

	/**
	 * Initialises memory based on a state
	 * 
	 * @param state
	 *            initial state
	 * @throws InvalidStrategyStateException
	 *             if the initial distribution function is undefined for the
	 *             given state
	 */
	public void init(int state) throws InvalidStrategyStateException;

	/**
	 * Updates memory
	 * 
	 * @param action
	 *            action taken in the previous states
	 * @param state
	 *            the current state
	 * @throws InvalidStrategyStateException
	 *             if memory update function is not defined for the given
	 *             action, state and the current strategy's memory state.
	 */
	public void updateMemory(int action, int state) throws InvalidStrategyStateException;

	/**
	 * Next move function
	 * 
	 * @param state
	 *            current state
	 * @return the distribution on actions prescribed by the strategy in a
	 *         state.
	 * @throws InvalidStrategyStateException
	 *             if next move function is undefined for the given state in
	 *             current strategy's memory state.
	 */
	public Distribution getNextMove(int state) throws InvalidStrategyStateException;

	/**
	 * Resets the strategy to uninitialised state
	 */
	public void reset();

	/**
	 * Exports adversary to a given file
	 * 
	 * @param file
	 *            file name to which adversary will be exported
	 */
	public void exportToFile(String file);

	/**
	 * Builds the product of the model and the strategy..
	 * 
	 * @param model
	 *            The model for which the strategy is defined.
	 * 
	 */
	public Model buildProduct(Model model);

	/**
	 * Get textual description of the strategy
	 *
	 * @return the textual description of the strategy
	 */
	public String getInfo();

	/**
	 * Set textual description of the strategy
	 *
	 * @param info strategy information
	 */
	public void setInfo(String info);

	/**
	 * Returns the size of memory of the strategy.
	 * @return size of memory
	 */
	public int getMemorySize();

	/**
	 * Returns strategy type
	 *
	 * @return type of the strategy
	 */
	public String getType();

	/**
	 * Returns the current memory element that fully describes state of the strategy
	 * @return
	 */
	public Object getCurrentMemoryElement();

	/**
	 * Updates the strategy's state to the one provided
	 * 
	 * @param memory memory element representing the state of the strategy
	 * @throws InvalidStrategyStateException if the memory element is not recognised by the strategy
	 */
	public void setMemory(Object memory) throws InvalidStrategyStateException;

	/**
	 * Returns the textual description of the current state of the strategy
	 * (ideally, human readable) 
	 *
	 * @return textual description of the current state of the strategy
	 */
	public String getStateDescription();

}
