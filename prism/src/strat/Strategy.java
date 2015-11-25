package strat;

import prism.PrismException;
import prism.PrismLog;
import explicit.Distribution;
import explicit.Model;
import parser.Values;

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
	// Types of info stored for each choice
	public enum Choice {
		INDEX, ACTION, UNKNOWN, ARBITRARY, UNREACHABLE;
	};
	
	
	/**
	 * Store constants used to build model for this strategy.
	 * 
	 * @param lastConstants
	 */
	public void setConstants(Values lastConstants);
	
	/**
	 * Retrieve constants used to build model for this strategy.
	 */
	public Values getConstants();
	
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
	 * Builds the product of the model and the strategy. The product is built by
	 * adding extra integer variable to the state to represent the memory state
	 * of the strategy. The initial states are the first N states of the product
	 * where N is the size of the original model.
	 * 
	 * @param model
	 *            The model for which the strategy is defined.
	 * @throws PrismException
	 * 
	 */
	public Model buildProduct(Model model) throws PrismException;

	/**
	 * Get textual description of the strategy
	 * 
	 * @return the textual description of the strategy
	 */
	public String getInfo();

	/**
	 * Set textual description of the strategy
	 * 
	 * @param info
	 *            strategy information
	 */
	public void setInfo(String info);

	/**
	 * Returns the size of memory of the strategy.
	 * 
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
	 * Returns the current memory element that fully describes state of the
	 * strategy
	 * 
	 * @return
	 */
	public Object getCurrentMemoryElement();

	/**
	 * Updates the strategy's state to the one provided
	 * 
	 * @param memory
	 *            memory element representing the state of the strategy
	 * @throws InvalidStrategyStateException
	 *             if the memory element is not recognised by the strategy
	 */
	public void setMemory(Object memory) throws InvalidStrategyStateException;

	/**
	 * Returns a textual description giving key information about the strategy (ideally, human readable).
	 */
	public String getDescription();

	/**
	 * Returns the initial memory state of the strategy for the state. It is
	 * required by the model checker to determine which states can be treated as
	 * initial in the product.
	 * 
	 * @param s
	 *            the state for which to return initial memory element
	 * 
	 * @return non negative integer or -1 if product does not contain extra
	 *         variables
	 */
	public int getInitialStateOfTheProduct(int s);

	
	//	/**
	//	 * Retrieve the expected value that this strategy will achieve from it's
	//	 * current state
	//	 * @return the expect value of the function, return -1 if exp values are not defined
	//	 */
	//	public double getExpectedValue();
	//	
	//	/**
	//	 * Get expected value if a given action was taken and given state turned out to be a successor
	//	 * @param action action
	//	 * @param state state
	//	 * @return expectation
	//	 */
	//	public double getExpectedValue(int action, int state);

	// NEW METHODS:
	
	/**
	 * Export the strategy to a PrismLog, displaying strategy choices as action names.
	 */
	public void exportActions(PrismLog out);
	
	/**
	 * Export the strategy to a PrismLog, displaying strategy choices as indices.
	 */
	public void exportIndices(PrismLog out);
	
	/**
	 * Export the model induced by this strategy to a PrismLog.
	 */
	public void exportInducedModel(PrismLog out);
	
	/**
	 * Export the strategy to a dot file (of the model showing the strategy).
	 */
	public void exportDotFile(PrismLog out);
	
	/**
	 * Initialise the strategy, based on an initial model state.
	 * @param s Initial state of the model
	 */
	public void initialise(int s);

	/**
	 * Update the strategy, based on the next step in a model's history.
	 * @param action The action taken in the previous state of the model
	 * @param s The new state of the model
	 */
	public void update(Object action, int s);
	
	/**
	 * Get the action chosen by the strategy in the current state (assuming it is deterministic). 
	 */
	public Object getChoiceAction();
	
	/**
	 * Clear storage of the strategy.
	 */
	public void clear();
}
