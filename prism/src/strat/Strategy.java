//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <d.a.parker@cs.bham.ac.uk> (University of Birmingham/Oxford)
//	* Aistis Simaitis <aistis.aimaitis@cs.ox.ac.uk> (University of Oxford)
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

package strat;

import java.io.File;
import java.util.HashMap;

import explicit.Distribution;
import explicit.Model;
import prism.Prism.StrategyExportType;
import prism.PrismException;
import prism.PrismLog;

/**
 * Interface for classes to store strategies (for MDPs, games, etc.)
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
 */
public interface Strategy
{
	// Types of info stored for each choice
	public enum Choice {
		INDEX, ACTION, UNKNOWN, ARBITRARY, UNREACHABLE;
	};

	// Getters for basic strategy info
	
	/**
	 * Get a textual description of the strategy.
	 */
	public String getInfo();

	/**
	 * Get the size of the memory of the strategy.
	 */
	public int getMemorySize();

	/**
	 * Get the type of the strategy, as a string.
	 */
	public String getType();

	/**
	 * Returns a textual description giving key information about the strategy (ideally, human readable).
	 */
	public String getDescription();

	// Setters for basic strategy info
	
	/**
	 * Set a textual description for the strategy.
	 * @param info String describing the strategy.
	 */
	public void setInfo(String info);

	// Methods for interactive querying of strategy
	
	/**
	 * Initialise the strategy (set its memory) based on an initial model state.
	 * @param state Initial model state
	 * @throws InvalidStrategyStateException if the initial distribution function is undefined for the given state
	 */
	public void init(int state) throws InvalidStrategyStateException;

	/**
	 * Update the strategy (update its memory) based on an action and next state from the model.
	 * @param action Action taken in the previous model state
	 * @param state The next model state
	 * @throws InvalidStrategyStateException if memory update function is not defined for the given action, state and the current strategy's memory state
	 */
	public void updateMemory(int action, int state) throws InvalidStrategyStateException;

	/**
	 * Get the move of the strategy in model state {@code state}
	 * as a distribution over indices of choices that are available in that state,
	 * i.e., which specifies the probability with which each choice should be taken.  
	 * @param state Current model state
	 * @throws InvalidStrategyStateException if next move function is undefined for the given state in current strategy's memory state
	 */
	public Distribution getNextMove(int state) throws InvalidStrategyStateException;

	/**
	 * Get the move of the strategy in model state {@code state}
	 * as a distribution over indices of actions that are available in that state,
	 * i.e., which specifies the probability with which each action should be taken.  
	 * @param state Current model state
	 * @throws InvalidStrategyStateException if next move function is undefined for the given state in current strategy's memory state
	 */
	public HashMap<String,Double> getNextAction(int state) throws InvalidStrategyStateException;
	
	/**
	 * Get the current memory status of the strategy.
	 */
	public Object getCurrentMemoryElement();

	/**
	 * Set the current memory status of the strategy to the one provided.
	 * @param memory Memory element representing the state of the strategy
	 * @throws InvalidStrategyStateException if the memory element is not recognised by the strategy
	 */
	public void setMemory(Object memory) throws InvalidStrategyStateException;

	/**
	 * Reset the strategy to an uninitialised state.
	 */
	public void reset();

	// Product methods
	
	/**
	 * Build the product the strategy and a model. The product is built by
	 * adding an extra integer variable to the state to represent the memory state
	 * of the strategy. The initial states are the first N states of the product
	 * where N is the size of the original model.
	 * @param model The model for which the strategy is defined.
	 */
	public Model buildProduct(Model model) throws PrismException;

	/**
	 * Get the initial memory state of the strategy for a state.
	 * It is required by the model checker to determine which states can be treated as initial in the product.
	 * Returns -1 if the product does not contain extra variables  
	 * @param s The state for which to return initial memory element
	 */
	public int getInitialStateOfTheProduct(int s);

	// Export methods
	
	/**
	 * Export the strategy to a file.
	 */
	public void exportToFile(String file);
	
	// New export methods

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

	// Other new methods
	
	/**
	 * Export the strategy to a file.
	 */
	public void exportStratToFile(File file, StrategyExportType exportType);
	
	/**
	 * Restrict the strategy to the states that are reachable under the strategy.  
	 * @throws PrismException 
	 */
	public void restrictStrategyToReachableStates() throws PrismException;
	
	/**
	 * Clear storage of the strategy.
	 */
	public void clear();
}
