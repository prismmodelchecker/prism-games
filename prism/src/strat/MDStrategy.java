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
import prism.PrismException;
import prism.PrismLog;
import prism.Prism.StrategyExportType;

/**
 * Classes to store memoryless deterministic (MD) strategies.
 */
public abstract class MDStrategy implements Strategy
{
	/**
	 * Current state of model
	 */
	protected int currentState = -1;

	// Getters specifically for MD strategies
	
	/**
	 * Get the number of states of the model associated with this strategy. 
	 */
	public abstract int getNumStates();

	/**
	 * Is choice information stored for state s?
	 */
	public abstract boolean isChoiceDefined(int s);

	/**
	 * Get the type of choice information stored for state s.
	 */
	public abstract Strategy.Choice getChoice(int s);

	/**
	 * Get the index of the choice taken in state s.
	 * The index is defined with respect to a particular model, stored locally.
	 * Other possible values: -1 (unknown), -2 (arbitrary), -3 (unreachable)
	 */
	public abstract int getChoiceIndex(int s);

	/**
	 * Get the action taken in state s.
	 */
	public abstract Object getChoiceAction(int s);

	// Methods for Strategy
	
	@Override
	public String getInfo()
	{
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public int getMemorySize()
	{
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public String getType()
	{
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String getDescription()
	{
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void setInfo(String info)
	{
		// TODO Auto-generated method stub
	}

	@Override
	public void init(int state) throws InvalidStrategyStateException
	{
		currentState = state;
	}

	@Override
	public void updateMemory(int action, int state) throws InvalidStrategyStateException
	{
		currentState = state;
	}

	@Override
	public Distribution getNextMove(int state) throws InvalidStrategyStateException
	{
		// TODO Auto-generated method stub
		return null;
	}

	public abstract HashMap<String,Double> getNextAction(int state) throws InvalidStrategyStateException;;
	
	@Override
	public Object getCurrentMemoryElement()
	{
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void setMemory(Object memory) throws InvalidStrategyStateException
	{
		// TODO Auto-generated method stub
		
	}

	@Override
	public void reset()
	{
		// TODO Auto-generated method stub
		
	}

	@Override
	public Model buildProduct(Model model) throws PrismException
	{
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public int getInitialStateOfTheProduct(int s)
	{
		// TODO Auto-generated method stub
		return 0;
	}
	
	@Override
	public void exportToFile(String file)
	{
		// TODO Auto-generated method stub
	}
	
	@Override
	public void exportActions(PrismLog out)
	{
		int n = getNumStates();
		for (int s = 0; s < n; s++) {
			if (isChoiceDefined(s))
				out.println(s + ":" + getChoiceAction(s));
		}
	}

	@Override
	public void exportIndices(PrismLog out)
	{
		int n = getNumStates();
		for (int s = 0; s < n; s++) {
			if (isChoiceDefined(s))
				out.println(s + ":" + getChoiceIndex(s));
		}
	}
	
	@Override
	public abstract void exportInducedModel(PrismLog out);
	
	@Override
	public abstract void exportDotFile(PrismLog out);
	
	@Override
	public abstract void exportStratToFile(File file, StrategyExportType exportType);
	
	@Override
	public abstract void restrictStrategyToReachableStates() throws PrismException;
	
	@Override
	public abstract void clear();
}
