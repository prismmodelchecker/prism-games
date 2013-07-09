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

import explicit.Distribution;
import explicit.Model;
import prism.PrismException;
import prism.PrismLog;

/**
 * Classes to store memoryless deterministic (MD) strategies.
 */
public abstract class MDStrategy implements Strategy
{
	public abstract int getNumStates();
	public abstract int getChoice(int i);
	public abstract Object getChoiceAction(int i);
	
	public void export(PrismLog out)
	{
		int n = getNumStates();
		for (int i = 0; i < n; i++) {
			out.println(i + ":" + getChoice(i));
		}
	}

	// Temp stubs

	@Override
	public void init(int state) throws InvalidStrategyStateException
	{
		// TODO Auto-generated method stub
		
	}

	@Override
	public void updateMemory(int action, int state) throws InvalidStrategyStateException
	{
		// TODO Auto-generated method stub
		
	}

	@Override
	public Distribution getNextMove(int state) throws InvalidStrategyStateException
	{
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void reset()
	{
		// TODO Auto-generated method stub
		
	}

	@Override
	public void exportToFile(String file)
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
	public String getInfo()
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
	public String getStateDescription()
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
}
