//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <d.a.parker@cs.bham.ac.uk> (University of Birmingham/Oxford)
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
import java.util.List;

import dv.IntegerVector;
import explicit.Distribution;
import prism.Model;
import prism.Prism.StrategyExportType;
import prism.PrismException;
import prism.PrismLog;

/**
 * Class to store a memoryless deterministic (MD) strategy, as an IntegerVector (i.e. stored natively as an array).
 */
public class MDStrategyIV extends MDStrategy
{
	// Model associated with the strategy
	private Model model;
	// Other model info
	private int numStates;
	private List<String> actions;
	// Array storing MD strategy: *action* index (not choice index) for each state
	private IntegerVector iv;
	
	/**
	 * Create an MDStrategyIV from an IntegerVector.
	 */
	public MDStrategyIV(Model model, IntegerVector iv)
	{
		this.model = model;
		numStates = (int) model.getNumStates();
		actions = model.getSynchs();
		this.iv = iv;
	}
	
	// Methods for MDStrategy
	
	@Override
	public int getNumStates()
	{
		return numStates;
	}
	
	@Override
	public boolean isChoiceDefined(int s)
	{
		return iv.getElement(s) >= 0;
	}

	@Override
	public Strategy.Choice getChoice(int s)
	{
		int c = iv.getElement(s);
		switch (c) {
		case -1:
			return Choice.UNKNOWN;
		case -2:
			return Choice.ARBITRARY;
		case -3:
			return Choice.UNREACHABLE;
		default:
			return Choice.INDEX;
		}
	}
	
	@Override
	public int getChoiceIndex(int s)
	{
		return iv.getElement(s);
	}
	
	@Override
	public Object getChoiceAction(int s)
	{
		int c = iv.getElement(s);
		return c >= 0 ? actions.get(c) : c == -1 ? "?" : c == -2 ? "*" : "-";
	}
	
	// Methods for Strategy
	
	@Override
	public Distribution getNextMove(int state) throws InvalidStrategyStateException
	{
		// TODO Auto-generated method stub
		return null;
	}
	
	@Override
	public HashMap<String, Double> getNextAction(int state) throws InvalidStrategyStateException
	{
		if (iv == null || state >= iv.getSize() || state < 0)
			throw new InvalidStrategyStateException("Strategy not defined for state " + state + ".");
		
		if(iv.getElement(state)<0) return null;
		else{
			HashMap<String, Double> actionProbs = new HashMap<String, Double>();
			actionProbs.put(actions.get(iv.getElement(state)), 1.0);
			return actionProbs;
		}
	}
	
	@Override
	public void exportActions(PrismLog out)
	{
		/*try {
			model.exportStrategyToFile(iv, Prism.EXPORT_STRAT_ACTIONS, true, null);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (PrismException e) {
			e.printStackTrace();
		}*/
	}
	
	@Override
	public void exportIndices(PrismLog out)
	{
		/*try {
			model.exportStrategyToFile(iv, Prism.EXPORT_STRAT_INDICES, true, null);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (PrismException e) {
			e.printStackTrace();
		}*/
	}
	
	@Override
	public void exportInducedModel(PrismLog out)
	{
		/*try {
			model.exportStrategyToFile(iv, Prism.EXPORT_STRAT_INDUCED, true, null);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (PrismException e) {
			e.printStackTrace();
		}*/
	}

	@Override
	public void exportDotFile(PrismLog out)
	{
		/*try {
			model.exportStrategyToFile(iv, Prism.EXPORT_STRAT_DOT, true, null);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (PrismException e) {
			e.printStackTrace();
		}*/
	}
	
	@Override
	public void exportStratToFile(File file, StrategyExportType exportType)
	{
		/*try {
			switch (exportType) {
			case ACTIONS:
				model.exportStrategyToFile(iv, Prism.EXPORT_STRAT_ACTIONS, true, file);
				break;
			case INDICES:
				model.exportStrategyToFile(iv, Prism.EXPORT_STRAT_INDICES, true, file);
				break;
			case INDUCED_MODEL:
				model.exportStrategyToFile(iv, Prism.EXPORT_STRAT_INDUCED, true, file);
				break;
			case DOT_FILE:
				model.exportStrategyToFile(iv, Prism.EXPORT_STRAT_DOT, true, file);
				break;
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (PrismException e) {
			e.printStackTrace();
		}*/
	}

	@Override
	public void restrictStrategyToReachableStates() throws PrismException
	{
		// Nothing to do here. It is already done in 'sparse/PS_NondetUntil.cc'
	}

	@Override
	public void clear()
	{
		iv.clear();
		iv = null;
	}
}
