//==============================================================================
//	
//	Copyright (c) 2014-
//	Authors:
//	* Clemens Wiltsche <clemens.wiltsche@cs.ox.ac.uk> (University of Oxford)
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

import java.io.FileWriter;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

import parser.Values;
import prism.PrismException;
import prism.PrismLog;
import explicit.Distribution;
import explicit.Model;
import explicit.SMG;

public class StochasticUpdateStrategyProduct implements Strategy
{
	// model info
	protected Values lastConstants;
	
	// strategy info
    protected String info = "No information available.";
  
    protected List<StochasticUpdateStrategy> strategies;

    protected List<Integer> controlled_by;
    protected List<List<Integer>> local_state; // maps (global state, component) -> local state
    protected List<List<List<Integer>>> l_action; // maps (global state, global action, component) -> local action

    protected int lastState;

    protected int memorySize;

    public int getLocalState(int state, int component)
    {
	return local_state.get(state).get(component);
    }

    @Override
    public void init(int state) throws InvalidStrategyStateException
    {
	for(int i = 0; i < strategies.size(); i++) {
	    strategies.get(i).init(local_state.get(state).get(i));
	}

	lastState = state;
    }

    @Override
    public Distribution getNextMove(int state) throws InvalidStrategyStateException
    {
	Distribution result = new Distribution();
	int cb = controlled_by.get(state); // who controls the state
	if(cb<0) return new Distribution(); // not controlled by P1: return empty Distribution
	int ls = local_state.get(state).get(cb); // which state is the controlling component at
	Distribution local;
	try {
	    local = strategies.get(cb).getNextMove(ls);
	} catch (InvalidStrategyStateException e) {
	    throw new InvalidStrategyStateException(String.format("Inconsistent strategy state"));
	}

	// need to find global action for each local action in support
	l_action.get(state);
	for(Integer ls_succ : local.getSupport()) {
	    int ga = 0;
	    search_for_ga:
	    for(ga = 0; ga < l_action.get(state).size(); ga++) { // search through global actions
		if(l_action.get(state).get(ga).get(cb) == ls_succ) {
		    result.set(ga, local.get(ls_succ));
		    break search_for_ga;
		}
	    }
	}
	return result;
    }


    @Override
    public void updateMemory(int action, int state) throws InvalidStrategyStateException
    {
	for(int i = 0; i < strategies.size(); i++) {
	    int ls = local_state.get(state).get(i);
	    int la = l_action.get(lastState).get(action).get(i);
	    if(la>=0) { // move present locally
		strategies.get(i).updateMemory(la, ls);
	    }
	}
	lastState = state;
    }

    public List<Distribution> memoryUpdate(int action, int state) throws InvalidStrategyStateException
    {
	List<Distribution> result = new ArrayList<Distribution>(strategies.size());
	for(int i = 0; i < strategies.size(); i++) {
	    int ls = local_state.get(state).get(i);
	    int la = l_action.get(lastState).get(action).get(i);
	    if(la>=0) { // move present locally
		result.add(strategies.get(i).memoryUpdate(la, ls));
	    } else {
		result.add(null);
	    }
	}
	return result;
    }

    public String memoryUpdateString(int state, int choice, int next, NumberFormat df) throws InvalidStrategyStateException
    {
		Distribution dist = getNextMove(state);
		List<Distribution> mus = memoryUpdate(choice, next);
		String label = df.format(dist.get(choice)) + " [";
		boolean first1 = true;
		for (int d = 0; d < mus.size(); d++) {
			Distribution mu = mus.get(d);
			int l_next = getLocalState(next, d); // local next state
			if (first1)
				first1 = false;
			else
				label += ", ";
			if (mu == null) {
				label += "\u2013";
			} else {
				label += String.format("mu(%d): {", d);
				boolean first2 = true;
				for (Integer m : mu.getSupport()) {
					if (first2)
						first2 = false;
					else
						label += ", ";
					label += String.format("(%d, %d)=", l_next, m) + df.format(mu.get(m));
				}
				label += "}";
			}
		}
		return label + "]";
    }

    @Override
    public void reset()
    {
	for(StochasticUpdateStrategy strat : strategies) {
	    strat.reset();
	}
    }

    @Override
    public void exportToFile(String filename)
    {
	FileWriter out = null;
	try {
	    out = new FileWriter(filename);
	    out.write(this.toString());
	    out.flush();
	} catch (IOException e) {
	    e.printStackTrace();
	} finally {
	    if (out != null)
		try {
		    out.close();
		} catch (IOException e) {
		}
	}
    }

    @Override
    public String toString()
    {
	String strats = Strategies.FORMAT_STRING_SU_STRAT_COMP + "\n";
	strats += "// Composed Stochastic Update Strategy\n";
	for(int i = 0; i < strategies.size(); i++)
	    strats += String.format("%s\n", strategies.get(i).toString());
	return strats;
    }

    @Override
    public Model buildProduct(Model model) throws PrismException
    {
	if(!model.getClass().equals(SMG.class)) {
	    throw new PrismException("Unsupported model type");
	}
	throw new PrismException("Not supported for stochastic memory update strategy");
    }


    @Override
    public void setInfo(String info)
    {
	this.info = info;
    }

    @Override
    public String getInfo()
    {
	return info;
    }
    
    @Override
    public int getMemorySize()
    {
	return memorySize;
    }
    
    @Override
    public String getType()
    {
	return Strategies.FORMAT_STRING_SU_STRAT_COMP;
    }
    
    @Override
    public Object getCurrentMemoryElement()
    {
	List<Object> memory = new ArrayList<Object>();
	for(StochasticUpdateStrategy strat : strategies) {
	    memory.add(strat.getCurrentMemoryElement());
	}
	return memory;
    }
    
    @Override
    public void setMemory(Object memory) throws InvalidStrategyStateException
    {
	if(memory instanceof List && ((List)memory).size() == strategies.size()) {
	    for(int i = 0; i < strategies.size(); i++) {
		strategies.get(i).setMemory(((List)memory).get(i));
	    }
	    // set lastState
	    List<Integer> lss = new ArrayList<Integer>(strategies.size());
	    for(int i = 0; i < strategies.size(); i++) {
		lss.add(strategies.get(i).lastState);
	    }
	    lastState = local_state.indexOf(lss);
	} else {
	    throw new InvalidStrategyStateException("Memory has incompatible format");
	}
    }

    @Override
    public String getDescription()
    {
	String desc = "";
	desc += "Product of stochastic update strategies\n";
	desc += "Size of memory: " + getMemorySize() + "\n";
	return desc;
    }

    @Override
    public int getInitialStateOfTheProduct(int s)
    {
	return -1; // not available for SU strategies
    }

    private void construct(List<Strategy> strats) throws PrismException
    {
	this.strategies = new ArrayList<StochasticUpdateStrategy>();
	for(Strategy strat : strats) {
	    if(strat == null)
		throw new PrismException("At least one component does not have a strategy");
	    if(!(strat instanceof StochasticUpdateStrategy))
		throw new PrismException("Product construction only applicable to SU strategies");
	    strategies.add((StochasticUpdateStrategy)strat);
	}

	memorySize = 1;
	for(int i = 0; i < strats.size(); i++) {
	    this.memorySize *= strategies.get(i).getMemorySize(); // product of memory sizes
	}
    }

    public StochasticUpdateStrategyProduct(List<Strategy> strats) throws PrismException
    {	
	construct(strats);
    }

    // construct strategy from file
    public StochasticUpdateStrategyProduct(Scanner scan) throws PrismException
    {
	String nextLine;
	List<Strategy> strats = new ArrayList<Strategy>();

	nextLine = scan.nextLine(); // first line in file
	while(scan.hasNext()) {
	    if(nextLine.startsWith("start strategy")) {
		Strategy strat = new StochasticUpdateStrategy(scan);
		strats.add(strat);
	    } else {
		nextLine = scan.nextLine();
	    }
	}

	    construct(strats);
    }

    public void setComposition(SMG smg)
    {
	this.controlled_by = new ArrayList<Integer>(smg.getControlledBy());
	this.local_state = new ArrayList<List<Integer>>(smg.getLocalState());
	this.l_action = new ArrayList<List<List<Integer>>>(smg.getLAction());
    }

	@Override
	public void exportActions(PrismLog out)
	{
	    out.print(this.toString());
	}


	@Override
	public void exportIndices(PrismLog out)
	{
		// TODO Auto-generated method stub
		
	}

	@Override
	public void exportInducedModel(PrismLog out)
	{
		// TODO Auto-generated method stub
		
	}

	@Override
	public void initialise(int s)
	{
		// TODO Auto-generated method stub
		
	}

	@Override
	public void update(Object action, int s)
	{
		// TODO Auto-generated method stub
		
	}

	@Override
	public Object getChoiceAction()
	{
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void clear()
	{
		// TODO Auto-generated method stub
		
	}

	@Override
	public void exportDotFile(PrismLog out)
	{
		// TODO Auto-generated method stub
		
	}
	
	@Override
	public void setConstants(Values lastConstants) {
		// TODO Auto-generated method stub
		this.lastConstants = lastConstants;
	}

	@Override
	public Values getConstants() {
		// TODO Auto-generated method stub
		return lastConstants;
	}

}
