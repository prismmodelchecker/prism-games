//==============================================================================
//	
//	Copyright (c) 2002-
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

package simulator;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;
import java.util.Vector;

import parser.State;
import parser.VarList;
import parser.ast.ModulesFile;
import parser.type.TypeBool;
import parser.type.TypeDouble;
import parser.type.TypeInt;
import prism.ModelType;
import prism.PrismException;
import prism.PrismLangException;

public class ChoiceExplicit implements Choice
{
	private explicit.MDP model;
	private int state;
	private int choice;
	private String action;
	private int moduleOrActionIndex;

	private List<State> successors;
	private List<Double> probabilities;

	private VarList varList;

	public ChoiceExplicit(ModulesFile modulesFile, explicit.MDP model, int state, int choice) throws PrismLangException
	{
		this.model = model;
		this.state = state;
		this.choice = choice;
		this.action = (String) model.getAction(state, choice);
		try {
			String[] mns = modulesFile.getModuleNames();
			for (int i = 0; i < mns.length; i++) {
				if (mns[i].equals(action)) {
					moduleOrActionIndex = -i - 1;
					throw new Exception(); // break both loops
				}
			}
			Vector<String> ans = modulesFile.getSynchs();
			for (int i = 0; i < ans.size(); i++) {
				if (ans.get(i).equals(action)) {
					moduleOrActionIndex = i + 1;
					throw new Exception(); // break loop
				}
			}
		} catch (Exception e) {
		}

		successors = new ArrayList<State>();
		probabilities = new ArrayList<Double>();
		Iterator<Entry<Integer, Double>> transitions = model.getTransitionsIterator(state, choice);
		while (transitions.hasNext()) {
			Entry<Integer, Double> transition = transitions.next();
			successors.add(model.getStatesList().get(transition.getKey()));
			probabilities.add(transition.getValue());
		}

		this.varList = modulesFile.createVarList(); // for showing the updates

	}

	@Override
	public void scaleProbabilitiesBy(double d)
	{
		for (int i = 0; i < probabilities.size(); i++) {
			probabilities.set(i, probabilities.get(i) * d);
		}
	}

	@Override
	public int getModuleOrActionIndex()
	{
		return moduleOrActionIndex;
	}

	@Override
	public String getModuleOrAction()
	{
		return action;
	}

	@Override
	public int size()
	{
		return model.getNumTransitions(state, choice);
	}

	@Override
	public String getUpdateString(int i, State currentState) throws PrismLangException
	{
		String s = "";
		boolean first = true;
		// assume the variables are in order
		for (int j = 0; j < varList.getNumVars(); j++) {
			Object c_vj = currentState.varValues[j];
			Object s_vj = successors.get(i).varValues[j];
			if (!(s_vj.equals(c_vj))) { // value changed
				if (first)
					first = false;
				else
					s += ", ";
				s += String.format("%s'=%s", varList.getName(j), s_vj);
			}
		}
		return s;
	}

	@Override
	public String getUpdateStringFull(int i)
	{
		// not done quite right because we do not have access to the source model
		String s = "";
		boolean first = true;
		// assume the variables are in order
		for (int j = 0; j < varList.getNumVars(); j++) {
			Object c_vj = model.getStatesList().get(state).varValues[j];
			Object s_vj = successors.get(i).varValues[j];
			if (!(s_vj.equals(c_vj))) { // value changed
				// append string
				if (first)
					first = false;
				else
					s += ", ";
				s += String.format("%1$s'=%2$s", varList.getName(j), s_vj);
			}
		}
		return s;
	}

	@Override
	public State computeTarget(int i, State currentState) throws PrismLangException
	{
		if (i < 0 || i >= size())
			throw new PrismLangException("Choice does not have an element " + i);
		if (indexOf(model.getStatesList(), currentState) != state) {
			throw new PrismLangException("current state inconsistent");
		}
		return successors.get(i);
	}

	@Override
	public void computeTarget(int i, State currentState, State newState) throws PrismLangException
	{
		newState.varValues = computeTarget(i, currentState).varValues;
	}

	@Override
	public double getProbability(int i)
	{
		return probabilities.get(i);
	}

	@Override
	public double getProbabilitySum()
	{
		double sum = 0.0;
		for (int i = 0; i < probabilities.size(); i++) {
			sum += probabilities.get(i);
		}
		return sum;
	}

	@Override
	public int getIndexByProbabilitySum(double x)
	{
		int i, n;
		double d;
		n = size();
		i = 0;
		d = 0.0;
		for (i = 0; x >= d && i < n; i++) {
			d += probabilities.get(i);
		}
		return i - 1;
	}

	public void checkValid(ModelType modelType) throws PrismException
	{
		// TODO
	}

	@Override
	public void checkForErrors(State currentState, VarList varList) throws PrismException
	{
		if (indexOf(model.getStatesList(), currentState) != state) {
			throw new PrismLangException("current state inconsistent");
		}
	}

	@Override
	public String toString()
	{
		int i, n;
		boolean first = true;
		String s = "";
		n = size();
		for (i = 0; i < n; i++) {
			if (first)
				first = false;
			else
				s += " + ";
			s += getProbability(i) + ":" + successors.get(i);
		}
		return s;
	}

	private <T> int indexOf(List<T> list, T o)
	{ // indexOf based on reference, not equals(o) function
		for (int i = 0; i < list.size(); i++) {
			if (list.get(i) == o) {
				return i;
			}
		}
		return -1; // not in list
	}

}
