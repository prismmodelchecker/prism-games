//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <d.a.parker@cs.bham.ac.uk> (University of Birmingham/Oxford)
//	* Aistis Simaitis <aistis.aimaitis@cs.ox.ac.uk> (University of Oxford)
//	* Ganindu Prabhashana <ganindu88@gmail.com> (University of Freiburg)
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
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Scanner;

import explicit.Distribution;
import explicit.MDP;
import explicit.MDPSimple;
import explicit.MDPSparse;
import explicit.Model;
import explicit.NondetModel;
import explicit.SMG;
import explicit.STPGExplicit;
import prism.Prism.StrategyExportType;
import prism.PrismLog;

/**
 * Class to store a memoryless deterministic (MD) strategy, as a (Java) array of choice indices.
 */
public class MDStrategyArray extends MDStrategy
{
	// Model associated with the strategy
	private NondetModel model;
	// Index of choice taken in each state (wrt model above) 
	// Other possible values: -1 (unknown), -2 (arbitrary), -3 (unreachable)
	private int choices[];

	/**
	 * Create an MDStrategyArray from an integer array of choices.
	 * The array may later be modified/delete - take a copy if you want to keep it.
	 */
	public MDStrategyArray(NondetModel model, int choices[])
	{
		this.model = model;
		this.choices = choices;
	}
	
	/**
	 * Create an MDStrategyArray from an input stream provided by a scanner.
	 */
	public MDStrategyArray(Scanner scan)
	{
		// ignoring "Adv:" line
		scan.nextLine();
		HashMap<Integer, Integer> choicesMap = new HashMap<Integer, Integer>();
		int s, c;
		while (scan.hasNext()) {
			s = scan.nextInt();
			c = scan.nextInt();
			choicesMap.put(s, c);
		}
		this.choices = new int[choicesMap.size()];
		Iterator<Entry<Integer, Integer>> choicesIter = choicesMap.entrySet().iterator();
		while (choicesIter.hasNext()) {
			Entry<Integer, Integer> pair = choicesIter.next();
	        this.choices[pair.getKey()] = pair.getValue();
	    }
	}
	
	// Methods for MDStrategy

	@Override
	public int getNumStates()
	{
		return model.getNumStates();
	}

	@Override
	public boolean isChoiceDefined(int s)
	{
		return choices[s] >= 0;
	}

	@Override
	public Strategy.Choice getChoice(int s)
	{
		switch (choices[s]) {
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
		return choices[s];
	}

	@Override
	public Object getChoiceAction(int s)
	{
		int c = choices[s];
		return c >= 0 ? model.getAction(s, c) : c == -1 ? "?" : c == -2 ? "*" : "-";
	}

	// Methods for Strategy
	
	@Override
	public Distribution getNextMove(int state) throws InvalidStrategyStateException
	{
		if (choices == null || state >= choices.length || state < 0)
			throw new InvalidStrategyStateException("Strategy not defined for state " + state + ".");
		
		if(choices[state] >= 0){
			Distribution dist = new Distribution();
			dist.add(choices[state], 1);
			return dist;
		}
		else return null;		
	}
	
	@Override
	public HashMap<String, Double> getNextAction(int state) throws InvalidStrategyStateException
	{
		// TODO Auto-generated method stub
		return null;
	}
	
	@Override
	public Model buildProduct(Model model)
	{
		// checking for supported model types
		if (model.getClass().equals(MDPSimple.class)) {
			return this.buildProductMDPSimple((MDPSimple) model);
		}
		if (model.getClass().equals(MDPSparse.class)) {
			return this.buildProductMDPSparse((MDPSparse) model);
		}
		if (model.getClass().equals(STPGExplicit.class)) {
			return this.buildProductSTPGExplicit((STPGExplicit) model);
		}
		if (model.getClass().equals(SMG.class)) {
			return this.buildProductSMG((SMG) model);
		}
		throw new UnsupportedOperationException("The product building is not supported for this class of models");
	}
	
	/**
	 * Builds product of MDPSimple and a given strategy
	 * 
	 * @param model
	 *            model
	 */
	private MDPSimple buildProductMDPSimple(MDPSimple model)
	{
		MDPSimple mdp = new MDPSimple(model);
		int n = mdp.getNumStates();
		int c;
		Distribution distr;

		for (int s = 0; s < n; s++) {

			c = getChoiceIndex(s);

			// if for adversary choice is undefined, taking the first one as default
			if (c < 0)
				c = 0;

			// replacing the choices with the one prescribed by the strategy
			distr = mdp.getChoice(s, c);
			mdp.clearState(s);
			mdp.addChoice(s, distr);
		}
		
		return mdp;
	}

	/**
	 * Builds product of MDPSparse and a given strategy
	 * 
	 * @param model
	 *            model
	 */
	private MDPSparse buildProductMDPSparse(MDPSparse model)
	{
		return new MDPSparse(buildProductMDPSimple(new MDPSimple(model)));
	}

	/**
	 * Builds product between the given two player game and the strategy
	 * Implements strategy for player 1 only, thus returning MDP
	 * 
	 * @param model
	 *            the model
	 * @return strategy
	 */
	private Model buildProductSTPGExplicit(STPGExplicit model)
	{
		STPGExplicit stpg = new STPGExplicit(model);
		int n = stpg.getNumStates();
		int c;
		Distribution distr;

		for (int s = 0; s < n; s++) {
			// checking if the state belong to player 1
			if (stpg.getPlayer(s) != 1)
				// if not then doing nothing
				continue;

			c = getChoiceIndex(s);

			// if for adversary choice is undefined, taking the first one as
			// default
			if (c < 0)
				c = 0;

			// replacing the choices with the one prescribed by the strategy
			distr = stpg.getChoice(s, c);
			stpg.clearState(s);
			stpg.addChoice(s, distr);
		}
		return stpg;
	}

	/**
	 * 
	 * @param model
	 * @return
	 */
	private Model buildProductSMG(SMG model)
	{
		SMG smg = new SMG(model);
		int n = smg.getNumStates();
		int c;
		Distribution distr;

		for (int s = 0; s < n; s++) {
			// checking if the state belong to player 1
			if (smg.getPlayer(s) != 1)
				// if not then doing nothing
				continue;

			c = getChoiceIndex(s);

			// if for adversary choice is undefined, taking the first one as
			// default
			if (c < 0)
				c = 0;

			// replacing the choices with the one prescribed by the strategy
			distr = smg.getChoice(s, c);
			smg.clearState(s);
			smg.addChoice(s, distr);
		}

		return smg;
	}

	@Override
	public void exportToFile(String file)
	{
		// Print adversary
		//PrismLog out = new PrismFileLog(file);
		FileWriter out=null;
		try {
			out = new FileWriter(new File(file));
			out.write(Strategies.FORMAT_STRING_MD_STRAT);
			out.write("\n");
			out.write("Adv:");
			out.write("\n");
			for (int i = 0; i < choices.length; i++) {
				out.write(i + " " + choices[i]);
				out.write("\n");
			}
			out.flush();
		
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		finally{
			try {
				if(out!=null)
				out.close();
			} catch (IOException e) {
				// nothing we can do
			}
		}
	}

	@Override
	public void exportInducedModel(PrismLog out)
	{
		Model dtmcInd = model.constructInducedModel(this);
		dtmcInd.exportToPrismExplicitTra(out);
	}

	@Override
	public void exportDotFile(PrismLog out)
	{
		model.exportToDotFileWithStrat(out, null, choices);
	}

	@Override
	public void exportStratToFile(File file, StrategyExportType exportType)
	{
		// TODO Auto-generated method stub
		
	}

	@Override
	public void restrictStrategyToReachableStates()
	{
		MDP mdp = null;
		if (model instanceof MDP)
			mdp = (MDP) model;
		else
			return;
		BitSet restrict = new BitSet();
		BitSet explore = new BitSet();
		// Get initial states
		for (int is : mdp.getInitialStates()) {
			restrict.set(is);
			explore.set(is);
		}
		// Compute reachable states (store in 'restrict') 
		boolean foundMore = true;
		while (foundMore) {
			foundMore = false;
			for (int s = explore.nextSetBit(0); s >= 0; s = explore.nextSetBit(s + 1)) {
				explore.set(s, false);
				int choiceIndex = getChoiceIndex(s);
				if (choiceIndex >= 0) {
					Iterator<Map.Entry<Integer, Double>> iter = mdp.getTransitionsIterator(s, choiceIndex);
					while (iter.hasNext()) {
						Map.Entry<Integer, Double> e = iter.next();
						int dest = e.getKey();
						if (!restrict.get(dest)) {
							foundMore = true;
							restrict.set(dest);
							explore.set(dest);
						}
					}
				}
			}
		}
		// Set strategy choice for non-reachable state to -3
		int n = mdp.getNumStates();
		for (int s = restrict.nextClearBit(0); s < n; s = restrict.nextClearBit(s + 1)) {
			choices[s] = -3;
		}
	}
	
	@Override
	public void clear()
	{
		choices = null;
	}
	
	@Override
	public String toString()
	{
		return Arrays.toString(choices);
	}
}
