//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Aistis Simaitis <aistis.aimaitis@cs.ox.ac.uk> (University of Oxford)
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
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Scanner;

import explicit.Distribution;
import explicit.MDPSimple;
import explicit.MDPSparse;
import explicit.Model;
import explicit.SMG;
import explicit.STPGExplicit;
import prism.Prism.StrategyExportType;
import prism.PrismException;
import prism.PrismLog;

/**
 * Implementation of the memoryless deterministic strategy
 */
public class MemorylessDeterministicStrategy implements Strategy
{
	private Distribution[] choices;
	private String info = "No information available.";

	public MemorylessDeterministicStrategy(int[] choices)
	{
		this.choices = new Distribution[choices.length];
		Distribution dist;
		for (int i = 0; i < choices.length; i++) {
			dist = new Distribution();
			dist.add(choices[i] < 0 ? 0 : choices[i], 1);
			this.choices[i] = dist;
		}
	}

	/**
	 * Creates a MemorylessDeterministicStrategy from the input stream provided by the scanner.
	 *
	 * @param scan
	 */
	public MemorylessDeterministicStrategy(Scanner scan)
	{
		// ignoring "Adv:" line
		scan.nextLine();
		List<Distribution> dists = new LinkedList<Distribution>();
		Distribution dist;
		int s, c;
		while (scan.hasNext()) {
			s = scan.nextInt();
			c = scan.nextInt();
			dist = new Distribution();
			dist.add(c, 1);
			dists.add(dist);
		}
		this.choices = dists.toArray(new Distribution[] {});
	}

	@Override
	public String getInfo()
	{
		return info;
	}

	@Override
	public int getMemorySize()
	{
		return 0;
	}

	@Override
	public String getType()
	{
		return "Memoryless deterministic";
	}

	@Override
	public String getDescription()
	{
		String desc = "";
		desc += "Memoryless deterministic strategy\n";
		desc += "Size of memory: 0\n";
		desc += "Size of next move function: " + choices.length + "\n";
		return desc;
	}

	@Override
	public void setInfo(String info)
	{
		this.info = info;
	}

	@Override
	public void init(int state) throws InvalidStrategyStateException
	{
		// do nothing
	}

	@Override
	public void updateMemory(int action, int state) throws InvalidStrategyStateException
	{
		//currentState = state;
	}

	@Override
	public Distribution getNextMove(int state) throws InvalidStrategyStateException
	{

		if (choices == null || state >= choices.length || state < 0)
			throw new InvalidStrategyStateException("Strategy not defined for state " + state + ".");

		return choices[state];
	}

	@Override
	public HashMap<String, Double> getNextAction(int state) throws InvalidStrategyStateException
	{
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Object getCurrentMemoryElement()
	{
		// System.out.println("Memory element requested");
		return null;
	}

	@Override
	public void setMemory(Object memory) throws InvalidStrategyStateException
	{
		// do nothing
		// System.out.println("Set memory element");
	}

	@Override
	public void reset()
	{
		// do nothing
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
	 */
	public MDPSimple buildProductMDPSimple(MDPSimple model)
	{
		MDPSimple mdp = new MDPSimple(model);
		int n = mdp.getNumStates();
		int c;
		Distribution distr;

		for (int s = 0; s < n; s++) {

			// getting the choice of a strategy for this state
			try {
				c = this.getNextMove(s).keySet().iterator().next();

				// if for adversary choice is undefined, taking the first one as
				// default
				if (c < 0)
					c = 0;
			} catch (InvalidStrategyStateException e) {
				// strategy undefined for this state -- keep calm and carry
				// on
				continue;
			}

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
	public MDPSparse buildProductMDPSparse(MDPSparse model)
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

			// getting the choice of a strategy for this state
			try {
				c = this.getNextMove(s).keySet().iterator().next();

				// if for adversary choice is undefined, taking the first one as
				// default
				if (c < 0)
					c = 0;
			} catch (InvalidStrategyStateException e) {
				// strategy undefined for this state -- keep calm and carry
				// on
				continue;
			}

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

			// getting the choice of a strategy for this state
			try {
				c = this.getNextMove(s).keySet().iterator().next();

				// if for adversary choice is undefined, taking the first one as
				// default
				if (c < 0)
					c = 0;
			} catch (InvalidStrategyStateException e) {
				// strategy undefined for this state -- keep calm and carry
				// on
				continue;
			}

			// replacing the choices with the one prescribed by the strategy
			distr = smg.getChoice(s, c);
			smg.clearState(s);
			smg.addChoice(s, distr);
		}

		return smg;
	}

	@Override
	public int getInitialStateOfTheProduct(int s)
	{
		return -1;
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
			out.write(i + " " + choices[i].keySet().iterator().next());
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
	public void exportActions(PrismLog out)
	{
		// TODO Auto-generated method stub
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
	public void exportDotFile(PrismLog out)
	{
		// TODO Auto-generated method stub
		
	};
	
	@Override
	public void exportStratToFile(File file, StrategyExportType exportType)
	{
		// TODO Auto-generated method stub
		
	}

	@Override
	public void restrictStrategyToReachableStates() throws PrismException
	{
		// TODO Auto-generated method stub
		throw new PrismException("Reach option is not supported for this strategy type");
	}

	@Override
	public void clear()
	{
		// TODO Auto-generated method stub
	}

	@Override
	public String toString()
	{
		return Arrays.toString(choices);
	}
}
