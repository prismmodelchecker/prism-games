//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford)
//  * Gabriel Santos <gabriel.santos@cs.ox.ac.uk> (University of Oxford)
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
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import explicit.CSG;
import explicit.DTMCSimple;
import explicit.Distribution;
import explicit.MDPSimple;
import explicit.Model;
import explicit.ModelCheckerResult;
import parser.State;
import parser.VarList;
import parser.ast.Coalition;
import parser.ast.Declaration;
import parser.ast.DeclarationInt;
import parser.ast.Expression;
import prism.Prism.StrategyExportType;
import prism.PrismComponent;
import prism.PrismException;
import prism.PrismLog;
import prism.PrismPrintStreamLog;

public class CSGStrategy extends PrismComponent implements Strategy {

	protected CSG model;
	protected List<List<List<Map<BitSet, Double>>>> csgchoices; // player -> iteration -> state -> indexes -> value
	protected ModelCheckerResult[] prechoices;
	protected BitSet[] targets;
	protected Map<BitSet, BitSet> subgames;
	protected CSGStrategyType type;
	protected BitSet no;
	protected BitSet yes;
	protected BitSet inf;
	protected int numCoalitions;
	
	public enum CSGStrategyType {
		ZERO_SUM, EQUILIBRIA_M, EQUILIBRIA_P, EQUILIBRIA_R;
	}

	public CSGStrategy(CSG model, List<List<List<Map<BitSet, Double>>>> csgchoices, Map<BitSet, BitSet> subgames, int numCoalitions, CSGStrategyType type) {
		this.model = model;
		this.csgchoices = csgchoices;
		this.subgames = subgames;
		this.numCoalitions = numCoalitions;
		this.type = type;
	}
	
	public CSGStrategy(CSG model, List<List<List<Map<BitSet, Double>>>> csgchoices, ModelCheckerResult[] prechoices, BitSet[] targets, CSGStrategyType type) {
		this.model = model;
		this.csgchoices = csgchoices;
		this.prechoices = prechoices;
		this.targets = targets;
		this.type = type;
	}

	public CSGStrategy(CSG model, List<List<List<Map<BitSet, Double>>>> csgchoices, BitSet no, BitSet yes, BitSet inf, CSGStrategyType type) {
		this.model = model;
		this.csgchoices = csgchoices;
		this.prechoices = null;
		this.targets = null;
		this.no = no;
		this.yes = yes;
		this.inf = inf;
		this.type = type;
	}
	
	@Override
	public String getInfo() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public int getMemorySize() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public String getType() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String getDescription() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void setInfo(String info) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void init(int state) throws InvalidStrategyStateException {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void updateMemory(int action, int state) throws InvalidStrategyStateException {
		// TODO Auto-generated method stub
		
	}

	@Override
	public Distribution getNextMove(int state) throws InvalidStrategyStateException {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public HashMap<String, Double> getNextAction(int state) throws InvalidStrategyStateException {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Object getCurrentMemoryElement() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void setMemory(Object memory) throws InvalidStrategyStateException {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void reset() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public Model buildProduct(Model model) throws PrismException {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public int getInitialStateOfTheProduct(int s) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public void exportToFile(String file) {
		// TODO Auto-generated method stub
	}

	@Override
	public void exportActions(PrismLog out) {
		try {
			switch(type) {
				case ZERO_SUM:
					exportZeroSumStrategy(out);
					break;
				case EQUILIBRIA_M:
					exportMultiEquilibriaStrategy(out);
					break;
				case EQUILIBRIA_P:
				case EQUILIBRIA_R:
					exportEquilibriaStrategy(out);
					break;
			}
		}
		catch (Exception e) {
			e.printStackTrace();
		}
	}

	@Override
	public void exportIndices(PrismLog out) {
		// TODO Auto-generated method stub
	}

	@Override
	public void exportInducedModel(PrismLog out) {
		// TODO Auto-generated method stub
	}
	
	@Override
	public void exportDotFile(PrismLog out) {
		// TODO Auto-generated method stub
	}
	
	@Override
	public void exportStratToFile(File file, StrategyExportType exportType) {
		// TODO Auto-generated method stub
	}
	
	public void exportZeroSumStrategy(PrismLog out) throws PrismException {
		MDPSimple mdp = new MDPSimple();
		Map<Integer, Integer> onmap = new HashMap<Integer, Integer>();
		List<State> statelist = new ArrayList<State>();
		VarList varlist = null;
		State initial;
		BitSet explored =  new BitSet();
		int n, s;
		s = model.getFirstInitialState();
		initial = model.getStatesList().get(s);
		varlist = model.getVarList();
		n = mdp.addState();
		mdp.addInitialState(n);
		mdp.setVarList(varlist);
		statelist.add(n, initial);
		onmap.put(s, n);
		generateMDPZeroSum(mdp, onmap, statelist, explored, 0, 0, s);
		mdp.setStatesList(statelist);
		mdp.exportToDotFile(out, null, true);
		out.print("\n/*");
		out.print("\n -- Transitions --  \n");
		mdp.exportToPrismExplicitTra(out);	
		out.print("\n -- States --  \n");
		mdp.exportStates(0, mdp.getVarList(), out);
		out.print("*/\n");
		mainLog.println("Additional info on transitions and states added to file.");
	}
	
	public void exportEquilibriaStrategy(PrismLog out) throws PrismException, InvalidStrategyStateException {
		MDPSimple mdp = new MDPSimple();
		Map<Integer, Integer> onmap = new HashMap<Integer, Integer>();
		List<State> statelist = new ArrayList<State>();
		VarList varlist = null;
		State initial;
		BitSet explored =  new BitSet();
		int i, n, p, s;
		s = model.getFirstInitialState();
		initial = model.getStatesList().get(s);
		varlist = model.getVarList();
		n = mdp.addState();
		mdp.addInitialState(n);		
		mdp.setVarList(varlist);
		statelist.add(n, initial);
		onmap.put(s, n);
		if (targets[0].isEmpty() && targets[1].isEmpty()) {
			Distribution d = new Distribution();
			d.add(n, 1.0);
			mdp.addActionLabelledChoice(n, d, "CSG: Unsat(0) -- Unsat(1)");
		}
		else {
			BitSet[] reach = new BitSet[targets.length];
			for (p = 0; p < 2; p++) {
				reach[p] = new BitSet();
				for (i = 0; i < model.getNumStates(); i++) {
					if (type == CSGStrategyType.EQUILIBRIA_P)
						if (prechoices[p].soln[i] > 0)
							reach[p].set(i);
					if  (type == CSGStrategyType.EQUILIBRIA_R)
						if (prechoices[p].soln[i] < Double.POSITIVE_INFINITY)
							reach[p].set(i);
				}
			}
			generateMDPEquilibria(mdp, onmap, statelist, reach, explored, 0, s); 
			addPrecompStrategies(mdp, onmap, statelist, reach);		
		}
		mdp.setStatesList(statelist);
		mdp.exportToDotFile(out, null, true);
		out.print("\n/*");
		out.print("\n -- Transitions --  \n");
		mdp.exportToPrismExplicitTra(out);	
		out.print("\n -- States --  \n");
		mdp.exportStates(0, mdp.getVarList(), out);
		out.print("*/\n");
		mainLog.println("Additional info on transitions and states added to file.");
	}
	
	public void exportMultiEquilibriaStrategy(PrismLog out) throws PrismException {
		MDPSimple mdp = new MDPSimple();
		Map<Integer, Integer> onmap = new HashMap<Integer, Integer>();
		List<State> statelist = new ArrayList<State>();
		State initial;
		BitSet explored =  new BitSet();
		int n, s;
		s = model.getFirstInitialState();
		initial = model.getStatesList().get(s);
		n = mdp.addState();
		mdp.addInitialState(n);		
		mdp.setVarList(model.getVarList());
		statelist.add(n, initial);
		onmap.put(s, n);
		generateMDPMultiEquilibria(mdp, onmap, statelist, explored, 0, s);
		mdp.setStatesList(statelist);
		mdp.exportToDotFile(out, null, true);
		out.print("\n/*");
		out.print("\n -- Transitions --  \n");
		mdp.exportToPrismExplicitTra(out);	
		out.print("\n -- States --  \n");
		mdp.exportStates(0, mdp.getVarList(), out);
		out.print("*/\n");
		mainLog.println("Additional info on transitions and states added to file.");
	}
	
	public void addPrecompStrategies(MDPSimple mdp, Map<Integer, Integer> onmap, List<State> statelist, BitSet[] reach) throws PrismException, InvalidStrategyStateException {
		BitSet[] minus = new BitSet[targets.length];
		BitSet[] goals = new BitSet[targets.length];
		BitSet explored = new BitSet();
		int i, p;
		for (p = 0; p < 2; p++) {
			minus[p] = new BitSet();
			minus[p].or(targets[p]);	
			minus[p].andNot(targets[(p + 1) % 2]);
			goals[p] = new BitSet();
			goals[p].or(targets[p]);
		}
		for (p = 0; p < 2; p++) {
			for (i = minus[p].nextSetBit(0); i >= 0; i = minus[p].nextSetBit(i + 1)) {
				explored = new BitSet();
				if (statelist.contains(model.getStatesList().get(i)) && !(goals[0].get(i) && goals[1].get(i))) {
					addPrecompStrategies(mdp, onmap, statelist, goals, reach, explored, p, i);
				}
			}
		}
		
	}
	
	public void addPrecompStrategies(MDPSimple mdp, Map<Integer, Integer> onmap, List<State> statelist, BitSet[] goals, BitSet[] reach, BitSet explored, int p, int s) throws PrismException, InvalidStrategyStateException {
		Distribution d;
		String label = "MDP: ";
		String joint = "";
		double v;
		int c, i, m, n;
		if (!onmap.containsKey(s)) {
			n = mdp.addState();
			onmap.put(s, n);
			statelist.add(n, model.getStatesList().get(s));
		}
		else {
			n = onmap.get(s);
		}
		for (int t = 0; t < model.getNumChoices(s); t++) {
			for (int u : model.getChoice(s, t).getSupport()) {
				if (goals[0].get(s))
					goals[0].set(u);
				if (goals[1].get(s))
					goals[1].set(u);
			}
		}
		explored.set(s);
		if (goals[p].get(s) && goals[(p + 1) % 2].get(s)) {
			d = new Distribution();
			d.add(n, 1.0);
			mdp.addActionLabelledChoice(n, d, "CSG: Sat(0) -- Sat(1)");
		}
		else if (goals[p].get(s) && !goals[(p +1) % 2].get(s) && !reach[(p +1) % 2].get(s)) {
			d = new Distribution();
			d.add(n, 1.0);
			mdp.addActionLabelledChoice(n, d, "CSG: Sat(" + p + ") -- Unsat(" + (p + 1) % 2 + ")");
		}
		else if (!goals[(p +1) % 2].get(s)) {
			d = new Distribution();
			c = prechoices[(p + 1) % 2].strat.getNextMove(s).getSupport().size();
			for (int t : prechoices[(p + 1) % 2].strat.getNextMove(s).getSupport()) {
				v = prechoices[(p + 1) % 2].strat.getNextMove(s).get(t);
				for (int u : model.getChoice(s, t).getSupport()) {
					if (!onmap.containsKey(u)) {
						m = mdp.addState();
						onmap.put(u, m);
						statelist.add(m, model.getStatesList().get(u));
					}
					else {
						m = onmap.get(u);
					}
					if (!explored.get(u))
						addPrecompStrategies(mdp, onmap, statelist, goals, reach, explored, p, u);							
					d.add(m, v * model.getChoice(s, t).get(u));
				}
				for (i = 0; i < model.getActions(s, t).length; i++) {
					joint += "[" + model.getActions(s, t)[i] + "]";
				}
				c--;
				label += v + ": " + joint + ((c > 0)? " + " : "");
			}
			mdp.addActionLabelledChoice(n, d, label);
		}
	}
	
	public void localMixedProduct(Map<BitSet, Double> prods, BitSet prod, double v, int k, int p, int s) {
		if (p < csgchoices.size() - 1) {
			for (BitSet strat : csgchoices.get(p).get(k).get(s).keySet()) {
				BitSet newprod = new BitSet();
				double newv = v * csgchoices.get(p).get(k).get(s).get(strat);
				newprod.or(prod);
				newprod.or(strat);
				localMixedProduct(prods, newprod, newv, k, p + 1, s);
			}
		}
		else {
			for (BitSet strat : csgchoices.get(p).get(k).get(s).keySet()) {
				BitSet newprod = new BitSet();
				newprod.or(prod);
				newprod.or(strat);
				double newv = v * csgchoices.get(p).get(k).get(s).get(strat);
				prods.put(newprod, newv);
			}	
		}
	}
	
	public void generateMDPMultiEquilibria(MDPSimple mdp, Map<Integer, Integer> onmap, List<State> statelist, BitSet explored, int k, int s) {
		Distribution d;
		Map<BitSet, Double> prods = new HashMap<BitSet, Double>();
		BitSet tmp = new BitSet();
		BitSet sat = new BitSet();
		String[] action = new String[csgchoices.size()];
		String joint = null;
		String label = null;
		String lsubg = "";
		int c, i, m, n, q, p, t;
		boolean chck = true;
		boolean loop = false;
		explored.set(s);
		n = onmap.get(s);
		//System.out.println(subgames);
		for (BitSet subgame : subgames.keySet()) {
			if (subgames.get(subgame).get(s)) {
				sat.or(subgame);
			}
		}
		for (p = 0; p < numCoalitions; p++) {
			if (sat.get(p)) {
				lsubg += "Sat(" + p + ")";
			}
			else {
				lsubg += "Unsat(" + p + ")";
			}
			if (p < numCoalitions - 1)
				lsubg += " -- ";
		}
		if (sat.cardinality() == numCoalitions) {
			d = new Distribution();
			d.add(n, 1.0);
			mdp.addActionLabelledChoice(n, d, lsubg);	
		}
		else {
			for (p = 0; p < numCoalitions; p++) {
				chck = chck && csgchoices.get(p).get(k).get(s) != null;
				action[p] = "";
			}
			if (chck) {
				localMixedProduct(prods, new BitSet(), 1.0, 0, 0, s);
				d = new Distribution();
				for (t = 0; t < model.getNumChoices(s); t++) {
					tmp.clear();
					for (q = 0; q < model.getIndexes(s, t).length; q++) {
						i = model.getIndexes(s, t)[q];						
						tmp.set((i > 0)? i : model.getIdles()[q]);
					}
					if (prods.containsKey(tmp)) {						
						for (int u : model.getChoice(s, t).getSupport()) {
							if (!onmap.containsKey(u)) {
								m = mdp.addState();
								onmap.put(u, m);
								statelist.add(m, model.getStatesList().get(u));
								if (!explored.get(u))
									generateMDPMultiEquilibria(mdp, onmap, statelist, explored, k, u);
							}
							else {
								m = onmap.get(u);
								if (m == n && model.getChoice(s, t).getSupport().size() == 1 && model.getNumChoices(s) == 1)
									loop = true;
							}
							d.add(m, model.getChoice(s, t).get(u) * prods.get(tmp));
						}
					}
				}
				if (loop) {
					mdp.addActionLabelledChoice(n, d, lsubg);	
				}
				else if (!d.isEmpty()) {
					label = "CSG: ";
					for (p = 0; p < numCoalitions; p++) {
						c = csgchoices.get(p).get(k).get(s).keySet().size();
						for (BitSet act : csgchoices.get(p).get(k).get(s).keySet()) {
							joint = "";
							for (i = act.nextSetBit(0); i >= 0; i = act.nextSetBit(i + 1)) {
								joint += "[" + model.getActions().get(i - 1) + "]";
							}
							c--;
							action[p] += csgchoices.get(p).get(k).get(s).get(act) +": " + joint + ((c > 0)? " + " : ""); 
						}
						label += (p + 1 < csgchoices.size())? action[p] + " -- " : action[p];
					}
					mdp.addActionLabelledChoice(n, d, label);
				}
			}
		}
	}
	
	public void generateMDPEquilibria(MDPSimple mdp, Map<Integer, Integer> onmap, List<State> statelist, BitSet[] reach, BitSet explored, int k, int s) {
		Distribution d;
		String[] action = new String[csgchoices.size()];
		String label = null;
		String joint = null;
		BitSet tmp = new BitSet();
		Map<BitSet, Double> prods = new HashMap<BitSet, Double>();
		int c, i, m, n, p, q, t;
		boolean chck = true;
		n = onmap.get(s);
		explored.set(s);
		if (targets[0].get(s) && targets[1].get(s)) {
			d = new Distribution();
			d.add(n, 1.0);
			mdp.addActionLabelledChoice(n, d, "CSG: Sat(0) -- Sat(1)");
		}
		else if (targets[0].get(s) && !targets[1].get(s) && !reach[0].get(s)) {
			d = new Distribution();
			d.add(n, 1.0);
			mdp.addActionLabelledChoice(n, d, "CSG: Sat(0) -- Unsat(1)");
		}
		else if (!targets[0].get(s) && targets[1].get(s) && !reach[1].get(s)) {
			d = new Distribution();
			d.add(n, 1.0);
			mdp.addActionLabelledChoice(n, d, "CSG: Unsat(0) -- Sat(1)");
		}
		else if (!targets[0].get(s) && !targets[1].get(s) && !reach[0].get(s) && !reach[1].get(s)) {
			d = new Distribution();
			d.add(n, 1.0);
			mdp.addActionLabelledChoice(n, d, "CSG: Unsat(0) -- Unsat(1)");
		}
		else {
			for (p = 0; p < 2; p++) {
				chck = chck && csgchoices.get(p).get(k).get(s) != null;
				action[p] = "";
			}
			if (chck) {
				localMixedProduct(prods, new BitSet(), 1.0, 0, 0, s);
				d = new Distribution();
				for (t = 0; t < model.getNumChoices(s); t++) {
					tmp.clear();
					for (q = 0; q < model.getIndexes(s, t).length; q++) {
						i = model.getIndexes(s, t)[q];						
						tmp.set((i > 0)? i : model.getIdles()[q]);
					}
					if (prods.containsKey(tmp)) {
						for (int u : model.getChoice(s, t).getSupport()) {
							if (!onmap.containsKey(u)) {
								m = mdp.addState();
								onmap.put(u, m);
								statelist.add(m, model.getStatesList().get(u));
								if (!explored.get(u))
									generateMDPEquilibria(mdp, onmap, statelist, reach, explored, k, u);  // should check for explored?
							}
							else {
								m = onmap.get(u);
							}
							d.add(m, model.getChoice(s, t).get(u) * prods.get(tmp));
						}
					}
				}
				if (!d.isEmpty()) {
					label = "CSG: ";
					for (p = 0; p < 2; p++) {
						c = csgchoices.get(p).get(k).get(s).keySet().size();
						for (BitSet act : csgchoices.get(p).get(k).get(s).keySet()) {
							joint = "";
							for (i = act.nextSetBit(0); i >= 0; i = act.nextSetBit(i + 1)) {
								joint += "[" + model.getActions().get(i - 1) + "]";
							}
							c--;
							action[p] += csgchoices.get(p).get(k).get(s).get(act) +": " + joint + ((c > 0)? " + " : ""); 
						}
						label += (p + 1 < csgchoices.size())? action[p] + " -- " : action[p];
					}
					mdp.addActionLabelledChoice(n, d, label);
				}
			}
		}
	}
	
	public void generateMDPZeroSum(MDPSimple mdp, Map<Integer, Integer> onmap, List<State> statelist, BitSet explored, int k, int p, int s) {
		Distribution d;
		BitSet tmp1 = new BitSet();
		BitSet tmp2 = new BitSet();
		String act1 = null;
		String act2 = null;
		int i, m, n, q, t;
		n = onmap.get(s);
		explored.set(s);
		if (yes.get(s)) {
			d = new Distribution();
			d.add(n, 1.0);
			mdp.addActionLabelledChoice(n, d, "Sat");
		}
		else if (no.get(s)) {
			d = new Distribution();
			d.add(n, 1.0);
			mdp.addActionLabelledChoice(n, d, "Unsat");
		}
		else if (inf.get(s)) {
			d = new Distribution();
			d.add(n, 1.0);
			mdp.addActionLabelledChoice(n, d, "Infinity");
		}
		else if (csgchoices.get(p).get(k).get(s) != null) {
			for (t = 0; t < model.getNumChoices(s); t++) { // goes through the transitions of the original model
				tmp1.clear();
				for (q = 0; q < model.getIndexes(s, t).length; q++) {
					i = model.getIndexes(s, t)[q];						
					tmp1.set((i > 0)? i : model.getIdles()[q]); // indexes of a transition in the original model
				}
				for (BitSet act : csgchoices.get(p).get(k).get(s).keySet()) {
					d = null;
					act1 = "";
					act2 = "";
					tmp2.clear();
					tmp2.or(act);
					tmp2.andNot(tmp1);
					if (tmp2.isEmpty()) {
						d = new Distribution();
						for (int u : model.getChoice(s, t).getSupport()) {
							if (!onmap.containsKey(u)) {
								m = mdp.addState();
								onmap.put(u, m);
								statelist.add(m, model.getStatesList().get(u));
								if (!explored.get(u))
									generateMDPZeroSum(mdp, onmap, statelist, explored, k, p, u); // should check for explored?
							}
							else {
								m = onmap.get(u);
							}
							d.add(m, model.getChoice(s, t).get(u) * csgchoices.get(p).get(k).get(s).get(act));
						}
						for (i = tmp1.nextSetBit(0); i >= 0; i = tmp1.nextSetBit(i + 1)) {
							if (act.get(i))
								act1 += "[" + model.getActions().get(i - 1) + "]";
							else
								act2 += "[" + model.getActions().get(i - 1) + "]";
						}
						act1 = csgchoices.get(p).get(k).get(s).get(act) + ": " + act1;
					}
					if (d != null) {
						mdp.addActionLabelledChoice(n, d, act1 + "--" + act2);
					}
				}
			}
		}
	}

	@Override
	public void restrictStrategyToReachableStates() throws PrismException {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void clear() {
		// TODO Auto-generated method stub
		
	}

}
