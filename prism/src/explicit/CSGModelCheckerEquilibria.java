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

package explicit;

import java.math.BigDecimal;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.math3.util.Precision;

import explicit.CSGModelCheckerEquilibria.CSGResultStatus;
import explicit.ProbModelChecker.TermCrit;
import explicit.rewards.CSGRewards;
import explicit.rewards.CSGRewardsSimple;
import explicit.rewards.MDPRewards;
import parser.State;
import parser.VarList;
import parser.ast.Coalition;
import parser.ast.Declaration;
import parser.ast.DeclarationBool;
import parser.ast.DeclarationInt;
import parser.ast.Expression;
import parser.ast.ExpressionTemporal;
import prism.Pair;
import prism.PrismComponent;
import prism.PrismException;
import prism.PrismFileLog;
import prism.PrismLangException;
import prism.PrismNotSupportedException;
import prism.PrismSettings;
import prism.PrismUtils;
import strat.CSGStrategy;
import strat.CSGStrategy.CSGStrategyType;
import strat.Strategy;

public class CSGModelCheckerEquilibria extends CSGModelChecker {
	
	protected MDPModelChecker mdpmc;
	
	//{player i -> action j (index) -> [<[product of ids of actions], value for joint action]>] 
	private HashMap<Integer, HashMap<Integer, ArrayList<Pair<BitSet, Double>>>> assertions;
	//[player i][action j (index)][<product of ids of actions, value for joint action>]
	private ArrayList<ArrayList<ArrayList<Pair<BitSet, Double>>>> products;
	//[id payoff player i]
	private ArrayList<Integer> payoffs;	
	
	private HashMap<Integer, HashMap<Integer, ArrayList<Pair<BitSet, Double>>>> gradient;
	//{action j -> [player i, action j (index)]}
	private HashMap<Integer, int[]> mapActionIndex;
	
	private ArrayList<ArrayList<BitSet>> supports;
	private ArrayList<BitSet> allSupports;
	
	protected BitSet[] dominated;
	protected BitSet[] dominating;
	protected BitSet players;
	
	//smt solver for labelled polytopes 
	protected CSGLabeledPolytopes eqslps;
	//smt solver for support enumeration
	protected CSGSupportEnumeration eqssmt;
	
	protected String smtSolver;
		
	public enum CSGResultStatus {
		SAT, UNKNOWN, UNSAT;
	}
	
	public CSGModelCheckerEquilibria(PrismComponent parent) throws PrismException {
		super(parent);
		players = new BitSet();
		supports = new ArrayList<ArrayList<BitSet>>();
		allSupports =  new ArrayList<BitSet>();
				
		mapActionIndex = new HashMap<Integer, int[]>();
		products = new ArrayList<ArrayList<ArrayList<Pair<BitSet, Double>>>>();
		assertions = new HashMap<Integer, HashMap<Integer, ArrayList<Pair<BitSet, Double>>>>();
		gradient = new HashMap<Integer, HashMap<Integer, ArrayList<Pair<BitSet, Double>>>>();
		payoffs = new ArrayList<Integer>();
		
		mdpmc = new MDPModelChecker(parent);
		mdpmc.setVerbosity(0);
		mdpmc.setSilentPrecomputations(true);
		
		smtSolver = getSettings().getString(PrismSettings.PRISM_SMT_SOLVER);
		switch (smtSolver) {
			case "Z3":
				break;
			case "Yices":
				break;
			default:
				throw new PrismException("Unknown SMT solver \"" + smtSolver + "\"");
		}
	}
	
	public void setSolver() throws PrismException {
		switch (smtSolver) {
			case "Z3":
				eqslps = new CSGLabeledPolytopesZ3Stack(maxRows, maxCols);
				break;
			case "Yices":
				eqslps = new CSGLabeledPolytopesYicesStack();
		}
	}

	public void buildCoalitions(CSG csg, List<Coalition> coalitions) throws PrismException {
		if (coalitions == null || coalitions.isEmpty())
			throw new PrismException("Coalitions must not be empty");
		int c, p, all;
		all = 0;
		numPlayers = csg.getNumPlayers();
		numCoalitions = coalitions.size();
		coalitionIndexes = new BitSet[coalitions.size()];
		actionIndexes = new BitSet[coalitions.size()];
		Map<Integer, String> pmap = new HashMap<Integer, String>();
		for (p = 0; p < numPlayers; p++) {
			pmap.put(p + 1, csg.getPlayerName(p));
		}
		for (c = 0; c < coalitions.size(); c++) {
			coalitionIndexes[c] = new BitSet();
			actionIndexes[c] = new BitSet();
			for (p = 0; p < numPlayers; p++) {
				if (!coalitionIndexes[c].get(p)) {					
					if (coalitions.get(c).isPlayerIndexInCoalition(p + 1, pmap)) {
						coalitionIndexes[c].set(p);
						actionIndexes[c].or(csg.getIndexes()[p]);
						actionIndexes[c].set(csg.getIdleForPlayer(p));
					}
				}
				else {
					throw new PrismLangException("Repeated player in coalition " + coalitions.get(c));
				}
			}
			all += coalitionIndexes[c].cardinality();
		}
		if (all != numPlayers)
			throw new PrismLangException("All players must be in a coalition");
		players.clear();
		players.set(0, numCoalitions);
	}
	
	public void buildAllSupports() {
		BitSet support;
		for (int p = 0; p < numCoalitions; p++) {
			if (!dominating[p].isEmpty()) {
				supports.get(p).add(dominating[p]);
			}
			else { 
				buildSupportsPlayer(new BitSet(), p, 0);
			}
		}
		for (BitSet s : supports.get(0)) {
			support = new BitSet();
			support.or(s);
			buildAllSupportsAux(support, 1);
		}
	}
	
	public void buildAllSupportsAux(BitSet supp, int p) {
		for (BitSet s : supports.get(p)) {
			BitSet curr = (BitSet) supp.clone();
			curr.or(s);
			if(p == numCoalitions - 1) {
				if (!allSupports.contains(curr))
					allSupports.add(curr);
			}
			else {
				buildAllSupportsAux(curr, p + 1);
			}
		}		
	}

	public void buildSupportsPlayer(BitSet supp, int p, int a) {
		BitSet gt0 = (BitSet) supp.clone(); 
		if (!dominated[p].get(strategies.get(p).get(a))) {
			gt0.set(strategies.get(p).get(a));
		}
		BitSet eq0 = (BitSet) supp.clone();
		if (a == strategies.get(p).size() - 1) {
			if (!eq0.isEmpty()) {
				supports.get(p).add(eq0);
			}
			if (!gt0.isEmpty()) {
				supports.get(p).add(gt0);
			}
		}
		else {
			buildSupportsPlayer(eq0, p, a + 1);
			buildSupportsPlayer(gt0, p, a + 1);
		}	
	}
	
	public BitSet findDominated(int p) throws PrismException {
		Pair<BitSet, Double> pair1, pair2;
		BitSet domi = new BitSet();
		boolean domb;
		if (assertions.get(p).keySet().size() == 1) {
			return domi;
		}
		else {
			for (int act1 : assertions.get(p).keySet()) {
				for (int act2 : assertions.get(p).keySet()) {
					if (act1 != act2) {
						domb = true;
						for (int prod = 0; prod < assertions.get(p).get(act1).size(); prod++) {
							pair1 = assertions.get(p).get(act1).get(prod);
							pair2 = assertions.get(p).get(act2).get(prod);
							if (pair1.first.equals(pair2.first)) {
								domb = domb && Double.compare(pair1.second, pair2.second) < 0;
							}
							else {
								throw new PrismException("Error when comparing indexes");
							}
						}
						if (domb) {
							domi.set(strategies.get(p).get(act1));
							break;
						}
					}
				}
			}
		}
		return domi;
	}
	
	public BitSet findDominating(int p) throws PrismException {
		Pair<BitSet, Double> pair1, pair2;
		BitSet domi = new BitSet();
		boolean domb;
		if (assertions.get(p).keySet().size() == 1) {
			return domi;
		}
		else {
			for (int act1 : assertions.get(p).keySet()) {
				domb = true;
				for (int act2 : assertions.get(p).keySet()) {
					if (act1 != act2) {
						for (int prod = 0; prod < assertions.get(p).get(act1).size(); prod++) {
							pair1 = assertions.get(p).get(act1).get(prod);
							pair2 = assertions.get(p).get(act2).get(prod);
							if (pair1.first.equals(pair2.first)) {
								domb = domb && Double.compare(pair1.second, pair2.second) > 0;
							}
							else {
								throw new PrismException("Error when comparing indexes");
							}
						}
					}
				}
				if (domb) {
					domi.set(strategies.get(p).get(act1));
				}
			}
		}
		return domi;
	}
	
	public int[] findMaxIndexes(double[][] a) {
		int result[] = new int[2];
		result[0] = 0;
		result[1] = 0;	
		double value = 0.0;
		for(int r = 0; r < a.length; r++) {
			for(int c = 0; c < a[r].length; c++) {
				if(Double.compare(a[r][c], value) > 0) {
				//if(a[r][c] > value) {
					value = a[r][c];
					result[0] = r;
					result[1] = c;	
				}
			}
		}		
		return result;
	}
	
	public void buildStepGame(CSG csg, List<CSGRewards> rewards, List<Map<Integer, BitSet>> mmap, double[][] val, int s, boolean min) throws PrismException {
		Map<BitSet, Integer> imap = new HashMap<BitSet, Integer>();
		BitSet jidx;
		BitSet indexes = new BitSet();
		BitSet tmp = new BitSet();
		String act;
		double v;
		int c, i, p, t;
		int[] joint;
		int[] idle = new int[numPlayers];
		actions.clear();
		strategies.clear();
		utilities.clear();
		varIndex = 0;
		Arrays.fill(idle, -1);
		for (c = 0; c < numCoalitions; c++) {
			actions.add(c, new ArrayList<String>());
			strategies.add(c, new ArrayList<Integer>());
		}
		for (t = 0; t < csg.getNumChoices(s); t++) {
			jidx = new BitSet();
			joint = csg.getIndexes(s, t);
			indexes.clear();
			for (p = 0; p < numPlayers; p++) {
				if (joint[p] != -1)
					indexes.set(joint[p]);
				else 
					indexes.set(csg.getIdleForPlayer(p));
			}
			for (c = 0; c < numCoalitions; c++) {
				v = 0.0;
				tmp.clear();
				tmp.or(actionIndexes[c]);
				tmp.and(indexes);
				if (tmp.cardinality() != coalitionIndexes[c].cardinality()) {
					throw new PrismException("Error in coalition");					
				}
				else {
					if(!imap.keySet().contains(tmp)) {
						act = "";
						strategies.get(c).add(varIndex);
				    	if (mmap != null) 
				    		mmap.get(c).put(strategies.get(c).size() - 1, (BitSet) tmp.clone());
						for (i = tmp.nextSetBit(0); i >= 0; i = tmp.nextSetBit(i + 1)) {
							act += "[" + csg.getActions().get(i - 1) + "]";
						}
						actions.get(c).add(act);
						jidx.set(varIndex);
						imap.put((BitSet) tmp.clone(), varIndex);
						varIndex++;
					}
					else {
						jidx.set(imap.get(tmp));
					}
				}	
			}
			utilities.put(jidx, new ArrayList<Double>());
			for (c = 0; c < numCoalitions; c++) {
				v = 0.0;
				for (int d : csg.getChoice(s, t).getSupport()) {
					v += csg.getChoice(s, t).get(d) * val[c][d];
				}
				if (rewards != null) {
					if (rewards.get(c) != null)
						v += rewards.get(c).getTransitionReward(s, t);		
				}
				v = Precision.round(v, 12, BigDecimal.ROUND_FLOOR);
				utilities.get(jidx).add(c, (min)? -1.0 * v : v); // might have to add min (v, 1.0) due to assertions for probabilistic
			}
		}	
		//System.out.println("-- imap " + imap);
		//System.out.println("\n## state " + s);
		//System.out.println("-- strategies " + strategies);
		//System.out.println("-- actions " + actions);
		//System.out.println("-- utilities " + utilities);
	}
	
	public ArrayList<ArrayList<ArrayList<Double>>> buildBiMatrixGame(CSG csg, CSGRewards r1, CSGRewards r2, List<Map<Integer, BitSet>> mmap,  List<ArrayList<Integer>> nmap, double[][] val, int s, boolean min) throws PrismException {
		ArrayList<ArrayList<ArrayList<Double>>> bmgame = new ArrayList<ArrayList<ArrayList<Double>>>();
		ArrayList<CSGRewards> rewards = null;
		BitSet action = new BitSet();
		int col, p, row, irow, icol;
		if (numCoalitions > 2) 
			throw new PrismLangException("Multiplayer game not supported by this method");
		if (r1 != null || r2 != null) {
			rewards = new ArrayList<CSGRewards>();
			rewards.add(0, r1);
			rewards.add(1, r2);
		}
		buildStepGame(csg, rewards, mmap, val, s, min);
		//System.out.println("-- utilities " + utilities);
		//System.out.println("-- strategies " + strategies);
		//System.out.println("-- mmap " + mmap);
		clear();
		computeAssertions();
		//System.out.println("-- assertions " + assertions);
		//System.out.println("-- gradient " + gradient);
		for (p = 0; p < numCoalitions; p++) {
			dominated[p] = findDominated(p);
			dominating[p] = findDominating(p);
			//System.out.println("-- dominated " + p + ": " + dominated[p]);
			//System.out.println("-- dominating " + p + ": " + dominating[p]);
		}
		buildAllSupports();
		//System.out.println("-- supports " + allSupports);
		for (p = 0; p < 2; p++) {
			bmgame.add(p, new ArrayList<ArrayList<Double>>());
			irow = 0;
			for (row = 0; row < strategies.get(0).size(); row++) {
				if (!dominated[0].get(strategies.get(0).get(row))) {
					bmgame.get(p).add(irow, new ArrayList<Double>());
					action.clear();
					action.set(strategies.get(0).get(row));
					icol = 0;
					for (col = 0; col < strategies.get(1).size(); col++) {
						if (!dominated[1].get(strategies.get(1).get(col))) {
							action.set(strategies.get(1).get(col));
							if (utilities.containsKey(action))
								bmgame.get(p).get(irow).add(icol, utilities.get(action).get(p));
							else 
								throw new PrismException("Error in building bimatrix game");
							action.clear(strategies.get(1).get(col));
							if (p == 0 && irow == 0)
								nmap.get(1).add(icol, col);
							icol++;
						}
					}
					if (p == 0)
						nmap.get(0).add(irow, row);
					irow++;
				}
			}
		}
		//System.out.println("-- nmap " + nmap);
		return bmgame;
	}
	
	public void clear() {
		supports.clear();
		allSupports.clear();
		products.clear();
		assertions.clear();
		gradient.clear();
		payoffs.clear();
		mapActionIndex.clear();
		for(int c = 0; c < numCoalitions; c++) {
			supports.add(c, new ArrayList<BitSet>());
			products.add(c, new ArrayList<ArrayList<Pair<BitSet, Double>>>());
			assertions.put(c, new HashMap<Integer, ArrayList<Pair<BitSet, Double>>>());
			gradient.put(c, new HashMap<Integer, ArrayList<Pair<BitSet, Double>>>());
			for(int a = 0; a < strategies.get(c).size(); a++) {
				products.get(c).add(a, new ArrayList<Pair<BitSet, Double>>());
				mapActionIndex.put(strategies.get(c).get(a), new int[] {c, a});
			}
		}
		for(int c = 0; c < numCoalitions; c++) { // strategies allocated first
			payoffs.add(c, varIndex);
			varIndex++;
		}
	}
	
	public void computeAssertions() throws PrismException {
		int c, q;
		BitSet ps;
		BitSet acts = new BitSet();
		for (c = 0; c < numCoalitions; c++) {
			//System.out.println("\n-- player " + c + " ");
        	ps = (BitSet) players.clone();
        	ps.clear(c);
        	for (q = 0; q < strategies.get(c).size(); q++) {
            	//System.out.println("-- action " + q);
        		acts.clear();
        		acts.set(strategies.get(c).get(q));
        		assertions.get(c).put(q, prodAction(acts, ps, q, c));
        	}
        }
		/*
		System.out.println("-- actions ");
		System.out.println(actions);
		System.out.println("-- strategies ");
		System.out.println(strategies);		
		
		
		System.out.println("-- assertions ");
		for (c = 0; c < numCoalitions; c++) {
			System.out.println("--- player " + c);
        	for (q = 0; q < strategies.get(c).size(); q++)
        		System.out.println("---- action " + q + " " + assertions.get(c).get(q));
		}
		
		/*
		//System.out.println("-- gradient ");
		//System.out.println(gradient);
		System.out.println("-- map ");
		for (int i : map.keySet()) {
			System.out.println(i + "= " + Arrays.toString(map.get(i)));
		}
		*/
	}
	
	public ArrayList<Pair<BitSet, Double>> prodAction(BitSet prod, BitSet sp, int act, int p) throws PrismException {
		prodAction(new Pair<BitSet, Double>(new BitSet(), 0.0), prod, sp, act, p);
		ArrayList<Pair<BitSet, Double>> sum = new ArrayList<Pair<BitSet, Double>>();		
		for(int j = 0; j < products.get(p).get(act).size(); j++) {
			sum.add(products.get(p).get(act).get(j));
		}
		products.get(p).get(act).clear();
		return sum;
	}
	
	public void prodAction(Pair<BitSet, Double> expr, BitSet prod, BitSet sp, int act, int p) throws PrismException {
		BitSet set;
		BitSet curr = (BitSet) sp.clone();
		Pair<BitSet, Double> nexpr;
		if(products.get(p) == null) 
			products.add(p, new ArrayList<ArrayList<Pair<BitSet, Double>>>());
		else if(products.get(p).get(act) == null)
			products.get(p).add(act, new ArrayList<Pair<BitSet, Double>>());
		for(int cp = sp.nextSetBit(0); cp < sp.size() && cp != -1; cp = sp.nextSetBit(cp + 1)) {
			curr.clear(cp);
			for(int a = 0; a < actions.get(cp).size(); a++) {
				set = new BitSet();
				set.or(prod);
				set.set(strategies.get(cp).get(a));		
				nexpr = new Pair<BitSet, Double>(new BitSet(), 0.0);
				nexpr.first.or(expr.first);
				nexpr.first.set(strategies.get(cp).get(a));
				prodAction(nexpr, set, curr, act, p);
				if(sp.cardinality() == 1 && set.cardinality() == numCoalitions) { // should have to check for set size?	
					nexpr.first.or(expr.first);
					nexpr.first.set(strategies.get(cp).get(a));
					nexpr.second = utilities.get(set).get(p);			
					products.get(p).get(act).add(nexpr);
					if(gradient.get(p).get(act) == null)
						gradient.get(p).put(act, new ArrayList<Pair<BitSet, Double>>());	
					if(nexpr.second != 0.0)
						gradient.get(p).get(act).add(nexpr);
					BitSet der = new BitSet();
					for(int i = nexpr.first.nextSetBit(0); i >= 0; i = nexpr.first.nextSetBit(i + 1)) {
						der.or(nexpr.first);
						der.set(strategies.get(p).get(act));
						der.clear(i);
						if(gradient.get(mapActionIndex.get(i)[0]).get(mapActionIndex.get(i)[1]) == null)
							gradient.get(mapActionIndex.get(i)[0]).put(mapActionIndex.get(i)[1], new ArrayList<Pair<BitSet, Double>>());
						if(nexpr.second != 0.0)
							gradient.get(mapActionIndex.get(i)[0]).get(mapActionIndex.get(i)[1]).add(new Pair<BitSet, Double>((BitSet) der.clone(), nexpr.second));
						der.clear();
					}
				}
			}
		}
	}
	
	public double[][] computeBoundedReachProbs(CSG csg, BitSet target, int n) {
		double[][] sol = new double[n][csg.getNumStates()];
		double[] sol1 = new double[csg.getNumStates()];
		double[] sol2 = new double[csg.getNumStates()];
		double v, sum;
		int i, s, t;	
		for (s = 0; s < csg.getNumStates(); s++) {
			sol2[s] = sol1[s] = target.get(s)? 1.0 : 0.0;
		}
		for (i = 0; i < n; i++) {
			for (s = 0; s < csg.getNumStates(); s++) {
				v = 0.0;
				for (t = 0; t < csg.getNumChoices(s); t++) {
					sum = 0.0;
					for (int u : csg.getChoice(s, t).keySet()) {
						sum += csg.getChoice(s, t).get(u) * sol2[u];
					}
					v = (sum > v)? sum : v;
				}
				if (!target.get(s))
					sol1[s] = v;
			}
			sol2 =  Arrays.copyOf(sol1, sol1.length);
			sol[i] = sol1;
		}
		return sol;
	} 
	
	public ModelCheckerResult computeBoundedEquilibria(CSG csg, List<Coalition> coalitions, List<CSGRewards> rewards, List<ExpressionTemporal> exprs, BitSet[] targets, BitSet[] remain, int[] bounds, boolean min) throws PrismException {
		if (generateStrategy || exportAdv) {
			throw new PrismException("Strategy synthesis for bounded properties is not supported yet.");
		}
		ModelCheckerResult res = new ModelCheckerResult();
		List<CSGRewards> newRewards= null;
		BitSet[] only = new BitSet[coalitions.size()];
		BitSet[] phi1 = new BitSet[3];
		BitSet cpy =  new BitSet();
		double[][] sol = new double[coalitions.size()][csg.getNumStates()];
		double[][] tmp = new double[coalitions.size()][csg.getNumStates()];
		double[][] val = new double[coalitions.size()][csg.getNumStates()];
		double[][] eq;
		double[] sw;
		double[] r = new double[csg.getNumStates()];
		int i, j, n1, n2, k, s;
		boolean rew;
		long currentTime, timePrecomp;		
		
		rew = rewards != null;
		
		buildCoalitions(csg, coalitions);
		findMaxRowsCols(csg);
		setSolver();
		mainLog.println("Starting bounded equilibria computation (solver=" + eqslps.getSolverName() + ")...");
		dominated = new BitSet[numCoalitions];
		dominating = new BitSet[numCoalitions];
		
		// Case next
		if ((exprs.get(0).getOperator() == ExpressionTemporal.P_X) || (exprs.get(1).getOperator() == ExpressionTemporal.P_X)) {
			for (i = 0; i < 2; i++) {
				if (exprs.get(i).getOperator() == ExpressionTemporal.P_X) {
					for (s = 0; s < csg.getNumStates(); s++) {
						sol[i][s] = targets[i].get(s)? 1.0 : 0.0;
					}
				}
				else {
					sol[i] = mdpmc.computeBoundedUntilProbs(csg, remain[i], targets[i], bounds[i]-1, min).soln;
				}
			}
			for (s = 0; s < csg.getNumStates(); s++) {
				eq = stepEquilibria(csg, null, null, null, eqslps, null, sol, s, min);
				sw = swne(eq, null, min);
				tmp[0][s] = sw[1];
				tmp[1][s] = sw[2];
				r[s] = sw[1] + sw[2];
			}
			mainLog.println("\nCoalition results (initial state): (" + tmp[0][csg.getFirstInitialState()] + "," + tmp[1][csg.getFirstInitialState()] + ")");
			res.soln = r;
			res.numIters = 1;
			return res;		
		}
		
		if (!rew) {
			for (i = 0; i < 2; i++) {
				phi1[i] = new BitSet();
				if (remain[i] == null) 
					phi1[i].set(0, csg.getNumStates());
			else
				phi1[i].or(remain[i]);
			}
		}
		
		if (targets == null) {
			targets = new BitSet[coalitions.size()];
			for (i = 0; i < coalitions.size(); i++) { // Case for cumulative rewards
				targets[i] = new BitSet();
			}
		}
		
		for (i = 0; i < coalitions.size(); i++) {
			only[i] = new BitSet();
			only[i].or(targets[i]);
			for (j = 0; j < coalitions.size(); j++) {
				if (i != j)
					only[i].andNot(targets[j]);
			}
		}
		
		k = Math.abs(bounds[0] - bounds[1]);
		n1 = (bounds[0] > bounds[1])? k : 0;
		n2 = (bounds[1] > bounds[0])? k : 0;
		
		if (!rew) {
			phi1[2] = new BitSet();
			phi1[2].or(phi1[0]); 
			phi1[2].and(phi1[1]); // intersection of phi1(1) and phi1(2)
			cpy.clear();
			cpy.or(phi1[0]);
			phi1[0].andNot(phi1[1]); // phi1(1) minus phi1(2)
			phi1[1].andNot(cpy); // phi1(2) minus phi1(1)
			cpy.clear();
		}
		else {
			newRewards = new ArrayList<CSGRewards>();
		}
		
		//System.out.println("for bounds[0]");
		//double[][] pre0 = computeBoundedReachProbs(csg, targets[0], bounds[0]); 
		//System.out.println("for bounds[1]");
		//double[][] pre1 = computeBoundedReachProbs(csg, targets[1], bounds[1]); 
		
		timePrecomp = System.currentTimeMillis();
		if (rew) {
			if (bounds[0] > bounds[1]) {
				if (exprs.get(0).getOperator() == ExpressionTemporal.R_C)
					val[0] = mdpmc.computeCumulativeRewards(csg, rewards.get(0), n1, min).soln;
				else
					val[0] = mdpmc.computeInstantaneousRewards(csg, rewards.get(0), n1, min).soln;
			}
			if (bounds[1] > bounds[0]) {
				if (exprs.get(1).getOperator() == ExpressionTemporal.R_C)
					val[1] = mdpmc.computeCumulativeRewards(csg, rewards.get(1), n2, min).soln;
				else
					val[1] = mdpmc.computeInstantaneousRewards(csg, rewards.get(1), n2, min).soln;
			}
		}	
		timePrecomp = System.currentTimeMillis() - timePrecomp;
		
		while (true) {
			currentTime = System.currentTimeMillis();
			if (!rew) {
				if (n1 > 0) {
					if (remain[0] == null) 
						val[0] = mdpmc.computeBoundedReachProbs(csg, targets[0], n1, min).soln;
					else
						val[0] = mdpmc.computeBoundedUntilProbs(csg, remain[0], targets[0], n1, min).soln;
				}
				if (n2 > 0) {
					if (remain[1] == null) 
						val[1] = mdpmc.computeBoundedReachProbs(csg, targets[1], n2, min).soln;
					else 
						val[1] = mdpmc.computeBoundedUntilProbs(csg, remain[1], targets[1], n2, min).soln;
				}
			}
			timePrecomp += System.currentTimeMillis() - currentTime;
			if (Math.min(n1, n2) > 0) {
				for (s = 0; s < csg.getNumStates(); s++) {
					if (rew) {
						newRewards.clear();
						for (i = 0; i < 2; i++) {
							newRewards.add(i, rewards.get(i));
							if (!(exprs.get(i).getOperator() == ExpressionTemporal.R_C))
								newRewards.set(i, null);
						}
						eq = stepEquilibria(csg, newRewards.get(0), newRewards.get(1), null, eqslps, null, sol, s, min);	
						sw = swne(eq, null, min);
						tmp[0][s] = sw[1];
						tmp[1][s] = sw[2];
					} 
					else {
						if (targets[0].get(s) && targets[1].get(s)) {
							tmp[0][s] = 1.0;
							tmp[1][s] = 1.0;					
						}
						else if (only[0].get(s)) {
							tmp[0][s] = 1.0;
							tmp[1][s] = val[1][s];
						}
						else if (only[1].get(s)) {
							tmp[0][s] = val[0][s];
							tmp[1][s] = 1.0;		
						}
						else if(phi1[0].get(s)) {
							tmp[0][s] = val[0][s];
							tmp[1][s] = 0.0;	
						}
						else if(phi1[1].get(s)) {
							tmp[0][s] = 0.0;
							tmp[1][s] = val[1][s];
						}
						else if(!phi1[2].get(s)) {
							tmp[0][s] = 0.0;
							tmp[1][s] = 0.0;
						}
						else {
							eq = stepEquilibria(csg, null, null, null, eqslps, null, sol, s, min);
							sw = swne(eq, null, min);
							tmp[0][s] = sw[1];
							tmp[1][s] = sw[2];
						}
					}
				}
				for (s = 0; s < csg.getNumStates(); s++) {
					sol[0][s] = tmp[0][s];
					sol[1][s] = tmp[1][s];
					r[s] = sol[0][s] + sol[1][s];
				}
				/*
				String sols;
				sols = "(";
				for (p = 0; p < numCoalitions; p++) {
					if (p < numCoalitions - 1)
						sols += sol[p][csg.getFirstInitialState()] + ",";
					else
						sols += sol[p][csg.getFirstInitialState()] + ")";
				}
				mainLog.println(k + ": " + sols);
				*/
			}
			else {
				for (s = 0; s < csg.getNumStates(); s++) {
					if (rew) {
						if (n1 == 0 && n2 == 0) {
							sol[0][s] = (exprs.get(0).getOperator() == ExpressionTemporal.R_C)? 0.0 : rewards.get(0).getStateReward(s);
							sol[1][s] = (exprs.get(1).getOperator() == ExpressionTemporal.R_C)? 0.0 : rewards.get(1).getStateReward(s);
						}
						else if (n1 == 0) {
							sol[0][s] = (exprs.get(0).getOperator() == ExpressionTemporal.R_C)? 0.0 : rewards.get(0).getStateReward(s);
							sol[1][s] = val[1][s];
						}
						else {
							sol[0][s] = val[0][s];
							sol[1][s] = (exprs.get(1).getOperator() == ExpressionTemporal.R_C)? 0.0 : rewards.get(1).getStateReward(s);
						}
					}
					else {
						if (n1 == 0 && n2 == 0) {
							sol[0][s] = targets[0].get(s)? 1.0 : 0.0;
							sol[1][s] = targets[1].get(s)? 1.0 : 0.0;
						}
						else if (n1 == 0) {
							sol[0][s] = targets[0].get(s)? 1.0 : 0.0;
							sol[1][s] = targets[1].get(s)? 1.0 : val[1][s];
						}
						else {
							sol[0][s] = targets[0].get(s)? 1.0 : val[0][s];
							sol[1][s] = targets[1].get(s)? 1.0 : 0.0;
						}
					}
				}
			}
			if (k == Math.max(bounds[0], bounds[1])) {
				break;
			}
			k++;
			n1 = Math.min(n1 + 1, bounds[0]);
			n2 = Math.min(n2 + 1, bounds[1]);
		}
		mainLog.println("\nPrecomputation took " + timePrecomp / 1000.0 + " seconds.");
		mainLog.println("Coalition results (initial state): (" + sol[0][csg.getFirstInitialState()] + "," + sol[1][csg.getFirstInitialState()] + ")");
		res.soln = r;
		res.numIters = k;
		return res;		
	}

	public ModelCheckerResult computeMultiReachEquilibria(CSG csg, List<Coalition> coalitions, List<CSGRewards> rewards, BitSet[] targets, BitSet[] remain, boolean min) throws PrismException {
		ModelCheckerResult res = new ModelCheckerResult();
		Declaration[] goals = new Declaration[coalitions.size()];
		//Map<BitSet, BitSet> common = new HashMap<BitSet, BitSet>();
		List<List<List<Map<BitSet, Double>>>> lstrat = null;
		List<List<Map<BitSet, Double>>> sstrat = null;
		List<Map<Integer, BitSet>> mmap = null;
		Set<BitSet> games = new HashSet<BitSet>();
		BitSet unsat = new BitSet();
		BitSet ps;
		double[][] sol;
		double[][] val;
		double[][] tmp;
		double[][] eq;
		double[] r;
		double[] sw;
		int i, j, k, n, p, s, t, v;
		boolean done, rew;
		
		rew = rewards != null;
		
		buildCoalitions(csg, coalitions);
		dominated = new BitSet[numCoalitions];
		dominating = new BitSet[numCoalitions];
		
		VarList newvarlist = (VarList) csg.getVarList().clone();
		
		for (p = 0; p < coalitions.size(); p++) {
			ps = (BitSet) players.clone();
			ps.clear(p);
			goals[p] = new Declaration("_" + p, new DeclarationBool());
			buildSubGames(games, ps, p);
			newvarlist.addVar(goals[p].getName(), goals[p].getDeclType(), 0, csg.getConstantValues());
		}
		games.add((BitSet) players.clone());
		unsat.set(0, csg.getNumStates());		
		for (BitSet subgame : games) {
			//BitSet intersec = new BitSet();
			//intersec.set(0, csg.getNumStates()); // intersec only seems to be used to define common, which later does not seem to be used
			for (p = 0; p < coalitions.size(); p++) {
				unsat.andNot(targets[p]);
			//	if (subgame.get(p))
			//		intersec.and(targets[p]); 
			//	else
			//		intersec.andNot(targets[p]);
			}
			//common.put(subgame, intersec); // common does not seem to be used afterwards
		} 
		
		//System.out.println("-- common " + common);
		
		CSG newmodel = new CSG();
		newmodel.setActions(csg.getActions());
		newmodel.setPlayers(csg.getPlayers());
		newmodel.setIndexes(csg.getIndexes());
		newmodel.setIdles(csg.getIdles());
		
		List<State> newstatelist = new ArrayList<State>();
		BitSet explored = new BitSet();
		BitSet isubgame = new BitSet();
		
		s = csg.getFirstInitialState();
		
		State initial = new State(csg.getVarList().getNumVars() + numCoalitions);
		for (v = 0; v < csg.getVarList().getNumVars(); v++) {
			initial.setValue(v, csg.getStatesList().get(s).varValues[v]);
		}
		
		Map<Integer, Integer> extmap = new HashMap<Integer, Integer>();
				
		if (unsat.get(s)) {
			for (p = 0; p < numCoalitions; p++) {
				initial.setValue(newvarlist.getIndex(goals[p].getName()), false);
			}	
		}
		else {	
			for (p = 0; p < numCoalitions; p++) {
				if (targets[p].get(s)) {
					initial.setValue(newvarlist.getIndex(goals[p].getName()), true);
					isubgame.set(p);
				}
				else {
					initial.setValue(newvarlist.getIndex(goals[p].getName()), false);
				}
			}
		}
		
		//System.out.println("-- isubgame " + isubgame);
		
		mainLog.println("\nBuilding extended model...");

		Map<Integer, Integer> nomap = new HashMap<Integer, Integer>();	
		Map<BitSet, BitSet> subgames = new HashMap<BitSet, BitSet>();
		
		n = newmodel.addState();
		newstatelist.add(n , initial);
		newmodel.addInitialState(n);
		
		newmodel.setVarList(newvarlist);
		newmodel.setStatesList(newstatelist);

		nomap.put(n, s);
		extmap.put(initial.hashCode(), n);
		explored.set(n);
		exploreSat(csg, newmodel, nomap, extmap, subgames, newstatelist, newvarlist, goals, targets, explored, isubgame, n);
		
		//System.out.println("-- nomap " + nomap);
		
		/*
		for (BitSet subset : common.keySet()) {
			for (s = common.get(subset).nextSetBit(0); s >= 0; s = common.get(subset).nextSetBit(s + 1)) {
				subgames.put(subset, new BitSet());
				subgames.get(subset).set(onmap.get(s));
				exploreSat(csg, newmodel, subgames, onmap, newstatelist, newvarlist, goals, targets, subset, explored, onmap.get(s));
			}
		}
		*/
		
		System.out.println("-- subgames " + subgames);

		Path currentRelativePath = Paths.get("");
		String path = currentRelativePath.toAbsolutePath().toString();		
		PrismFileLog pflog = new PrismFileLog(path + "/newmodel.dot");
		newmodel.exportToDotFile(pflog, null, true);
	
		pflog = new PrismFileLog(path + "/model.dot");
		csg.exportToDotFile(pflog, null, true);
		
		mainLog.print("\n" + newmodel.infoStringTable() + "\n");
				
		sol = new double[coalitions.size()][newmodel.getNumStates()];
		val = new double[coalitions.size()][newmodel.getNumStates()];
		tmp = new double[coalitions.size()][newmodel.getNumStates()];
		r = new double[newmodel.getNumStates()];
				
		List<CSGRewards> newrewards = null;
		if (rew) {
			newrewards = new ArrayList<CSGRewards>();
			for (p = 0; p < numCoalitions; p++) {
				newrewards.add(p, new CSGRewardsSimple(newmodel.getNumStates()));
			}
			for (s = 0;  s < newmodel.getNumStates(); s++) {
				for (p = 0; p < numCoalitions; p++) {
					((CSGRewardsSimple) newrewards.get(p)).setStateReward(s, rewards.get(p).getStateReward(nomap.get(s)));
					for (t = 0; t < newmodel.getNumChoices(s); t++) {
						((CSGRewardsSimple) newrewards.get(p)).setTransitionReward(s, t, rewards.get(p).getTransitionReward(nomap.get(s), t));
					}
				}
			}
			for (BitSet subgame : subgames.keySet()) {
				for (p = 0; p < numCoalitions; p++) {
					if (subgame.get(p)) {
						for (s = subgames.get(subgame).nextSetBit(0); s >= 0; s = subgames.get(subgame).nextSetBit(s + 1)) {
							((CSGRewardsSimple) newrewards.get(p)).setStateReward(s, 0.0);
							for (t = 0; t < newmodel.getNumChoices(s); t++) {
								((CSGRewardsSimple) newrewards.get(p)).setTransitionReward(s, t, 0.0);
							}						
						}
					}
				}
			}
		}
		else {
			for (BitSet subgame : subgames.keySet()) {
				for (p = 0; p < numCoalitions; p++) {
					if (subgame.get(p)) {
						for (s = subgames.get(subgame).nextSetBit(0); s >= 0; s = subgames.get(subgame).nextSetBit(s + 1)) {
							sol[p][s] = 1.0;
						}
					}
				}
			}
		}
		/*
		for (s = 0; s < newmodel.getNumStates(); s++) {
			System.out.println("## old " + nomap.get(s) + " " + csg.getStatesList().get(nomap.get(s)));
			System.out.println("## new " + s + " " + newmodel.getStatesList().get(s));
			for (p = 0; p < numCoalitions; p++) {
				System.out.println("-- player " + p);
				System.out.println("-- state rewards old " + rewards.get(p).getStateReward(nomap.get(s)));
				System.out.println("-- state rewards new " + newrewards.get(p).getStateReward(s));
				for (t = 0; t < newmodel.getNumChoices(s); t++) {
					System.out.println("-- trans rewards old " + csg.getAction(nomap.get(s), t) + " " + rewards.get(p).getTransitionReward(nomap.get(s), t));
					System.out.println("-- trans rewards new " + newmodel.getAction(s, t) + " " + newrewards.get(p).getTransitionReward(s, t));
				}
			}
		}
		
		System.exit(1);
		*/
		if (generateStrategy || exportAdv) {
			mmap = new ArrayList<Map<Integer, BitSet>>();
			sstrat = new ArrayList<List<Map<BitSet, Double>>>();
			lstrat = new ArrayList<List<List<Map<BitSet, Double>>>>();
			for (i = 0; i < coalitions.size(); i++) {
        		mmap.add(i, new HashMap<Integer, BitSet>());
				lstrat.add(i, new ArrayList<List<Map<BitSet, Double>>>());
				lstrat.get(i).add(0, new ArrayList<Map<BitSet, Double>>());
				for (j = 0; j < newmodel.getNumStates(); j++) {	
					lstrat.get(i).get(0).add(j, null);
				}
			} 
		}
				
		for (p = 0; p < numCoalitions; p++) {
			Arrays.fill(tmp[p], 0.0);
			Arrays.fill(val[p], 0.0);
		}
		String sols;
		done = true;
		k = 0;
		while (true) {
			//System.out.print(k + ": ");
			for (s = 0; s < newmodel.getNumStates(); s++) {
				//System.out.println("\n## s " + s + " " + newmodel.getStatesList().get(s));		
				if (generateStrategy || exportAdv) {
					sstrat = new ArrayList<List<Map<BitSet, Double>>>();
					for (i = 0; i < coalitions.size(); i++) {
		        		mmap.get(i).clear();
		        	}
				}
				if (rew) {
					eq = stepEquilibria(newmodel, newrewards, mmap, sstrat, sol, s, min);
					addStateRewards(eq, newrewards, s, min);
				}
				else {
					eq = stepEquilibria(newmodel, null, mmap, sstrat, sol, s, min);
				}
				sw = swne(eq, sstrat, min);
				for (p = 0; p < numCoalitions; p++) {
					val[p][s] = sw[p + 1];
				}
				if (generateStrategy || exportAdv) {			
					for (p = 0; p < coalitions.size(); p++) {
						if (lstrat.get(p).get(0).get(s) == null) {
							lstrat.get(p).get(0).set(s, sstrat.get(0).get(p));
						}
						else if (!lstrat.get(0).get(0).get(s).equals(sstrat.get(0).get(p)) && checkEquilibriumChange(sol, sw, s)) { // will have to change for min
							lstrat.get(p).get(0).set(s, sstrat.get(0).get(p));
						}
					}		
				}
			}
			for (s = 0; s < newmodel.getNumStates(); s++) {
				for (p = 0; p < numCoalitions; p++) {
					sol[p][s] = val[p][s];
				}
				r[s] = 0.0;
				for (p = 0; p < numCoalitions; p++) {
					r[s] += sol[p][s];
				}
			}
			//System.out.println();
			sols = "(";
			for (p = 0; p < numCoalitions; p++) {
				if (p < numCoalitions - 1)
					sols += sol[p][newmodel.getFirstInitialState()] + ",";
				else
					sols += sol[p][newmodel.getFirstInitialState()] + ")";
				//System.out.println(p + " " + sol[p][newmodel.getFirstInitialState()]);
				done = done & PrismUtils.doublesAreClose(sol[p], tmp[p], 1e-6, true);
			}
			System.out.println(k + ": " + sols);
			k++;
			if (done) {
				break;
			}
			else if (!done && k == maxIters) {
				throw new PrismException("Could not converge after " + k + " iterations");
			}
			else {
				done = true;
				for (p = 0; p < numCoalitions; p++) {
					//System.out.println(Arrays.toString(sol[p]));
					tmp[p] = Arrays.copyOf(sol[p], sol[p].length);
				}
			}
		}
		mainLog.println();
		for (p = 0; p < numCoalitions; p++) {
			mainLog.println("Result for coalition " + coalitions.get(p) + ": " + sol[p][newmodel.getFirstInitialState()] + " (value in the intial state).");
		}
		r = new double[csg.getNumStates()];
		for (s = 0; s < newmodel.getNumStates(); s++) {
			r[nomap.get(s)] = 0.0;
			for (p = 0; p < numCoalitions; p++) {
				r[nomap.get(s)] += sol[p][s];
			}
		}
		if (generateStrategy || exportAdv)
			res.strat = new CSGStrategy(newmodel, lstrat, subgames, numCoalitions, CSGStrategyType.EQUILIBRIA_M);
		res.soln = r;
		res.numIters = k;
		return res;
	}
	
	public boolean checkEquilibriumChange(double[][] sol, double[] eq, int s) {
		int p;
		boolean result = true;
		for (p = 0; p < numCoalitions; p++) {
			result = result && Double.compare(sol[p][s], eq[p + 1]) == 0;
			if (!result)
				return true;
		}
		return false;
	}
	
	public void exploreSat(CSG csg, CSG newmodel, Map<Integer, Integer> nomap, Map<Integer, Integer> extmap, Map<BitSet, BitSet> subgames, List<State> newstatelist, 
									VarList newvarlist, Declaration[] goals, BitSet[] targets, BitSet explored, BitSet subgame, int s) {
		Distribution d;
		int l, m, p, t, v;
		explored.set(s);
		if (!subgame.isEmpty()) { 
			if (!subgames.containsKey(subgame))
				subgames.put(subgame, new BitSet());
			subgames.get(subgame).set(s);
		}
		for (t = 0; t < csg.getNumChoices(nomap.get(s)); t++) {
			d = new Distribution();
			for (int u : csg.getChoice(nomap.get(s), t).getSupport()) {
				State ustate = new State(newstatelist.get(s));
				BitSet usubgame = new BitSet();
				usubgame.or(subgame);
				for (v = 0; v < csg.getVarList().getNumVars(); v++) {
					ustate.setValue(v, csg.getStatesList().get(u).varValues[v]);
				}
				for (p = 0; p < numCoalitions; p++) {
					if (targets[p].get(u)) {
						ustate.setValue(newvarlist.getIndex(goals[p].getName()), true);
						usubgame.set(p);
					}
				}
				if (!extmap.containsKey(ustate.hashCode())) {
					m = newmodel.addState();
					nomap.put(m, u);
					extmap.put(ustate.hashCode(), m);
					newstatelist.add(m, ustate);
				}
				else {
					m = extmap.get(ustate.hashCode());
				}
				if (!explored.get(m)) 
					exploreSat(csg, newmodel, nomap, extmap, subgames, newstatelist, newvarlist, goals, targets, explored, usubgame, m);
				d.add(m, csg.getChoice(nomap.get(s), t).get(u));
			}	
			l = newmodel.addActionLabelledChoice(s, d, csg.getAction(nomap.get(s), t));
			newmodel.setIndexes(s, l, csg.getIndexes(nomap.get(s), t));
		}
	}
	
	public void buildSubGames(Set<BitSet> games, BitSet sp, int p) {
		BitSet prod = new BitSet();
		prod.set(p);
		games.add((BitSet) prod.clone());
		for(int cp = sp.nextSetBit(0); cp >= 0; cp = sp.nextSetBit(cp + 1)) {
			BitSet newprod = new BitSet();
			newprod.or(prod);
			newprod.set(cp);
			games.add(newprod);
		}
	}
	
	public ModelCheckerResult computeReachEquilibria(CSG csg, List<Coalition> coalitions, List<CSGRewards> rewards, BitSet[] targets, BitSet[] remain, boolean min) throws PrismException {
		ModelCheckerResult[] obj = new ModelCheckerResult[coalitions.size()];
		ModelCheckerResult res = new ModelCheckerResult();
		List<List<List<Map<BitSet, Double>>>> lstrat = null;
		List<List<Map<BitSet, Double>>> sstrat = null;
		List<Map<Integer, BitSet>> mmap = null;
		BitSet[] only = new BitSet[targets.length];
		BitSet[] phi1 = new BitSet[3];
		BitSet known = new BitSet();
		BitSet temp =  new BitSet();
		double[][] sol = new double[coalitions.size()][csg.getNumStates()];
		double[][] val = new double[coalitions.size()][csg.getNumStates()];
		double[][] tmp = new double[coalitions.size()][csg.getNumStates()];
		double[][] eq;
		double[] r = new double[csg.getNumStates()];
		double[] sw;
		int i, j, k, p, s;
		boolean done, rew;
		long timePrecomp;
				
		// player -> iteration -> state -> indexes -> value
		if (generateStrategy || exportAdv) {
			//mdpmc.setGenerateStrategy(true);
			mdpmc.setGenStrat(true);
			mmap = new ArrayList<Map<Integer, BitSet>>();
			sstrat = new ArrayList<List<Map<BitSet, Double>>>();
			lstrat = new ArrayList<List<List<Map<BitSet, Double>>>>();
			for (i = 0; i < coalitions.size(); i++) {
        		mmap.add(i, new HashMap<Integer, BitSet>());
				lstrat.add(i, new ArrayList<List<Map<BitSet, Double>>>());
				lstrat.get(i).add(0, new ArrayList<Map<BitSet, Double>>());
				for (j = 0; j < csg.getNumStates(); j++) {	
					lstrat.get(i).get(0).add(j, null);
				}
			} 
		}
		rew = rewards != null;
		for (i = 0; i < targets.length; i++) {
			only[i] = new BitSet();
			only[i].or(targets[i]);
			for (j = 0; j < targets.length; j++) {
				if (i != j)
					only[i].andNot(targets[j]);
			}
			known.or(targets[i]);
		}		
		if (!rew) {
			for (i = 0; i < 2; i++) {
				phi1[i] = new BitSet();
				if (remain[i] == null) 
					phi1[i].set(0, csg.getNumStates());
				else
					phi1[i].or(remain[i]);
			}
			phi1[2] = new BitSet();
			phi1[2].or(phi1[0]); 
			phi1[2].and(phi1[1]); // intersection of phi1(1) and phi1(2)
			temp.clear();
			temp.or(phi1[0]);
			phi1[0].andNot(phi1[1]); // phi1(1) minus phi1(2)
			phi1[1].andNot(temp); // phi1(2) minus phi1(1)
			known.or(phi1[0]);
			known.or(phi1[1]);
			temp.clear();
			temp.set(0, csg.getNumStates());
			temp.andNot(phi1[2]);
			known.or(temp);
		}		
		buildCoalitions(csg, coalitions);
		dominated = new BitSet[numCoalitions];
		dominating = new BitSet[numCoalitions];
		mainLog.println();
		findMaxRowsCols(csg);
		setSolver();
		mainLog.println("Starting equilibria computation (solver=" + eqslps.getSolverName() + ")...");
		
		// Add an exception somewhere for reachability rewards and states with positive rewards that do not reach the targets
		k = 0;
		if (rew) {
			// Precompuation for rewards
			timePrecomp = System.currentTimeMillis();
			for (i = 0; i < targets.length; i++) {
				obj[i] = mdpmc.computeReachRewards((MDP) csg, (MDPRewards) rewards.get(i), targets[i], min);
				val[i] = obj[i].soln;
			}
			timePrecomp = System.currentTimeMillis() - timePrecomp;
			for (s = 0; s < csg.getNumStates(); s++) {
				if (targets[0].get(s) && targets[1].get(s)) {
					sol[0][s] = 0.0;
					sol[1][s] = 0.0;
				}
				else if (only[0].get(s)) {
					sol[0][s] = 0.0;
					sol[1][s] = val[1][s];
				}
				else if (only[1].get(s)) {
					sol[0][s] = val[0][s];
					sol[1][s] = 0.0;
				}
			}
		}
		else {
			// Precomputation for probabilistic
			timePrecomp = System.currentTimeMillis();
			for (i = 0; i < targets.length; i++) {		
				if (remain[i] != null)
					obj[i] = mdpmc.computeUntilProbs(csg, remain[i], targets[i], min);
				else 
					obj[i] = mdpmc.computeReachProbs((MDP) csg, targets[i], min);
				val[i] = obj[i].soln;
			}
			timePrecomp = System.currentTimeMillis() - timePrecomp;
			for (s = 0; s < csg.getNumStates(); s++) {
				if (targets[0].get(s) && targets[1].get(s)) {
					sol[0][s] = 1.0;
					sol[1][s] = 1.0;
				}
				else if (only[0].get(s)) {
					sol[0][s] = 1.0;
					sol[1][s] = val[1][s];
				}
				else if (only[1].get(s)) {
					sol[0][s] = val[0][s];
					sol[1][s] = 1.0;
				}
				else if (phi1[0].get(s)) {
					sol[0][s] = val[0][s];
					sol[1][s] = 0.0;
				}
				else if (phi1[1].get(s)) {
					sol[0][s] = 0.0;
					sol[1][s] = val[1][s];
				}
				else if (!phi1[2].get(s)) {
					sol[0][s] = 0.0;
					sol[1][s] = 0.0;
				}
			}	
		}
		mainLog.println();
		done = true;
		dominated = new BitSet[numCoalitions];
		dominating = new BitSet[numCoalitions];
		while (true) {
			for (s = 0; s < csg.getNumStates(); s++) {
				if (!known.get(s)) {
					if (generateStrategy || exportAdv) {
						sstrat = new ArrayList<List<Map<BitSet, Double>>>();
						mmap.clear();
			    		for (p = 0; p < 2; p++) {
			        		mmap.add(p, new HashMap<Integer, BitSet>());
			        	}
					}
					if (rew) {
						eq = stepEquilibria(csg, rewards.get(0), rewards.get(1), mmap, eqslps, sstrat, sol, s, min);
					}
					else {
						eq = stepEquilibria(csg, null, null, mmap, eqslps, sstrat, sol, s, min);
					}
					sw = swne(eq, sstrat, min);
					val[0][s] = sw[1];
					val[1][s] = sw[2];
					// player -> iteration -> state -> indexes -> value
					if (generateStrategy || exportAdv) {			
						for (p = 0; p < coalitions.size(); p++) {
							if (lstrat.get(p).get(0).get(s) == null) {
								lstrat.get(p).get(0).set(s, sstrat.get(0).get(p));
							}
							else if (!lstrat.get(0).get(0).get(s).equals(sstrat.get(0).get(p)) && checkEquilibriumChange(sol, sw, s)) {
								lstrat.get(p).get(0).set(s, sstrat.get(0).get(p));
							}
						}		
					}					
				}
				// loop over states
			}
			for (s = 0; s < csg.getNumStates(); s++) {
				if (!known.get(s)) {
					sol[0][s] = val[0][s];
					sol[1][s] = val[1][s];
				}
				r[s] = sol[0][s] + sol[1][s];
			}
			/*
			String sols;
			sols = "(";
			for (p = 0; p < numCoalitions; p++) {
				if (p < numCoalitions - 1)
					sols += sol[p][csg.getFirstInitialState()] + ",";
				else
					sols += sol[p][csg.getFirstInitialState()] + ")";
			}
			mainLog.println(k + ": " + sols);
			*/
			done = done & PrismUtils.doublesAreClose(sol[0], tmp[0], termCritParam, termCrit == TermCrit.ABSOLUTE);
			done = done & PrismUtils.doublesAreClose(sol[1], tmp[1], termCritParam, termCrit == TermCrit.ABSOLUTE);
			if (done) {
				break;
			}
			else if (!done && k == maxIters) {
				throw new PrismException("Could not converge after " + k + " iterations");
			}
			else {
				done = true;
				tmp[0] = Arrays.copyOf(sol[0], sol[0].length);
				tmp[1] = Arrays.copyOf(sol[1], sol[1].length);
			}
			k++;
		}
		mainLog.println("\nValue iteration converged after " + k + " iterations.");
		mainLog.println("\nPrecomputation took " + timePrecomp / 1000.0 + " seconds.");
		mainLog.println("Coalition results (initial state): (" + sol[0][csg.getFirstInitialState()] + "," + sol[1][csg.getFirstInitialState()] + ")");
		res.soln = r;
		if (generateStrategy || exportAdv) 	{
			if (rew)
				res.strat = new CSGStrategy(csg, lstrat, obj, targets, CSGStrategyType.EQUILIBRIA_R);
			else
				res.strat = new CSGStrategy(csg, lstrat, obj, targets, CSGStrategyType.EQUILIBRIA_P);
		}
		res.numIters = k;
		return res;		
	}
	
	public double[] swne(double[][] eqs, List<List<Map<BitSet, Double>>> strats, boolean min) {
		List<Map<BitSet, Double>> dist = null;
		BitSet isw = new BitSet();
		double[] eq = new double[numCoalitions+1];
		double sum;
		int p;
		eq[0] = 0.0;
		for (p = 0; p < numCoalitions; p++) { // gets first equilibrium
			eq[0] += eqs[0][p];
			eq[p+1] = eqs[0][p];
		}
		isw.set(0);
		if (strats != null) { 
			dist = strats.get(0);
		}
		for (int i = 1; i < eqs.length; i++) { // if there are more than one
			sum = 0.0;
			for (p = 0; p < numCoalitions; p++) { // computes the sum
				sum += eqs[i][p];
			}
			if (Double.compare(sum, eq[0]) > 0) { // selects equilibrium if it has a higher sum
				eq[0] = 0.0;
				for (p = 0; p < numCoalitions; p++) {
					eq[0] += eqs[i][p];
					eq[p+1] = eqs[i][p];
				}
				isw.clear();
				isw.set(i);
				if(strats != null)
					dist = strats.get(i);
			}
			else if (Double.compare(sum, eq[0]) == 0) {
				isw.set(i);
			}
		}
		if (isw.cardinality() != 1) {
			int idx = findMaxEqIndexes(isw, eqs, eq);
			sum = 0.0;
			for (p = 0; p < numCoalitions; p++) { // computes the sum
				sum += eqs[idx][p];
				eq[p+1] = eqs[idx][p];
			}
			eq[0] = sum;
			if(strats != null)
				dist = strats.get(idx);
		}
		final List<Map<BitSet, Double>> strat = dist;
		if (strats != null)
			strats.removeIf((List<Map<BitSet, Double>> e) -> !e.equals(strat));
		if (min) {
			for (int i = 0; i < eq.length; i++)
				eq[i] = -1.0 * eq[i];
		}
		return eq;
	}
	
	public int findMaxEqIndexes(BitSet indexes, double[][] eqs, double eq[]) {
		int idx;
		BitSet tmp = new BitSet();
		BitSet[] maxindexes = new BitSet[numCoalitions];
		double max; 
		for (int p = 0; p < numCoalitions; p++) {
			maxindexes[p] = new BitSet();
			max = eq[p+1];
			for (int i = indexes.nextSetBit(0); i >= 0; i = indexes.nextSetBit(i+1)) {
				if (Double.compare(eqs[i][p], max) > 0) {
					maxindexes[p].clear();
					maxindexes[p].set(i);
					max = eqs[i][p];
				}
				else if (Double.compare(eqs[i][p], max) == 0) {
					maxindexes[p].set(i);
				}
			}
		}
		if (maxindexes[0].cardinality() == 1) {
			idx =  maxindexes[0].nextSetBit(0);
			return idx;
		}
		else {
			tmp.or(maxindexes[0]);
			for (int p = 1; p < numCoalitions; p++) {
				tmp.and(maxindexes[p]);
				if (tmp.cardinality() == 1) {
					idx = tmp.nextSetBit(0);
					return idx;
				}
			}
			// if this part of the code is reached, all players get the same payoff for all equilibria in tmp
			idx = maxindexes[0].nextSetBit(0);
			return idx;
		}
	}
	
	public void extractStrategyFromSupport(ArrayList<Map<BitSet, Double>> eqstrat, List<Map<Integer, BitSet>> mmap, BitSet support) {
		BitSet indx;
		int a, i, p;
		for (p = 0; p < numCoalitions; p++) {
			indx = new BitSet();
			for (a = 0; a < strategies.get(p).size(); a++) {
				indx.set(strategies.get(p).get(a));
			}
			indx.and(support);
			i = indx.nextSetBit(0);
			eqstrat.add(p, new HashMap<BitSet, Double>());
			eqstrat.get(p).put(mmap.get(p).get(strategies.get(p).indexOf(i)), 1.0); // indexOf should be changed
		}
	}
	
	public void addStateRewards(double[][] eqs, CSGRewards csgRewards1, CSGRewards csgRewards2, int s, boolean min) {
		for (int e = 0; e < eqs.length; e++) {
			if (csgRewards1 != null)
				eqs[e][0] += ((min)? -1 * csgRewards1.getStateReward(s) : csgRewards1.getStateReward(s));
			if (csgRewards2 != null)
				eqs[e][1] += ((min)? -1 * csgRewards2.getStateReward(s) : csgRewards2.getStateReward(s));
		}
	}
	
	public void addStateRewards(double[][] eqs, List<CSGRewards> rewards, int s, boolean min) {
		int e, p;
		for (e = 0; e < eqs.length; e++) {
			for (p = 0; p < numCoalitions; p++) {
				if (rewards.get(p) != null)
					eqs[e][p] +=  ((min)? -1 * rewards.get(p).getStateReward(s) : rewards.get(p).getStateReward(s));
			}
		}
	}
	
	public double[][] stepEquilibria(CSG csg, List<CSGRewards> rewards, List<Map<Integer, BitSet>> mmap, List<List<Map<BitSet, Double>>> strats, double[][] val, int s, boolean min) throws PrismException {
		Pair<CSGResultStatus, ArrayList<Double>> eqsresult;
		ArrayList<Double> eq;
		ArrayList<ArrayList<Double>> eqs = new ArrayList<ArrayList<Double>>();
		ArrayList<Map<BitSet, Double>> eqstrat = null;
		BitSet active;
		BitSet support = null;
		long time;
		double[][] result;
		double sumt, sumv, v;
		int n, p, p1, p2;
			
		buildStepGame(csg, rewards, mmap, val, s, min);
		active = csg.getConcurrentPlayers(s);

		time = System.currentTimeMillis();
		if (active.cardinality() == 1) { // in case just one player has a choice
			result = new double[1][numCoalitions];
			p1 = active.nextSetBit(0);
			v = Double.NEGATIVE_INFINITY;
			sumv = Double.NEGATIVE_INFINITY;
			sumt = Double.NEGATIVE_INFINITY;
			for (BitSet entry : utilities.keySet()) {
				sumv = 0.0;
				for (p2 = 0; p2 < numCoalitions; p2++) {
					sumv += utilities.get(entry).get(p2); // computes sum of utilities
				}
				if (utilities.get(entry).get(p1) > v) { // maximizes for player who has a choice
					support = entry;
					sumt = 0.0;
					v = utilities.get(entry).get(p1);
					for (p2 = 0; p2 < numCoalitions; p2++) {
						result[0][p2] = utilities.get(entry).get(p2);
						sumt += utilities.get(entry).get(p2); // sum of the utilities for the selected entry
					}
				}
				else if (utilities.get(entry).get(p1) == v && sumv > sumt) { // case utility for player is the same but sum is higher
					support = entry;
					sumt = 0.0;
					for (p2 = 0; p2 < numCoalitions; p2++) {
						result[0][p2] = utilities.get(entry).get(p2);
						sumt += utilities.get(entry).get(p2);
					}
				}
			}
			if (generateStrategy || exportAdv) {
				eqstrat = new ArrayList<Map<BitSet, Double>>();
				extractStrategyFromSupport(eqstrat, mmap, support);
				strats.add(eqstrat);
			}
			return result;
		} 
		else {	
			boolean allzero = true;
			for (BitSet entry : utilities.keySet()) {
				for (p = 0; p < numCoalitions; p++) {
					allzero = allzero && Double.compare(utilities.get(entry).get(p), 0.0) == 0;
					if (!allzero)
						break;
				}
				if (!allzero)
					break;
			}
			if (allzero) {
				result = new double[1][numCoalitions];
				Arrays.fill(result[0], 0.0);
				if (generateStrategy || exportAdv) {
					eqstrat = new ArrayList<Map<BitSet, Double>>();
					extractStrategyFromSupport(eqstrat, mmap, (BitSet) utilities.keySet().toArray()[0]);
					strats.add(eqstrat);
				}
				return result;
			}
			if (utilities.keySet().size() == 1) {
				result = new double[1][numCoalitions];
				for (BitSet entry : utilities.keySet()) {
					for (p = 0; p < numCoalitions; p++) {
						result[0][p] = utilities.get(entry).get(p);
					}
					if (generateStrategy || exportAdv) {
						eqstrat = new ArrayList<Map<BitSet, Double>>();
						extractStrategyFromSupport(eqstrat, mmap, entry);
						strats.add(eqstrat);
					}
				}
				//System.out.println("-- actions " + actions);
				//System.out.println("-- strategies " + strategies);
				//System.out.println(utilities);
				return result;
			}
			clear();
			computeAssertions();
			for (p = 0; p < numCoalitions; p++) {
				dominated[p] = findDominated(p);
				dominating[p] = findDominating(p);
				//dominated[p] = new BitSet();
				//dominating[p] = new BitSet();
				//if (!dominated[p].isEmpty()) 
					//System.out.println("dominated " + p + " " + dominated[p]);
				//if (!dominating[p].isEmpty()) 
					//System.out.println("dominating " + p + " " + dominating[p]);
			}
			//System.out.println("## " + s);
			buildAllSupports();
			
			eqssmt.setIndexes(strategies);
			eqssmt.setNumPlayers(numCoalitions);
			eqssmt.init();
			eqssmt.translateAssertions(assertions, mapActionIndex);
			
			/*
			eqss.setNumPlayers(numCoalitions);
			eqss.setIndexes(strategies);
			eqss.setAssertions(assertions);
			eqss.setGradient(gradient);
			eqss.init();
			*/
			
			/*
			//if (s == ?) {
				System.out.println("\n-- assertions ");
				for (int c = 0; c < numCoalitions; c++) {
					System.out.println("--- player " + c);
		        	for (int q = 0; q < strategies.get(c).size(); q++)
		        		System.out.println("---- action " + actions.get(c).get(q) + " is " + strategies.get(c).get(q)  + " " + assertions.get(c).get(q));
				}
			//}
			*/
			
			/*
			if (generateStrategy || exportAdv) {
				eqstrat = new ArrayList<Map<BitSet, Double>>();
				for (int p = 0; p < 2; p++) {
					eqstrat.add(p, new HashMap<BitSet, Double>());
					for (int t : nash.getStrat().get(e).get(p).getSupport()) {
						eqstrat.get(p).put(mmap.get(p).get(t), nash.getStrat().get(e).get(p).get(t));
					}								
				} 
				strats.add(e, eqstrat);
			}
			*/
			
			if (allSupports.size() == 1) {
				result = new double[1][numCoalitions];
				for (BitSet entry : allSupports) {
					for (p = 0; p < numCoalitions; p++) {
						result[0][p] = utilities.get(entry).get(p);
					}
					/*
					String act = "";
					if (s == ?) {
						//for (int i = entry.nextSetBit(0); i >= 0; i = entry.nextSetBit(i + 1)) {	
						//	act += "[" +  csg.getActions().get(i) + "]";
						//}
						System.out.println("## Equilibrium (unique entry or dominating) " + entry + " " + act);
					}
					*/
					if (generateStrategy || exportAdv) {
						eqstrat = new ArrayList<Map<BitSet, Double>>();
						extractStrategyFromSupport(eqstrat, mmap, entry);
						strats.add(eqstrat);
					}
				}
				return result;
			}
			else {
				//System.out.println("Solver call");
				//System.out.println(utilities);
				for (int supp = 0; supp < allSupports.size(); supp++) {
					//if (allSupports.get(supp).cardinality() == 25)
					//System.out.println("\n ## Support " + allSupports.get(supp));
					//System.out.println("## Support " + allSupports.get(supp));
					eqsresult = eqssmt.computeEquilibria(allSupports.get(supp), mapActionIndex, s);
					//eqsresult = new Pair<CSGResultStatus, ArrayList<Double>>(CSGResultStatus.UNKNOWN, new ArrayList<Double>());
					if (eqsresult.first == CSGResultStatus.SAT) {
						/*
						String act = "";
						if (s == 0) {
							//System.out.println(csg.getActions());
							//for (int i = allSupports.get(supp).nextSetBit(0); i >= 0; i = allSupports.get(supp).nextSetBit(i + 1)) {	
							//	act += "[" +  csg.getActions().get(i) + "]";
							//}
							System.out.println("## Equilibrium " + allSupports.get(supp) + " " + act);
						}
						*/
						eq = eqsresult.second;
						eqs.add(eq);
						//System.out.println("-- mmap " + mmap);
						//System.out.println("-- strat " + eqssmt.getStrat());
						if (generateStrategy || exportAdv) {
							eqstrat = new ArrayList<Map<BitSet, Double>>();
							for (p = 0; p < numCoalitions; p++) {
								eqstrat.add(p, new HashMap<BitSet, Double>());
								for (int t : eqssmt.getStrat().get(p).getSupport()) {
									eqstrat.get(p).put(mmap.get(p).get(t), eqssmt.getStrat().get(p).get(t));
								}
								//System.out.println("-- strat " + p + " " + eqstrat.get(p));
							}
							strats.add(eqstrat);
							//System.out.println("-- strats " + strats);
						}
						
						time = System.currentTimeMillis() - time;
						//mainLog.println("\nFound an equilibirum in  " + time / 1000.0 + " seconds.");
						//System.out.println(eq);
						//System.exit(1);
					}
					/*
					else if (eqsresult.first == CSGResultStatus.UNKNOWN) {// && allSupports.get(supp).cardinality() == 18) {
						//System.out.println("UNKNOWN " + allSupports.get(supp));
						eqss.setNumPlayers(numCoalitions);
						eqss.setMap(map);
						eqss.setIndexes(strategies);
						eqss.setAssertions(assertions);
						eqss.setGradient(gradient);
						eqss.init();
						//System.out.println("## Support " + allSupports.get(supp));
						eqsresult = eqss.computeEquilibria(allSupports.get(supp), map, s);
						if (eqsresult.first == CSGResultStatus.SAT) {
							//System.out.println("## Equilibrium found " + allSupports.get(supp));
							eq = eqsresult.second;
							eqs.add(eq);
							
							time = System.currentTimeMillis() - time;
							mainLog.println("\nFound an equilibirum in  " + time / 1000.0 + " seconds.");
							System.out.println(eq);
							System.exit(1);	
						}
					}
					*/
				}
			}
		
			//System.out.println(eqs);
						
			result = new double[eqs.size()][numCoalitions];
		
			for (n = 0; n < eqs.size(); n++) {
				for (p = 0; p < numCoalitions; p++) {
					result[n][p] = eqs.get(n).get(p);
				}
			}
		}
		return result;
	}
	
	public double[][] stepEquilibria(CSG csg, CSGRewards csgRewards1, CSGRewards csgRewards2, List<Map<Integer, BitSet>> mmap, CSGLabeledPolytopes nash,
			List<List<Map<BitSet, Double>>> strats, double[][] val, int s, boolean min) throws PrismException {
		Map<BitSet, Double> d1 = null;
		Map<BitSet, Double> d2 = null;
		ArrayList<Map<BitSet, Double>> eqstrat;
		ArrayList<ArrayList<Integer>> nmap;
		ArrayList<ArrayList<ArrayList<Double>>> bmgame;
		double[][] val1s, val2s, result;
		double val1, val2, ent1, ent2;
		int[] mIndxs;
		int nrows, ncols, mrow, mcol;
		boolean equalA, equalB;

		mmap = new ArrayList<Map<Integer, BitSet>>();
		nmap = new ArrayList<ArrayList<Integer>>();
		for (int p = 0; p < 2; p++) {
			mmap.add(p, new HashMap<Integer, BitSet>());
			nmap.add(p, new ArrayList<Integer>());
		}
		bmgame = buildBiMatrixGame(csg, csgRewards1, csgRewards2, mmap, nmap, val, s, min);	
		nrows = bmgame.get(0).size();
		ncols = bmgame.get(0).get(0).size();
		val1s = new double[nrows][ncols];
		val2s = new double[nrows][ncols];

		/*
		// --- Uncomment to print matrices ---
		//if (s == csg.getFirstInitialState()) {
			System.out.println("\n-- matrices for state " + s + " " + csg.getStatesList().get(s));
			for (int p = 0; p < 2; p++) {
				System.out.println("-- player " + p);
				for (int r = 0; r < nrows; r++) {
					System.out.println("-- row " + r + " " + bmgame.get(p).get(r));
				}
			}
			System.out.println(actions);
			System.out.println(strategies);			
			System.out.println(mmap);
			for (Map<Integer, BitSet> lmap : mmap) {	
				for (int i : lmap.keySet()) {
					System.out.println("-- " + i);
					System.out.println(lmap.get(i));
					for (int j = lmap.get(i).nextSetBit(0); j >= 0; j = lmap.get(i).nextSetBit(j+1)) {
						System.out.print(csg.getActions().get(j-1) + " ");
					}
					System.out.println();
				}
			}
			System.out.println();
		//} 
		*/

		if (nrows > 1 && ncols > 1) { // both players have choices
			equalA = true;
			equalB = true;
			ent1 = bmgame.get(0).get(0).get(0);
			ent2 = bmgame.get(1).get(0).get(0);
			for (int r = 0; r < nrows; r++) {
				for (int c = 0; c < ncols; c++) {
					val1 = bmgame.get(0).get(r).get(c);
					val2 = bmgame.get(1).get(r).get(c);
					equalA = equalA && Double.compare(ent1, val1) == 0;
					equalB = equalB && Double.compare(ent2, val2) == 0;
					val1s[r][c] = val1;
					val2s[r][c] = val2;		
				}
			}
			if (!(equalA && equalB)) { // at least one has different entries
				if(equalA || equalB) { // if all entries of one of them are the same
					result = new double[1][2];
					if (equalA) { 
						mIndxs = findMaxIndexes(val2s);
						mrow = mIndxs[0];
						mcol = mIndxs[1];
					}
					else {
						mIndxs = findMaxIndexes(val1s);
						mrow = mIndxs[0];
						mcol = mIndxs[1];
					}
					result[0][0] = val1s[mrow][mcol];
					result[0][1] = val2s[mrow][mcol];
					if (generateStrategy || exportAdv) {
						eqstrat = new ArrayList<Map<BitSet, Double>>();
						eqstrat.add(0, new HashMap<BitSet, Double>());
						eqstrat.get(0).put(mmap.get(0).get(nmap.get(0).get(mrow)), 1.0);
						eqstrat.add(1, new HashMap<BitSet, Double>());
						eqstrat.get(1).put(mmap.get(1).get(nmap.get(1).get(mcol)), 1.0);
						strats.add(0, eqstrat);
					}
					addStateRewards(result, csgRewards1, csgRewards2, s, min);
				}
				else { // both players have choices and matrices are not trivial, call solver
					nash.update(nrows, ncols, val1s, val2s);
					nash.compEq();
					nash.compPayoffs();
					result = new double[nash.getNeq()][2];
					for (int e = 0; e < nash.getNeq(); e++) {
						result[e][0] = nash.getP1p()[e];
						result[e][1] = nash.getP2p()[e];
						if (generateStrategy || exportAdv) {
							eqstrat = new ArrayList<Map<BitSet, Double>>();
							for (int p = 0; p < 2; p++) {
								eqstrat.add(p, new HashMap<BitSet, Double>());
								//System.out.println("-- strat from solver " + nash.getStrat().get(e).get(p).getSupport());
								for (int t : nash.getStrat().get(e).get(p).getSupport()) {
									eqstrat.get(p).put(mmap.get(p).get(nmap.get(p).get(t)), nash.getStrat().get(e).get(p).get(t));
								}
							} 
							strats.add(e, eqstrat);
						}
					}
					addStateRewards(result, csgRewards1, csgRewards2, s, min);
				}
			}
			else { // all entries in both are the same
				result = new double[1][2];
				result[0][0] = ent1;
				result[0][1] = ent2;
				if (generateStrategy || exportAdv) {
					eqstrat = new ArrayList<Map<BitSet, Double>>();
					eqstrat.add(0, new HashMap<BitSet, Double>());
					eqstrat.get(0).put(mmap.get(0).get(nmap.get(0).get(0)), 1.0);
					eqstrat.add(1, new HashMap<BitSet, Double>());
					eqstrat.get(1).put(mmap.get(1).get(nmap.get(0).get(0)), 1.0);
					strats.add(0, eqstrat);
				}
				addStateRewards(result, csgRewards1, csgRewards2, s, min);
			}
		} 
		else { // just one of the players have choices
			result = new double[1][2];
			double vt1, vt2, sumv, sumt;
			if(generateStrategy || exportAdv) {
				d1 = new HashMap<BitSet, Double>();
				d2 = new HashMap<BitSet, Double>();
			}
			val1 = Double.NEGATIVE_INFINITY; 
			val2 = Double.NEGATIVE_INFINITY;
			sumv = Double.NEGATIVE_INFINITY;
			if (nrows > 1 && ncols == 1) {
				for (int r = 0; r < nrows; r++) {
					vt1 = bmgame.get(0).get(r).get(0);
					vt2 = bmgame.get(1).get(r).get(0);
					sumt = vt1 + vt2;
					if (Double.compare(vt1, val1) > 0 || (Double.compare(vt1, val1) == 0 && Double.compare(sumt, sumv) > 0)) {
						if(generateStrategy || exportAdv) {
							d1.clear();
							d1.put(mmap.get(0).get(nmap.get(0).get(r)), 1.0);
						}
						val2 = vt2;
						val1 = vt1;
						sumv = val1 + val2;
					}
				}
				if(generateStrategy || exportAdv) 
					d2.put(mmap.get(1).get(nmap.get(1).get(0)), 1.0);
			} 
			else if (nrows == 1 && ncols > 1) {
				for (int c = 0; c < ncols; c++) {
					vt1 = bmgame.get(0).get(0).get(c);
					vt2 = bmgame.get(1).get(0).get(c);
					sumt = vt1 + vt2;
					if (Double.compare(vt2, val2) > 0 || (Double.compare(vt2, val2) == 0 && Double.compare(sumt, sumv) > 0)) {
						if(generateStrategy || exportAdv) {
							d2.clear();
							d2.put(mmap.get(1).get(nmap.get(1).get(c)), 1.0);
						}
						val2 = vt2;
						val1 = vt1;
						sumv = val1 + val2;
					}
				}
				if(generateStrategy || exportAdv) 
					d1.put(mmap.get(0).get(nmap.get(0).get(0)), 1.0);
			} 
			else if (nrows == 1 && ncols == 1) {
				val1 = bmgame.get(0).get(0).get(0);
				val2 = bmgame.get(1).get(0).get(0);
				if(generateStrategy || exportAdv) {
					d1.put(mmap.get(0).get(nmap.get(0).get(0)), 1.0);
					d2.put(mmap.get(1).get(nmap.get(1).get(0)), 1.0);
				}
			} 
			else {
				throw new PrismException("Error with matrix rank");
			}
			if (generateStrategy || exportAdv) {
				eqstrat = new ArrayList<Map<BitSet, Double>>();
				eqstrat.add(0, d1);
				eqstrat.add(1, d2);	 
				strats.add(0, eqstrat);
			} 
			result[0][0] = val1;
			result[0][1] = val2;
			addStateRewards(result, csgRewards1, csgRewards2, s, min);
		}
		return result;
	}
	
}
