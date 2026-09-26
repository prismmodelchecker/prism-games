//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <david.parker@cs.ox.ac.uk> (University of Oxford)
//	* Gabriel Santos <gabriel.santos@cs.ox.ac.uk> (University of Oxford)
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
import java.util.ArrayList;
import soplex.SoPlex;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.apache.commons.math3.util.Precision;

import explicit.CSGModelCheckerEquilibria.CSGResultStatus;
import explicit.rewards.CSGRewards;
import explicit.rewards.MDPRewards;
import parser.ast.Coalition;
import parser.ast.ExpressionTemporal;
import prism.PrismComponent;
import prism.PrismDevNullLog;
import prism.PrismException;
import prism.PrismLangException;
import prism.PrismNotSupportedException;
import prism.PrismSettings;
import prism.PrismUtils;
import strat.CSGEquilibriumStrategy;
import strat.Strategy;

public class CSGModelCheckerEquilibria extends CSGModelChecker
{
	protected MDPModelChecker mdpmc;
	
	/** Pure supports for a normal form game (indexed by coalition) */
	private ArrayList<BitSet> psupports;

	/** Correlated equilibria: index of each joint action (variable of the joint distribution) */
	private HashMap<BitSet, Integer> ceVarMap;

	/** Dominated actions */
	protected BitSet[] dominated;
	
	/** Optimal two-player Nash equilibria (SW/SF) over the LCP encoding */
	protected CSGNashLCP nashSolver;
	/** Stage-game solver for more than two coalitions (created on first use) */
	protected StageGameSolverScip multiSolver = null;
	/** Solver for correlated equilibria */
	protected CSGCorrelated ceSolver;
	/** Name of the SMT solver */
	protected String smtSolver;
	/** Whether to check for the assumption for equilibria model checking */	

	protected boolean assumptionCheck = false;

	/** Types and criteria for equilibria */
	public static final int NASH = 1;
	public static final int CORR = 2;
	public static final int SWEQ = 3;
	public static final int FAIR = 4;


	/** Different status for SMT equilibria computation */
	public enum CSGResultStatus {
		SAT, UNKNOWN, UNSAT;
	}
	
	/**
	 * Create a new CSGModelCheckerEquilibria, inherit basic state from parent (unless null).
	 */
	public CSGModelCheckerEquilibria(PrismComponent parent) throws PrismException {
		super(parent);
		psupports = new ArrayList<BitSet>();
		ceVarMap = new HashMap<BitSet, Integer>();
		mdpmc = new MDPModelChecker(parent);
		mdpmc.setVerbosity(0);
		mdpmc.setSilentPrecomputations(true);		
		assumptionCheck = getSettings().getBoolean(PrismSettings.PRISM_EQ_ASSUMPTION_CHECK);
		smtSolver = getSettings().getString(PrismSettings.PRISM_SMT_SOLVER);
		switch (smtSolver) {
			case "Z3":
				break;
			case "Yices":
				break;
			case "SCIP":
				break;
			default:
				throw new PrismException("Unknown SMT solver \"" + smtSolver + "\"");
		}
	}
	

	/**
	 * Sets the solver according to the settings and the equilibrium type.
	 * 
	 * @param eqType Correlated/Nash
	 * @return Name
	 * @throws PrismException
	 */
	public String setSolver(int eqType) throws PrismException {
		String name = null;
		switch (eqType) {
			case CORR:
				switch (lpSolver) {
					case "Z3":
						ceSolver = new CSGCorrelatedZ3(maxRows * maxCols, numCoalitions);
						name = ceSolver.getSolverName();
						break;
					case "SoPlex": {
						CSGCorrelatedSoPlex spx = new CSGCorrelatedSoPlex(maxRows * maxCols, numCoalitions);
						spx.setScaler(SoPlex.scalerFromName(soplexScaling));
						ceSolver = spx;
						name = ceSolver.getSolverName();
						break;
					}
					default: throw new PrismException("Unsupported solver for correlated equilibria computation");
				}
				break;
			default: {
				switch (smtSolver) {
					case "Z3":
						nashSolver = new CSGNashLCPZ3();
						name = nashSolver.getSolverName();
						break;
					case "Yices":
						nashSolver = new CSGNashLCPYices();
						name = nashSolver.getSolverName();
						break;
					case "SCIP":
						nashSolver = new CSGNashLCPScip();
						name = nashSolver.getSolverName();
				}
			}
		}
		return name;
	}

	/**
	 * Compute and store information about coalitions (for a nonzero-sum problem):
	 * 
	 * @param csg The CSG
	 * @param coalitions The list of coalitions
	 * @throws PrismException
	 */
	public void buildCoalitions(CSG<Double> csg, List<Coalition> coalitions) throws PrismException {
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
					if (coalitions.get(c).isPlayerIndexInCoalition(p, pmap)) {
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
	}
	
	/**
	 * Finds maximum and average number of actions for all coalitions.
	 * 
	 * @param csg The CSG
	 */
	public void findMaxAvgAct(CSG<Double> csg) {
		String max = "(";
		String avg = "(";
		int c, p, n, s;
		maxRows = 0;
		maxCols = 0;
		avgNumActions = new double[numCoalitions];
		Arrays.fill(avgNumActions, 0.0);
		maxNumActions  = new int[numCoalitions];
		Arrays.fill(maxNumActions, 0);
		for (s = 0; s < csg.getNumStates(); s++) {
			for (c = 0; c < numCoalitions; c++) {
				n = 1;
				for (p = coalitionIndexes[c].nextSetBit(0); p >= 0; p = coalitionIndexes[c].nextSetBit(p + 1)) {
					n *= csg.getIndexesForPlayer(s, p).cardinality();
				}
				maxNumActions[c] = (maxNumActions[c] < n)? n : maxNumActions[c];
				avgNumActions[c] += n;
			}
		}
		for (c = 0; c < numCoalitions; c++) {
			avgNumActions[c] /= csg.getNumStates();
			max += (c < numCoalitions -1)? maxNumActions[c] + "," : maxNumActions[c] + ")";
			avg += (c < numCoalitions -1)? PrismUtils.formatDouble2dp(avgNumActions[c]) + "," : PrismUtils.formatDouble2dp(avgNumActions[c]) + ")";
		}
		mainLog.println("Max/avg (actions): " + max + "/" + avg);
	}
	
	/**
	 * Payoff tables of the stage game built by buildStepGame: for coalition c and action position q (index in
	 * strategies.get(c)), a map from the other coalitions' actions (a BitSet of their action ids, one per coalition)
	 * to the payoff of c when it plays q against them. Used for dominance and for the correlated equilibria constraints.
	 */
	public ArrayList<ArrayList<HashMap<BitSet, Double>>> buildPayoffTables() throws PrismException {
		ArrayList<ArrayList<HashMap<BitSet, Double>>> tables = new ArrayList<>();
		// position of each action id in strategies.get(c)
		HashMap<Integer, Integer> position = new HashMap<>();
		for (int c = 0; c < numCoalitions; c++) {
			ArrayList<HashMap<BitSet, Double>> tc = new ArrayList<>();
			for (int q = 0; q < strategies.get(c).size(); q++) {
				tc.add(new LinkedHashMap<BitSet, Double>());
				position.put(strategies.get(c).get(q), q);
			}
			tables.add(tc);
		}
		BitSet own = new BitSet();
		for (Entry<BitSet, ArrayList<Double>> e : utilities.entrySet()) {
			for (int c = 0; c < numCoalitions; c++) {
				own.clear();
				own.or(psupports.get(c));
				own.and(e.getKey());
				if (own.cardinality() != 1)
					throw new PrismException("Joint action " + e.getKey() + " does not have exactly one action of coalition " + c);
				int id = own.nextSetBit(0);
				BitSet others = (BitSet) e.getKey().clone();
				others.clear(id);
				tables.get(c).get(position.get(id)).put(others, e.getValue().get(c));
			}
		}
		return tables;
	}

	/**
	 * Actions of coalition p (as action ids) strictly dominated by another of its actions, from the payoff tables.
	 */
	public BitSet findDominated(int p, ArrayList<ArrayList<HashMap<BitSet, Double>>> tables) throws PrismException {
		BitSet domi = new BitSet();
		ArrayList<HashMap<BitSet, Double>> tp = tables.get(p);
		for (int a1 = 0; a1 < tp.size(); a1++) {
			for (int a2 = 0; a2 < tp.size(); a2++) {
				if (a1 == a2)
					continue;
				boolean dominated = true;
				for (Entry<BitSet, Double> e : tp.get(a1).entrySet()) {
					Double v2 = tp.get(a2).get(e.getKey());
					if (v2 == null)
						throw new PrismException("Incomplete payoff table for coalition " + p);
					if (!(e.getValue() < v2)) {
						dominated = false;
						break;
					}
				}
				if (dominated) {
					domi.set(strategies.get(p).get(a1));
					break;
				}
			}
		}
		return domi;
	}

	/**
	 * Finds row and column indexes for the maximum entry in a matrix.
	 * 
	 * @param a Matrix
	 * @return
	 */
	public int[] findMaxIndexes(double[][] a) {
		int result[] = new int[2];
		result[0] = 0;
		result[1] = 0;	
		double value = Double.NEGATIVE_INFINITY;
		for(int r = 0; r < a.length; r++) {
			for(int c = 0; c < a[r].length; c++) {
				if(Double.compare(a[r][c], value) > 0) {
					value = a[r][c];
					result[0] = r;
					result[1] = c;	
				}
			}
		}		
		return result;
	}
	
	/**
	 * 
	 * 
	 * @param eqs
	 * @param csgRewards1
	 * @param csgRewards2
	 * @param s
	 * @param min
	 */
	public void addStateRewards(double[][] eqs, CSGRewards<Double> csgRewards1, CSGRewards<Double> csgRewards2, int s, boolean min) {
		for (int e = 0; e < eqs.length; e++) {
			if (csgRewards1 != null)
				eqs[e][0] += ((min)? -1 * csgRewards1.getStateReward(s) : csgRewards1.getStateReward(s));
			if (csgRewards2 != null)
				eqs[e][1] += ((min)? -1 * csgRewards2.getStateReward(s) : csgRewards2.getStateReward(s));
		}
	}
	
	/**
	 * 
	 * 
	 * @param eqs
	 * @param rewards
	 * @param s
	 * @param min
	 */
	public void addStateRewards(double[][] eqs, List<CSGRewards<Double>> rewards, int s, boolean min) {
		int e, p;
		for (e = 0; e < eqs.length; e++) {
			for (p = 0; p < numCoalitions; p++) {
				if (rewards.get(p) != null)
					eqs[e][p] +=  ((min)? -1 * rewards.get(p).getStateReward(s) : rewards.get(p).getStateReward(s));
			}
		}
	}
	
	/**
	 * 
	 * 
	 * @param eqs
	 * @param rewards
	 * @param s
	 * @param min
	 */
	public void addStateRewards(double[] eqs, List<CSGRewards<Double>> rewards, int s, boolean min) {
		for (int p = 0; p < numCoalitions; p++) {
			if (rewards.get(p) != null)
				eqs[p+1] +=  ((min)? -1.0 * rewards.get(p).getStateReward(s) : rewards.get(p).getStateReward(s));
		}
	}
	
	/**
	 * Finds the SWNE for when just one player has choices. 
	 * 
	 * @param mmap Index map
	 * @param strats Overall strategy
	 * @param eqstrat Strategy for the current state
	 * @param active Active player
	 * @return
	 */
	public double[][] findSWNEOnePlayer(List<Map<Integer, BitSet>> mmap, List<List<Map<BitSet, Double>>> strats, List<Map<BitSet, Double>> eqstrat, BitSet active) {
		BitSet support = null;
		double[][] result;
		double sumt, sumv, v;
		int p1, p2;
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
			if (Double.compare(utilities.get(entry).get(p1), v) > 0) { // maximizes for player who has a choice
				support = entry;
				sumt = 0.0;
				v = utilities.get(entry).get(p1);
				for (p2 = 0; p2 < numCoalitions; p2++) {
					result[0][p2] = utilities.get(entry).get(p2);
					sumt += utilities.get(entry).get(p2); // sum of the utilities for the selected entry
				}
			}
			else if (Double.compare(utilities.get(entry).get(p1), v) == 0 && Double.compare(sumv, sumt) > 0) { // case utility for player is the same but sum is higher
				support = entry;
				sumt = 0.0;
				for (p2 = 0; p2 < numCoalitions; p2++) {
					result[0][p2] = utilities.get(entry).get(p2);
					sumt += utilities.get(entry).get(p2);
				}
			}
		}
		if (genStrat) {
			eqstrat = new ArrayList<Map<BitSet, Double>>();
			extractStrategyFromSupport(mmap, eqstrat, support);
			strats.add(eqstrat);
		}
		return result;
	}
	
	/**
	 * Extracts the strategy for a given support.
	 * 
	 * @param mmap Index map
	 * @param eqstrat Strategy for the current state
	 * @param support Support
	 */
	public void extractStrategyFromSupport(List<Map<Integer, BitSet>> mmap, List<Map<BitSet, Double>> eqstrat, BitSet support) {
		BitSet indx = new BitSet();
		int a, i, p;
		for (p = 0; p < numCoalitions; p++) {
			indx.clear();
			for (a = 0; a < strategies.get(p).size(); a++) {
				indx.set(strategies.get(p).get(a));
			}
			indx.and(support);
			i = indx.nextSetBit(0);
			eqstrat.add(p, new HashMap<BitSet, Double>());
			eqstrat.get(p).put(mmap.get(p).get(strategies.get(p).indexOf(i)), 1.0);
		}
	}
	
	/**
	 * Build info needed for the utility table to solve a CSG state s. 
	 * 
	 * @param csg The CSG
	 * @param rewards List of rewards
	 * @param mmap Index map
	 * @param val Current values for each state
	 * @param s State index
	 * @param min Whether minimising/maximising
	 * @throws PrismException
	 */
	public void buildStepGame(CSG<Double> csg, List<CSGRewards<Double>> rewards, List<Map<Integer, BitSet>> mmap, BitSet D, BitSet E,
							  double[][] val, int s, boolean min) throws PrismException {
		Map<BitSet, Integer> imap = new HashMap<BitSet, Integer>();
		BitSet jidx;
		BitSet indexes = new BitSet();
		BitSet tmp = new BitSet();
		String act;
		double v;
		int c, i, p, t;
		int[] joint;
		int[] idle = new int[numPlayers];
		ceVarMap.clear();
		actions.clear();
		psupports.clear();
		strategies.clear();
		utilities.clear();
		varIndex = 0;
		Arrays.fill(idle, -1);
		for (c = 0; c < numCoalitions; c++) {
			actions.add(c, new ArrayList<String>());
			psupports.add(c, new BitSet());
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
						psupports.get(c).set(varIndex);
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
			ceVarMap.put(jidx, utilities.keySet().size() - 1);
			for (c = 0; c < numCoalitions; c++) {
				v = 0.0;
				if (D != null && D.get(c)) {
					if (rewards == null) 
						v = 1.0;
				}
				else if (E != null && E.get(c)) {
					v = 0.0; // failed (until): value 0, no incentives
				}
				else {
					for (int d : csg.getChoice(s, t).getSupport()) {
						if (!Double.isNaN(val[c][d])) {
							v += csg.getChoice(s, t).get(d) * val[c][d];
						}
						else {
							mainLog.println("val[c][d]: " + val[c][d]);
							mainLog.println("\n## state " + s);
							mainLog.println("-- strategies " + strategies);
							mainLog.println("-- actions " + actions);
							mainLog.println("-- utilities " + utilities);
							throw new PrismException("Error in building game for state " + s);
						} 
					} 
					if (rewards != null) {
						if (rewards.get(c) != null)
							v += (Double) rewards.get(c).getTransitionReward(s, t);		
					}
					v = Precision.round(v, 12, BigDecimal.ROUND_HALF_EVEN);
				}
				utilities.get(jidx).add(c, (min)? -1.0 * v : v); // might have to add min (v, 1.0) due to assertions for probabilistic
			}
		}	
		//System.out.println("-- imap " + imap);
		//if (s == csg.getFirstInitialState()) {
			//System.out.println("\n## state " + s);
			//System.out.println("-- strategies " + strategies);
			//System.out.println("-- actions " + actions);
			//System.out.println("-- utilities " + utilities);
			//System.out.println("-- mmap " + mmap);
		//}
	}
	
	/**
	 * Builds a bimatrix game (two-player case). 
	 * 
	 * @param csg The CSG
	 * @param r1 Rewards for the first coalition
	 * @param r2 Rewards for the second coalition
	 * @param mmap Index map
	 * @param nmap Reduced index map
	 * @param val Current values for each state 
	 * @param s State index
	 * @param min Whether minimising/maximising
	 * @return
	 * @throws PrismException
	 */
	public ArrayList<ArrayList<ArrayList<Double>>> buildBimatrixGame(CSG<Double> csg, CSGRewards<Double> r1, CSGRewards<Double> r2, List<Map<Integer, BitSet>> mmap,  List<ArrayList<Integer>> nmap, double[][] val, int s, boolean min) throws PrismException {
		ArrayList<ArrayList<ArrayList<Double>>> bmgame = new ArrayList<ArrayList<ArrayList<Double>>>();
		ArrayList<CSGRewards<Double>> rewards = null;
		BitSet action = new BitSet();
		int col, p, row, irow, icol;
		if (numCoalitions > 2) 
			throw new PrismLangException("Multiplayer game not supported by this method");
		if (r1 != null || r2 != null) {
			rewards = new ArrayList<>();
			rewards.add(0, r1);
			rewards.add(1, r2);
		}
		buildStepGame(csg, rewards, mmap, null, null, val, s, min);
		//System.out.println("-- utilities " + utilities);
		//System.out.println("-- strategies " + strategies);
		//System.out.println("-- mmap " + mmap);
		ArrayList<ArrayList<HashMap<BitSet, Double>>> tables = buildPayoffTables();
		for (p = 0; p < numCoalitions; p++) {
			dominated[p] = findDominated(p, tables);
		}
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
	
	/**
	 * Deal with two-player bounded equilibria.
	 * 
	 * @param csg
	 * @param coalitions
	 * @param rewards
	 * @param exprs
	 * @param targets
	 * @param remain
	 * @param bounds
	 * @param eqType
	 * @param crit
	 * @param min
	 * @return
	 * @throws PrismException
	 */
	public ModelCheckerResult computeBoundedEquilibria(CSG<Double> csg, List<Coalition> coalitions, List<CSGRewards<Double>> rewards, List<ExpressionTemporal> exprs, BitSet[] targets, BitSet[] remain, int[] bounds, int eqType, int crit, boolean min) throws PrismException {
		if (genStrat) {
			throw new PrismNotSupportedException("Strategy synthesis for bounded properties is not supported yet");
		}
		ModelCheckerResult res = new ModelCheckerResult();
		List<CSGRewards<Double>> newRewards = null;
		BitSet[] only = new BitSet[coalitions.size()];
		BitSet[] phi1 = new BitSet[3];
		BitSet cpy =  new BitSet();
		double[][] sol = new double[coalitions.size()][csg.getNumStates()];
		double[][] tmp = new double[coalitions.size()][csg.getNumStates()];
		double[][] val = new double[coalitions.size()][csg.getNumStates()];
		double[] eq;
		double[] r = new double[csg.getNumStates()];
		int i, j, n1, n2, k, s;
		boolean rew;
		long currentTime, timePrecomp;		
		
		rew = rewards != null;
		
		buildCoalitions(csg, coalitions);
		findMaxRowsCols(csg);
		mainLog.println("Starting bounded equilibria computation (solver=" + setSolver(eqType) + ")...");
		dominated = new BitSet[numCoalitions];
		
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
				eq = stepEquilibriaTwoPlayer(csg, null, null, null, sol, s, eqType, crit, rew, min);
				tmp[0][s] = eq[1];
				tmp[1][s] = eq[2];
				r[s] = eq[1] + eq[2];
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
			newRewards = new ArrayList<>();
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
						eq = stepEquilibriaTwoPlayer(csg, newRewards, null, null, sol, s, eqType, crit, rew, min);	
						tmp[0][s] = eq[1];
						tmp[1][s] = eq[2];
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
							eq = stepEquilibriaTwoPlayer(csg, null, null, null, sol, s, eqType, crit, rew, min);
							tmp[0][s] = eq[1];
							tmp[1][s] = eq[2];
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

	
	/**
	 * Deal with multi-player infinite-horizon equilibria (unfinished).
	 * 
	 * @param csg
	 * @param coalitions
	 * @param rewards
	 * @param targets
	 * @param remain
	 * @param eqType
	 * @param crit
	 * @param min
	 * @return
	 * @throws PrismException
	 */
	public ModelCheckerResult computeMultiEquilibria(CSG<Double> csg, List<Coalition> coalitions, List<CSGRewards<Double>> rewards, List<ExpressionTemporal> exprs,
			BitSet bounded, BitSet[] targets, BitSet[] remain, int[] bounds, int eqType, int crit, boolean min) throws PrismException {
		mainLog.println("\n# Running multi-player equilibria...\n");
		ModelCheckerResult res = new ModelCheckerResult();
		double[][] sol;
		double[] r;
		long timeTaken;
		int c, s;
		
		buildCoalitions(csg, coalitions);
		dominated = new BitSet[numCoalitions];
		findMaxAvgAct(csg);
		if (eqType == CORR) {
			int maxSize = 1;
			for (c = 0; c < numCoalitions; c++)
				maxSize = maxSize * maxNumActions[c];
			switch (lpSolver) {
				case "Z3":
					ceSolver = new CSGCorrelatedZ3(maxSize, numCoalitions);
					break;
				case "SoPlex": {
					CSGCorrelatedSoPlex spx = new CSGCorrelatedSoPlex(maxSize, numCoalitions);
					spx.setScaler(SoPlex.scalerFromName(soplexScaling));
					ceSolver = spx;
					break;
				}
				default:
					throw new PrismException("Unsupported solver for correlated equilibria computation");
			}
		}

		// Kind of objective and bound of each coalition
		int[] kind = new int[numCoalitions];
		int[] bound = new int[numCoalitions];
		for (c = 0; c < numCoalitions; c++) {
			int op = exprs.get(c).getOperator();
			bound[c] = bounded.get(c) ? bounds[c] : -1;
			if (!bounded.get(c))
				kind[c] = CSGMultiObjectives.OBJ_UNBOUNDED;
			else if (op == ExpressionTemporal.P_X)
				kind[c] = CSGMultiObjectives.OBJ_NEXT;
			else if (op == ExpressionTemporal.R_C)
				kind[c] = CSGMultiObjectives.OBJ_CUMULATIVE;
			else if (op == ExpressionTemporal.R_I)
				kind[c] = CSGMultiObjectives.OBJ_INSTANTANEOUS;
			else
				kind[c] = CSGMultiObjectives.OBJ_BOUNDED;
		}

		BitSet unboundedObjs = new BitSet();
		unboundedObjs.set(0, numCoalitions);
		unboundedObjs.andNot(bounded);
		checkStopping(csg, targets, remain, unboundedObjs, rewards != null);

		timeTaken = System.currentTimeMillis();
		CSGMultiObjectives objectives = new CSGMultiObjectives(numCoalitions, rewards != null, kind, bound, targets, remain);
		MultiSubgames subgames = new MultiSubgames(csg, rewards, objectives, eqType, crit, min);
		sol = bounded.isEmpty() ? subgames.solve(new BitSet(), new BitSet(), true) : subgames.timed(new BitSet(), new BitSet(), 0);
		timeTaken = System.currentTimeMillis() - timeTaken;
		if (genStrat) {
			CSGEquilibriumStrategy strat = new CSGEquilibriumStrategy(csg, objectives, eqType == CORR, subgames.localStrategies);
			if (rewards == null)
				strat.setHopeless(hopeless(csg, targets, remain));
			res.strat = strat;
		}

		mainLog.println();
		if (!bounded.isEmpty())
			mainLog.println("Horizon of the bounded objectives: " + objectives.horizon + " steps (" + subgames.numTimed + " subgame steps computed)");
		mainLog.println("Unbounded subgames solved: " + subgames.memo.size() + " (value iteration: " + subgames.totalIters + " iterations in total)");
		for (c = 0; c < numCoalitions; c++) {
			mainLog.println("Result for coalition " + coalitions.get(c) + ": " + sol[c][csg.getFirstInitialState()] + " (value in the initial state).");
		}
		if (genStrat) {
			DTMCModelChecker dtmcmc = new DTMCModelChecker(this);
			dtmcmc.inheritSettings(this);
			dtmcmc.setSilentPrecomputations(true);
			dtmcmc.setLog(new PrismDevNullLog());
			double[] computed = new double[numCoalitions];
			for (c = 0; c < numCoalitions; c++)
				computed[c] = sol[c][csg.getFirstInitialState()];
			checkStrategyValues(coalitions, computed, ((CSGEquilibriumStrategy) res.strat).achievedValues(rewards, dtmcmc));
		}
		r = new double[csg.getNumStates()];
		for (s = 0; s < csg.getNumStates(); s++) {
			r[s] = 0.0;
			for (c = 0; c < numCoalitions; c++) {
				r[s] += sol[c][s];
			}
		}
		res.soln = r;
		res.numIters = subgames.totalIters;
		res.timeTaken = timeTaken / 1000.0;
		return res;
	}

	/**
	 * Value iteration for multi-player (unbounded) equilibria, over the subgames (D, E): D are the coalitions that
	 * have reached their targets (value 1 for probabilities; no more rewards), E those that have failed their until
	 * objective (value 0). Coalitions in D and E have no incentives (free) and keep acting. In a state where further
	 * coalitions reach their targets (or fail), the values are those of the subgame (D + newD, E + newE), which has
	 * more coalitions done and is solved first (recursively) and only once (memoised).
	 */
	private class MultiSubgames
	{
		final CSG<Double> csg;
		final List<CSGRewards<Double>> rewards;
		/** Rewards accumulated step by step (as rewards, but null for instantaneous objectives) */
		final List<CSGRewards<Double>> stepRewards;
		/** The objectives, and when each is decided */
		final CSGMultiObjectives obj;
		final int eqType, crit;
		final boolean min, rew;
		/** Solved (unbounded) subgames, keyed by D (bits 0..n-1) and E (bits n..2n-1) */
		final Map<BitSet, double[][]> memo = new HashMap<>();
		/** Values of the subgames at each step t < horizon (while some bounded objective is undecided) */
		final List<Map<BitSet, double[][]>> memoTimed = new ArrayList<>();
		int totalIters = 0, numTimed = 0;

		/**
		 * Local strategies (if generating strategies), per memory (the subgame and, while bounded objectives are
		 * undecided, the step: see memoryKey) and state. As for two coalitions: Nash, one map per coalition from its
		 * actions (player action indices) to probabilities; correlated, one map (index 0) from joint actions to
		 * probabilities. Null for states where the play moves to a larger subgame (and so is decided there).
		 */
		final Map<BitSet, List<Map<BitSet, Double>>[]> localStrategies = new HashMap<>();
		/** Local strategy of the last stage game solved (by stageValues) */
		private List<Map<BitSet, Double>> lastStrategy;

		/** Memory of a strategy: D (bits 0..n-1), E (bits n..2n-1) and, for t >= 0, the step (bit 2n + 1 + t) */
		BitSet memoryKey(BitSet D, BitSet E, int t)
		{
			BitSet k = key(D, E);
			if (t >= 0)
				k.set(2 * numCoalitions + 1 + t);
			return k;
		}

		@SuppressWarnings("unchecked")
		private List<Map<BitSet, Double>>[] strategiesFor(BitSet memory)
		{
			return localStrategies.computeIfAbsent(memory, m -> (List<Map<BitSet, Double>>[]) new List[csg.getNumStates()]);
		}

		MultiSubgames(CSG<Double> csg, List<CSGRewards<Double>> rewards, CSGMultiObjectives obj, int eqType, int crit, boolean min)
		{
			this.csg = csg;
			this.rewards = rewards;
			this.obj = obj;
			this.eqType = eqType;
			this.crit = crit;
			this.min = min;
			this.rew = rewards != null;
			for (int t = 0; t <= obj.horizon; t++)
				memoTimed.add(new HashMap<>());
			if (rew) {
				stepRewards = new ArrayList<>();
				for (int c = 0; c < numCoalitions; c++)
					stepRewards.add(obj.kind[c] == CSGMultiObjectives.OBJ_INSTANTANEOUS ? null : rewards.get(c));
			} else {
				stepRewards = null;
			}
		}

		/**
		 * Values of the subgame (D, E) after t steps. Once all bounded objectives are decided (their coalitions are
		 * in D or E), these are the values of the unbounded subgame, whatever t. Otherwise (t < horizon), one step of
		 * backward induction from the values after t + 1 steps; in states where coalitions become done or fail,
		 * the values are those of the larger subgame at the same step (plus the reward of instantaneous objectives
		 * that end there).
		 */
		double[][] timed(BitSet D, BitSet E, int t) throws PrismException
		{
			if (obj.boundedDecided(D, E))
				return solve(D, E, false);
			BitSet k = key(D, E);
			double[][] sol = memoTimed.get(t).get(k);
			if (sol != null)
				return sol;

			int n = csg.getNumStates();
			sol = new double[numCoalitions][n];
			double[][] next = null;
			for (int s = 0; s < n; s++) {
				BitSet uniD = (BitSet) D.clone(), uniE = (BitSet) E.clone();
				obj.classify(s, t, uniD, uniE);
				if (!uniD.equals(D) || !uniE.equals(E)) {
					double[][] sub = timed(uniD, uniE, t);
					for (int c = 0; c < numCoalitions; c++) {
						sol[c][s] = sub[c][s];
						if (uniD.get(c) && !D.get(c) && obj.kind[c] == CSGMultiObjectives.OBJ_INSTANTANEOUS)
							sol[c][s] += rewards.get(c).getStateReward(s);
					}
					continue;
				}
				if (t >= obj.horizon)
					throw new PrismException("Bounded objective undecided after " + t + " steps in state " + s);
				if (next == null)
					next = timed(D, E, t + 1);
				double[] stage = stageValues(D, E, next, s, stepRewards);
				for (int c = 0; c < numCoalitions; c++) {
					if (D.get(c))
						sol[c][s] = rew ? 0.0 : 1.0;
					else if (E.get(c))
						sol[c][s] = 0.0;
					else
						sol[c][s] = stage[c] + (rew ? stateReward(c, s, stepRewards) : 0.0);
				}
				if (genStrat)
					strategiesFor(memoryKey(D, E, t))[s] = lastStrategy;
			}
			memoTimed.get(t).put(k, sol);
			numTimed++;
			return sol;
		}

		BitSet key(BitSet D, BitSet E)
		{
			BitSet k = (BitSet) D.clone();
			for (int c = E.nextSetBit(0); c >= 0; c = E.nextSetBit(c + 1))
				k.set(numCoalitions + c);
			return k;
		}

		/** Values of the subgame (D, E) for all coalitions and states (from the memo if already solved) */
		double[][] solve(BitSet D, BitSet E, boolean main) throws PrismException
		{
			BitSet k = key(D, E);
			double[][] sol = memo.get(k);
			if (sol == null) {
				sol = iterate(D, E, main);
				memo.put(k, sol);
			}
			return sol;
		}

		private double[][] iterate(BitSet D, BitSet E, boolean main) throws PrismException
		{
			int n = csg.getNumStates();
			double[][] sol = new double[numCoalitions][n];
			double[][] prev = new double[numCoalitions][n];
			double[] stage;
			int c, s, iters;
			boolean done;

			// All coalitions done: fixed values
			if (D.cardinality() + E.cardinality() == numCoalitions) {
				for (c = 0; c < numCoalitions; c++)
					Arrays.fill(sol[c], (!rew && D.get(c)) ? 1.0 : 0.0);
				return sol;
			}

			// States where further coalitions become done (or fail): values from the larger subgames
			double[][][] fromSub = new double[n][][];
			for (s = 0; s < n; s++) {
				// (only unbounded objectives are undecided here, so the step does not matter)
				BitSet uniD = (BitSet) D.clone(), uniE = (BitSet) E.clone();
				obj.classify(s, -1, uniD, uniE);
				if (!uniD.equals(D) || !uniE.equals(E))
					fromSub[s] = solve(uniD, uniE, false);
			}

			// Initial values: done coalitions fixed, others 0
			for (c = 0; c < numCoalitions; c++)
				Arrays.fill(sol[c], (!rew && D.get(c)) ? 1.0 : 0.0);
			iters = 0;
			while (true) {
				for (c = 0; c < numCoalitions; c++)
					prev[c] = Arrays.copyOf(sol[c], n);
				for (s = 0; s < n; s++) {
					if (fromSub[s] != null) {
						for (c = 0; c < numCoalitions; c++)
							sol[c][s] = fromSub[s][c][s];
						continue;
					}
					stage = stageValues(D, E, prev, s, rewards);
					for (c = 0; c < numCoalitions; c++) {
						if (D.get(c))
							sol[c][s] = rew ? 0.0 : 1.0;
						else if (E.get(c))
							sol[c][s] = 0.0;
						else
							sol[c][s] = stage[c] + (rew ? stateReward(c, s, rewards) : 0.0);
					}
					if (genStrat) {
						// as for two coalitions: keep the first strategy, replace it only when the values change
						List<Map<BitSet, Double>>[] strat = strategiesFor(memoryKey(D, E, -1));
						boolean changed = false;
						for (c = 0; c < numCoalitions; c++)
							changed = changed || Double.compare(sol[c][s], prev[c][s]) != 0;
						if (strat[s] == null || (changed && !strat[s].equals(lastStrategy)))
							strat[s] = lastStrategy;
					}
				}
				iters++;
				done = true;
				for (c = 0; c < numCoalitions; c++)
					done = done & PrismUtils.doublesAreClose(sol[c], prev[c], termCritParam, termCrit == TermCrit.ABSOLUTE);
				if (main) {
					StringBuilder sb = new StringBuilder("(");
					for (c = 0; c < numCoalitions; c++)
						sb.append(c > 0 ? "," : "").append(sol[c][csg.getFirstInitialState()]);
					mainLog.println(iters + ": " + sb + ")");
				}
				if (done && iters > 1)
					break;
				if (iters == maxIters) {
					String msg = "Value iteration did not converge within " + iters + " iterations (subgame D=" + D + ", E=" + E + ")";
					if (errorOnNonConverge)
						throw new PrismException(msg);
					mainLog.printWarning(msg + "; the values are those of the last iteration");
					break;
				}
			}
			totalIters += iters;
			return sol;
		}

		/** State reward of coalition c in s (values here are true values: for min, the stage games are negated and
		 *  their results negated back) */
		private double stateReward(int c, int s, List<CSGRewards<Double>> rs)
		{
			CSGRewards<Double> r = rs.get(c);
			return r == null ? 0.0 : r.getStateReward(s);
		}

		/** Equilibrium values (without state rewards) of the stage game at s, given the values val of the successors */
		private double[] stageValues(BitSet D, BitSet E, double[][] val, int s, List<CSGRewards<Double>> rs) throws PrismException
		{
			double[] v = new double[numCoalitions];
			double[] eq;
			// index maps and strategy (if generating strategies): mmap.get(c) maps coalition c's action positions
			// to its players' action indices, which is how the local strategies are expressed
			List<Map<Integer, BitSet>> mmap = null;
			List<List<Map<BitSet, Double>>> strats = null;
			if (genStrat) {
				mmap = new ArrayList<>();
				for (int c = 0; c < numCoalitions; c++)
					mmap.add(new HashMap<Integer, BitSet>());
				strats = new ArrayList<>();
			}
			if (eqType == CORR) {
				// stepCorrelatedEquilibria adds the state rewards of all coalitions: removed here (added by the caller when due)
				eq = stepCorrelatedEquilibria(csg, rs, mmap, strats, D, E, val, s, min, crit);
				for (int c = 0; c < numCoalitions; c++)
					v[c] = eq[c + 1] - (rew ? stateReward(c, s, rs) : 0.0);
			} else {
				eq = toEquilibrium(stepEquilibriaMulti(csg, rs, mmap, strats, D, E, val, s, min, crit), strats, min);
				for (int c = 0; c < numCoalitions; c++)
					v[c] = eq[c + 1];
			}
			lastStrategy = genStrat ? strats.get(0) : null;
			return v;
		}
	}

	/**
	 * All subsets of {0, ..., n-1} with k elements.
	 */
	public static List<BitSet> subsetsOfSize(int n, int k) {
		List<BitSet> result = new ArrayList<BitSet>();
		subsetsOfSize(n, k, 0, new BitSet(), result);
		return result;
	}

	private static void subsetsOfSize(int n, int k, int from, BitSet current, List<BitSet> result) {
		if (current.cardinality() == k) {
			result.add((BitSet) current.clone());
			return;
		}
		for (int i = from; i < n; i++) {
			current.set(i);
			subsetsOfSize(n, k, i + 1, current, result);
			current.clear(i);
		}
	}
	
	
	
	/**
	 * 
	 * 
	 * @param sol
	 * @param eq
	 * @param s
	 * @return
	 */
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
	
	/**
	 * The strategy computed for two coalitions (unbounded objectives), as a CSGEquilibriumStrategy: with neither
	 * coalition decided, the local strategies of the stage games (lstrat); once one is done or has failed, the
	 * other's optimal strategy in the MDP where the coalitions choose actions jointly (from the precomputation obj),
	 * as a (deterministic) joint action.
	 */
	@SuppressWarnings("unchecked")
	private CSGEquilibriumStrategy twoPlayerStrategy(CSG<Double> csg, List<List<List<Map<BitSet, Double>>>> lstrat, ModelCheckerResult[] obj,
			BitSet[] targets, BitSet[] remain, boolean rew, boolean corr)
	{
		int n = csg.getNumStates();
		int[] kind = { CSGMultiObjectives.OBJ_UNBOUNDED, CSGMultiObjectives.OBJ_UNBOUNDED };
		int[] bound = { -1, -1 };
		CSGMultiObjectives objectives = new CSGMultiObjectives(2, rew, kind, bound, targets, rew ? null : remain);
		Map<BitSet, List<Map<BitSet, Double>>[]> local = new HashMap<>();
		// neither decided: stage game strategies
		List<Map<BitSet, Double>>[] none = (List<Map<BitSet, Double>>[]) new List[n];
		for (int s = 0; s < n; s++) {
			List<Map<BitSet, Double>> ls = new ArrayList<>();
			for (int c = 0; c < (corr ? 1 : 2); c++) {
				Map<BitSet, Double> m = lstrat.get(c).get(0).get(s);
				if (m == null) {
					ls = null;
					break;
				}
				ls.add(m);
			}
			none[s] = ls;
		}
		local.put(new BitSet(), none);
		// one decided (done or failed): the other's MDP strategy
		for (int p = 0; p < 2; p++) {
			int q = 1 - p;
			List<Map<BitSet, Double>>[] mdp = (List<Map<BitSet, Double>>[]) new List[n];
			for (int s = 0; s < n; s++) {
				int t = obj[q].strat == null ? -1 : obj[q].strat.getChoiceIndex(s, -1);
				if (t < 0)
					continue;
				BitSet joint = new BitSet();
				int[] indexes = csg.getIndexes(s, t);
				for (int i = 0; i < indexes.length; i++)
					joint.set(indexes[i] > 0 ? indexes[i] : csg.getIdles()[i]);
				List<Map<BitSet, Double>> ls = new ArrayList<>();
				if (corr) {
					ls.add(Collections.singletonMap(joint, 1.0));
				} else {
					for (int c = 0; c < 2; c++) {
						BitSet b = (BitSet) joint.clone();
						b.and(actionIndexes[c]);
						ls.add(Collections.singletonMap(b, 1.0));
					}
				}
				mdp[s] = ls;
			}
			BitSet done = new BitSet(), failed = new BitSet();
			done.set(p);
			failed.set(2 + p);
			local.put(done, mdp);
			local.put(failed, mdp);
		}
		CSGEquilibriumStrategy strat = new CSGEquilibriumStrategy(csg, objectives, corr, local);
		if (!rew)
			strat.setHopeless(hopeless(csg, targets, remain));
		return strat;
	}

	/** Per coalition, the states from which its (probabilistic) objective cannot be satisfied under any profile */
	private BitSet[] hopeless(CSG<Double> csg, BitSet[] targets, BitSet[] remain)
	{
		BitSet[] h = new BitSet[targets.length];
		for (int c = 0; c < targets.length; c++)
			h[c] = mdpmc.prob0((MDP<Double>) csg, remain == null ? null : remain[c], targets[c], false, null);
		return h;
	}

	/**
	 * Compares the values achieved by a synthesised strategy (in the initial state) with the computed ones, warning if
	 * they differ: the strategy then does not realise the equilibrium that was computed (e.g. in non-stopping games).
	 */
	protected void checkStrategyValues(List<Coalition> coalitions, double[] computed, double[] achieved)
	{
		double tol = Math.max(1e-6, 10 * termCritParam);
		StringBuilder bad = new StringBuilder();
		double maxDiff = 0.0;
		mainLog.println("\nChecking the synthesised strategy (values achieved in the initial state):");
		for (int c = 0; c < computed.length; c++) {
			mainLog.println("Coalition " + coalitions.get(c) + ": achieved " + achieved[c] + " (computed " + computed[c] + ")");
			double diff = (Double.isInfinite(computed[c]) || Double.isInfinite(achieved[c]))
					? (computed[c] == achieved[c] ? 0.0 : Double.POSITIVE_INFINITY) : Math.abs(achieved[c] - computed[c]);
			if (diff > tol * Math.max(1.0, Math.abs(computed[c]))) {
				bad.append(bad.length() > 0 ? ", " : "").append(coalitions.get(c));
				maxDiff = Math.max(maxDiff, diff);
			}
		}
		if (bad.length() > 0) {
			if (Double.isInfinite(maxDiff))
				mainLog.printWarning("The synthesised strategy does not achieve the computed values for coalition(s) " + bad
						+ ", so it may not be an equilibrium strategy");
			else
				mainLog.printWarning("The values achieved by the synthesised strategy differ from the computed ones for coalition(s) " + bad
						+ " (by up to " + maxDiff + "): either value iteration stopped before converging (try a smaller -epsilon)"
						+ " or the strategy does not realise the computed equilibrium");
		}
	}

	/**
	 * Stopping assumption for the unbounded objectives of an equilibria-based property: from every state, under all
	 * profiles, each unbounded objective is decided (its target reached, or a state reached from which it can no longer
	 * be satisfied) with probability 1. Bounded objectives are always decided by their bound. The check is optional
	 * (-eqassumptioncheck) since it can be expensive and is stronger than needed (some non-stopping games converge);
	 * in either case only warnings are issued.
	 * @param unbounded The indices of the unbounded objectives
	 */
	protected void checkStopping(CSG<Double> csg, BitSet[] targets, BitSet[] remain, BitSet unbounded, boolean rew) throws PrismException
	{
		if (unbounded.isEmpty())
			return;
		if (!assumptionCheck) {
			mainLog.printWarning("Equilibria computation assumes the game is stopping for this property (not checked; use -eqassumptioncheck)");
			return;
		}
		mainLog.println("Checking whether the game is stopping for this property...");
		int n = csg.getNumStates();
		StringBuilder failed = new StringBuilder();
		for (int i = unbounded.nextSetBit(0); i >= 0; i = unbounded.nextSetBit(i + 1)) {
			BitSet rem = (!rew && remain != null) ? remain[i] : null;
			// decided: target reached, or the objective can no longer be satisfied under any profile
			BitSet decided = mdpmc.prob0((MDP<Double>) csg, rem, targets[i], false, null);
			decided.or(targets[i]);
			if (mdpmc.prob1((MDP<Double>) csg, null, decided, true, null).cardinality() != n)
				failed.append(failed.length() > 0 ? ", " : "").append(i + 1);
		}
		if (failed.length() > 0)
			mainLog.printWarning("The game is not stopping for this property: objective(s) " + failed
					+ " not decided with probability 1 under all profiles, so the result may not correspond to an equilibrium");
	}

	/**
	 * 
	 * 
	 * @param csg
	 * @param coalitions
	 * @param rewards
	 * @param targets
	 * @param remain
	 * @param eqType
	 * @param crit
	 * @param min
	 * @return
	 * @throws PrismException
	 */
	public ModelCheckerResult computeReachEquilibria(CSG<Double> csg, List<Coalition> coalitions, List<CSGRewards<Double>> rewards, BitSet[] targets, BitSet[] remain, int eqType, int crit, boolean min) throws PrismException {
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
		double[] eq;
		double[] r = new double[csg.getNumStates()];
		//double[] sw;
		int i, j, k, p, s;
		boolean done, rew;
		long timePrecomp;
				
		// player -> iteration -> state -> indexes -> value
		if (genStrat) {
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
		mainLog.println();
		findMaxRowsCols(csg);
		
		mainLog.println("Starting equilibria computation (solver=" + setSolver(eqType) + ")...");
		BitSet unboundedObjs = new BitSet();
		unboundedObjs.set(0, targets.length);
		checkStopping(csg, targets, remain, unboundedObjs, rew);

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
		while (true) {
			for (s = 0; s < csg.getNumStates(); s++) {
				if (!known.get(s)) {
					if (genStrat) {
						sstrat = new ArrayList<List<Map<BitSet, Double>>>();
						mmap.clear();
			    		for (p = 0; p < 2; p++) {
			        		mmap.add(p, new HashMap<Integer, BitSet>());
			        	}
					}
					eq = stepEquilibriaTwoPlayer(csg, rewards, mmap, sstrat, sol, s, eqType, crit, rew, min);
					val[0][s] = eq[1];
					val[1][s] = eq[2];
					// player -> iteration -> state -> indexes -> value
					if (genStrat) {
						switch (eqType) {
							case CORR: {
								if (lstrat.get(0).get(0).get(s) == null) {
									lstrat.get(0).get(0).set(s, sstrat.get(0).get(0));
								}
								else if (!lstrat.get(0).get(0).get(s).equals(sstrat.get(0).get(0)) && checkEquilibriumChange(sol, eq, s)) {
									lstrat.get(0).get(0).set(s, sstrat.get(0).get(0));
								}
								break;
							}
							default: {
								for (p = 0; p < coalitions.size(); p++) {
									if (lstrat.get(p).get(0).get(s) == null) {
										lstrat.get(p).get(0).set(s, sstrat.get(0).get(p));
									}
									else if (!lstrat.get(p).get(0).get(s).equals(sstrat.get(0).get(p)) && checkEquilibriumChange(sol, eq, s)) {
										lstrat.get(p).get(0).set(s, sstrat.get(0).get(p));
									}
								}
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
			
			String sols;
			sols = "(";
			for (p = 0; p < numCoalitions; p++) {
				if (p < numCoalitions - 1)
					sols += sol[p][csg.getFirstInitialState()] + ",";
				else
					sols += sol[p][csg.getFirstInitialState()] + ")";
			}
			mainLog.println(k + ": " + sols);
			
			done = done & PrismUtils.doublesAreClose(sol[0], tmp[0], termCritParam, termCrit == TermCrit.ABSOLUTE);
			done = done & PrismUtils.doublesAreClose(sol[1], tmp[1], termCritParam, termCrit == TermCrit.ABSOLUTE);
			if (done) {
				break;
			}
			else if (!done && k == maxIters) {
				String msg = "Value iteration did not converge within " + k + " iterations";
				if (errorOnNonConverge)
					throw new PrismException(msg);
				mainLog.printWarning(msg + "; the values are those of the last iteration");
				break;
			}
			else {
				done = true;
				tmp[0] = Arrays.copyOf(sol[0], sol[0].length);
				tmp[1] = Arrays.copyOf(sol[1], sol[1].length);
			}
			k++;
		}
		if (done)
			mainLog.println("\nValue iteration converged after " + k + " iterations.");
		mainLog.println("\nPrecomputation took " + timePrecomp / 1000.0 + " seconds.");
		mainLog.println("Coalition results (initial state): (" + sol[0][csg.getFirstInitialState()] + "," + sol[1][csg.getFirstInitialState()] + ")");
		res.soln = r;

		if (genStrat) {
			CSGEquilibriumStrategy strat = twoPlayerStrategy(csg, lstrat, obj, targets, remain, rew, eqType == CORR);
			res.strat = strat;
			DTMCModelChecker dtmcmc = new DTMCModelChecker(this);
			dtmcmc.inheritSettings(this);
			dtmcmc.setSilentPrecomputations(true);
			dtmcmc.setLog(new PrismDevNullLog());
			double[] computed = { sol[0][csg.getFirstInitialState()], sol[1][csg.getFirstInitialState()] };
			checkStrategyValues(coalitions, computed, strat.achievedValues(rewards, dtmcmc));
		}
		res.numIters = k;
		return res;		
	}
	
	/**
	 * The equilibrium computed for a stage game, as returned by the step methods: {sum, payoff of coalition 0, ...},
	 * with the signs restored for min (the stage games are solved with negated payoffs). The stage solvers return
	 * the selected (SW or SF) equilibrium only, so eqs has a single row (and strats, if given, a single strategy).
	 * 
	 * @param eqs The equilibrium payoffs (one row)
	 * @param strats Strategy computed for the stage game (or null)
	 * @param min If minimising
	 */
	public double[] toEquilibrium(double[][] eqs, List<List<Map<BitSet, Double>>> strats, boolean min) throws PrismException {
		if (eqs.length != 1)
			throw new PrismException("Expected a single equilibrium from the stage game (got " + eqs.length + ")");
		double[] eq = new double[numCoalitions + 1];
		for (int c = 0; c < numCoalitions; c++) {
			eq[c + 1] = eqs[0][c];
			eq[0] += eqs[0][c];
		}
		if (strats != null && strats.size() > 1)
			strats.subList(1, strats.size()).clear();
		if (min) {
			for (int i = 0; i < eq.length; i++)
				eq[i] = -1.0 * eq[i];
		}
		return eq;
	}

	/**
	 *
	 *
	 * @param csg
	 * @param rewards
	 * @param mmap
	 * @param strats
	 * @param val
	 * @param s
	 * @param min
	 * @param crit
	 * @return
	 * @throws PrismException
	 */
	public double[] stepCorrelatedEquilibria(CSG<Double> csg, List<CSGRewards<Double>> rewards, List<Map<Integer, BitSet>> mmap, List<List<Map<BitSet, Double>>> strats, 
											 BitSet D, BitSet E, double[][] val, int s, boolean min, int crit) throws PrismException {
		EquilibriumResult result;
		ArrayList<Map<BitSet, Double>> eqstrat = null;
		BitSet idx = null, tmp = null;
		double[] eqs = new double[numCoalitions+1];
		int c, i;
		buildStepGame(csg, rewards, mmap, D, E, val, s, min);
		// for coalition c and action position q: the others' actions -> payoff of c (the CE incentive constraints)
		ArrayList<ArrayList<HashMap<BitSet, Double>>> ceConstraints = buildPayoffTables();
		if (genStrat) {
			eqstrat = new ArrayList<Map<BitSet, Double>>();
			eqstrat.add(new HashMap<BitSet, Double>());
			tmp = new BitSet();
		}
		if (utilities.size() == 1) {
			for (BitSet e : utilities.keySet()) {
				for (c = 0; c < numCoalitions; c++) {
					eqs[0] += utilities.get(e).get(c);
					eqs[c+1] = utilities.get(e).get(c);
				}
				if (genStrat) {
					idx = new BitSet();
					for (c = 0; c < numCoalitions; c++) {
						tmp.clear();
						tmp.or(psupports.get(c));
						tmp.and(e);
						i = tmp.nextSetBit(0);
						idx.or(mmap.get(c).get(strategies.get(c).indexOf(i)));
					}
					eqstrat.get(0).put(idx, 1.0);
				}
			}
			if (genStrat) {
				strats.add(0, eqstrat);
			}
		}
		else {
			// coalitions that are done (D) or have failed (E) are not part of the objectives
			BitSet active = new BitSet();
			active.set(0, numCoalitions);
			if (D != null)
				active.andNot(D);
			if (E != null)
				active.andNot(E);
			ceSolver.setActiveCoalitions(active);
			result = ceSolver.computeEquilibrium(utilities, ceConstraints, strategies, ceVarMap, crit);
			if (result.getStatus() == CSGResultStatus.SAT) {
				eqs[0] = 0.0;
				for (Double d : result.getPayoffVector()) {
					eqs[0] += d;
				}
				for (c = 0; c < numCoalitions; c++) {
					eqs[c+1] = result.getPayoffVector().get(c);
				}
				for (BitSet e : ceVarMap.keySet()) {
					if (genStrat) {
						idx = new BitSet();
						for (c = 0; c < numCoalitions; c++) {
							tmp.clear();
							tmp.or(psupports.get(c));
							tmp.and(e);
							i = tmp.nextSetBit(0);
							idx.or(mmap.get(c).get(strategies.get(c).indexOf(i)));
						}
						if (Double.compare(result.getStrategy().get(0).get(ceVarMap.get(e)), 0.0) > 0)
							eqstrat.get(0).put(idx, result.getStrategy().get(0).get(ceVarMap.get(e)));
					}
				}
				if (genStrat) {
					strats.add(0, eqstrat);
				}
			}
			else {
				throw new PrismException(ceSolver.getSolverName() + " could not find an optimal solution for state " + s);
			}
		}
		if (rewards != null) {
			addStateRewards(eqs, rewards, s, min);
			// keep eqs[0] = sum of the payoffs
			eqs[0] = 0.0;
			for (i = 1; i < eqs.length; i++)
				eqs[0] += eqs[i];
		}
		if (min) {
			for (i = 0; i < eqs.length; i++)
				eqs[i] = -1.0 * eqs[i];
		}
		return eqs;
	}
	
	
	/**
	 * 
	 * 
	 * @param csg
	 * @param rewards
	 * @param mmap
	 * @param strats
	 * @param val
	 * @param s
	 * @param min
	 * @return
	 * @throws PrismException
	 */
	public double[][] stepEquilibriaMulti(CSG<Double> csg, List<CSGRewards<Double>> rewards, List<Map<Integer, BitSet>> mmap, List<List<Map<BitSet, Double>>> strats, 
									 	  BitSet D, BitSet E, double[][] val, int s, boolean min, int crit) throws PrismException {
		double[][] result;
		BitSet active;
		int c, q;

		buildStepGame(csg, rewards, mmap, D, E, val, s, min);
		active = csg.getConcurrentPlayers(s);

		// Only one joint action
		if (utilities.size() == 1) {
			BitSet joint = utilities.keySet().iterator().next();
			result = new double[1][numCoalitions];
			for (c = 0; c < numCoalitions; c++)
				result[0][c] = utilities.get(joint).get(c);
			if (genStrat) {
				ArrayList<Map<BitSet, Double>> eqstrat = new ArrayList<Map<BitSet, Double>>();
				extractStrategyFromSupport(mmap, eqstrat, joint);
				strats.add(eqstrat);
			}
			return result;
		}
		// Only one player has a choice (SW: best for that player, ties broken by the sum; SF goes to the solver)
		if (active.cardinality() == 1 && crit != FAIR) {
			return findSWNEOnePlayer(mmap, strats, null, active);
		}

		// Stage game: coalition actions indexed as in strategies/mmap; coalitions that are done (D) or have failed (E) are free
		int[] numActions = new int[numCoalitions];
		for (c = 0; c < numCoalitions; c++)
			numActions[c] = strategies.get(c).size();
		StageGame<Double> game = new StageGame<>(numActions);
		int[] joint = new int[numCoalitions];
		for (Entry<BitSet, ArrayList<Double>> e : utilities.entrySet()) {
			for (c = 0; c < numCoalitions; c++) {
				joint[c] = -1;
				for (q = 0; q < numActions[c]; q++) {
					if (e.getKey().get(strategies.get(c).get(q))) {
						joint[c] = q;
						break;
					}
				}
				if (joint[c] == -1)
					throw new PrismException("Error in building the stage game for state " + s);
			}
			for (c = 0; c < numCoalitions; c++)
				game.setPayoff(c, joint, e.getValue().get(c));
		}
		for (c = 0; c < numCoalitions; c++)
			game.setFree(c, (D != null && D.get(c)) || (E != null && E.get(c)));

		if (multiSolver == null) {
			try {
				multiSolver = new StageGameSolverScip();
			} catch (PrismException e) {
				throw new PrismException("SCIP is required for Nash equilibria with more than two coalitions. " + e.getMessage());
			}
			multiSolver.setFairnessWelfareTieBreak(true); // SF: gap, then sum, then each coalition (as for two-player Nash)
		}
		StageGameResult<Double> res = multiSolver.solve(game, StageGameSolver.Concept.NASH,
				crit == FAIR ? StageGameSolver.Criterion.SOCIAL_FAIRNESS : StageGameSolver.Criterion.SOCIAL_WELFARE, false);
		if (res.getStatus() != StageGameResult.Status.OPTIMAL)
			throw new PrismException(multiSolver.getSolverName() + " could not find an optimal equilibrium for state " + s + " (" + res.getStatus()
					+ (res.getMessage() != null ? ": " + res.getMessage() : "") + ")");
		result = new double[1][numCoalitions];
		for (c = 0; c < numCoalitions; c++)
			result[0][c] = res.getPayoffs().get(c);
		if (genStrat) {
			ArrayList<Map<BitSet, Double>> eqstrat = new ArrayList<Map<BitSet, Double>>();
			for (c = 0; c < numCoalitions; c++) {
				eqstrat.add(c, new HashMap<BitSet, Double>());
				List<Double> x = res.getStrategies().get(c);
				for (q = 0; q < x.size(); q++)
					if (x.get(q) > 0.0)
						eqstrat.get(c).put(mmap.get(c).get(q), x.get(q));
			}
			strats.add(eqstrat);
		}
		return result;
	}
	
	/**
	 * Returns the equilibrium (array of values) and updates strategies for the two-player case.
	 * 
	 * @param csg
	 * @param rewards
	 * @param mmap
	 * @param strats
	 * @param val
	 * @param s
	 * @param eqType
	 * @param crit
	 * @param rew
	 * @param min
	 * @return
	 * @throws PrismException
	 */
	public double[] stepEquilibriaTwoPlayer(CSG<Double> csg, List<CSGRewards<Double>> rewards, List<Map<Integer, BitSet>> mmap, List<List<Map<BitSet, Double>>> strats,
			 								double[][] val, int s, int eqType, int crit, boolean rew, boolean min) throws PrismException {
		double[][] equilibria;
		double[] equilibrium;
		
		switch (eqType) {
			case CORR : {
				if (rew) {
					equilibrium = stepCorrelatedEquilibria(csg, rewards, mmap, strats, null, null, val, s, min, crit);
				}
				else 
					equilibrium = stepCorrelatedEquilibria(csg, null, mmap, strats, null, null, val, s, min, crit);
				break;
			}
			default : {

				if (rew) {
					equilibria = stepNashEquilibria(csg, rewards.get(0), rewards.get(1), mmap, strats, val, s, min, crit);
				}
				else {
					equilibria = stepNashEquilibria(csg, null, null, mmap, strats, val, s, min, crit);
				}
				equilibrium = toEquilibrium(equilibria, strats, min);
			}
		}
		return equilibrium;
	}
	
	/**
	 * Computes Nash equilibria for a bimatrix game.
	 * 
	 * @param csg
	 * @param csgRewards1
	 * @param csgRewards2
	 * @param mmap
	 * @param strats
	 * @param val
	 * @param s
	 * @param min
	 * @return
	 * @throws PrismException
	 */
	public double[][] stepNashEquilibria(CSG<Double> csg, CSGRewards<Double> csgRewards1, CSGRewards<Double> csgRewards2, List<Map<Integer, BitSet>> mmap,
									 	 List<List<Map<BitSet, Double>>> strats, double[][] val, int s, boolean min, int crit) throws PrismException {
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
		bmgame = buildBimatrixGame(csg, csgRewards1, csgRewards2, mmap, nmap, val, s, min);	
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
				if((equalA || equalB) && crit != FAIR) { // if all entries of one of them are the same (SW: pure, max of the other)
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
					if (genStrat) {
						eqstrat = new ArrayList<Map<BitSet, Double>>();
						eqstrat.add(0, new HashMap<BitSet, Double>());
						eqstrat.get(0).put(mmap.get(0).get(nmap.get(0).get(mrow)), 1.0);
						eqstrat.add(1, new HashMap<BitSet, Double>());
						eqstrat.get(1).put(mmap.get(1).get(nmap.get(1).get(mcol)), 1.0);
						strats.add(0, eqstrat);
					}
					addStateRewards(result, csgRewards1, csgRewards2, s, min);
				}
				else { // both players have choices and matrices are not trivial (or SF), call solver
					result = solveNash(nrows, ncols, val1s, val2s, crit, mmap, nmap, strats);
					addStateRewards(result, csgRewards1, csgRewards2, s, min);
				}
			}
			else { // all entries in both are the same
				result = new double[1][2];
				result[0][0] = ent1;
				result[0][1] = ent2;
				if (genStrat) {
					eqstrat = new ArrayList<Map<BitSet, Double>>();
					eqstrat.add(0, new HashMap<BitSet, Double>());
					eqstrat.get(0).put(mmap.get(0).get(nmap.get(0).get(0)), 1.0);
					eqstrat.add(1, new HashMap<BitSet, Double>());
					eqstrat.get(1).put(mmap.get(1).get(nmap.get(1).get(0)), 1.0);
					strats.add(0, eqstrat);
				}
				addStateRewards(result, csgRewards1, csgRewards2, s, min);
			}
		} 
		else if (crit == FAIR && nrows * ncols > 1) { // just one of the players has choices, SF: the other may be mixed
			for (int r = 0; r < nrows; r++) {
				for (int c = 0; c < ncols; c++) {
					val1s[r][c] = bmgame.get(0).get(r).get(c);
					val2s[r][c] = bmgame.get(1).get(r).get(c);
				}
			}
			result = solveNash(nrows, ncols, val1s, val2s, crit, mmap, nmap, strats);
			addStateRewards(result, csgRewards1, csgRewards2, s, min);
		}
		else { // just one of the players has choices
			result = new double[1][2];
			double vt1, vt2, sumv, sumt;
			if(genStrat) {
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
						if(genStrat) {
							d1.clear();
							d1.put(mmap.get(0).get(nmap.get(0).get(r)), 1.0);
						}
						val2 = vt2;
						val1 = vt1;
						sumv = val1 + val2;
					}
				}
				if(genStrat)
					d2.put(mmap.get(1).get(nmap.get(1).get(0)), 1.0);
			} 
			else if (nrows == 1 && ncols > 1) {
				for (int c = 0; c < ncols; c++) {
					vt1 = bmgame.get(0).get(0).get(c);
					vt2 = bmgame.get(1).get(0).get(c);
					sumt = vt1 + vt2;
					if (Double.compare(vt2, val2) > 0 || (Double.compare(vt2, val2) == 0 && Double.compare(sumt, sumv) > 0)) {
						if(genStrat) {
							d2.clear();
							d2.put(mmap.get(1).get(nmap.get(1).get(c)), 1.0);
						}
						val2 = vt2;
						val1 = vt1;
						sumv = val1 + val2;
					}
				}
				if(genStrat)
					d1.put(mmap.get(0).get(nmap.get(0).get(0)), 1.0);
			} 
			else if (nrows == 1 && ncols == 1) {
				val1 = bmgame.get(0).get(0).get(0);
				val2 = bmgame.get(1).get(0).get(0);
				if(genStrat) {
					d1.put(mmap.get(0).get(nmap.get(0).get(0)), 1.0);
					d2.put(mmap.get(1).get(nmap.get(1).get(0)), 1.0);
				}
			} 
			else {
				throw new PrismException("Error with matrix rank");
			}
			if (genStrat) {
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

	/**
	 * Optimal (SW or SF) Nash equilibrium of the bimatrix game (a, b) with the SMT solver; returns its values as a
	 * one-row array (as expected by toEquilibrium()) and, if strategies are generated, adds its strategy to strats.
	 */
	private double[][] solveNash(int nrows, int ncols, double[][] a, double[][] b, int crit, List<Map<Integer, BitSet>> mmap,
								 ArrayList<ArrayList<Integer>> nmap, List<List<Map<BitSet, Double>>> strats) throws PrismException {
		nashSolver.compute(nrows, ncols, a, b, crit);
		double[][] result = new double[1][2];
		result[0][0] = nashSolver.getPayoffs()[0];
		result[0][1] = nashSolver.getPayoffs()[1];
		if (genStrat) {
			ArrayList<Map<BitSet, Double>> eqstrat = new ArrayList<Map<BitSet, Double>>();
			for (int p = 0; p < 2; p++) {
				eqstrat.add(p, new HashMap<BitSet, Double>());
				Distribution<Double> d = nashSolver.getStrategies().get(p);
				for (int t : d.getSupport()) {
					eqstrat.get(p).put(mmap.get(p).get(nmap.get(p).get(t)), d.get(t));
				}
			}
			strats.add(0, eqstrat);
		}
		return result;
	}
}
