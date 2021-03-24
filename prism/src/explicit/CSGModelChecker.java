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

import java.io.FileOutputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.math.BigDecimal;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Vector;
import java.util.stream.Collectors;

import org.apache.commons.math3.util.Precision;

import acceptance.AcceptanceOmega;
import acceptance.AcceptanceOmega.LiftBitSet;
import acceptance.AcceptanceRabin.RabinPair;
import acceptance.AcceptanceRabin;
import acceptance.AcceptanceReach;
import acceptance.AcceptanceType;
import automata.DA;
import automata.DASimplifyAcceptance;
import automata.LTL2DA;
import automata.LTL2NBA;
import automata.LTL2RabinLibrary;
import explicit.LTLModelChecker.LTLProduct;
import explicit.ProbModelChecker.TermCrit;
import explicit.rewards.CSGRewards;
import explicit.rewards.CSGRewardsSimple;
import explicit.rewards.MDPRewards;
import explicit.rewards.MDPRewardsSimple;
import explicit.rewards.StateRewardsConstant;
import lpsolve.LpSolve;
import lpsolve.LpSolveException;
import parser.State;
import parser.VarList;
import parser.ast.Coalition;
import parser.ast.Declaration;
import parser.ast.DeclarationInt;
import parser.ast.Expression;
import parser.ast.ExpressionTemporal;
import prism.IntegerBound;
import prism.PrismComponent;
import prism.PrismException;
import prism.PrismFileLog;
import prism.PrismLangException;
import prism.PrismLog;
import prism.PrismUtils;
import prism.PrismSettings;
import strat.CSGStrategy;
import strat.CSGStrategy.CSGStrategyType;
import strat.MDStrategyArray;
import strat.MemorylessDeterministicStrategy;
import strat.Strategy;

public class CSGModelChecker extends ProbModelChecker {

	public static final int R_INFINITY = 0;
	public static final int R_CUMULATIVE = 1;
	public static final int R_ZERO = 2;
	
	protected HashMap<BitSet, ArrayList<Double>> utilities;
	protected HashMap<BitSet, ArrayList<Distribution>> probabilities;
	protected ArrayList<ArrayList<Integer>> strategies; 
	protected ArrayList<ArrayList<String>> actions;

	protected HashMap<Integer, Distribution> adv; // probably shouldn't be global or be here at all
	
	protected BitSet[] coalitionIndexes;
	protected BitSet[] actionIndexes;
	
	protected double[] avgNumActions;
	
	protected double minEntry; 
	protected double scaleFactor = getSettings().getDouble(PrismSettings.PRISM_ZS_LP_SCALE_FACTOR);
	
	protected int numCoalitions;
	protected int numPlayers;
	protected int varIndex;
	
	protected int maxRows;
	protected int maxCols;
	protected long timerVal;	 
	
	protected boolean allEqual;

	public CSGModelChecker(PrismComponent parent) throws PrismException {
		super(parent);
		utilities = new HashMap<BitSet, ArrayList<Double>>();
		probabilities = new HashMap<BitSet, ArrayList<Distribution>>();
		actions = new ArrayList<ArrayList<String>>();
		strategies = new ArrayList<ArrayList<Integer>>();
	}
		
	public ModelCheckerResult computeReachRewards(CSG csg, Coalition coalition, CSGRewards rewards, BitSet target, boolean min1, boolean min2, int unreachingSemantics) throws PrismException {
		switch (unreachingSemantics) {
		case R_INFINITY:
			return computeReachRewardsInfinity(csg, coalition, rewards, target, min1, min2);
		case R_CUMULATIVE:
			return computeReachRewardsCumulative(csg, coalition, rewards, target, min1, min2, false);
		case R_ZERO:
			throw new PrismException("F0 is not yet supported for CSGs.");
		default:
			throw new PrismException("Unknown semantics for runs unreaching the target in CSGModelChecker: " + unreachingSemantics);
		}
	}
	
	public BitSet A(ArrayList<ArrayList<Distribution>> mdist, BitSet v2, BitSet y) {
		BitSet result = new BitSet();
		int col, row;
		boolean b;
		for (row = 0; row < mdist.size(); row++) {
			b = true;
			for (col = 0; col < mdist.get(row).size(); col++) {
				b = b && (mdist.get(row).get(col).isSubsetOf(y) || v2.get(col));
			}
			if (b) {
				result.set(row);
			}
		}
		return result;
	}
	
	public BitSet B(ArrayList<ArrayList<Distribution>> mdist, BitSet v1, BitSet x) {
		BitSet result = new BitSet();
		int col, row;
		boolean b;
		for (col = 0; col < mdist.get(0).size(); col++) {
			b = false;
			for (row = v1.nextSetBit(0); row >= 0; row = v1.nextSetBit(row + 1)) {
				b = mdist.get(row).get(col).containsOneOf(x);
				if (b) {
					result.set(col);
					break;
				}
			}
		}
		return result;
	}
	
	public BitSet apreXY(CSG csg, BitSet x, BitSet y) throws PrismException {
		ArrayList<ArrayList<Distribution>> mdist;
		BitSet result = new BitSet();
		int s;
		for (s = 0; s < csg.getNumStates(); s++) {
			mdist = buildMatrixDist(csg, s);
			result.set(s, B(mdist, A(mdist, new BitSet(), y), x).cardinality() == mdist.get(0).size());
		}
		return result;
	}
	
	public BitSet AF(CSG csg, BitSet b) throws PrismException {
		int n = csg.getNumStates();
		BitSet x, y, sol1;
		x = new BitSet();
		y = new BitSet();
		y.set(0, n);
		sol1 = new BitSet();
		boolean done_x, done_y;
		done_y = false;
		while (!done_y) {
			done_x = false;
			x.clear();
			while (!done_x) {
				sol1.clear();
				sol1.or(apreXY(csg, x, y));
				sol1.or(b);
				done_x = x.equals(sol1);
				x.clear();
				x.or(sol1);
			}
			done_y = y.equals(x);
			y.clear();
			y.or(x);
		}
		return y;
	}

	public boolean apreXYZ(ArrayList<ArrayList<Distribution>> mdist, BitSet x, BitSet y, BitSet z) {
		BitSet a, ab, e, v;
		e = new BitSet();
		a = new BitSet();
		ab = new BitSet();
		v = new BitSet();
		v.set(0, mdist.size());
		boolean done_v = false;
		while (!done_v) {
			a.clear();
			ab.clear();
			a.or(A(mdist, e, z));
			ab.or(A(mdist, B(mdist, v, x), y));
			a.and(ab);
			done_v = v.equals(a);
			v.clear();
			v.or(a);
		}
		return !(v.isEmpty());
	}
	
	public BitSet apreXYZ(CSG csg, BitSet x, BitSet y, BitSet z) throws PrismException {
		ArrayList<ArrayList<Distribution>> mdist;
		BitSet result = new BitSet();
		for (int s = 0; s < csg.getNumStates(); s++) {
			mdist = buildMatrixDist(csg, s);
			result.set(s, apreXYZ(mdist, x, y ,z));
		}
		return result;
	}
	
	public BitSet AFG(CSG csg, BitSet b) throws PrismException {
		int n = csg.getNumStates();
		BitSet x, y, z, sol1, sol2;
		x = new BitSet();
		y = new BitSet();
		z = new BitSet();
		sol1 = new BitSet();
		sol2 = new BitSet();
		boolean done_x, done_y, done_z;
		done_z = false;
		z.set(0, n);
		while (!done_z) {
			done_x = false;
			x.clear();
			while (!done_x) {
				done_y = false;
				y.clear();
				y.set(0, n);
				while (!done_y) {
					sol1.clear();
					sol2.clear();
					sol1 = apreXYZ(csg, x, y, z);
					sol1.and(b);
					sol2.or(b);
					sol2.flip(0, n);
					sol2.and(apreXY(csg, x, z));
					sol1.or(sol2);
					done_y = y.equals(sol1);
					y.clear();
					y.or(sol1);
				}
				done_x = x.equals(y);
				x.clear();
				x.or(y);
			}
			done_z = z.equals(x);
			z.clear();
			z.or(x);
		}
		return z;
	}
	
	public boolean pre1(ArrayList<ArrayList<Distribution>> mdist, BitSet x) {
		int row, col;
		boolean b = false;
		for (row = 0; row < mdist.size(); row++) {
			b = true;
			for (col = 0; col < mdist.get(row).size(); col++) {
				b = b && mdist.get(row).get(col).isSubsetOf(x);
			}
			if (b) 
				return true;
		}
		return b;
	}
	
	public void pre1(CSG csg, BitSet x, BitSet sol) throws PrismException {
		ArrayList<ArrayList<Distribution>> mdist;
		for (int s = 0; s < csg.getNumStates(); s++) {
			mdist = buildMatrixDist(csg, s);
			sol.set(s, pre1(mdist, x));
		}
	}
	
	public BitSet G(CSG csg, BitSet b) throws PrismException {
		int n = csg.getNumStates();
		BitSet sol1, x;
		sol1 = new BitSet();
		x = new BitSet();
		boolean done_x = false;
		x.set(0, n);
		while (!done_x) {
			sol1.clear();
			pre1(csg, x, sol1);
			sol1.and(b);
			done_x = x.equals(sol1);
			x.clear();
			x.or(sol1);
		}
		return x;
	}
	
	public void buildCoalitions(CSG csg, Coalition coalition, boolean min) throws PrismException {
		if (coalition == null || coalition.isEmpty())
			throw new PrismException("Coalitions must not be empty");
		int p;
		numPlayers = csg.getNumPlayers();
		coalitionIndexes = new BitSet[2];
		actionIndexes = new BitSet[2];
		Map<Integer, String> pmap = new HashMap<Integer, String>();
		for (p = 0; p < numPlayers; p++) {
			pmap.put(p + 1, csg.getPlayerName(p));
		}
		for (p = 0; p < 2; p++) {
			coalitionIndexes[p] = new BitSet();
			actionIndexes[p]= new BitSet();
		}
		for (p = 0; p < numPlayers; p++) {
			if (min) {
				if(coalition.isPlayerIndexInCoalition(p + 1, pmap)) {
					coalitionIndexes[1].set(p);
					actionIndexes[1].or(csg.getIndexes()[p]);
					actionIndexes[1].set(csg.getIdleForPlayer(p));
				}
				else {
					coalitionIndexes[0].set(p);
					actionIndexes[0].or(csg.getIndexes()[p]);
					actionIndexes[0].set(csg.getIdleForPlayer(p));
				}
			}
			else {
				if(coalition.isPlayerIndexInCoalition(p + 1, pmap)) {
					coalitionIndexes[0].set(p);
					actionIndexes[0].or(csg.getIndexes()[p]);
					actionIndexes[0].set(csg.getIdleForPlayer(p));
				}
				else {
					coalitionIndexes[1].set(p);
					actionIndexes[1].or(csg.getIndexes()[p]);
					actionIndexes[1].set(csg.getIdleForPlayer(p));
				}
			}
		}
		numCoalitions = 2;
		findMaxRowsCols(csg);
	}

	public void findMaxRowsCols(CSG csg) {
		int p, mc, mr, s;
		maxRows = 0;
		maxCols = 0;
		avgNumActions = new double[numCoalitions];
		Arrays.fill(avgNumActions, 0.0);
		for (s = 0; s < csg.getNumStates(); s++) {
			mc = 1;
			mr = 1;
			for (p = coalitionIndexes[0].nextSetBit(0); p >= 0; p = coalitionIndexes[0].nextSetBit(p + 1)) {
				mr *= csg.getIndexesForPlayer(s, p).cardinality();
			}
			avgNumActions[0] += mr;
			for (p = coalitionIndexes[1].nextSetBit(0); p >= 0; p = coalitionIndexes[1].nextSetBit(p + 1)) {
				mc *= csg.getIndexesForPlayer(s, p).cardinality();		
			}
			avgNumActions[1] += mc;
			maxRows = (maxRows < mr)? mr : maxRows;
			maxCols = (maxCols < mc)? mc : maxCols;
		}
		avgNumActions[0] /= csg.getNumStates();
		avgNumActions[1] /= csg.getNumStates();
		mainLog.println("Max/avg (actions): " + "(" + maxRows + ","  + maxCols + ")/(" + PrismUtils.formatDouble2dp(avgNumActions[0]) + "," + PrismUtils.formatDouble2dp(avgNumActions[1]) + ")");
	}
	
	public void buildStepGame(CSG csg, List<CSGRewards> rewards, Map<BitSet, Integer> imap, double[] val, int s) throws PrismException {
		BitSet jidx;
		BitSet indexes = new BitSet();
		BitSet tmp = new BitSet();
		String act;
		double u, v;
		int c, i, p, t;
		int[] joint;
		actions.clear();
		strategies.clear();
		utilities.clear();
		probabilities.clear();
		varIndex = 0;
		minEntry = Double.POSITIVE_INFINITY;
		allEqual = true;
		u = Double.NaN;
		if (imap != null)
			imap.clear();
		else
			imap = new HashMap<BitSet, Integer>();
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
			for (c = 0; c < 2; c++) {
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
			probabilities.put(jidx, new ArrayList<Distribution>());
			v = 0.0;
			if (val != null) {
				for (int d : csg.getChoice(s, t).getSupport()) 
					v += csg.getChoice(s, t).get(d) * val[d];
			}
			if (rewards != null)
				v += rewards.get(0).getTransitionReward(s, t);
			if (u != Double.NaN)
				allEqual = allEqual && Double.compare(u, v) == 0;
			utilities.get(jidx).add(0, v);
			probabilities.get(jidx).add(0, csg.getChoice(s, t));
			minEntry = (minEntry > v)? v : minEntry;
			u = v;
		}		
		//System.out.println("## state " + csg.getStatesList().get(s));
		//System.out.println("-- actions " + actions);
		//System.out.println("-- strategies " + strategies);		
		//System.out.println("-- utilities " + utilities);
		//System.out.println("-- imap " + imap);
	}
	
	public ArrayList<ArrayList<Double>> buildMatrixGame(CSG csg, CSGRewards r, Map<Integer, BitSet> mmap, double[] val,  int s, boolean min) throws PrismException {
		ArrayList<ArrayList<Double>> mgame = new ArrayList<ArrayList<Double>>();
		ArrayList<CSGRewards> rewards = null;
		Map<BitSet, Integer> imap = new HashMap<BitSet, Integer>();
		Map<Integer, BitSet> rmap;
		BitSet action = new BitSet();
		int col, row;
		if (r != null) {
			rewards = new ArrayList<CSGRewards>();
			rewards.add(0, r);
		}
        buildStepGame(csg, rewards, imap, val, s);
        rmap = imap.entrySet().stream().collect(Collectors.toMap(Map.Entry::getValue, Map.Entry::getKey));
    	if (mmap != null) {
    		if (min) {
        		for (col = 0; col < strategies.get(1).size(); col++) {
            		mmap.put(col, rmap.get(strategies.get(1).get(col)));
        		}
    		}
    		else {
    			for (row = 0; row < strategies.get(0).size(); row++) {
            		mmap.put(row, rmap.get(strategies.get(0).get(row)));
    			}        
    		}
    	}
        for (row = 0; row < strategies.get(0).size(); row++) {
        	mgame.add(row, new ArrayList<Double>());
        	action.clear();
        	action.set(strategies.get(0).get(row));
        	for (col = 0; col < strategies.get(1).size(); col++) {
            	action.set(strategies.get(1).get(col));
    			if (utilities.containsKey(action)) {
    				mgame.get(row).add(col, utilities.get(action).get(0));
    			}
    			else 
    				throw new PrismException("Error in building matrix game");
    			action.clear(strategies.get(1).get(col));
            }
        }
		return mgame;
	}
	
	public ArrayList<ArrayList<Distribution>> buildMatrixDist(CSG csg, int s) throws PrismException {
		ArrayList<ArrayList<Distribution>> mdist = new ArrayList<ArrayList<Distribution>>();       
		BitSet action = new BitSet();
		int col, row;
        buildStepGame(csg, null, null, null, s);
        for (row = 0; row < strategies.get(0).size(); row++) {
        	mdist.add(row, new ArrayList<Distribution>());
        	action.clear();
        	action.set(strategies.get(0).get(row));
        	for (col = 0; col < strategies.get(1).size(); col++) {
            	action.set(strategies.get(1).get(col));
    			if (probabilities.containsKey(action))
    				mdist.get(row).add(col, probabilities.get(action).get(0));
    			else 
    				throw new PrismException("Error in building distribution matrix");
    			action.clear(strategies.get(1).get(col));
            }
        }
		return mdist;
	}
	
	public void buildLPLpsolve(LpSolve lp, ArrayList<ArrayList<Double>> mgame, boolean rew, boolean min) throws LpSolveException {	
		int nrows = (min)? mgame.get(0).size() : mgame.size(); // Number of rows
		int ncols = (min)? mgame.size() : mgame.get(0).size(); // Number of columns
		int[] vari = new int[nrows+1]; // Indexes of variables, should be m + 1 for an m x n matrix
		double[] row = new double[nrows+1];
		lp.setColName(1, "v"); 
		// Sets bounds for each p variable
		for (int i = 2; i <= nrows + 1; i++) {
			if (min)
				lp.setColName(i, actions.get(1).get(i - 2));
			else 
				lp.setColName(i, actions.get(0).get(i - 2));
			lp.setBounds(i, 0, 1.0); 
		}
		// Rewards mode
		if (rew) {
			lp.setBounds(1, -1.0 * lp.getInfinite(), lp.getInfinite());
		}
		else {
			lp.setBounds(1, 0, 1.0);
		}
		lp.setAddRowmode(true); 
		int k = 0;
		// Builds each column considering the n + 1 constraints as a matrix
		for (int j = 0; j < ncols; j++) {
			vari[k] = 1;
			row[k] = scaleFactor;
			for(int i = 0; i < nrows; i++) {
				k++;
				vari[k] = k+1;
				if (min)
					row[k] = -1.0 * scaleFactor * mgame.get(j).get(i); 
				else
					row[k] = -1.0 * scaleFactor * mgame.get(i).get(j); 
			}
			if (min)
				lp.addConstraintex(nrows+1, row, vari, LpSolve.GE, 0.0);
			else
				lp.addConstraintex(nrows+1, row, vari, LpSolve.LE, 0.0);
			k = 0;
		}
		// Sets name for each row (over ncols as they represent the constraints)
		for (k = 0; k < ncols; k++) {
			if (min)
				lp.setRowName(k + 1, actions.get(0).get(k));
			else
				lp.setRowName(k + 1, actions.get(1).get(k));
		
		}
		for (k = 0; k < nrows + 1; k++) {
			vari[k] = k + 1;
			row[k] = (k > 0)? 1 : 0;
		}
		lp.addConstraintex(k, row, vari, LpSolve.EQ, 1.0);
		lp.setAddRowmode(false);
		for (k = 0; k < nrows + 1; k++) {
			vari[k] = k + 1;
			row[k] = (k == 0)? 1 : 0;
		}
		lp.setObjFnex(k, row, vari);
		if (min)
			lp.setMinim();
		else
			lp.setMaxim();
		//lp.printLp();
	}
	
	public int valInfinity(ArrayList<ArrayList<Double>> mgame) {
		BitSet hasInf = new BitSet();
		int row, col;
		boolean infRow;
		for (row  = 0; row < mgame.size(); row++) {
			infRow = true;
			for (col = 0; col < mgame.get(0).size(); col++) {
				if (mgame.get(row).get(col) == Double.POSITIVE_INFINITY) {
					hasInf.set(col);
				}
				else {
					infRow = false;
				}
			}
			if (infRow) 
				return row;
		}
		if (!hasInf.isEmpty()) {
			ArrayList<ArrayList<Double>> mcopy = new ArrayList<ArrayList<Double>>(mgame);
			int ncol, nrow;
			nrow = 0;
			mgame.clear();
			for (row = 0; row < mcopy.size(); row++) {
				mgame.add(nrow, new ArrayList<Double>());
				ncol = 0;
				for (col = 0; col < mcopy.get(0).size(); col++) {
					if (!hasInf.get(col)) { 
						mgame.get(nrow).add(ncol, mcopy.get(row).get(col));
						ncol++;
					}
				}
				nrow++;
			}
		}
		return -1;
	}
	
	public double val(LpSolve lp, ArrayList<ArrayList<Double>> mgame, List<Map<BitSet, Double>> strat, Map<Integer, BitSet> rmap, int s, boolean rew, boolean min) throws PrismException {
		long timer = System.currentTimeMillis();
		int nrows = mgame.size(); // Number of rows
		int ncols = mgame.get(0).size(); // Number of columns
		double res = Double.NaN;
		Map<BitSet, Double> d = new HashMap<BitSet, Double> ();		
		if (allEqual) {
			if(generateStrategy || exportAdv) {
				d.put(rmap.get(0), 1.0);
				strat.set(s, d);
			}
			return minEntry;
		}
		else if (nrows == 1) {
			int srow = 0;
			res = Double.POSITIVE_INFINITY;
			for (int col = 0; col < ncols; col++) {
				if (res > mgame.get(0).get(col)) {
					res = mgame.get(0).get(col);
					srow = (min)? col : 0;
				}
			}
			if (generateStrategy || exportAdv) {
				d.put(rmap.get(srow), 1.0); // in case of min, rmap is actually "cmap"
				strat.set(s, d);
			}
			return res;
		}
		else if (ncols == 1) {
			int scol = 0;
			res = Double.NEGATIVE_INFINITY;
			for (int row = 0; row < nrows; row++) {
				if (res < mgame.get(row).get(0)) {
					res = mgame.get(row).get(0);
					scol = (min)? 0 : row;
				}
			}
			if(generateStrategy || exportAdv) {
				d.put(rmap.get(scol), 1.0);
				strat.set(s, d);
			}
			return res;
		}
		else {
			// System.out.println("$$ state " + s + " is concurrent");
			// Should add check for trivial games
			int infty;
			infty = valInfinity(mgame);
			if (infty != -1) {
				res = Double.POSITIVE_INFINITY;
				if (generateStrategy || exportAdv) {
					d.put(rmap.get(infty), 1.0);
				}
				return res;
			}
			else {
				try {
					buildLPLpsolve(lp, mgame, rew, min);
				} 
				catch (LpSolveException e1) {
					throw new PrismException("Exception raised by lpSolve when building linear program for state  " + s);
				}
				//lp.unscale();
				//lp.setPresolve(LpSolve.PRESOLVE_ROWS, lp.getPresolveloops());
				//lp.setPresolve(LpSolve.PRESOLVE_COLS, lp.getPresolveloops());
				//lp.setPresolve(LpSolve.PRESOLVE_ROWDOMINATE, lp.getPresolveloops());
				//lp.setPresolve(LpSolve.PRESOLVE_COLDOMINATE, lp.getPresolveloops());
				//lp.setPresolve(LpSolve.PRESOLVE_BOUNDS, lp.getPresolveloops());
				//lp.setPresolve(LpSolve.PRESOLVE_REDUCEGCD, lp.getPresolveloops());
				//lp.setScaling(LpSolve.SCALE_GEOMETRIC);
				//lp.setScaling(LpSolve.SCALE_POWER2);
				//lp.setScaling(LpSolve.SCALE_EQUILIBRATE);
				try {
					int status = lp.solve();
					if (status == LpSolve.OPTIMAL) {
						res = lp.getObjective();
						double[] values = (min)? new double[maxCols + 1] : new double[maxRows + 1];
						lp.getVariables(values);
						if (generateStrategy || exportAdv) {
							for (int row = 1; row <= nrows; row++) {
								if (values[row] > 0)
									d.put(rmap.get(row - 1), values[row]);
							}
							strat.set(s, d);
						}
					}
					else {
						throw new PrismException("lpSolve could not find an optimal solution for state " + s);
					}
				}
				catch (Exception e) {
					//e.printStackTrace();
					//lp.printLp();
					mainLog.println("Exception raised by lpSolve when computing value for state " + s + ". lpSolve status: " + lp.getStatustext(lp.getStatus()));
					mainLog.println("Rounding up entries...");
					ArrayList<ArrayList<Double>> mcopy = new ArrayList<ArrayList<Double>>(mgame);
					mgame.clear();
					for (int row = 0; row < nrows; row++) {
						mgame.add(row, new ArrayList<Double>());
						for (int col = 0; col < ncols; col++) {
							mgame.get(row).add(col, Precision.round(mcopy.get(row).get(col), 9, BigDecimal.ROUND_FLOOR));
						}	
					}
					mcopy.clear();
					try {
						if (min)
							lp.resizeLp(0, maxCols + 1);
						else
							lp.resizeLp(0, maxRows + 1);
						buildLPLpsolve(lp, mgame, rew, min);
						int status = lp.solve();
						if (status == LpSolve.OPTIMAL) {
							res = lp.getObjective();
							double[] values = (min)? new double[maxCols + 1] : new double[maxRows + 1];
							lp.getVariables(values);
							if (generateStrategy || exportAdv) {
								for (int row = 1; row <= nrows; row++) {
									if (values[row] > 0)
										d.put(rmap.get(row - 1), values[row]);
								}
								strat.set(s, d);
							}
						}
						else {
							throw new PrismException("Rounding up failed for state " + s + ". Failed to compute solution");	
						}
					} 
					catch (LpSolveException e2) {
						throw new PrismException("Rounding up failed for state " + s + ". Failed to compute solution");
					}
					return  res;
				}
			}
		}
		timer = System.currentTimeMillis() - timer;
		timerVal += timer;
		return res;
	}
	
	public ModelCheckerResult computeReachProbsValIter(CSG csg, BitSet no, BitSet yes, int limit, boolean bounded, boolean min) throws PrismException {
		if ((generateStrategy || exportAdv) && bounded) {
			throw new PrismException("Strategy synthesis for bounded properties is not supported yet.");
		}
		ModelCheckerResult res = new ModelCheckerResult();
		LpSolve lp;
		ArrayList<ArrayList<Double>> mgame = new ArrayList<ArrayList<Double>>();
		List<List<List<Map<BitSet, Double>>>> lstrat = null;
		List<Map<BitSet, Double>> kstrat = null;
		Map<Integer, BitSet> mmap = null;
		BitSet known = new BitSet();
		double[] nsol = new double[csg.getNumStates()];
		double[] ntmp = new double[csg.getNumStates()];
		long timer;
		int i, k, s;
		boolean done = false;
		if (generateStrategy || exportAdv) {
			mmap = new HashMap<Integer, BitSet>();
			kstrat = new ArrayList<Map<BitSet, Double>>();
			lstrat = new ArrayList<List<List<Map<BitSet, Double>>>>();
			lstrat.add(0, new ArrayList<List<Map<BitSet, Double>>>());
			lstrat.get(0).add(0, new ArrayList<Map<BitSet, Double>>());
			for (i = 0; i < csg.getNumStates(); i++) {	
				lstrat.get(0).get(0).add(i, null);
				kstrat.add(i, null);
			}
			if (bounded) {
				for (k = 1; k < limit; k++) {
					lstrat.get(0).add(k, new ArrayList<Map<BitSet, Double>>());
					for (i = 0; i < csg.getNumStates(); i++) {	
						lstrat.get(0).get(k).add(i, null);
					}
				}
			}
		}
		timer = System.currentTimeMillis();
		if (verbosity >= 1)
			mainLog.println("\nStarting value iteration...");
		try {
			lp = LpSolve.makeLp(maxCols + 1, maxRows + 1);
			lp.setVerbose(LpSolve.CRITICAL);	
		}
		catch(Exception e) {
			e.printStackTrace();
			throw new PrismException(e.toString());
		}
		known.or(no);
		known.or(yes);
		for (s = 0; s < csg.getNumStates(); s++) {
			nsol[s] = ntmp[s] = no.get(s)? 0.0 : yes.get(s)? 1.0 : 0.0;
		}
		k = 0;
		while(!done) {
			for(s = 0; s < csg.getNumStates(); s++) {
				if (!known.get(s)) {
					mgame = buildMatrixGame(csg, null, mmap, ntmp, s, min);
					try {
						if (min)
							lp.resizeLp(0, maxCols + 1);
						else
							lp.resizeLp(0, maxRows + 1);
					
					} 
					catch (LpSolveException e) {
						throw new PrismException("Exception raised by lpSolve when resizing linear program for state " +  s + " at iteration " + k);
					}
					nsol[s] = val(lp, mgame, kstrat, mmap, s, false, min);
					// player -> iteration -> state -> indexes -> value
					if (bounded && (generateStrategy || exportAdv)) {
						if (lstrat.get(0).get(k).get(s) == null || !lstrat.get(0).get(k - 1).get(s).equals(kstrat.get(s))) {
							lstrat.get(0).get(k).set(s, kstrat.get(s));
						}
						else {
							lstrat.get(0).get(k).set(s, lstrat.get(0).get(k - 1).get(s));
						}
					}
					else if (generateStrategy || exportAdv) {
						if (lstrat.get(0).get(0).get(s) == null) {
							lstrat.get(0).get(0).set(s, kstrat.get(s));
						}
						else if (!lstrat.get(0).get(0).get(s).equals(kstrat.get(s))) {
							lstrat.get(0).get(0).set(s, kstrat.get(s));
						}
					}
				}
				else if (generateStrategy || exportAdv) {
					lstrat.get(0).get(0).add(s, null);
				}
			}
			k++;
			done = PrismUtils.doublesAreClose(nsol, ntmp, termCritParam, termCrit == TermCrit.RELATIVE);
			if (!done && k == maxIters) {
				throw new PrismException("Could not converge after " + maxIters + " iterations");
			}
			else if(k == limit) {
				done = true;
			}
			else {
				ntmp = Arrays.copyOf(nsol, nsol.length);
			}
		}		
		mainLog.println("\nValue iteration converged after " + k + " iterations.");
		timer = System.currentTimeMillis() - timer;
		res.soln = nsol;
		res.numIters = k;
		if (generateStrategy || exportAdv)
			res.strat = new CSGStrategy(csg, lstrat, no, yes, new BitSet(), CSGStrategyType.ZERO_SUM);
		res.timeTaken = timer / 1000.0;	
		return res;
	}
	
	public ModelCheckerResult computeReachRewardsValIter(CSG csg, CSGRewards rewards, BitSet target, BitSet known, BitSet inf, double init[], int limit, boolean bounded, boolean min) throws PrismException {
		if ((generateStrategy || exportAdv) && bounded) {
			throw new PrismException("Strategy synthesis for bounded properties is not supported yet.");
		}
		ModelCheckerResult res = new ModelCheckerResult();
		LpSolve lp;
		ArrayList<ArrayList<Double>> mgame = new ArrayList<ArrayList<Double>>();
		List<List<List<Map<BitSet, Double>>>> lstrat = null;
		List<Map<BitSet, Double>> kstrat = null;
		Map<Integer, BitSet> mmap = null;
		BitSet unknown = new BitSet();
		double[] nsol = new double[csg.getNumStates()];
		double[] ntmp = new double[csg.getNumStates()];
		long timer;
		int i, k, s;
		boolean done = false;
		if (generateStrategy || exportAdv) {
			mmap = new HashMap<Integer, BitSet>();
			kstrat = new ArrayList<Map<BitSet, Double>>();
			lstrat = new ArrayList<List<List<Map<BitSet, Double>>>>(1);
			lstrat.add(0, new ArrayList<List<Map<BitSet, Double>>>());
			lstrat.get(0).add(0, new ArrayList<Map<BitSet, Double>>());
			for (i = 0; i < csg.getNumStates(); i++) {	
				lstrat.get(0).get(0).add(i, null);
				kstrat.add(i, null);
			}
			if (bounded) {
				for (k = 1; k < limit; k++) {
					lstrat.get(0).add(k, new ArrayList<Map<BitSet, Double>>());
					for (i = 0; i < csg.getNumStates(); i++) {	
						lstrat.get(0).get(k).add(i, null);
					}
				}
			}
		}
		timer = System.currentTimeMillis();
		if (verbosity >= 1)
			mainLog.println("\nStarting value iteration...");
		try {
			lp = LpSolve.makeLp(maxCols + 1, maxRows + 1);
			lp.setVerbose(LpSolve.CRITICAL);	
		}
		catch(Exception e) {
			e.printStackTrace();
			throw new PrismException(e.toString());
		}
		if (init != null) {
			if (known !=null) {
				for (i = 0; i < csg.getNumStates(); i++) 
					nsol[i] = ntmp[i] = known.get(i)? init[i] : target.get(i)? 0.0 : inf.get(i)? Double.POSITIVE_INFINITY : init[i];
			}
			else {
				for (i = 0; i < csg.getNumStates(); i++) 
					nsol[i] = ntmp[i] = target.get(i)? 0.0 : inf.get(i)? Double.POSITIVE_INFINITY : init[i];
			}
		}
		else {
			for (i = 0; i < csg.getNumStates(); i++)
				nsol[i] = ntmp[i] = target.get(i)? 0.0 : inf.get(i)? Double.POSITIVE_INFINITY : 0.0;
		}
		unknown.set(0, csg.getNumStates());
		unknown.andNot(target);
		unknown.andNot(inf);
		k = 0;
		while(!done) {
			for(s = 0; s < csg.getNumStates(); s++) {
				if (unknown.get(s)) {
					mgame = buildMatrixGame(csg, rewards, mmap, ntmp, s, min);
					try {
						if (min)
							lp.resizeLp(0, maxCols + 1);
						else
							lp.resizeLp(0, maxRows + 1);
					} 
					catch (LpSolveException e) {
						throw new PrismException("Exception raised by lpSolve when resizing linear program for state " +  s + " at iteration " + k);
					}
					nsol[s] = val(lp, mgame, kstrat, mmap, s, true, min);
					nsol[s] += rewards.getStateReward(s);
					if (bounded && (generateStrategy || exportAdv)) {
						// player -> iteration -> state -> indexes -> value
						if (lstrat.get(0).get(k).get(s) == null || !lstrat.get(0).get(k - 1).get(s).equals(kstrat.get(s))) {
							lstrat.get(0).get(k).set(s, kstrat.get(s));
						}
						else {
							lstrat.get(0).get(k).set(s, lstrat.get(0).get(k - 1).get(s));
						}
					}
					else if (generateStrategy || exportAdv) {
						if (lstrat.get(0).get(0).get(s) == null) {
							lstrat.get(0).get(0).set(s, kstrat.get(s));
						}
						else if (!lstrat.get(0).get(0).get(s).equals(kstrat.get(s))) {
							lstrat.get(0).get(0).set(s, kstrat.get(s));
						}
					}
				}
			}
			k++;
			done = PrismUtils.doublesAreClose(nsol, ntmp, termCritParam, termCrit == TermCrit.RELATIVE);
			if (!done && k == maxIters) {
				throw new PrismException("Could not converge after " + maxIters + " iterations");
			}
			else if(k == limit) {
				done = true;
			}
			else {
				ntmp = Arrays.copyOf(nsol, nsol.length);
			}
		}
		mainLog.println("\nValue iteration converged after " + k + " iterations.");
		timer = System.currentTimeMillis() - timer;
		res.soln = nsol;
		res.numIters = k;
		if (generateStrategy || exportAdv)
			res.strat = new CSGStrategy(csg, lstrat,  new BitSet(), target, inf, CSGStrategyType.ZERO_SUM);
		res.timeTaken = timer / 1000.0;	
		return res;
	}
	
	public ModelCheckerResult computeNextProbs(CSG csg, BitSet target, boolean min1, boolean min2, Coalition coalition) throws PrismException {		
		ModelCheckerResult res = null;		
		if (verbosity >= 1)
			mainLog.println("\nStarting next probabilistic reachability for CSGs...");
		buildCoalitions(csg, coalition, min1);
		switch (solnMethod) {
			case VALUE_ITERATION:
				res = computeReachProbsValIter(csg, new BitSet(), new BitSet(), 1, true, min1);
				break;
			default:
				throw new PrismException("Unknown CSG solution method " + solnMethod);
		}
		return res;
	}
	
	public ModelCheckerResult computeTotalRewards(CSG csg, CSGRewards csgRewards, boolean min1, boolean min2, Coalition coalition) throws PrismException {
		return computeTotalRewards(csg, csgRewards, coalition, min1, min2);
	}
	
	public ModelCheckerResult computeReachRewards(CSG csg, CSGRewards rewards, BitSet target, int unreachingSemantics, boolean min1, boolean min2, Coalition coalition) throws PrismException {		
		return computeReachRewards(csg, coalition, rewards, target, min1, min2, unreachingSemantics);
	}

	public ModelCheckerResult computeUntilProbs(CSG csg, BitSet remain, BitSet target, boolean min1, boolean min2, Coalition coalition) throws PrismException {
		return computeUntilProbs(csg, remain, target, maxIters, min1, min2, coalition);
	}
	
	public ModelCheckerResult computeBoundedUntilProbs(CSG csg, BitSet remain, BitSet target, int k, boolean min1, boolean min2, Coalition coalition) throws PrismException{
		return computeUntilProbs(csg, remain, target, k, min1, min2, coalition);
	}
	
	public ModelCheckerResult computeUntilProbs(CSG csg, BitSet remain, BitSet target, int bound, boolean min1, boolean min2, Coalition coalition) throws PrismException {
		ModelCheckerResult res = null;
		BitSet no, tmp, yes;
		int n, numYes, numNo;
		long timerProb0;
		if (verbosity >= 1)
			mainLog.println("\nStarting probabilistic reachability (until)...");
		n = csg.getNumStates();
		BitSet all = new BitSet();
		all.set(0, csg.getNumStates());
		no = new BitSet();
		tmp = new BitSet();
		yes = new BitSet();
		tmp = (BitSet) all.clone();
		tmp.andNot(target);
		tmp.andNot(remain);		
				
		// If <<C>>Pmax=?[F phi], sets coaltiion to compute <<N\C>>Pmax>=1 [G ¬phi] to compute the set "no", that is, the 
		// set of states from which N\C can ensure that C will not to reach phi
		// If <<C>>Pmin=?[F phi], sets coalition to compute <<C>>Pmax>=1 [G ¬phi] to compute the set "no", that is, the 
		// set of states from which C can ensure not to reach phi
		
		buildCoalitions(csg, coalition, !min1);
		timerProb0 = System.currentTimeMillis();
		if (precomp && prob0) { 
			no.or(target);
			no.flip(0, n);
			no = G(csg, no);
		}
		timerProb0 = System.currentTimeMillis() - timerProb0;
		
		no.or(tmp);
		tmp.clear();
		yes.or(target);				
		numYes = yes.cardinality();
		numNo = no.cardinality();
	
		buildCoalitions(csg, coalition, min1);
		
		if (verbosity >= 1)
			mainLog.println("target=" + target.cardinality() + ", yes=" + numYes + ", no=" + numNo + ", maybe=" + (n - (numYes + numNo)));
		switch (solnMethod) {
			case VALUE_ITERATION:
				res = computeReachProbsValIter(csg, no, yes, bound, false, min1);
				break;
			default:
				throw new PrismException("Unknown CSG solution method " + solnMethod);
		}
		
		res.timeProb0 = timerProb0 / 1000.0;
		if (verbosity >= 1)
			mainLog.println("Precomputation took " + timerProb0 / 1000.0 + " seconds.");
		return res;		
	}
	
	public ModelCheckerResult computeBoundedReachProbs(CSG csg, BitSet remain, BitSet target, int k, boolean min1, boolean min2, Coalition coalition, boolean genAdv) throws PrismException {
		ModelCheckerResult res = null;
		BitSet no = new BitSet();
		BitSet yes = new BitSet();
		if (verbosity >= 1)
			mainLog.println("\nStarting bounded probabilistic reachability...");
		buildCoalitions(csg, coalition, min1);
		no.set(0, csg.getNumStates());		
		no.andNot(target);
		no.andNot(remain);
		yes.or(target);
		switch (solnMethod) {
			case VALUE_ITERATION:
				res = computeReachProbsValIter(csg, no, yes, k, true, min1);
				break;
			default:
				throw new PrismException("Unknown CSG solution method " + solnMethod);
		}
		return res;
	}
	
	public ModelCheckerResult computeReachProbs(CSG csg, BitSet remain, BitSet target, boolean min1, boolean min2, int bound, Coalition coalition, boolean genAdv) throws PrismException {
		ModelCheckerResult res = null;
		BitSet no, yes;
		int n, numYes, numNo;
		long timerProb0, timerProb1;
		if (verbosity >= 1)
			mainLog.println("\nStarting probabilistic reachability...");
		n = csg.getNumStates();
		
		BitSet all = new BitSet();
		all.set(0, csg.getNumStates());
		no = new BitSet();
		yes = new BitSet();
		
		yes.or(target);
				
		// If <<C>>Pmax=?[F phi], sets coaltiion to compute <<N\C>>Pmax>=1 [G ¬phi] to compute the set "no", that is, the 
		// set of states from which N\C can ensure that C will not to reach phi
		// If <<C>>Pmin=?[F phi], sets coalition to compute <<C>>Pmax>=1 [G ¬phi] to compute the set "no", that is, the 
		// set of states from which C can ensure not to reach phi
		
		buildCoalitions(csg, coalition, !min1);
		timerProb0 = System.currentTimeMillis();
		if (precomp && prob0) { 
			no.or(target);
			no.flip(0, n);
			no = G(csg, no);
		}
		timerProb0 = System.currentTimeMillis() - timerProb0;
			
		// If <<C>>Pmax=?[F phi], sets coalition to compute <<C>>Pmax>=1 [F phi] to compute the set "yes"", that is, the
		// set of states from which C can reach phi with proability 1
		// If <<C>>Pmin=? [F phi], sets coalition to compute <<N\C>>Pmax>=1 [F phi] to compute the set "yes", that is, the
		// set of state from which N\C can force C to reach phi with probability 1
		
		buildCoalitions(csg, coalition, min1);
		timerProb1 = System.currentTimeMillis();
		if (precomp && prob1) {
			yes.or(AF(csg, target));
		}
		timerProb1 = System.currentTimeMillis() - timerProb1;
		
		numYes = yes.cardinality();
		numNo = no.cardinality();

		if (verbosity >= 1)
			mainLog.println("target=" + target.cardinality() + ", yes=" + numYes + ", no=" + numNo + ", maybe=" + (n - (numYes + numNo)));

		switch (solnMethod) {
			case VALUE_ITERATION:
				res = computeReachProbsValIter(csg, no, yes, bound, false, min1);
				break;
			default:
				throw new PrismException("Unknown CSG solution method " + solnMethod);
		}
		
		res.timeProb0 = timerProb0 / 1000.0;
		res.timePre = (timerProb0 + timerProb1) / 1000.0;
		if (verbosity >= 1)
			mainLog.println("Precomputation took " + res.timePre / 1000.0 + " seconds.");
		return res;
	}
	
	public ModelCheckerResult computeTotalRewards(CSG csg, CSGRewards rewards, Coalition coalition, boolean min1, boolean min2) throws PrismException {
		ModelCheckerResult res = new ModelCheckerResult();
		BitSet inf;
		double rew;
		long timer, timerPrecomp;
		int nonzero, infty, n, s, t;
		
		if (verbosity >= 1)
			mainLog.println("\nStarting total expected reward...");
		timer = System.currentTimeMillis();

		n = csg.getNumStates();
		
		BitSet zero = new BitSet();
		BitSet all = new BitSet();
		all.set(0, n);
		
		for(s = 0; s < n; s++) {
			if (rewards.getStateReward(s) != 0) {
				zero.set(s);
			}
			nonzero = 0;
			infty = 0;
			for (t = 0; t < csg.getNumChoices(s); t++) {
				rew = rewards.getTransitionReward(s, t);
				if (rew > 0.0 && rew != Double.POSITIVE_INFINITY && rew != Double.NEGATIVE_INFINITY) {
					nonzero++;
					zero.set(s);
				} 
				else if (rew == Double.POSITIVE_INFINITY || rew == Double.NEGATIVE_INFINITY)
					infty++;
			}
			if (nonzero != 0 && nonzero != csg.getNumChoices(s) - infty) {
				mainLog.println(csg.getStatesList().get(s));
				for (int c = 0; c < csg.getNumChoices(s); c++) {
					mainLog.println(csg.getAction(s, c) + ": " + rewards.getTransitionReward(s, c));
				}
				throw new PrismException("If a transition reward is nonzero, all reward transitions going from the state must be");
			}
		}
		
		zero.flip(0, n); // states with zero rewards
				
		// If <<C>>Rmax=?[C], builds coalition to compute the set of zero rewards states that N/C can force C to reach and stay (FG zero),
		// all other states are then states in which C can accumulate rewards indefinitely
		
		// If <<C>>Rmin=?[C], builds coalition to compute the set of zero rewards states that C can reach and stay without accumulating rewards,
		// all other states are then states in which N/C can force C to accumulate rewards indefinitely
	
		timerPrecomp = System.currentTimeMillis();
		buildCoalitions(csg, coalition, !min1);
		inf = AFG(csg, zero); 
		inf.flip(0, n);
		timerPrecomp = System.currentTimeMillis() - timerPrecomp;
		buildCoalitions(csg, coalition, min1);
		
		// Standard max reward calculation, but with empty target set
		switch (solnMethod) {
			case VALUE_ITERATION:
				res = computeReachRewardsValIter(csg, rewards, new BitSet(), null, inf, null, maxIters, false, min1);
				break;
			default:
				throw new PrismException("Unknown CSG solution method " + solnMethod);
		}

		timer = System.currentTimeMillis() - timer;
		res.timeTaken = timer;
		if (verbosity >= 1)
			mainLog.println("Precomputation took " + timerPrecomp / 1000.0 + " seconds.");
		mainLog.println("Expected total reward took " + timer / 1000.0 + " seconds.");
		return res;
	}
	
	public ModelCheckerResult computeInstantaneousRewards(CSG csg, CSGRewards csgRewards, Coalition coalition, int k, boolean min1, boolean min2) throws PrismException {
		LpSolve lp;
		ModelCheckerResult res = new ModelCheckerResult();
		ArrayList<ArrayList<Double>> mgame = new ArrayList<ArrayList<Double>>();
		List<Map<BitSet, Double>> kstrat = (generateStrategy || exportAdv)? new ArrayList<Map<BitSet, Double>>() : null;
		double nsol[], nsoln2[], ntmp[];
		long timer;
		int i, n;
		mainLog.println("\nStarting backwards instantaneous rewards computation...");
		n = csg.getNumStates();
		timer = System.currentTimeMillis();
		nsol = new double[n];
		nsoln2 = new double[n];
		for (i = 0; i < n; i++)
			nsol[i] = csgRewards.getStateReward(i);
		
		buildCoalitions(csg, coalition, min1);
		try {
			lp = LpSolve.makeLp(maxCols + 1, maxRows + 1);
			lp.setVerbose(LpSolve.IMPORTANT);
		}
		catch(Exception e) {
			e.printStackTrace();
			throw new PrismException(e.toString());
		}
		for (i = 0; i < k; i++) {
			for(int s = 0; s < csg.getNumStates(); s++) {				
				mgame.clear();
				mgame = buildMatrixGame(csg, null, null, nsol, s, min1);
				try {
					if (min1)
						lp.resizeLp(0, maxCols + 1);
					else
						lp.resizeLp(0, maxRows + 1);
					nsoln2[s] = val(lp, mgame, kstrat, null, s, true, min1);
				}
				catch(Exception e) {
					e.printStackTrace();
					throw new PrismException(e.toString());
				}
			}
			ntmp = nsol;
			nsol = nsoln2;
			nsoln2 = ntmp;
		}		
		timer = System.currentTimeMillis() - timer;
		mainLog.println("Backwards transient instantaneous rewards computation took " + i + " iters and " + timer / 1000.0 + " seconds.");
		res.soln = nsol;
		res.lastSoln = nsoln2;
		res.numIters = i;
		res.timeTaken = timer / 1000.0;
		res.timePre = 0.0;
		return res;
	}
	
	public ModelCheckerResult computeReachRewardsCumulative(CSG csg, Coalition coalition, CSGRewards rewards, BitSet target, boolean min1, boolean min2, boolean genAdv) throws PrismException {
		ModelCheckerResult res = new ModelCheckerResult();
		BitSet inf;
		double rew;
		long timer, timerPrecomp;
		int nonzero, infty, n, s, t, numTarget, numInf;;
		
		if (verbosity >= 1)
			mainLog.println("\nStarting expected reachability...");
		timer = System.currentTimeMillis();
		
		n = csg.getNumStates();
		
		BitSet zero = new BitSet();
		BitSet all = new BitSet();
		all.set(0, n);
		
		for(s = 0; s < n; s++) {
			if (target.get(s))
				continue;
			if (rewards.getStateReward(s) != 0) {
				zero.set(s);
			}
			nonzero = 0;
			infty = 0;
			for (t = 0; t < csg.getNumChoices(s); t++) {
				rew = rewards.getTransitionReward(s, t);
				if (rew > 0.0 && rew != Double.POSITIVE_INFINITY && rew != Double.NEGATIVE_INFINITY) {
					nonzero++;
					zero.set(s);
				} 
				else if (rew == Double.POSITIVE_INFINITY || rew == Double.NEGATIVE_INFINITY)
					infty++;
			}
			if (nonzero != 0 && nonzero != csg.getNumChoices(s) - infty) {
				mainLog.println(csg.getStatesList().get(s));
				for (int c = 0; c < csg.getNumChoices(s); c++) {
					mainLog.println(csg.getAction(s, c) + ": " + rewards.getTransitionReward(s, c));
				}
				throw new PrismException("If a transition reward is nonzero, all reward transitions going from the state must be");
			}
			
		}
		zero.flip(0, n); // States with zero state rewards
	
		// If <<C>>Rmax=?[C], builds coalition to compute the set of zero rewards states that N/C can force C to reach and stay (FG zero),
		// all other states are then states in which C can accumulate rewards indefinitely
		
		// If <<C>>Rmin=?[C], builds coalition to compute the set of zero rewards states that C can reach and stay without accumulating rewards,
		// all other states are then states in which N/C can force C to accumulate rewards indefinitely

		timerPrecomp = System.currentTimeMillis();
		buildCoalitions(csg, coalition, !min1);
		inf = AFG(csg, zero); 
		inf.flip(0, n);
		timerPrecomp = System.currentTimeMillis() - timerPrecomp;
		
		numTarget = target.cardinality();
		numInf = inf.cardinality();
		if (verbosity >= 1)
			mainLog.println("target=" + numTarget + ", inf=" + numInf + ", rest=" + (n - (numTarget + numInf)));
		
		buildCoalitions(csg, coalition, min1);
			
		switch (solnMethod) {
			case VALUE_ITERATION:
				res = computeReachRewardsValIter(csg, rewards, target, null, inf, null, maxIters, false, min1);
				break;
			default:
				throw new PrismException("Unknown CSG solution method " + solnMethod);
		}
		
		timer = System.currentTimeMillis() - timer;
		if (verbosity >= 1)
			mainLog.println("Precomputation took " + timerPrecomp / 1000.0 + " seconds.");
		mainLog.println("Expected total reward took " + timer / 1000.0 + " seconds.");
		res.timeTaken = timer;
		return res;
	}	

	public ModelCheckerResult computeReachRewardsInfinity(CSG csg, Coalition coalition, CSGRewards rewards, BitSet target, boolean min1, boolean min2)
															throws PrismException {
		ModelCheckerResult res = new ModelCheckerResult();
		BitSet inf;
		double[] init;
		long timer, timerProb1, timerApprox;
		int n, numTarget, numInf;

		if (verbosity >= 1)
			mainLog.println("\nStarting expected reachability...");

		timerApprox = 0;
		timer = System.currentTimeMillis();
		n = csg.getNumStates();

		BitSet all = new BitSet();
		all.set(0, n);
		
		init = new double[n];

		// If <<C>>Rmax=?[F t], builds coalition to compute <<N\C>>P=1[F t] that is, the set of states from which C can be forced to reach t with prob 1. 
		// The complement of this set is then <<C>>P<1[F t] as to maximize rewards C will try not reach the target and accumulate rewards indefinitely. 
		
		// If <<C>>Rmin=?[F t], builds coalition to compute <<C>>P=1[F t] that is, the set of states from which C can reach the target with prob 1 and thus
		// minimize rewards. The complement of this set is the set of states from which C can potentially be forced to accumulate rewards indefinitely. 
		
		timerProb1 = System.currentTimeMillis();
		buildCoalitions(csg, coalition, !min1);		
		inf = AF(csg, target);
		inf.flip(0, n);
		timerProb1 = System.currentTimeMillis() - timerProb1;
		buildCoalitions(csg, coalition, min1);

		numTarget = target.cardinality();
		numInf = inf.cardinality();
		if (verbosity >= 1)
			mainLog.println("target=" + numTarget + ", inf=" + numInf + ", rest=" + (n - (numTarget + numInf)));
		
		// Computes rewards with epsilon instead of zero. This is used to get the
		// over-approximation of the real result, which deals with the problem of staying 
		// in zero components for free when infinity should be gained.

		// First, get the minimum nonzero reward and maximal reward, will be
		// used as a basis for epsilon.
		// Also, check if by any chance all rewards are nonzero, then we don't
		// need to precompute.
		double minimumReward = Double.POSITIVE_INFINITY;
		double maximumReward = 0.0;
		boolean allNonzero = true;
		double r;
		for (int i = 0; i < n; i++) {
			r = rewards.getStateReward(i);
			if (r > 0.0 && r < minimumReward)
				minimumReward = r;
			if (r > maximumReward)
				maximumReward = r;
			allNonzero = allNonzero && r > 0;

			for (int j = 0; j < csg.getNumChoices(i); j++) {
				r = rewards.getTransitionReward(i, j);
				if (r > 0.0 && r < minimumReward)
					minimumReward = r;
				if (r > maximumReward)
					maximumReward = r;
				allNonzero = allNonzero && rewards.getTransitionReward(i, j) > 0;

				for (int k = 0; k < csg.getNumNestedChoices(i, j); k++) {
					r = rewards.getNestedTransitionReward(i, j, k);
					if (r > 0.0 && r < minimumReward)
						minimumReward = r;
					if (r > maximumReward)
						maximumReward = r;
					allNonzero = allNonzero && r > 0;
				}
			}
		}
		
		if (!allNonzero && !(rewards instanceof StateRewardsConstant)) {
			timerApprox = System.currentTimeMillis();
			// A simple heuristic that gives small epsilon, but still is
			// hopefully safe floating-point-wise
			double epsilon = Math.min(minimumReward, maximumReward * 0.01);

			if (verbosity >= 1) {
				mainLog.println("Computing the upper bound where " + epsilon + " is used instead of 0.0");
			}

			// Modifies the rewards
			double origZeroReplacement;
			if (rewards instanceof CSGRewardsSimple) {
				origZeroReplacement = ((CSGRewardsSimple) rewards).getZeroReplacement();
				((CSGRewardsSimple) rewards).setZeroReplacement(epsilon);
			} else {
				throw new PrismException("To compute expected reward I need to modify the reward structure. But I don't know how to modify" + rewards.getClass().getName());
			}
			// Computes the value when rewards are nonzero
			switch (solnMethod) {
				case VALUE_ITERATION:
					init = computeReachRewardsValIter(csg, rewards, target, null, inf, null, maxIters, false, min1).soln;
					break;
				default:
					throw new PrismException("Unknown CSG solution method " + solnMethod);
			}
			
			// Set the value iteration result to be the initial solution for the
			// next part in which "proper" zero rewards are used

			// Returns the rewards to the original state
			if (rewards instanceof CSGRewardsSimple) {
				((CSGRewardsSimple) rewards).setZeroReplacement(origZeroReplacement);
			}
			timerApprox = System.currentTimeMillis() - timerApprox;

			if (verbosity >= 1) {
				mainLog.println("Computed an over-approximation of the solution (in " + timerApprox / 1000.0 + " seconds), this will now be used to get the solution");
			}
		}

		// Compute real rewards
		switch (solnMethod) {
			case VALUE_ITERATION:
				res = computeReachRewardsValIter(csg, rewards, target, null, inf, init, maxIters, false, min1);
				break;
			default:
				throw new PrismException("Unknown CSG solution method " + solnMethod);
		}

		timer = System.currentTimeMillis() - timer;
		if (verbosity >= 1)
			mainLog.println("Expected reachability took " + timer / 1000.0 + " seconds.");

		res.timeTaken = timer / 1000.0;
		res.timePre = (timerProb1 + timerApprox) / 1000.0;
		if (verbosity >= 1)
			mainLog.println("Precomputation took " + timerProb1 / 1000.0 + " seconds.");
		return res;
	}
	
	public ModelCheckerResult computeCumulativeRewards(CSG csg, CSGRewards csgRewards, Coalition coalition, int k, boolean min1, boolean min2, boolean genAdv) throws PrismException {
		ModelCheckerResult res = new ModelCheckerResult();
		buildCoalitions(csg, coalition, min1);
		switch (solnMethod) {
			case VALUE_ITERATION:
				res = computeReachRewardsValIter(csg, csgRewards, new BitSet(), new BitSet(), new BitSet(), null, k, true, min1);
				break;
			default:
				throw new PrismException("Unknown CSG solution method " + solnMethod);
		}
		return res;
	}
	
	public ModelCheckerResult computeProbBoundedEquilibria(CSG csg, List<Coalition> coalitions, List<ExpressionTemporal> exprs, BitSet[] targets, BitSet[] remain, int[] bounds, boolean min) throws PrismException {
		ModelCheckerResult res = new ModelCheckerResult();
		CSGModelCheckerEquilibria csgeq = new CSGModelCheckerEquilibria(this);
		res = csgeq.computeBoundedEquilibria(csg, coalitions, null, exprs, targets, remain, bounds, min);
		return res;		
	}
	
	public ModelCheckerResult computeProbReachEquilibria(CSG csg, List<Coalition> coalitions, BitSet[] targets, BitSet[] remain, boolean min) throws PrismException {
		ModelCheckerResult res = new ModelCheckerResult();
		CSGModelCheckerEquilibria csgeq = new CSGModelCheckerEquilibria(this);
		res = csgeq.computeReachEquilibria(csg, coalitions, null, targets, remain, min);
		return res;
	}
	
	public ModelCheckerResult computeRewBoundedEquilibria(CSG csg, List<Coalition> coalitions, List<CSGRewards> rewards, List<ExpressionTemporal> exprs, int[] bounds, boolean min) throws PrismException {
		ModelCheckerResult res = new ModelCheckerResult();
		CSGModelCheckerEquilibria csgeq = new CSGModelCheckerEquilibria(this);
		res = csgeq.computeBoundedEquilibria(csg, coalitions, rewards, exprs, null, null, bounds, min);
		return res;
	}
	
	public ModelCheckerResult computeRewReachEquilibria(CSG csg, List<Coalition> coalitions, List<CSGRewards> rewards, BitSet[] targets, boolean min) throws PrismException {
		ModelCheckerResult res = new ModelCheckerResult();
		CSGModelCheckerEquilibria csgeq = new CSGModelCheckerEquilibria(this);
		res = csgeq.computeReachEquilibria(csg, coalitions, rewards, targets, null, min);
		return res;
	}
	
	public ModelCheckerResult computeMultiProbReachEquilibria(CSG csg, List<Coalition> coalitions, BitSet[] targets, BitSet[] remain, boolean min) throws PrismException {
		ModelCheckerResult res = new ModelCheckerResult();
		CSGModelCheckerEquilibria csgeq = new CSGModelCheckerEquilibria(this);
		res = csgeq.computeMultiReachEquilibria(csg, coalitions, null, targets, remain, min);
		return res;
	}
		
	public ModelCheckerResult computeMultiRewReachEquilibria(CSG csg, List<Coalition> coalitions, List<CSGRewards> rewards, BitSet[] targets, boolean min) throws PrismException {
		ModelCheckerResult res = new ModelCheckerResult();
		CSGModelCheckerEquilibria csgeq = new CSGModelCheckerEquilibria(this);
		res = csgeq.computeMultiReachEquilibria(csg, coalitions, rewards, targets, null, min);
		return res;
	}

	public ModelCheckerResult computeMixedEquilibria(CSG csg, List<Coalition> coalitions, List<CSGRewards> rewards, List<ExpressionTemporal> exprs, BitSet bounded, BitSet[] targets, BitSet[] remain, int[] bounds, boolean min) throws PrismException {
		ModelCheckerResult res = new ModelCheckerResult();
		CSGModelCheckerEquilibria csgeq = new CSGModelCheckerEquilibria(this);
		LTLModelChecker ltlmc = new LTLModelChecker(this);
		LTLProduct<CSG> product;
		List<CSGRewards> newrewards = new ArrayList<CSGRewards>();
		BitSet newremain[];
		BitSet newtargets[];
		int index, i, s, t;
		boolean rew;
		
		/*
		Path currentRelativePath = Paths.get("");
		String path = currentRelativePath.toAbsolutePath().toString();
		*/
		
		rew = rewards != null;
		index = bounded.nextSetBit(0);
		
		newremain = new BitSet[remain.length];
		newtargets = new BitSet[targets.length];
		
		Arrays.fill(newremain, null);
		
		/*		
		LTL2DA ltl2da = new LTL2DA(this);
		DA<BitSet,? extends AcceptanceOmega> daex = ltl2da.convertLTLFormulaToDA(exprs.get(index), csg.getConstantValues(), AcceptanceType.RABIN);
		try(OutputStream out1 = 
				new FileOutputStream(path + "/dra.dot")) {
					try (PrintStream printStream = 
							new PrintStream(out1)) {
								da.printDot(printStream);
								printStream.close();
					}
		}
		catch(Exception e) {
			e.printStackTrace();
		}
		*/
		
		AcceptanceType[] allowedAcceptance = {
				AcceptanceType.RABIN,
				AcceptanceType.REACH,
				AcceptanceType.BUCHI,
				AcceptanceType.STREETT,
				AcceptanceType.GENERIC
		};

		//csg.exportToDotFile(path + "/model.dot");
		
		if (rew) {			
			BitSet all = new BitSet();
			all.set(0, csg.getNumStates());
			DA<BitSet, AcceptanceRabin> da = null;
			Vector<BitSet> labelBS = new Vector<BitSet>();
			labelBS.add(0, all);
			targets[index] = all;
			switch (exprs.get(index).getOperator()) {
				case ExpressionTemporal.R_I:
					da = constructDRAForInstant("L0", false, new IntegerBound(null, false, bounds[index] + 1, false));
					break;
				case ExpressionTemporal.R_C:
					da = constructDRAForInstant("L0", false, new IntegerBound(null, false, bounds[index], false));
					break;	
			}			
						
			DASimplifyAcceptance.simplifyAcceptance(this, da, AcceptanceType.REACH);
			
			/*
			try(OutputStream out1 = 
					new FileOutputStream(path + "/dra.dot")) {
						try (PrintStream printStream = 
								new PrintStream(out1)) {
									da.printDot(printStream);
									printStream.close();
						}
			}
			catch(Exception e) {
				e.printStackTrace();
			}
			*/
			
			mainLog.println("\nConstructing CSG-" + da.getAutomataType() + " product...");
			product = ltlmc.constructProductModel(da, csg, labelBS, null);
			mainLog.print("\n" + product.getProductModel().infoStringTable());
		}
		else {
			product = ltlmc.constructProductCSG(this, csg, exprs.get(index), null, allowedAcceptance);
		}
		
		product.productModel.clearInitialStates();
		product.productModel.addInitialState(product.getModelState(csg.getFirstInitialState()));

		//product.productModel.exportToDotFile(path + "/product.dot");
	
		/*
		try {
			PrismFileLog pflog = new PrismFileLog(path + "/product.dot");
			System.out.println(path + "/product.dot");
			product.productModel.exportToDotFile(pflog, null, true);
		}
		catch(Exception e) {
			e.printStackTrace();
		}
		*/
								
		if (product.getAcceptance() instanceof AcceptanceReach) {
			mainLog.println("\nSkipping BSCC computation since acceptance is defined via goal states...");
			newtargets[index] = ((AcceptanceReach) product.getAcceptance()).getGoalStates();
		} else {
			mainLog.println("\nFinding accepting BSCCs...");
			newtargets[index] = ltlmc.findAcceptingBSCCs(product.getProductModel(), product.getAcceptance());
		}
				
		for (i = 0; i < targets.length; i++) {
			if (i != index) 
				newtargets[i] = product.liftFromModel(targets[i]);
		}
		
		for (i = 0; i < remain.length; i++) {
			if (remain[i] != null) {
				newremain[i] = product.liftFromModel(remain[i]);
			}
		}
						
		if (rew) {			
			for (i = 0; i < coalitions.size(); i++) {
				CSGRewards reward = new CSGRewardsSimple(product.productModel.getNumStates());
				if (i != index) {
					for (s = 0; s < product.productModel.getNumStates(); s++){
						((CSGRewardsSimple) reward).addToStateReward(s, rewards.get(i).getStateReward(product.getModelState(s)));
						for (t = 0; t < product.productModel.getNumChoices(s); t++) {
							((CSGRewardsSimple) reward).addToTransitionReward(s, t, rewards.get(i).getTransitionReward(product.getModelState(s), t));
						}
					}
				}
				else {
					switch (exprs.get(index).getOperator()) {
						case ExpressionTemporal.R_I:
							for (s = 0; s < product.productModel.getNumStates(); s++){
								for (t = 0; t < product.productModel.getNumChoices(s); t++) {
									for (int u : product.productModel.getChoice(s, t).getSupport()) {
										if (newtargets[index].get(u)) {
											((CSGRewardsSimple) reward).addToStateReward(s, rewards.get(i).getStateReward(product.getModelState(s)));
										}
									}
									((CSGRewardsSimple) reward).addToTransitionReward(s, t, 0.0);
								}
							}
							break;
						case ExpressionTemporal.R_C:
							for (s = 0; s < product.productModel.getNumStates(); s++){
								if (!newtargets[index].get(s)) {
									((CSGRewardsSimple) reward).addToStateReward(s, rewards.get(i).getStateReward(product.getModelState(s)));
									for (t = 0; t < product.productModel.getNumChoices(s); t++) {
										((CSGRewardsSimple) reward).addToTransitionReward(s, t, rewards.get(i).getTransitionReward(product.getModelState(s), t));
									}
								}
								else {
									((CSGRewardsSimple) reward).addToStateReward(s, 0.0);
									for (t = 0; t < product.productModel.getNumChoices(s); t++) {
										((CSGRewardsSimple) reward).addToTransitionReward(s, t, 0.0);
									}
								}
							}
							break;
					}	
				}
				newrewards.add(i, reward);
			}
			/*** Optional filtering ***/
			/*
			CSG csg_rm = new CSG(csg.getPlayers());
			List<CSGRewards> csg_rew_rm = new ArrayList<CSGRewards>();
			map_state = new HashMap<Integer, Integer>();
			list_state = new ArrayList<State>();
			map_state.put(product.productModel.getFirstInitialState(), csg_rm.addState());
			csg_rm.addInitialState(map_state.get(product.productModel.getFirstInitialState()));
			for (i = 0; i < rewards.size(); i++) {
				csg_rew_rm.add(i, new CSGRewardsSimple(product.productModel.getNumStates()));
			}
			filterStates(product.productModel, csg_rm, newrewards, csg_rew_rm, product.productModel.getFirstInitialState());
			csg_rm.setVarList(csg.getVarList());
			csg_rm.setStatesList(list_state);
			csg_rm.setActions(csg.getActions());
			csg_rm.setPlayers(csg.getPlayers());
			csg_rm.setIndexes(csg.getIndexes());
			csg_rm.setIdles(csg.getIdles());
			csg_rm.exportToDotFile(path + "/filtered.dot");
			BitSet[] filtered_targets = new BitSet[targets.length];
			for (i = 0; i < targets.length; i++) {
				filtered_targets[i] = new BitSet();
				for (s = newtargets[i].nextSetBit(0); s >= 0; s = newtargets[i].nextSetBit(s + 1)) {
					if (map_state.get(s) != null)
						filtered_targets[i].set(map_state.get(s));
					System.out.println(map_state.get(s));
				}
			}
			res = csgeq.computeReachEquilibria(csg_rm, coalitions, csg_rew_rm, filtered_targets, null);			
			*/
			res = csgeq.computeReachEquilibria(product.productModel, coalitions, newrewards, newtargets, null, min);		
		}
		else {
			res = csgeq.computeReachEquilibria(product.productModel, coalitions, null, newtargets, newremain, min);		
		}
		return res;
	}

	public static DA<BitSet,AcceptanceRabin> constructDRAForInstant(String labelA, boolean negateA, IntegerBound bounds) {
		DA<BitSet,AcceptanceRabin> dra;
		List<String> apList = new ArrayList<String>();
		BitSet edge_no, edge_yes;
		BitSet accL, accK;

		int saturation = bounds.getMaximalInterestingValue();
		
		// [0,saturation] + yes
		int states = saturation + 2;
		
		apList.add(labelA);
		
		dra = new DA<BitSet,AcceptanceRabin>(states);
		dra.setAcceptance(new AcceptanceRabin());
		dra.setAPList(apList);
		dra.setStartState(0);
		
		// edge labeled with the target label
		edge_yes  = new BitSet();
		// edge not labeled with the target label
		edge_no = new BitSet();
		
		if (negateA) {
			edge_no.set(0);  // no = a, yes = !a
		} 
		else {
			edge_yes.set(0); // yes = a, no = !a
		}

		int yes_state = states - 1;
		int next_counter;
		
		for (int counter = 0; counter <= saturation; counter++) {
			next_counter = counter + 1;
			if (next_counter > saturation) next_counter = saturation;
			
			if (bounds.isInBounds(counter)) {
				dra.addEdge(counter, edge_no, next_counter);
				if (counter <= saturation - 2) 
					dra.addEdge(counter, edge_yes, next_counter);
				else {
					dra.addEdge(counter, edge_yes, yes_state);
				}
			} 
			else {
				dra.addEdge(counter, edge_no, next_counter);
				dra.addEdge(counter, edge_yes, next_counter);
			}
		}

		// yes state = true loop
		dra.addEdge(yes_state, edge_no, yes_state);
		dra.addEdge(yes_state, edge_yes, yes_state);

		// acceptance =
		// infinitely often yes_state,
		// not infinitely often saturation state
		// this allows complementing by switching L and K.
		accL = new BitSet();
		accL.set(saturation);
		accK = new BitSet();
		accK.set(yes_state);

		dra.getAcceptance().add(new RabinPair(accL, accK));

		return dra;
	}
	
	protected HashMap<Integer, Integer> map_state;
	protected List<State> list_state;
	 
	public void filterStates(CSG csg, CSG csg_rm, List<CSGRewards> rewards, List<CSGRewards> rew_rm, int s) {
		Distribution d;
		int i;
		list_state.add(csg.getStatesList().get(s));
   		for (int c = 0; c < csg.getNumChoices(s); c++) { // gets all choices
   			d = new Distribution();
			for (int t : csg.getChoice(s, c).getSupport()) { // gets all targets
				if(!map_state.keySet().contains(t)) { //if not yet explored
					map_state.put(t, csg_rm.addState());
					filterStates(csg, csg_rm, rewards, rew_rm, t);
				}
				d.add(map_state.get(t), csg.getChoice(s, c).get(t)); // adds target to distribution
			}
			i = csg_rm.addActionLabelledChoice(map_state.get(s), d, csg.getAction(s, c));
			csg_rm.setIndexes(map_state.get(s), i, csg.getIndexes(s, c));
			if (rewards != null && rew_rm != null) {	
				for (int r = 0; r < rewards.size(); r++) {
					if (rewards.get(r) != null && rew_rm.get(r) != null)
						((CSGRewardsSimple) rew_rm.get(r)).addToTransitionReward(map_state.get(s), i, rewards.get(r).getTransitionReward(s, c));
				}
			}
		}
		if (rewards != null && rew_rm != null) {	
			for (int r = 0; r < rewards.size(); r++) {
				if (rewards.get(r) != null && rew_rm.get(r) != null)
					((CSGRewardsSimple) rew_rm.get(r)).addToStateReward(map_state.get(s), rewards.get(r).getStateReward(s));
			}
		}
	}
	
}
