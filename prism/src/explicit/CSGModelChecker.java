//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <david.parker@cs.ox.ac.uk> (University of Oxford)
// 	* Gabriel Santos <gabriel.santos@cs.ox.ac.uk> (University of Oxford)
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
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Vector;
import java.util.stream.Collectors;

import explicit.rewards.STPGRewards;
import org.apache.commons.math3.util.Precision;

import acceptance.AcceptanceRabin;
import acceptance.AcceptanceRabin.RabinPair;
import acceptance.AcceptanceReach;
import acceptance.AcceptanceType;
import automata.DA;
import automata.DASimplifyAcceptance;
import explicit.LTLModelChecker.LTLProduct;
import explicit.rewards.CSGRewards;
import explicit.rewards.CSGRewardsSimple;
import explicit.rewards.StateRewardsConstant;
import lpsolve.LpSolve;
import lpsolve.LpSolveException;
import soplex.SoPlex;
import parser.State;
import parser.ast.Coalition;
import parser.ast.Expression;
import parser.ast.ExpressionBinaryOp;
import parser.ast.ExpressionStrategyQual;
import parser.ast.ExpressionUnaryOp;
import parser.type.TypeBool;
import prism.PrismNotSupportedException;
import parser.ast.ExpressionTemporal;
import prism.IntegerBound;
import prism.PrismComponent;
import prism.PrismDevNullLog;
import prism.PrismException;
import prism.PrismSettings;
import prism.PrismUtils;
import strat.CSGStrategy;

/**
 * Explicit-state model checker for concurrent stochastic games (CSGs).
 */
public class CSGModelChecker extends ProbModelChecker
{
	public static final int R_INFINITY = 0;
	public static final int R_CUMULATIVE = 1;
	public static final int R_ZERO = 2;

	protected double scaleFactor = getSettings().getDouble(PrismSettings.PRISM_ZS_LP_SCALE_FACTOR);

	// Info about the current coalitions for model checking
	// (here, there are two, and the first always maximises)
	
	/** Number of players */
	protected int numPlayers;
	/** Number of coalitions */
	protected int numCoalitions;
	/** For each coalition, which players are in it (BitSet over 0-indexed player indices) */
	protected BitSet[] coalitionIndexes;
	/** For each coalition, which actions are used by players in it
	 * (BitSet over 1-indexed action indices, and including "idle" actions) */
	protected BitSet[] actionIndexes;
	
	/** For the current coalition, the max number of rows in the matrix game across all CSG states */
	protected int maxRows;
	/** For the current coalition, the max number of columns in the matrix game across all CSG states */
	protected int maxCols;
	/** For each coalition, the max number of actions across all CSG states */
	protected int[] maxNumActions;
	/** For each coalition, the average number of actions across all CSG states */
	protected double[] avgNumActions;

	// Info about the current matrix game being built/solved for state s
	// (as above, here, we assume that the first/second coalition maximise/mimimise
	
	/** For each coalition, a description of each coalition action in s */
	protected ArrayList<ArrayList<String>> actions;
	/** For each coalition, the indices of all coalition actions in s */ 
	protected ArrayList<ArrayList<Integer>> strategies;
	/** Matrix game values, stored as mapping from pairs of coalition action indices,
	 *  as a BitSet, to the value, stored in a 1-element ArrayList */
	protected HashMap<BitSet, ArrayList<Double>> utilities;
	/** CSG distribution used to compute value, again as mapping from BitSet */
	protected HashMap<BitSet, ArrayList<Distribution<Double>>> probabilities;
	
	/** Total number of coalition actions in s */
	protected int varIndex;
	/** Minimum matrix value */ 
	protected double minEntry;
	/** True if all matrix values are equal */
	protected boolean allEqual;

	/** Which solver to use for linear programming */
	protected String lpSolver;
	/** SoPlex instance reused for all matrix games (created on first use when -lpsolver soplex) */
	protected SoPlex soplex = null;
	/** SoPlex scaling method (see PrismSettings.PRISM_CSG_SOPLEX_SCALING) */
	protected String soplexScaling;

	protected long timerVal;


	/**
	 * Create a new CSGModelChecker, inherit basic state from parent (unless null).
	 */
	public CSGModelChecker(PrismComponent parent) throws PrismException
	{
		super(parent);
		utilities = new HashMap<BitSet, ArrayList<Double>>();
		probabilities = new HashMap<BitSet, ArrayList<Distribution<Double>>>();
		actions = new ArrayList<ArrayList<String>>();
		strategies = new ArrayList<ArrayList<Integer>>();
		lpSolver = getSettings().getString(PrismSettings.PRISM_CSG_LPSOLVER);
		soplexScaling = getSettings().getString(PrismSettings.PRISM_CSG_SOPLEX_SCALING);
	}

	// Numerical computation functions

	/**
	 * Compute next-state probabilities.
	 * i.e. compute the probability of being in a state in {@code target} in the next step.
	 * 
	 * @param csg The CSG
	 * @param target Target states
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 * @param coalition The coalition of players which define player 1
	 */
	public ModelCheckerResult computeNextProbs(CSG<Double> csg, BitSet target, boolean min1, boolean min2, Coalition coalition) throws PrismException
	{
		ModelCheckerResult res = new ModelCheckerResult();
		LpSolve lp;
		ArrayList<ArrayList<Double>> mgame = new ArrayList<ArrayList<Double>>();
		Map<Integer, BitSet> mmap = null;
		List<List<List<Map<BitSet, Double>>>> lstrat = null;
		List<Map<BitSet, Double>> kstrat = null;
		double[] nsol = new double[csg.getNumStates()];
		long timer;
		int s, i;

		if (genStrat) {
			mmap = new HashMap<Integer, BitSet>();
			kstrat = new ArrayList<Map<BitSet, Double>>();
			lstrat = new ArrayList<List<List<Map<BitSet, Double>>>>(1);
			lstrat.add(0, new ArrayList<List<Map<BitSet, Double>>>());
			lstrat.get(0).add(0, new ArrayList<Map<BitSet, Double>>());
			for (i = 0; i < csg.getNumStates(); i++) {
				lstrat.get(0).get(0).add(i, null);
				kstrat.add(i, null);
			}
		}
		timer = System.currentTimeMillis();
		buildCoalitions(csg, coalition, min1);
		for (s = 0; s < csg.getNumStates(); s++) {
			nsol[s] = target.get(s) ? 1.0 : 0.0;
		}
		try {
			lp = LpSolve.makeLp(maxCols + 1, maxRows + 1);
			lp.setVerbose(LpSolve.CRITICAL);
		} catch (Exception e) {
			e.printStackTrace();
			throw new PrismException(e.toString());
		}
		for (s = 0; s < csg.getNumStates(); s++) {
			mgame = buildMatrixGame(csg, null, mmap, nsol, s, min1);
			nsol[s] = val(lp, mgame, kstrat, mmap, s, true, min1);
			if (genStrat) {
				updateStrategy(kstrat, lstrat, 0, s, false);
			}
		}
		timer = System.currentTimeMillis() - timer;
		res.soln = nsol;
		res.lastSoln = null;
		res.numIters = 1;
		res.timeTaken = timer / 1000.0;
		res.timePre = 0.0;

		if (genStrat) {
			CSGStrategy strat = new CSGStrategy(csg, lstrat, new BitSet(), target, new BitSet());
			strat.setCoalitionActions(actionIndexes[min1 ? 1 : 0]);
			strat.setComplemented(negatedPath);
			res.strat = strat;
		}
		return res;
	}

	/**
	 * Compute bounded reachability/until probabilities.
	 * i.e. compute the min/max probability of reaching a state in {@code target},
	 * within k steps, and while remaining in states in {@code remain}.
	 * 
	 * @param csg The CSG
	 * @param remain Remain in these states (optional: null means "all")
	 * @param target Target states
	 * @param k Bound
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 * @param coalition The coalition of players which define player 1
	 */
	public ModelCheckerResult computeBoundedReachProbs(CSG<Double> csg, BitSet remain, BitSet target, int k, boolean min1, boolean min2, Coalition coalition,
			boolean genAdv) throws PrismException
	{
		// TODO: confirm that the case min1==min2 is not handled  
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

	/**
	 * Compute until probabilities.
	 * i.e. compute the min/max probability of reaching a state in {@code target},
	 * and while remaining in states in {@code remain}.
	 * 
	 * @param csg The CSG
	 * @param remain Remain in these states (optional: null means "all")
	 * @param target Target states
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 * @param coalition The coalition of players which define player 1
	 */
	public ModelCheckerResult computeUntilProbs(CSG<Double> csg, BitSet remain, BitSet target, boolean min1, boolean min2, Coalition coalition) throws PrismException
	{
		return computeUntilProbs(csg, remain, target, maxIters, min1, min2, coalition);
	}

	/**
	 * Compute bounded until probabilities.
	 * i.e. compute the min/max probability of reaching a state in {@code target},
	 * within k steps, and while remaining in states in {@code remain}.
	 * 
	 * @param csg The CSG
	 * @param remain Remain in these states (optional: null means "all")
	 * @param target Target states
	 * @param k Bound
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 * @param coalition The coalition of players which define player 1
	 */
	public ModelCheckerResult computeBoundedUntilProbs(CSG<Double> csg, BitSet remain, BitSet target, int k, boolean min1, boolean min2, Coalition coalition)
			throws PrismException
	{
		ModelCheckerResult res;
		Boolean tmp;
		tmp = precomp;
		precomp = false;
		res = computeUntilProbs(csg, remain, target, k, min1, min2, coalition);
		precomp = tmp;
		return res;
	}

	/**
	 * Compute (possibly bounded) reachability probabilities.
	 * i.e. compute the min/max probability of reaching a state in {@code target}.
	 * 
	 * @param csg The CSG
	 * @param target Target states
	 * @param bound Bound
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 * @param coalition The coalition of players which define player 1
	 */
	public ModelCheckerResult computeReachProbs(CSG<Double> csg, BitSet target, boolean min1, boolean min2, int bound, Coalition coalition,
			boolean genAdv) throws PrismException
	{
		// TODO: confirm that the case min1==min2 is not handled  
		ModelCheckerResult res = null;
		BitSet no, yes;
		int n, numYes, numNo;
		long timerProb0, timerProb1;
		boolean bounded = (bound != maxIters);
		if (bounded && genStrat)
			throw new PrismNotSupportedException("Strategy synthesis for bounded properties is not supported yet");
		
		if (verbosity >= 1)
			mainLog.println("\nStarting probabilistic reachability...");
		n = csg.getNumStates();

		BitSet all = new BitSet();
		all.set(0, csg.getNumStates());
		no = new BitSet();
		yes = new BitSet();

		yes.or(target);

		// If <<C>>Pmax=?[F phi], sets coalition to compute <<N\C>>Pmax>=1 [G ¬phi] to compute the set "no", that is, the 
		// set of states from which N\C can ensure that C will not to reach phi
		// If <<C>>Pmin=?[F phi], sets coalition to compute <<C>>Pmax>=1 [G ¬phi] to compute the set "no", that is, the 
		// set of states from which C can ensure not to reach phi

		buildCoalitions(csg, coalition, !min1);
		timerProb0 = System.currentTimeMillis();
		if (!bounded && precomp && prob0) {
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
		
		if (!bounded && precomp && prob1) {
			if (!genStrat) {
				yes.or(AF(csg, target));
			} else {
				mainLog.println("Disabling Prob1 precomputation to allow strategy generation");
			}
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

	/**
	 * Compute (possibly bounded) reachability/until probabilities.
	 * i.e. compute the min/max probability of reaching a state in {@code target},
	 * and while remaining in states in {@code remain}.
	 * 
	 * @param csg The CSG
	 * @param remain Remain in these states (optional: null means "all")
	 * @param target Target states
	 * @param bound Bound
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 * @param coalition The coalition of players which define player 1
	 */
	public ModelCheckerResult computeUntilProbs(CSG<Double> csg, BitSet remain, BitSet target, int bound, boolean min1, boolean min2, Coalition coalition)
			throws PrismException
	{
		// TODO: confirm that the case min1==min2 is not handled  
		ModelCheckerResult res = null;
		BitSet no, tmp, yes;
		int n, numYes, numNo;
		long timerProb0, timerProb1;

		boolean bounded = (bound != maxIters);
		if (bounded && genStrat)
			throw new PrismNotSupportedException("Strategy synthesis for bounded properties is not supported yet");
		
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

		// If <<C>>Pmax=?[phi1 U phi2], sets coalition to compute <<N\C>>Pmax>=1 [G ¬phi2] to compute the set "no", that is, the 
		// set of states from which N\C can ensure that C will not to reach phi2
		// If <<C>>Pmin=?[phi1 U phi2], sets coalition to compute <<C>>Pmax>=1 [G ¬phi2] to compute the set "no", that is, the 
		// set of states from which C can ensure not to reach phi2

		buildCoalitions(csg, coalition, !min1);
		timerProb0 = System.currentTimeMillis();
		if (!bounded && precomp && prob0) {
			no.or(target);
			no.flip(0, n);
			no = G(csg, no);
		}
		timerProb0 = System.currentTimeMillis() - timerProb0;
		no.or(tmp);
		tmp.clear();
		tmp.set(0, n);
		tmp.andNot(remain);
		tmp.andNot(target);
		no.or(tmp);

		// If <<C>>Pmax=?[phi1 U ph2], sets coalition to compute <<C>>Pmax>=1 [F phi2] to compute the set "yes"", that is, the
		// set of states from which C can reach phi with proability 1 while remaining in phi1
		// If <<C>>Pmin=? [phi1 U phi2], sets coalition to compute <<N\C>>Pmax>=1 [F phi2] to compute the set "yes", that is, the
		// set of state from which N\C can force C to reach phi with probability 1 while remaining in phi1
		
		buildCoalitions(csg, coalition, min1);
		timerProb1 = System.currentTimeMillis();
		
		if (!bounded && precomp && prob1) {
			if (!genStrat) {
				yes.or(AF(csg, remain, target));
			} else {
				mainLog.println("Disabling Prob1 precomputation to allow strategy generation");
			}
		}
		timerProb1 = System.currentTimeMillis() - timerProb1;
		yes.or(target);
		
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
			mainLog.println("Precomputation took " + res.timePre + " seconds.");
		return res;
	}

	/**
	 * Compute reachability probabilities using value iteration.
	 * 
	 * @param csg The CSG
	 * @param no Probability 0 states
	 * @param yes Probability 1 states
	 * @param limit Bound
	 * @param bounded Is the problem (step) bounded?
	 * @param min Min or max probabilities for player 1 (true=min, false=max)
	 **/
	public ModelCheckerResult computeReachProbsValIter(CSG<Double> csg, BitSet no, BitSet yes, int limit, boolean bounded, boolean min) throws PrismException
	{
		if (genStrat && bounded) {
			throw new PrismNotSupportedException("Strategy synthesis for bounded properties is not supported yet");
		}
		LpSolve lp;
		ArrayList<ArrayList<Double>> mgame;
		List<List<List<Map<BitSet, Double>>>> lstrat = null;
		List<Map<BitSet, Double>> kstrat = null;
		Map<Integer, BitSet> mmap = null;
		BitSet known = new BitSet();
		double[] nsol = new double[csg.getNumStates()];
		double[] ntmp = new double[csg.getNumStates()];
		long timer;
		int i, k, s;
		boolean done = false;
		if (genStrat) {
			mmap = new HashMap<Integer, BitSet>();
			kstrat = new ArrayList<Map<BitSet, Double>>();
			// player -> iteration -> state -> indexes -> value
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
		} catch (Exception e) {
			e.printStackTrace();
			throw new PrismException(e.toString());
		}
		known.or(no);
		known.or(yes);
		for (s = 0; s < csg.getNumStates(); s++) {
			nsol[s] = ntmp[s] = no.get(s) ? 0.0 : yes.get(s) ? 1.0 : 0.0;
		}
		k = 0;
		while (!done) {
			for (s = 0; s < csg.getNumStates(); s++) {
				if (!known.get(s)) {
					mgame = buildMatrixGame(csg, null, mmap, ntmp, s, min);
					nsol[s] = val(lp, mgame, kstrat, mmap, s, false, min);
					// player -> iteration -> state -> indexes -> value
					if (genStrat) {
						updateStrategy(kstrat, lstrat, k, s, bounded);
					}
				} else if (genStrat) {
					lstrat.get(0).get(0).set(s, null);
				}
			}
			k++;
			done = PrismUtils.doublesAreClose(nsol, ntmp, termCritParam, termCrit == TermCrit.RELATIVE);
			if (!done && k == maxIters) {
				String msg = "Value iteration did not converge within " + k + " iterations";
				if (errorOnNonConverge)
					throw new PrismException(msg);
				mainLog.printWarning(msg + "; the values are those of the last iteration");
				done = true;
			} else if (k == limit) {
				done = true;
			} else {
				ntmp = Arrays.copyOf(nsol, nsol.length);
			}
		}
		if (k < maxIters)
			mainLog.println("\nValue iteration converged after " + k + " iterations.");
		timer = System.currentTimeMillis() - timer;
		ModelCheckerResult res = new ModelCheckerResult();
		res.soln = nsol;
		res.numIters = k;
		if (genStrat) {
			CSGStrategy strat = new CSGStrategy(csg, lstrat, no, yes, new BitSet());
			strat.setCoalitionActions(actionIndexes[min ? 1 : 0]);
			strat.setComplemented(negatedPath);
			res.strat = strat;
			if (limit == maxIters)
				checkZeroSumStrategy(strat, null, nsol[csg.getFirstInitialState()], min);
		}
		res.timeTaken = timer / 1000.0;
		return res;
	}

	/**
	 * Compute instantaneous expected rewards,
	 * i.e. compute the min/max expected reward of the states after {@code k} steps.
	 * 
	 * @param csg The CSG
	 * @param csgRewards The rewards
	 * @param coalition The coalition of players which define player 1
	 * @param k Bound
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 */
	public ModelCheckerResult computeInstantaneousRewards(CSG<Double> csg, CSGRewards<Double> csgRewards, Coalition coalition, int k, boolean min1, boolean min2)
			throws PrismException
	{
		// TODO: confirm that the case min1==min2 is not handled  
		LpSolve lp;
		ModelCheckerResult res = new ModelCheckerResult();
		ArrayList<ArrayList<Double>> mgame = new ArrayList<ArrayList<Double>>();
		List<Map<BitSet, Double>> kstrat = (genStrat) ? new ArrayList<Map<BitSet, Double>>() : null;
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
		} catch (Exception e) {
			e.printStackTrace();
			throw new PrismException(e.toString());
		}
		for (i = 0; i < k; i++) {
			for (int s = 0; s < csg.getNumStates(); s++) {
				mgame.clear();
				mgame = buildMatrixGame(csg, null, null, nsol, s, min1);
				try {
					nsoln2[s] = val(lp, mgame, kstrat, null, s, true, min1);
				} catch (Exception e) {
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

	/**
	 * Compute expected cumulative rewards,
	 * i.e. compute the min/max expected reward accumulated within {@code k} steps.
	 * 
	 * @param csg The CSG
	 * @param csgRewards The rewards
	 * @param coalition The coalition of players which define player 1
	 * @param k Bound
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 * @param genAdv Whether or not to generate a strategy
	 */
	public ModelCheckerResult computeCumulativeRewards(CSG<Double> csg, CSGRewards<Double> csgRewards, Coalition coalition, int k, boolean min1, boolean min2, boolean genAdv)
			throws PrismException
	{
		// TODO: confirm that the case min1==min2 is not handled  
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

	/**
	 * Compute expected total rewards,
	 * i.e. compute the min/max expected reward accumulated.
	 * 
	 * @param csg The CSG
	 * @param rewards The rewards
	 * @param coalition The coalition of players which define player 1
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 */
	public ModelCheckerResult computeTotalRewards(CSG<Double> csg, CSGRewards<Double> rewards, boolean min1, boolean min2, Coalition coalition) throws PrismException
	{
		// TODO: confirm that the case min1==min2 is not handled  
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

		for (s = 0; s < n; s++) {
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
				} else if (rew == Double.POSITIVE_INFINITY || rew == Double.NEGATIVE_INFINITY)
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

	/**
	 * Compute expected reachability rewards.
	 * i.e. compute the min/max reward accumulated to reach a state in {@code target}.
	 * 
	 * @param csg The CSG
	 * @param coalition The coalition of players which define player 1
	 * @param rewards The rewards
	 * @param target Target states
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 * @param unreachingSemantics How to handle paths not reaching target
	 */
	public ModelCheckerResult computeReachRewards(CSG<Double> csg, CSGRewards<Double> rewards, BitSet target, int unreachingSemantics, boolean min1, boolean min2,
			Coalition coalition) throws PrismException
	{
		// TODO: confirm that the case min1==min2 is not handled  
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

	/**
	 * Compute expected reachability rewards (infinite if not reaching target).
	 * i.e. compute the min/max reward accumulated to reach a state in {@code target}.
	 * 
	 * @param csg The CSG
	 * @param coalition The coalition of players which define player 1
	 * @param rewards The rewards
	 * @param target Target states
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 */
	public ModelCheckerResult computeReachRewardsInfinity(CSG<Double> csg, Coalition coalition, CSGRewards<Double> rewards, BitSet target, boolean min1, boolean min2)
			throws PrismException
	{
		// TODO: confirm that the case min1==min2 is not handled  
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

			// Computes the value when rewards are nonzero
			switch (solnMethod) {
			case VALUE_ITERATION:
				// (no strategy needed for the over-approximation)
				boolean genStratSave = genStrat;
				genStrat = false;
				try {
					init = computeReachRewardsValIter(csg, replaceZeroRewards(rewards, epsilon), target, null, inf, null, maxIters, false, min1).soln;
				} finally {
					genStrat = genStratSave;
				}
				break;
			default:
				throw new PrismException("Unknown CSG solution method " + solnMethod);
			}

			// Set the value iteration result to be the initial solution for the
			// next part in which "proper" zero rewards are used

			timerApprox = System.currentTimeMillis() - timerApprox;

			if (verbosity >= 1) {
				mainLog.println(
						"Computed an over-approximation of the solution (in " + timerApprox / 1000.0 + " seconds), this will now be used to get the solution");
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

	/**
	 * Compute expected reachability rewards (cumnulated value if not reaching target).
	 * i.e. compute the min/max reward accumulated to reach a state in {@code target}.
	 * 
	 * @param csg The CSG
	 * @param coalition The coalition of players which define player 1
	 * @param rewards The rewards
	 * @param target Target states
	 * @param min1 Min or max probabilities for player 1 (true=min, false=max)
	 * @param min2 Min or max probabilities for player 2 (true=min, false=max)
	 */
	public ModelCheckerResult computeReachRewardsCumulative(CSG<Double> csg, Coalition coalition, CSGRewards<Double> rewards, BitSet target, boolean min1, boolean min2,
			boolean genAdv) throws PrismException
	{
		// TODO: confirm that the case min1==min2 is not handled  
		ModelCheckerResult res = new ModelCheckerResult();
		BitSet inf;
		double rew;
		long timer, timerPrecomp;
		int nonzero, infty, n, s, t, numTarget, numInf;
		;

		if (verbosity >= 1)
			mainLog.println("\nStarting expected reachability...");
		timer = System.currentTimeMillis();

		n = csg.getNumStates();

		BitSet zero = new BitSet();
		BitSet all = new BitSet();
		all.set(0, n);

		for (s = 0; s < n; s++) {
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
				} else if (rew == Double.POSITIVE_INFINITY || rew == Double.NEGATIVE_INFINITY)
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

	/**
	 * Compute expected reachability rewards using value iteration.
	 * 
	 * @param csg The CSG
	 * @param rewards The rewards
	 * @param target Target states
	 * @param known States whose value is known
	 * @param inf States whose value is infinite
	 * @param limit Bound
	 * @param bounded Is the problem (step) bounded?
	 * @param min Min or max probabilities for player 1 (true=min, false=max)
	 **/
	public ModelCheckerResult computeReachRewardsValIter(CSG<Double> csg, CSGRewards<Double> rewards, BitSet target, BitSet known, BitSet inf, double init[], int limit,
			boolean bounded, boolean min) throws PrismException
	{
		if (genStrat && bounded) {
			throw new PrismNotSupportedException("Strategy synthesis for bounded properties is not supported yet");
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
		if (genStrat) {
			mmap = new HashMap<Integer, BitSet>();
			kstrat = new ArrayList<Map<BitSet, Double>>();
			// player -> iteration -> state -> indexes -> value
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
		} catch (Exception e) {
			e.printStackTrace();
			throw new PrismException(e.toString());
		}
		if (init != null) {
			if (known != null) {
				for (i = 0; i < csg.getNumStates(); i++)
					nsol[i] = ntmp[i] = known.get(i) ? init[i] : target.get(i) ? 0.0 : inf.get(i) ? Double.POSITIVE_INFINITY : init[i];
			} else {
				for (i = 0; i < csg.getNumStates(); i++)
					nsol[i] = ntmp[i] = target.get(i) ? 0.0 : inf.get(i) ? Double.POSITIVE_INFINITY : init[i];
			}
		} else {
			for (i = 0; i < csg.getNumStates(); i++)
				nsol[i] = ntmp[i] = target.get(i) ? 0.0 : inf.get(i) ? Double.POSITIVE_INFINITY : 0.0;
		}
		unknown.set(0, csg.getNumStates());
		unknown.andNot(target);
		unknown.andNot(inf);
		k = 0;
		while (!done) {
			for (s = 0; s < csg.getNumStates(); s++) {
				if (unknown.get(s)) {
					mgame = buildMatrixGame(csg, rewards, mmap, ntmp, s, min);
					nsol[s] = val(lp, mgame, kstrat, mmap, s, true, min);
					nsol[s] += rewards.getStateReward(s);
					if (genStrat) {
						// player -> iteration -> state -> indexes -> value
						updateStrategy(kstrat, lstrat, k, s, bounded);
					}
				}
			}
			k++;
			done = PrismUtils.doublesAreClose(nsol, ntmp, termCritParam, termCrit == TermCrit.RELATIVE);
			if (!done && k == maxIters) {
				String msg = "Value iteration did not converge within " + k + " iterations";
				if (errorOnNonConverge)
					throw new PrismException(msg);
				mainLog.printWarning(msg + "; the values are those of the last iteration");
				done = true;
			} else if (k == limit) {
				done = true;
			} else {
				ntmp = Arrays.copyOf(nsol, nsol.length);
			}
		}
		if (k < maxIters)
			mainLog.println("\nValue iteration converged after " + k + " iterations.");
		timer = System.currentTimeMillis() - timer;
		res.soln = nsol;
		res.numIters = k;
		if (genStrat) {
			CSGStrategy strat = new CSGStrategy(csg, lstrat, new BitSet(), target, inf);
			strat.setCoalitionActions(actionIndexes[min ? 1 : 0]);
			res.strat = strat;
			if (limit == maxIters)
				checkZeroSumStrategy(strat, rewards, nsol[csg.getFirstInitialState()], min);
		}
		res.timeTaken = timer / 1000.0;
		return res;
	}

	/**
	 * Deals with two-player bounded probabilistic reachability formulae
	 *
	 * @param csg The CSG
	 * @param coalitions A list of two coalitions
	 * @param exprs The list of objectives
	 * @param targets The list of sets of target states
	 * @param remain The list of sets of states we need to remain in (in case of until)
	 * @param bounds The list of the objectives' bounds (if applicable)
	 * @param min Whether we're minimising for the first coalition
	 * @return
	 * @throws PrismException
	 */
	public ModelCheckerResult computeProbBoundedEquilibria(CSG<Double> csg, List<Coalition> coalitions, List<ExpressionTemporal> exprs, BitSet[] targets,
			BitSet[] remain, int[] bounds, int eqType, int crit, boolean min) throws PrismException
	{
		ModelCheckerResult res = new ModelCheckerResult();
		CSGModelCheckerEquilibria csgeq = new CSGModelCheckerEquilibria(this);
		csgeq.inheritSettings(this);
		res = csgeq.computeBoundedEquilibria(csg, coalitions, null, exprs, targets, remain, bounds, eqType, crit, min);
		return res;
	}

	/**
	 * Deals with two-player probabilistic reachability formulae
	 * 
	 * @param csg The CSG
	 * @param coalitions A list of two coalitions
	 * @param targets The list of sets of target states
	 * @param remain The list of sets of states we need to remain in (in case of until)
	 * @param min Whether we're minimising for the first coalition
	 * @return
	 * @throws PrismException
	 */
	public ModelCheckerResult computeProbReachEquilibria(CSG<Double> csg, List<Coalition> coalitions, BitSet[] targets, BitSet[] remain, int eqType, int crit, boolean min)
			throws PrismException
	{
		ModelCheckerResult res = new ModelCheckerResult();
		CSGModelCheckerEquilibria csgeq = new CSGModelCheckerEquilibria(this);
		csgeq.inheritSettings(this);
		res = csgeq.computeReachEquilibria(csg, coalitions, null, targets, remain, eqType, crit, min);
		return res;
	}

	/**
	 * Deals with two-player bounded reachability rewards
	 * 
	 * @param csg The CSG
	 * @param coalitions A list of two coalitions
	 * @param rewards The list of reward structures
	 * @param exprs The list of objectives
	 * @param bounds The list of the objectives' bounds (if applicable)
	 * @param min Whether we're minimising for the first coalition
	 * @return
	 * @throws PrismException
	 */
	public ModelCheckerResult computeRewBoundedEquilibria(CSG<Double> csg, List<Coalition> coalitions, List<CSGRewards<Double>> rewards, List<ExpressionTemporal> exprs,
			int[] bounds, int eqType, int crit, boolean min) throws PrismException
	{
		ModelCheckerResult res = new ModelCheckerResult();
		CSGModelCheckerEquilibria csgeq = new CSGModelCheckerEquilibria(this);
		csgeq.inheritSettings(this);
		res = csgeq.computeBoundedEquilibria(csg, coalitions, rewards, exprs, null, null, bounds, eqType, crit, min);
		return res;
	}

	/**
	 * Deals with two-player reachability rewards formulae
	 * 
	 * @param csg The CSG
	 * @param coalitions A list of two coalitions
	 * @param rewards The list of reward structures
	 * @param targets The list of sets of target states
	 * @param min Whether we're minimising for the first coalition
	 * @return
	 * @throws PrismException
	 */
	public ModelCheckerResult computeRewReachEquilibria(CSG<Double> csg, List<Coalition> coalitions, List<CSGRewards<Double>> rewards, BitSet[] targets, int eqType, int crit, boolean min)
			throws PrismException
	{
		ModelCheckerResult res = new ModelCheckerResult();
		CSGModelCheckerEquilibria csgeq = new CSGModelCheckerEquilibria(this);
		csgeq.inheritSettings(this);
		res = csgeq.computeReachEquilibria(csg, coalitions, rewards, targets, null, eqType, crit, min);
		return res;
	}

	/**
	 * Equilibria for more than two coalitions: reachability, until and reachability rewards, bounded or not
	 * (including different bounds and mixtures of bounded and unbounded objectives).
	 * 
	 * @param csg The CSG
	 * @param coalitions The list of coalitions
	 * @param rewards The list of reward structures (null for probabilistic objectives)
	 * @param exprs The objectives (temporal operators)
	 * @param bounded The objectives which are bounded
	 * @param targets The target states of each objective (null for cumulative/instantaneous rewards)
	 * @param remain The states to remain in (until), or null
	 * @param bounds The bounds of the bounded objectives
	 * @param min Whether minimising
	 */
	public ModelCheckerResult computeMultiEquilibria(CSG<Double> csg, List<Coalition> coalitions, List<CSGRewards<Double>> rewards, List<ExpressionTemporal> exprs,
			BitSet bounded, BitSet[] targets, BitSet[] remain, int[] bounds, int eqType, int crit, boolean min) throws PrismException
	{
		CSGModelCheckerEquilibria csgeq = new CSGModelCheckerEquilibria(this);
		csgeq.inheritSettings(this);
		return csgeq.computeMultiEquilibria(csg, coalitions, rewards, exprs, bounded, targets, remain, bounds, eqType, crit, min);
	}

	/**
	 * Deals with computing mixed bounded and unbounded equilibria.
	 * 
	 * @param csg The CSG
	 * @param coalitions The list of coalitions
	 * @param rewards The list of reward structures
	 * @param exprs The list of objectives
	 * @param bounded Index of the objectives which are bounded
	 * @param targets The list of sets of target states
	 * @param remain The list of sets of states we need to remain in (in case of until)
	 * @param bounds The list of the objectives' bounds (if applicable)
	 * @param min Whether we're minimising for the first coalition
	 * @return
	 * @throws PrismException
	 */
	public ModelCheckerResult computeMixedEquilibria(CSG<Double> csg, List<Coalition> coalitions, List<CSGRewards<Double>> rewards, List<ExpressionTemporal> exprs,
			BitSet bounded, BitSet[] targets, BitSet[] remain, int[] bounds, int eqType, int crit, boolean min) throws PrismException
	{
		ModelCheckerResult res = new ModelCheckerResult();
		CSGModelCheckerEquilibria csgeq = new CSGModelCheckerEquilibria(this);
		csgeq.inheritSettings(this);
		LTLModelChecker ltlmc = new LTLModelChecker(this);
		LTLProduct<CSG<Double>> product;
		List<CSGRewards<Double>> newrewards = new ArrayList<>();
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

		AcceptanceType[] allowedAcceptance = { AcceptanceType.RABIN, AcceptanceType.REACH, AcceptanceType.BUCHI, AcceptanceType.STREETT,
				AcceptanceType.GENERIC };

		/*
		csg.exportToDotFile(path + "/model.dot");
		*/

		if (rew) {
			BitSet all = new BitSet();
			all.set(0, csg.getNumStates());
			DA<BitSet, AcceptanceRabin> da = null;
			Vector<BitSet> labelBS = new Vector<BitSet>();
			labelBS.add(0, all);
			targets[index] = all;
			switch (exprs.get(index).getOperator()) {
			case ExpressionTemporal.R_I:
				da = constructDRAForInstant("L0", new IntegerBound(null, false, bounds[index] + 1, false));
				break;
			case ExpressionTemporal.R_C:
				da = constructDRAForInstant("L0", new IntegerBound(null, false, bounds[index], false));
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
		} else {
			product = ltlmc.constructProductCSG(this, csg, exprs.get(index), null, allowedAcceptance);
		}

		((ModelExplicit<Double>) product.productModel).clearInitialStates();
		((ModelExplicit<Double>) product.productModel).addInitialState(product.getModelState(csg.getFirstInitialState()));

		/*
		product.productModel.exportToDotFile(path + "/product.dot");
		*/

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
				CSGRewards<Double> reward = new CSGRewardsSimple<>(product.productModel.getNumStates());
				if (i != index) {
					for (s = 0; s < product.productModel.getNumStates(); s++) {
						((CSGRewardsSimple<Double>) reward).addToStateReward(s, rewards.get(i).getStateReward(product.getModelState(s)));
						for (t = 0; t < product.productModel.getNumChoices(s); t++) {
							((CSGRewardsSimple<Double>) reward).addToTransitionReward(s, t, rewards.get(i).getTransitionReward(product.getModelState(s), t));
						}
					}
				} else {
					switch (exprs.get(index).getOperator()) {
					case ExpressionTemporal.R_I:
						for (s = 0; s < product.productModel.getNumStates(); s++) {
							for (t = 0; t < product.productModel.getNumChoices(s); t++) {
								for (Iterator<Integer> iter = product.productModel.getSuccessorsIterator(s, t); iter.hasNext();) {
									int u = iter.next();
									if (newtargets[index].get(u)) {
										((CSGRewardsSimple<Double>) reward).addToStateReward(s, rewards.get(i).getStateReward(product.getModelState(s)));
									}
								}
								((CSGRewardsSimple<Double>) reward).addToTransitionReward(s, t, 0.0);
							}
						}
						break;
					case ExpressionTemporal.R_C:
						for (s = 0; s < product.productModel.getNumStates(); s++) {
							if (!newtargets[index].get(s)) {
								((CSGRewardsSimple<Double>) reward).addToStateReward(s, rewards.get(i).getStateReward(product.getModelState(s)));
								for (t = 0; t < product.productModel.getNumChoices(s); t++) {
									((CSGRewardsSimple<Double>) reward).addToTransitionReward(s, t, rewards.get(i).getTransitionReward(product.getModelState(s), t));
								}
							} else {
								((CSGRewardsSimple<Double>) reward).addToStateReward(s, 0.0);
								for (t = 0; t < product.productModel.getNumChoices(s); t++) {
									((CSGRewardsSimple<Double>) reward).addToTransitionReward(s, t, 0.0);
								}
							}
						}
						break;
					}
				}
				newrewards.add(i, reward);
			}
			res = csgeq.computeReachEquilibria(product.productModel, coalitions, newrewards, newtargets, null, eqType, crit, min);
		} else {
			res = csgeq.computeReachEquilibria(product.productModel, coalitions, null, newtargets, newremain, eqType, crit, min);
		}
		return res;
	}

	// Precomputation methods (almost surely)
	
	/*
	 * Auxiliary method (as defined in L. Alfaro and T. Henzinger, Concurrent Omega-Regular Games)
	 */
	public BitSet A(ArrayList<ArrayList<Distribution<Double>>> mdist, BitSet v2, BitSet y)
	{
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

	/*
	 * Auxiliary method (as defined in L. Alfaro and T. Henzinger, Concurrent Omega-Regular Games)
	 */
	public BitSet B(ArrayList<ArrayList<Distribution<Double>>> mdist, BitSet v1, BitSet x)
	{
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

	/*
	 * Auxiliary method for AF (as defined in L. Alfaro and T. Henzinger, Concurrent Omega-Regular Games)
	 */
	public BitSet apreXY(CSG<Double> csg, BitSet x, BitSet y) throws PrismException
	{
		ArrayList<ArrayList<Distribution<Double>>> mdist;
		BitSet result = new BitSet();
		int s;
		for (s = 0; s < csg.getNumStates(); s++) {
			mdist = buildMatrixDist(csg, s);
			result.set(s, B(mdist, A(mdist, new BitSet(), y), x).cardinality() == mdist.get(0).size());
		}
		return result;
	}

	/*
	 * Eventually b
	 */
	public BitSet AF(CSG<Double> csg, BitSet b) throws PrismException
	{
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
	
	/*
	 * Eventually b, while remaining in a
	 */
	public BitSet AF(CSG<Double> csg, BitSet a, BitSet b) throws PrismException
	{
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
				sol1.and(a);
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

	/*
	 * Auxiliary method for AFG (as defined in L. Alfaro and T. Henzinger, Concurrent Omega-Regular Games)
	 */
	public boolean apreXYZ(ArrayList<ArrayList<Distribution<Double>>> mdist, BitSet x, BitSet y, BitSet z)
	{
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

	/*
	 * Auxiliary method for AFG (as defined in L. Alfaro and T. Henzinger, Concurrent Omega-Regular Games)
	 */
	public BitSet apreXYZ(CSG<Double> csg, BitSet x, BitSet y, BitSet z) throws PrismException
	{
		ArrayList<ArrayList<Distribution<Double>>> mdist;
		BitSet result = new BitSet();
		for (int s = 0; s < csg.getNumStates(); s++) {
			mdist = buildMatrixDist(csg, s);
			result.set(s, apreXYZ(mdist, x, y, z));
		}
		return result;
	}

	/*
	 * Eventually globally b
	 */
	public BitSet AFG(CSG<Double> csg, BitSet b) throws PrismException
	{
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


	/*
	 * Almost-sure Buchi (G F b): visit b infinitely often, forced almost-surely. Per
	 * LICS 2000 eq. 3 (S4.1), the strategy-complexity asymmetry reported for concurrent Buchi vs
	 * co-Buchi games (Buchi admits positional optimal strategies, co-Buchi may not) --
	 * unlike AFG (co-Buchi, eq. 4, S4.2), which needs the heavier 3-argument AFpre1/
	 * apreXYZ because "eventually always" has to commit to a single invariant, "infinitely
	 * often" only needs the same 2-argument Apre1/apreXY already used by AF, re-triggered
	 * every outer iteration:
	 *   Y = nu Y. mu X. [ (not B and Apre1(Y,X)) or (B and Pre1(Y)) ]
	 * The B-branch uses plain Pre1(Y) (not Apre1) to bank a completed visit and return to
	 * Y with certainty; the not-B branch uses apreXY (almost-surely progress toward the
	 * next visit while staying in Y).
	 */
	public BitSet AGF(CSG<Double> csg, BitSet b) throws PrismException
	{
		int n = csg.getNumStates();
		BitSet x, y, sol1, sol2, notb;
		x = new BitSet();
		y = new BitSet();
		sol1 = new BitSet();
		sol2 = new BitSet();
		notb = new BitSet();
		notb.set(0, n);
		notb.andNot(b);
		boolean done_x, done_y;
		done_y = false;
		y.set(0, n);
		while (!done_y) {
			done_x = false;
			x.clear();
			while (!done_x) {
				sol1.clear();
				sol1.or(apreXY(csg, x, y));
				sol1.and(notb);
				sol2.clear();
				pre1(csg, y, sol2);
				sol2.and(b);
				sol1.or(sol2);
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

	public boolean pre1(ArrayList<ArrayList<Distribution<Double>>> mdist, BitSet x)
	{
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

	public void pre1(CSG<Double> csg, BitSet x, BitSet sol) throws PrismException
	{
		ArrayList<ArrayList<Distribution<Double>>> mdist;
		for (int s = 0; s < csg.getNumStates(); s++) {
			mdist = buildMatrixDist(csg, s);
			sol.set(s, pre1(mdist, x));
		}
	}

	/*
	 * Globally b
	 */
	public BitSet G(CSG<Double> csg, BitSet b) throws PrismException
	{
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

	/*
	 * Sure-mode reachability: eventually b, forced with certainty (no probability
	 * involved, unlike AF/LF) -- mu X.(B or Pre1(X)), the plain-mu mirror of G's
	 * plain-nu safety fixpoint above. Reuses pre1 directly: unlike Apre1/Lpre1, Pre1
	 * already gives outright exists-forall certainty, so no apreXY/lpreXY-style
	 * extra-argument helper is needed for sure mode.
	 */
	public BitSet SF(CSG<Double> csg, BitSet b) throws PrismException
	{
		BitSet sol1, x;
		sol1 = new BitSet();
		x = new BitSet();
		boolean done_x = false;
		while (!done_x) {
			sol1.clear();
			pre1(csg, x, sol1);
			sol1.or(b);
			done_x = x.equals(sol1);
			x.clear();
			x.or(sol1);
		}
		return x;
	}

	/*
	 * Sure-mode co-Buchi: eventually globally b, forced with certainty. Unlike AFG/LFG,
	 * which need the dedicated apreXYZ/lpreXYZ-based 3-level fixpoint because almost/
	 * limit-sure winning has to re-accumulate probability across infinitely many attempts,
	 * the sure case has no such subtlety: F G B is exactly "sure-reach the maximal
	 * B-invariant" -- the standard automata/game-theory decomposition of co-Buchi as a
	 * reachability fixpoint (SF) composed with a safety fixpoint (G), computed once and
	 * fed in, not iteratively coupled the way genuine Buchi (infinitely-often) needs.
	 * (Soundness: any strategy witnessing F G B must, from the point where G B starts
	 * holding forever, exhibit a B-invariant set closed under the strategy's moves --
	 * every such set is a subset of the *maximal* B-invariant G(csg,b), so s must sure-
	 * reach G(csg,b); conversely sure-reaching G(csg,b) and then following its own
	 * safety-fixpoint witness clearly forces F G B.)
	 * Reuses G and SF directly -- no new fixpoint code needed.
	 */
	public BitSet SFG(CSG<Double> csg, BitSet b) throws PrismException
	{
		return SF(csg, G(csg, b));
	}

	/*
	 * Sure-mode Buchi (G F b): visit b infinitely often, forced with certainty. Unlike
	 * SFG (co-Buchi), which decomposes into two independent fixpoints computed once and
	 * composed (G then SF) -- because "eventually always" only ever needs to commit to a
	 * target invariant once -- genuine "infinitely often" needs the coupled nu Y. mu X.
	 * attractor-of-attractor construction: the classical Buchi-game fixpoint, sure mode's
	 * analogue of AGF/LGF's eq. 3 shape with Apre1/Lpre1 collapsing to plain Pre1 (no
	 * probability to accumulate, so the not-B branch needs only certain progress within
	 * the growing X attractor, not a separate "stay safe" argument the way Apre1(Y,X)
	 * tracks):
	 *   Y = nu Y. mu X. [ (B and Pre1(Y)) or (not B and Pre1(X)) ]
	 */
	public BitSet SGF(CSG<Double> csg, BitSet b) throws PrismException
	{
		int n = csg.getNumStates();
		BitSet x, y, sol1, sol2, notb;
		x = new BitSet();
		y = new BitSet();
		sol1 = new BitSet();
		sol2 = new BitSet();
		notb = new BitSet();
		notb.set(0, n);
		notb.andNot(b);
		boolean done_x, done_y;
		done_y = false;
		y.set(0, n);
		while (!done_y) {
			done_x = false;
			x.clear();
			while (!done_x) {
				sol1.clear();
				pre1(csg, x, sol1);
				sol1.and(notb);
				sol2.clear();
				pre1(csg, y, sol2);
				sol2.and(b);
				sol1.or(sol2);
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

	// Limit-mode operators (as defined in L. Alfaro and T. Henzinger, Concurrent Omega-Regular Games)
	
	public boolean lpreXY(ArrayList<ArrayList<Distribution<Double>>> mdist, BitSet x, BitSet y) 
	{
	    BitSet w = new BitSet();               
	    boolean done_w = false;
	    while (!done_w) {
	        BitSet next = A(mdist, B(mdist, w, x), y);
	        done_w = w.equals(next);
	        w.clear();
	        w.or(next);
	    }
	    return B(mdist, w, x).cardinality() == mdist.get(0).size();
	}
	
	public BitSet lpreXY(CSG<Double> csg, BitSet x, BitSet y) throws PrismException 
	{
	    ArrayList<ArrayList<Distribution<Double>>> mdist;
	    BitSet result = new BitSet();
	    for (int s = 0; s < csg.getNumStates(); s++) {
	        mdist = buildMatrixDist(csg, s);
	        result.set(s, lpreXY(mdist, x, y));
	    }
	    return result;
	}

	public BitSet LF(CSG<Double> csg, BitSet b) throws PrismException
	{
	    int n = csg.getNumStates();
	    BitSet x = new BitSet(), y = new BitSet(), sol1 = new BitSet();
	    y.set(0, n);
	    boolean done_x, done_y = false;
	    while (!done_y) {
	        done_x = false;
	        x.clear();
	        while (!done_x) {
	            sol1.clear();
	            sol1.or(lpreXY(csg, x, y));
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
	
	public boolean lpreXYZ(ArrayList<ArrayList<Distribution<Double>>> mdist, BitSet x, BitSet y, BitSet z) 
	{
	    BitSet v = new BitSet();
	    v.set(0, mdist.size());                
	    boolean done_v = false;
	    while (!done_v) {
	        BitSet ay = A(mdist, B(mdist, v, x), y);
	        BitSet w = new BitSet();
	        boolean done_w = false;
	        while (!done_w) {
	            BitSet az = A(mdist, B(mdist, w, x), z);
	            BitSet sol = (BitSet) az.clone();
	            sol.and(ay);
	            done_w = w.equals(sol);
	            w.clear();
	            w.or(sol);
	        }
	        done_v = v.equals(w);
	        v.clear();
	        v.or(w);
	    }
	    return !(v.isEmpty());
	}
	
	public BitSet lpreXYZ(CSG<Double> csg, BitSet x, BitSet y, BitSet z) throws PrismException 
	{
	    ArrayList<ArrayList<Distribution<Double>>> mdist;
	    BitSet result = new BitSet();
	    for (int s = 0; s < csg.getNumStates(); s++) {
	        mdist = buildMatrixDist(csg, s);
	        result.set(s, lpreXYZ(mdist, x, y, z));
	    }
	    return result;
	}
	
	public BitSet LFG(CSG<Double> csg, BitSet b) throws PrismException 
	{
	    int n = csg.getNumStates();
	    BitSet x = new BitSet(), y = new BitSet(), z = new BitSet();
	    BitSet sol1 = new BitSet(), sol2 = new BitSet();
	    boolean done_x, done_y, done_z = false;
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
	                sol1 = lpreXYZ(csg, x, y, z);
	                sol1.and(b);
	                sol2.or(b);
	                sol2.flip(0, n);
	                sol2.and(lpreXY(csg, x, z));
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
	
	/*
	 * Limit-sure Buchi (G F b): visit b infinitely often, forced in the limit. Mirrors
	 * AGF exactly, with lpreXY (Lpre1) substituted for apreXY (Apre1), as LF/LFG do
	 * relative to AF/AFG:
	 *   Y = nu Y. mu X. [ (not B and Lpre1(Y,X)) or (B and Pre1(Y)) ]
	 */
	public BitSet LGF(CSG<Double> csg, BitSet b) throws PrismException
	{
		int n = csg.getNumStates();
		BitSet x, y, sol1, sol2, notb;
		x = new BitSet();
		y = new BitSet();
		sol1 = new BitSet();
		sol2 = new BitSet();
		notb = new BitSet();
		notb.set(0, n);
		notb.andNot(b);
		boolean done_x, done_y;
		done_y = false;
		y.set(0, n);
		while (!done_y) {
			done_x = false;
			x.clear();
			while (!done_x) {
				sol1.clear();
				sol1.or(lpreXY(csg, x, y));
				sol1.and(notb);
				sol2.clear();
				pre1(csg, y, sol2);
				sol2.and(b);
				sol1.or(sol2);
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

	// General Rabin-chain operators (m pairs), as defined in L. Alfaro and T. Henzinger,
	// Concurrent Omega-Regular Games, S5.2 "General Rabin-chain games" and Appendix 6
	// "Computation of Predecessor Operators" -- generalizes SF/SFG/SGF, AF/AFG/AGF,
	// LF/LFG/LGF (all special/degenerate 1-pair cases of the construction below) to an
	// arbitrary number m of Rabin pairs. Pairs are given directly as BitSets over the CSG's
	// states (no LTL-to-automaton product).
	//
	// A Rabin-chain condition with m pairs is defined by a chain of state-sets
	// S = U_0 superseteq U_1 superseteq ... superseteq U_{2m-1} superseteq U_{2m} = emptyset,
	// with pair i (i=0..m-1) given by (L_i, K_i) = (U_{2i+1}, U_{2i}) -- i.e. Fin(L_i) and
	// Inf(K_i), matching AcceptanceRabin.RabinPair's (L,K) convention -- and winning
	// condition rho = OR_i (GF K_i and not GF L_i). Colors are B_k = U_k \ U_{k+1}.

	/** The three "modes" of the one-step predecessor combinator <alpha/beta/gamma, Y, X>
	 * from LICS 2000 Appendix 6. GAMMA and ALPHA only ever occur as the base (first,
	 * index-0) triple of a theta sequence (matching how ARpre_{2m-1,1}^m / LRpre_{2m-1,1}^m are 
	 * always the sole seed of the recursion) --
	 * GAMMA is Apre1's base (almost mode), ALPHA is Lpre1's base (limit mode); BETA is
	 * used for every triple appended afterwards, alternating with ALPHA as the chain
	 * recurses down through the remaining pairs. */
	private static final int RPRE_GAMMA = 0;
	private static final int RPRE_ALPHA = 1;
	private static final int RPRE_BETA = 2;

	/** One <mode, Y, X> triple of the theta-sequence combinator (Y, X are STATE sets --
	 * the current values of the enclosing state-level fixpoint variables). */
	private static class RabinPreTriple
	{
		int mode;
		BitSet Y;
		BitSet X;
		RabinPreTriple(int mode, BitSet Y, BitSet X)
		{
			this.mode = mode;
			this.Y = Y;
			this.X = X;
		}
	}

	/**
	 * This triple's own one-step term A_Y(B_X(var)) (or, for a GAMMA base triple with no
	 * bound variable, the constant A_Y(emptyset)) -- the elementary building block that
	 * apreXY/lpreXY/apreXYZ/lpreXYZ are themselves built from (see A/B above).
	 */
	private BitSet rpreTerm(ArrayList<ArrayList<Distribution<Double>>> mdist, RabinPreTriple t, BitSet var)
	{
		if (t.mode == RPRE_GAMMA) {
			return A(mdist, new BitSet(), t.Y);
		}
		return A(mdist, B(mdist, var, t.X), t.Y);
	}

	/**
	 * The single shared body C(theta): fold every triple's own term into one expression,
	 * left to right, via OR (if that triple's mode is ALPHA) or AND (if BETA) -- exactly
	 * C(theta o <alpha,Y,X>) = C(theta) or A_Y(B_X(W_X)) / C(theta o <beta,Y,X>) =
	 * C(theta) and A_Y(B_X(V_X)) from Appendix 6, unrolled over the whole sequence.
	 * {@code vars[i]} holds triple i's current bound-variable value (ignored for a GAMMA
	 * triple, which contributes only its constant term).
	 */
	private BitSet rpreCombine(ArrayList<ArrayList<Distribution<Double>>> mdist, RabinPreTriple[] theta, BitSet[] vars)
	{
		BitSet acc = rpreTerm(mdist, theta[0], vars[0]);
		for (int i = 1; i < theta.length; i++) {
			BitSet term = rpreTerm(mdist, theta[i], vars[i]);
			acc = (BitSet) acc.clone();
			if (theta[i].mode == RPRE_ALPHA) {
				acc.or(term);
			} else {
				acc.and(term);
			}
		}
		return acc;
	}

	/** Smallest triple index with an actual bound variable (mode != GAMMA); GAMMA can
	 * only be theta[0], so this is 0 (limit-mode base), 1 (almost-mode base with a longer
	 * chain), or -1 (the single-triple almost-mode base alone, i.e. plain Apre1). */
	private int rpreInnermostVarLevel(RabinPreTriple[] theta)
	{
		for (int i = 0; i < theta.length; i++) {
			if (theta[i].mode != RPRE_GAMMA) {
				return i;
			}
		}
		return -1;
	}

	/**
	 * Resolve triple `level`'s bound variable to its own fixpoint (mu/least if ALPHA,
	 * nu/greatest if BETA), given theta[level+1..] already fixed in `vars` by the
	 * enclosing (outer) call(s). At the innermost variable level, each candidate is
	 * obtained by re-evaluating the single shared body rpreCombine (fresh every
	 * iteration, folding in every other triple's *current* value); at every level above
	 * it, each candidate is simply the just-converged value returned by the recursive
	 * call one level down -- the "v <- converged w" pattern apreXYZ/lpreXYZ use for their
	 * 2-triple case, generalised to an arbitrary chain length.
	 */
	private BitSet rpreResolve(ArrayList<ArrayList<Distribution<Double>>> mdist, RabinPreTriple[] theta, int level, int innermost, BitSet[] vars)
	{
		RabinPreTriple t = theta[level];
		boolean isMu = (t.mode == RPRE_ALPHA);
		BitSet var = new BitSet();
		if (!isMu) {
			var.set(0, mdist.size());
		}
		vars[level] = var;
		boolean done = false;
		while (!done) {
			BitSet next;
			if (level == innermost) {
				next = rpreCombine(mdist, theta, vars);
			} else {
				next = rpreResolve(mdist, theta, level - 1, innermost, vars);
			}
			done = vars[level].equals(next);
			vars[level] = next;
		}
		return vars[level];
	}

	/**
	 * Evaluate, for one CSG state (given its one-step matrix game mdist), whether that
	 * state satisfies the predecessor predicate for the theta sequence -- i.e. whether it
	 * belongs to <theta> per Appendix 6's final membership check. If the outermost
	 * (last) triple is BETA, membership is "value nonempty" (matches apreXYZ/lpreXYZ's
	 * {@code !v.isEmpty()}); otherwise (GAMMA alone, or ALPHA outermost) membership is
	 * "B_X(value) == Gamma_2(s)" (matches apreXY/lpreXY's cardinality check).
	 */
	private boolean rpre(ArrayList<ArrayList<Distribution<Double>>> mdist, RabinPreTriple[] theta)
	{
		int n = theta.length;
		int innermost = rpreInnermostVarLevel(theta);
		RabinPreTriple outer = theta[n - 1];
		if (innermost == -1) {
			// Single GAMMA triple: plain Apre1(Y,X).
			BitSet value = A(mdist, new BitSet(), theta[0].Y);
			return B(mdist, value, theta[0].X).cardinality() == mdist.get(0).size();
		}
		BitSet[] vars = new BitSet[n];
		BitSet value = rpreResolve(mdist, theta, n - 1, innermost, vars);
		if (outer.mode == RPRE_BETA) {
			return !value.isEmpty();
		} else {
			return B(mdist, value, outer.X).cardinality() == mdist.get(0).size();
		}
	}

	/** State-set (BitSet over all CSG states) version of {@link #rpre}. */
	private BitSet rpre(CSG<Double> csg, RabinPreTriple[] theta) throws PrismException
	{
		ArrayList<ArrayList<Distribution<Double>>> mdist;
		BitSet result = new BitSet();
		for (int s = 0; s < csg.getNumStates(); s++) {
			mdist = buildMatrixDist(csg, s);
			result.set(s, rpre(mdist, theta));
		}
		return result;
	}

	/**
	 * Build the theta sequence for lambdaRpre_{k,1}^m (k=0..2m-1), per the recursive
	 * definitions of LICS 2000 S5.2: the top pair (i=m-1) seeds the sequence with
	 * <gamma/alpha, Y_2m, X_{2m-1}> (index 2m-1) then, for the co-index 2m-2, appends
	 * <beta, Y_{2m-2}, X_{2m-1}>; every subsequent (lower) pair i=m-2 downto 0 appends
	 * <alpha, Y_{2i+2}, X_{2i+1}> (for the odd index 2i+1) then <beta, Y_{2i}, X_{2i+1}>
	 * (for the even index 2i), stopping as soon as the requested k has been reached.
	 * `var[j]` is the CURRENT value of the state-level fixpoint variable at subscript j
	 * (Y_j for even j, X_j for odd j), with var[2*m] the extra outermost nu.
	 */
	private RabinPreTriple[] rpreBuildTheta(boolean almost, BitSet[] var, int m, int k)
	{
		int baseMode = almost ? RPRE_GAMMA : RPRE_ALPHA;
		ArrayList<RabinPreTriple> theta = new ArrayList<RabinPreTriple>();
		theta.add(new RabinPreTriple(baseMode, var[2 * m], var[2 * m - 1]));
		if (k == 2 * m - 1) {
			return theta.toArray(new RabinPreTriple[0]);
		}
		theta.add(new RabinPreTriple(RPRE_BETA, var[2 * m - 2], var[2 * m - 1]));
		if (k == 2 * m - 2) {
			return theta.toArray(new RabinPreTriple[0]);
		}
		for (int i = m - 2; i >= 0; i--) {
			theta.add(new RabinPreTriple(RPRE_ALPHA, var[2 * i + 2], var[2 * i + 1]));
			if (k == 2 * i + 1) {
				return theta.toArray(new RabinPreTriple[0]);
			}
			theta.add(new RabinPreTriple(RPRE_BETA, var[2 * i], var[2 * i + 1]));
			if (k == 2 * i) {
				return theta.toArray(new RabinPreTriple[0]);
			}
		}
		throw new IllegalArgumentException("Rabin-chain color index " + k + " out of range for m=" + m);
	}

	/**
	 * State-level general Rabin-chain fixpoint (formula 14): winning states for
	 * almost-sure (almost=true) or limit-sure (almost=false) satisfaction of the m-pair
	 * chain condition rho = OR_{i=0}^{m-1} (GF K_i and not GF L_i), given as parallel
	 * lists of L_i, K_i (Fin/Inf sets, matching AcceptanceRabin.RabinPair's convention;
	 * L.get(m-1)/K.get(m-1) is the *top* pair, i.e. i=m-1, and the chain runs
	 * S = U_0 = K.get(0) superseteq U_1 = L.get(0) superseteq ... superseteq
	 * U_{2m-1} = L.get(m-1) superseteq U_{2m} = emptyset -- callers must supply a
	 * genuinely nested chain, not an arbitrary Rabin condition). Structurally this is
	 * nu Y_{2m}. mu X_{2m-1}. nu Y_{2m-2}. ... mu X_1. nu Y_0. [ OR_k (B_k and
	 * lambdaRpre_{k,1}^m) ], computed via the same "innermost level does the real
	 * combine-and-check, every level above it just copies the converged value one level
	 * down" pattern used by AGF/AFG/LGF/LFG's nested while-loops, written as recursion
	 * so that it generalises to any chain depth.
	 */
	public BitSet rabinChainWin(CSG<Double> csg, List<BitSet> L, List<BitSet> K, boolean almost) throws PrismException
	{
		int m = L.size();
		if (K.size() != m || m < 1) {
			throw new PrismException("Rabin-chain L/K lists must be non-empty and of equal length");
		}
		int n = csg.getNumStates();
		// U_{2i} = K_i, U_{2i+1} = L_i, U_{2m} = empty; colors B_k = U_k \ U_{k+1}.
		BitSet[] U = new BitSet[2 * m + 1];
		for (int i = 0; i < m; i++) {
			U[2 * i] = K.get(i);
			U[2 * i + 1] = L.get(i);
		}
		U[2 * m] = new BitSet();
		final BitSet[] B = new BitSet[2 * m];
		for (int k = 0; k < 2 * m; k++) {
			B[k] = (BitSet) U[k].clone();
			B[k].andNot(U[k + 1]);
		}
		BitSet[] var = new BitSet[2 * m + 1];
		return rabinChainLevel(csg, n, almost, m, B, var, 2 * m);
	}

	private BitSet rabinChainLevel(CSG<Double> csg, int n, boolean almost, int m, BitSet[] B, BitSet[] var, int j) throws PrismException
	{
		boolean isNu = (j % 2 == 0);
		BitSet v = new BitSet();
		if (isNu) {
			v.set(0, n);
		}
		var[j] = v;
		boolean done = false;
		while (!done) {
			BitSet next;
			if (j == 0) {
				next = rabinChainBody(csg, almost, m, B, var);
			} else {
				next = rabinChainLevel(csg, n, almost, m, B, var, j - 1);
			}
			done = var[j].equals(next);
			var[j] = next;
		}
		return var[j];
	}

	private BitSet rabinChainBody(CSG<Double> csg, boolean almost, int m, BitSet[] B, BitSet[] var) throws PrismException
	{
		BitSet result = new BitSet();
		for (int k = 0; k < 2 * m; k++) {
			if (B[k].isEmpty()) {
				continue;
			}
			RabinPreTriple[] theta = rpreBuildTheta(almost, var, m, k);
			BitSet term = rpre(csg, theta);
			term.and(B[k]);
			result.or(term);
		}
		return result;
	}

	/**
	 * Sure-mode general Rabin-chain fixpoint: unlike almost/limit, sure winning needs no
	 * one-step matrix-game combinator at all (Pre1 already gives outright certainty), so
	 * the whole construction collapses to a direct analogue of SGF's coupled
	 * attractor-of-attractor, generalized from 1 to m pairs: mu X_{2m-1}. nu Y_{2m-2}.
	 * ... mu X_1. nu Y_0. [ OR_k (B_k and Pre1(var_k)) ], with no extra outer Y_{2m}
	 * and each color's
	 * predecessor being plain Pre1 of *that color's own* bound variable (Y_k if k even,
	 * X_k if k odd) rather than a cross-pair combinator.
	 */
	public BitSet rabinChainWinSure(CSG<Double> csg, List<BitSet> L, List<BitSet> K) throws PrismException
	{
		int m = L.size();
		if (K.size() != m || m < 1) {
			throw new PrismException("Rabin-chain L/K lists must be non-empty and of equal length");
		}
		int n = csg.getNumStates();
		BitSet[] U = new BitSet[2 * m + 1];
		for (int i = 0; i < m; i++) {
			U[2 * i] = K.get(i);
			U[2 * i + 1] = L.get(i);
		}
		U[2 * m] = new BitSet();
		final BitSet[] B = new BitSet[2 * m];
		for (int k = 0; k < 2 * m; k++) {
			B[k] = (BitSet) U[k].clone();
			B[k].andNot(U[k + 1]);
		}
		final BitSet[] var = new BitSet[2 * m];
		return rabinChainSureLevel(csg, n, B, var, 2 * m - 1);
	}

	private BitSet rabinChainSureLevel(CSG<Double> csg, int n, BitSet[] B, BitSet[] var, int j) throws PrismException
	{
		boolean isNu = (j % 2 == 0);
		BitSet v = new BitSet();
		if (isNu) {
			v.set(0, n);
		}
		var[j] = v;
		boolean done = false;
		while (!done) {
			BitSet next;
			if (j == 0) {
				next = rabinChainSureBody(csg, B, var);
			} else {
				next = rabinChainSureLevel(csg, n, B, var, j - 1);
			}
			done = var[j].equals(next);
			var[j] = next;
		}
		return var[j];
	}

	private BitSet rabinChainSureBody(CSG<Double> csg, BitSet[] B, BitSet[] var) throws PrismException
	{
		int n = csg.getNumStates();
		BitSet result = new BitSet();
		BitSet pre = new BitSet();
		for (int k = 0; k < B.length; k++) {
			if (B[k].isEmpty()) {
				continue;
			}
			pre.clear();
			pre1(csg, var[k], pre);
			BitSet term = (BitSet) pre.clone();
			term.and(B[k]);
			result.or(term);
		}
		return result;
	}

	/**
	 * Checks that the synthesised strategy of the coalition guarantees the computed value (in the initial state):
	 * the optimal value of the other players in the MDP induced by the strategy. Warns if they differ.
	 */
	protected void checkZeroSumStrategy(CSGStrategy strat, CSGRewards<Double> rewards, double computed, boolean min) throws PrismException
	{
		MDPModelChecker mc = new MDPModelChecker(this);
		mc.inheritSettings(this);
		mc.setSilentPrecomputations(true);
		mc.setLog(new PrismDevNullLog());
		double guaranteed = strat.guaranteedValue(mc, rewards, min);
		// for complemented path formulae (e.g. G a, computed as 1 - F !a), report the values of the one asked for
		boolean flip = rewards == null && negatedPath;
		double gShown = flip ? 1.0 - guaranteed : guaranteed, cShown = flip ? 1.0 - computed : computed;
		mainLog.println("\nChecking the synthesised strategy: it guarantees " + gShown + " (computed " + cShown + ")");
		double tol = Math.max(1e-6, 10 * termCritParam);
		double diff = (Double.isInfinite(computed) || Double.isInfinite(guaranteed)) ? (computed == guaranteed ? 0.0 : Double.POSITIVE_INFINITY)
				: Math.abs(guaranteed - computed);
		if (diff > tol * Math.max(1.0, Math.abs(computed))) {
			boolean worse = min ? guaranteed > computed : guaranteed < computed;
			if (worse)
				mainLog.printWarning("The synthesised strategy does not guarantee the computed value (it guarantees " + gShown
						+ "): value iteration may not have converged (try a smaller -epsilon), or optimal strategies may not exist (only near-optimal ones)");
			else
				mainLog.printWarning("The synthesised strategy guarantees more than the computed value (" + gShown
						+ "): value iteration may have stopped before converging (try a smaller -epsilon)");
		}
	}

	// Utility methods for CSG solving
	
	/**
	 * Compute and store information about coalitions (for a zero-sum problem):
	 * - numPlayers and numCoalitions
	 * - coalitionIndexes (players in each coalition)
	 * - actionIndexes (actions for each coalition)
	 * 
	 * This assumes that the first coalition (0) maximises and the second (1) minimises.
	 * 
	 * Also find the max size of the resulting matrix game needed across any CSG state.
	 * 
	 * @param csg The CSG
	 * @param coalition The coalition of players which aims to minimise/maximise
	 * @param min Min or max values for the coalition (true=min, false=max)
	 */
	public void buildCoalitions(CSG<Double> csg, Coalition coalition, boolean min) throws PrismException
	{
		if (coalition == null || coalition.isEmpty()) {
			throw new PrismException("Coalitions must not be empty");
		}
		numPlayers = csg.getNumPlayers();
		numCoalitions = 2;
		coalitionIndexes = new BitSet[2];
		actionIndexes = new BitSet[2];
		Map<Integer, String> pmap = new HashMap<Integer, String>();
		for (int p = 0; p < numPlayers; p++) {
			pmap.put(p + 1, csg.getPlayerName(p));
		}
		for (int c = 0; c < 2; c++) {
			coalitionIndexes[c] = new BitSet();
			actionIndexes[c] = new BitSet();
		}
		for (int p = 0; p < numPlayers; p++) {
			if (min) {
				if (coalition.isPlayerIndexInCoalition(p, pmap)) {
					coalitionIndexes[1].set(p);
					actionIndexes[1].or(csg.getIndexes()[p]);
					actionIndexes[1].set(csg.getIdleForPlayer(p));
				} else {
					coalitionIndexes[0].set(p);
					actionIndexes[0].or(csg.getIndexes()[p]);
					actionIndexes[0].set(csg.getIdleForPlayer(p));
				}
			} else {
				if (coalition.isPlayerIndexInCoalition(p, pmap)) {
					coalitionIndexes[0].set(p);
					actionIndexes[0].or(csg.getIndexes()[p]);
					actionIndexes[0].set(csg.getIdleForPlayer(p));
				} else {
					coalitionIndexes[1].set(p);
					actionIndexes[1].or(csg.getIndexes()[p]);
					actionIndexes[1].set(csg.getIdleForPlayer(p));
				}
			}
		}
		findMaxRowsCols(csg);
	}

	/**
	 * Qualitative (sure/almost/limit) strategy operator for CSGs -- {@code <<C>> sure/almost/
	 * limit [ phi ]}. Restricted for now to safety (G psi), reachability (F psi), Buchi
	 * (G F psi) and co-Buchi (F G psi), for a state (non-path) formula psi (no rPATL/LTL
	 * surface syntax for CSGs beyond this).
	 *
	 * Strictly zero-sum: the coalition C is the achieving (maximising) role for <<C>>
	 * ({@code forAll == false}), or the adversarial (minimising) role for [[C]]
	 * ({@code forAll == true}) -- buildCoalitions is called with min=forAll below,
	 * consistent with buildMatrixDist/buildMatrixGame's convention that rows correspond to
	 * the maximising coalition. This is a role flip only: C's own
	 * membership is untouched, and every downstream fixpoint operator (Apre1/Lpre1,
	 * rabinChainWin/rabinChainWinSure, ...) consumes coalitionIndexes[0]/[1] generically, so
	 * nothing else in this method changes between the two operators. [[C]] Op[phi] is
	 * definitionally "the complement of C plays the achieving role for phi instead" -- not a
	 * negation of phi -- so this needs no determinacy assumption (unlike, say, relating
	 * <<C>> phi to [[C]] via double negation, which would run into concurrent games' known
	 * sure-mode non-determinacy gaps).
	 *
	 * Dispatch by (mode, shape) covers all twelve combinations: G is Pre1/G for all three
	 * modes (safety provably collapses sure=almost=limit, LICS 2000 Theorem 3(1)); F is
	 * SF/AF/LF; F G (co-Buchi) is SFG/AFG/LFG; G F (Buchi) is SGF/AGF/LGF (see the dispatch
	 * block below).
	 */
	@Override
	protected StateValues checkExpressionStrategyQual(Model<?> model, ExpressionStrategyQual expr, boolean forAll, Coalition coalition, BitSet statesOfInterest) throws PrismException
	{
		CSG<Double> csg = (CSG<Double>) model;

		// Strip any outer parentheses
		Expression phi = expr.getExpression();
		while (Expression.isParenth(phi)) {
			phi = ((ExpressionUnaryOp) phi).getOperand();
		}

		ExpressionStrategyQual.QualMode mode = expr.getMode();

		// Classify the shape of phi. Four atomic shapes (G, F, G F, F G) extract a single
		// wrapped state formula psi and fall through to the (mode, shape) dispatch below.
		// A fifth, non-atomic shape is ALSO handled directly here, with no automaton
		// product: a single "G F psi1" (Buchi/recurrence) obligation combined, by
		// conjunction or disjunction, with a single "F G psi2" (co-Buchi/persistence)
		// obligation. This isn't an ad hoc special case -- for ARBITRARY psi1, psi2 (no
		// relationship between their target sets required), both
		//   G F psi1  &  F G psi2   (chain U0=S,U1=S,U2=A cup (not B),U3=not B,U4={}, pair 0
		//                            disabled -- F G B == not G F (not B), so the avoid-set
		//                            here is (not B), not B itself)
		//   G F psi1  |  F G psi2   (chain U0=S,U1=A cup (not B),U2=A,U3={},U4={})
		// are provably realisable by a FIXED per-state colouring -- i.e. memorylessly, on
		// the base CSG state space -- via rabinChainWin/rabinChainWinSure below (this is
		// the classical fact that a single Rabin pair is a [1,3]-parity condition). That
		// is exactly where it stops: two INDEPENDENT G F obligations conjoined
		// (generalised Buchi), or two independent F G obligations disjoined (its dual),
		// are provably NOT memoryless on the base state space -- no per-state colouring
		// realises them, regardless of implementation effort -- and genuinely need an
		// LTL-to-automaton product (constructProductCSG plus a chain-normalisation step
		// on whatever Rabin acceptance comes out, since that output is not naturally
		// chain-shaped either -- see LtlWiringProbe). That machinery is not implemented,
		// and the final else branch below says so explicitly instead of mis-handling or
		// silently under-approximating such formulas.
		String shape;
		Expression psi = null;
		Expression psi1 = null;
		Expression psi2 = null;
		int comboOp = -1;
		if (ExpressionTemporal.isGloballyFinally(phi)) {
			// G F psi (Buchi)
			shape = "G F";
			psi = ((ExpressionTemporal) ((ExpressionTemporal) phi).getOperand2()).getOperand2();
		} else if (ExpressionTemporal.isFinallyGlobally(phi)) {
			// F G psi (co-Buchi)
			shape = "F G";
			psi = ((ExpressionTemporal) ((ExpressionTemporal) phi).getOperand2()).getOperand2();
		} else if (ExpressionTemporal.isGlobally(phi)) {
			// G psi (safety)
			shape = "G";
			psi = ((ExpressionTemporal) phi).getOperand2();
		} else if (ExpressionTemporal.isFinally(phi)) {
			// F psi (reachability)
			shape = "F";
			psi = ((ExpressionTemporal) phi).getOperand2();
		} else if (phi instanceof ExpressionBinaryOp && matchGFandFGCombo((ExpressionBinaryOp) phi) != null) {
			// (G F psi1) & (F G psi2), or (G F psi1) | (F G psi2), in either operand order
			Object[] combo = matchGFandFGCombo((ExpressionBinaryOp) phi);
			comboOp = (Integer) combo[0];
			psi1 = (Expression) combo[1];
			psi2 = (Expression) combo[2];
			shape = (comboOp == ExpressionBinaryOp.AND) ? "G F & F G" : "G F | F G";
		} else {
			throw new PrismNotSupportedException("The " + mode + " operator currently supports G, F, G F, F G path formulas, and "
					+ "(G F ...) combined with (F G ...) by conjunction or disjunction; anything beyond that fragment (e.g. two "
					+ "independent G F obligations conjoined, or two independent F G obligations disjoined) is provably not "
					+ "memoryless on the game's state space and needs an LTL-to-automaton product construction, which is not "
					+ "implemented -- got " + phi);
		}

		// Coalition is the achieving role for <<C>>, the adversarial role for [[C]] -- see
		// method doc above
		buildCoalitions(csg, coalition, forAll);

		// The one non-atomic shape: dispatch straight to rabinChainWin/rabinChainWinSure
		// on the hand-built 2-pair chain described above, and return -- this shape has no
		// place in the (mode, shape) dispatch table below, which is only for the four
		// atomic shapes and their dedicated per-shape operators.
		if (comboOp != -1) {
			if (!(psi1.getType() instanceof TypeBool) || !(psi2.getType() instanceof TypeBool)) {
				throw new PrismException("The " + mode + " operator's path formula must wrap Boolean state formulas, not " + psi1 + ", " + psi2);
			}
			BitSet a = checkExpression(model, psi1, null).getBitSet();
			BitSet b = checkExpression(model, psi2, null).getBitSet();
			BitSet full = new BitSet();
			full.set(0, csg.getNumStates());
			List<BitSet> L;
			List<BitSet> K;
			if (comboOp == ExpressionBinaryOp.AND) {
				// G F a & F G b == G F a & not(G F (not b)) (since F G b == not G F (not b)),
				// so this reuses the "G F x & not G F y" chain with y := not b, NOT y := b:
				// U0=S, U1=S (pair 0 disabled), U2 = a cup (not b), U3 = not b, U4={}
				BitSet notB = (BitSet) b.clone();
				notB.flip(0, csg.getNumStates());
				BitSet aOrNotB = (BitSet) a.clone();
				aOrNotB.or(notB);
				L = Arrays.asList((BitSet) full.clone(), (BitSet) notB.clone());
				K = Arrays.asList((BitSet) full.clone(), aOrNotB);
			} else {
				// G F a | F G b: U0=S, U1=a cup (not b), U2=a, U3={}, U4={}
				BitSet notB = (BitSet) b.clone();
				notB.flip(0, csg.getNumStates());
				BitSet u1 = (BitSet) a.clone();
				u1.or(notB);
				L = Arrays.asList(u1, new BitSet());
				K = Arrays.asList((BitSet) full.clone(), (BitSet) a.clone());
			}
			BitSet comboWin;
			if (mode == ExpressionStrategyQual.QualMode.SURE) {
				comboWin = rabinChainWinSure(csg, L, K);
			} else if (mode == ExpressionStrategyQual.QualMode.ALMOST) {
				comboWin = rabinChainWin(csg, L, K, true);
			} else {
				comboWin = rabinChainWin(csg, L, K, false);
			}
			return StateValues.createFromBitSet(comboWin, csg);
		}

		if (!(psi.getType() instanceof TypeBool)) {
			throw new PrismException("The " + mode + " operator's path formula must wrap a Boolean state formula, not " + psi);
		}

		// Model check psi to get the target state set
		BitSet target = checkExpression(model, psi, null).getBitSet();

		// Dispatch on (mode, shape) to the operator that actually computes it. All twelve
		// (mode, shape) combinations are covered:
		//  - (SURE, G)      -> G/pre1      (safety, Pre1 fixpoint)
		//  - (ALMOST, G)    -> G/pre1      (safety collapses across modes -- see below)
		//  - (LIMIT, G)     -> G/pre1      (safety collapses across modes -- see below)
		//  - (SURE, F)      -> SF/pre1     (sure reachability, mu X.(B or Pre1(X)))
		//  - (SURE, F G)    -> SFG=SF(G)   (sure co-Buchi, reach-the-maximal-invariant)
		//  - (SURE, G F)    -> SGF/pre1    (sure Buchi, classical attractor-of-attractor)
		//  - (ALMOST, F)    -> AF/apreXY   (almost-sure reachability, eq. 2)
		//  - (LIMIT, F)     -> LF/lpreXY   (limit-sure reachability, eq. 2 w/ Lpre1)
		//  - (ALMOST, F G)  -> AFG/apreXYZ (almost-sure co-Buchi, eq. 4, S4.2)
		//  - (LIMIT, F G)   -> LFG/lpreXYZ (limit-sure co-Buchi, eq. 4 w/ Lpre1)
		//  - (ALMOST, G F)  -> AGF/apreXY  (almost-sure Buchi, eq. 3, S4.1)
		//  - (LIMIT, G F)   -> LGF/lpreXY  (limit-sure Buchi, eq. 3 w/ Lpre1)
		//
		// (ALMOST, G) and (LIMIT, G) reuse the same G/pre1 fixpoint as (SURE, G), rather
		// than needing Apre1/Lpre1-based operators of their own, because safety winning
		// regions provably coincide across all three qualitative modes, by
		// Theorem 3(1) of the LICS 2000 paper (S3.1): sure/almost/limit safety all reduce to
		// Pre1 alone, with Spre1/Apre1/Lpre1 collapsing to the same operator, and Table 1(a)
		// gives winning and spoiling strategies the same restricted class across all three
		// modes. Independently: any per-round risk of leaving B, however small, recurs at
		// every future visit to that state (safety is a per-step condition with no "spend
		// the risk budget once" escape the way reachability's LF/AF divergence exploits),
		// so a positive per-round breach probability compounds to certain eventual breach
		// over an infinite horizon -- ruling out any state where Apre1/Lpre1 could do better
		// than Pre1 for G alone.
		BitSet win;
		if (shape.equals("G") && mode == ExpressionStrategyQual.QualMode.SURE) {
			win = G(csg, target);
		} else if (shape.equals("G") && mode == ExpressionStrategyQual.QualMode.ALMOST) {
			win = G(csg, target);
		} else if (shape.equals("G") && mode == ExpressionStrategyQual.QualMode.LIMIT) {
			win = G(csg, target);
		} else if (shape.equals("F") && mode == ExpressionStrategyQual.QualMode.SURE) {
			win = SF(csg, target);
		} else if (shape.equals("F") && mode == ExpressionStrategyQual.QualMode.ALMOST) {
			win = AF(csg, target);
		} else if (shape.equals("F") && mode == ExpressionStrategyQual.QualMode.LIMIT) {
			win = LF(csg, target);
		} else if (shape.equals("F G") && mode == ExpressionStrategyQual.QualMode.SURE) {
			win = SFG(csg, target);
		} else if (shape.equals("F G") && mode == ExpressionStrategyQual.QualMode.ALMOST) {
			win = AFG(csg, target);
		} else if (shape.equals("F G") && mode == ExpressionStrategyQual.QualMode.LIMIT) {
			win = LFG(csg, target);
		} else if (shape.equals("G F") && mode == ExpressionStrategyQual.QualMode.SURE) {
			win = SGF(csg, target);
		} else if (shape.equals("G F") && mode == ExpressionStrategyQual.QualMode.ALMOST) {
			win = AGF(csg, target);
		} else if (shape.equals("G F") && mode == ExpressionStrategyQual.QualMode.LIMIT) {
			win = LGF(csg, target);
		} else {
			throw new PrismNotSupportedException("The " + mode + " [ " + shape + " ... ] operator is not implemented -- "
					+ "no " + mode + "-mode operator exists in this class for the " + shape + " shape");
		}

		return StateValues.createFromBitSet(win, csg);
	}

	/**
	 * If phi is "(G F psi1) op (F G psi2)" or "(F G psi2) op (G F psi1)" for op in
	 * {AND, OR} -- outer parentheses on either operand stripped, either operand order
	 * accepted -- returns {op (as an Integer, one of ExpressionBinaryOp.AND/OR), psi1,
	 * psi2}. Returns null for anything else, including AND/OR of any other shapes (in
	 * particular, this deliberately does NOT match G F/G F or F G/F G combinations --
	 * see the long comment in checkExpressionStrategyQual for why those need a real
	 * automaton product and are not handled here).
	 *
	 * @param phi A binary AND/OR expression (already known to be one, or not -- checked here)
	 */
	private static Object[] matchGFandFGCombo(Expression phi)
	{
		if (!(phi instanceof ExpressionBinaryOp)) {
			return null;
		}
		ExpressionBinaryOp bin = (ExpressionBinaryOp) phi;
		int op = bin.getOperator();
		if (op != ExpressionBinaryOp.AND && op != ExpressionBinaryOp.OR) {
			return null;
		}
		Expression op1 = bin.getOperand1();
		while (Expression.isParenth(op1)) {
			op1 = ((ExpressionUnaryOp) op1).getOperand();
		}
		Expression op2 = bin.getOperand2();
		while (Expression.isParenth(op2)) {
			op2 = ((ExpressionUnaryOp) op2).getOperand();
		}
		Expression gf;
		Expression fg;
		if (ExpressionTemporal.isGloballyFinally(op1) && ExpressionTemporal.isFinallyGlobally(op2)) {
			gf = op1;
			fg = op2;
		} else if (ExpressionTemporal.isGloballyFinally(op2) && ExpressionTemporal.isFinallyGlobally(op1)) {
			gf = op2;
			fg = op1;
		} else {
			return null;
		}
		Expression psi1 = ((ExpressionTemporal) ((ExpressionTemporal) gf).getOperand2()).getOperand2();
		Expression psi2 = ((ExpressionTemporal) ((ExpressionTemporal) fg).getOperand2()).getOperand2();
		return new Object[] { op, psi1, psi2 };
	}

	/**
	 * Find the max size of the matrix game needed across any CSG state,
	 * for the current coalition (as stored in coalitionIndexes).
	 * Also compute/report the average number of actions for each coalition across all CSG states.
	 * 
	 * @param csg The CSG
	 */
	public void findMaxRowsCols(CSG<Double> csg)
	{
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
			maxRows = (maxRows < mr) ? mr : maxRows;
			maxCols = (maxCols < mc) ? mc : maxCols;
		}
		avgNumActions[0] /= csg.getNumStates();
		avgNumActions[1] /= csg.getNumStates();
		mainLog.println("Max/avg (actions): " + "(" + maxRows + "," + maxCols + ")/(" + PrismUtils.formatDouble2dp(avgNumActions[0]) + ","
				+ PrismUtils.formatDouble2dp(avgNumActions[1]) + ")");
	}

	/**
	 * Build a matrix game with the transition distributions for state s. 
	 * 
	 * @param csg The CSG
	 * @param s Index of state to build matrix game for 
	 * @return
	 * @throws PrismException
	 */
	public ArrayList<ArrayList<Distribution<Double>>> buildMatrixDist(CSG<Double> csg, int s) throws PrismException
	{
		ArrayList<ArrayList<Distribution<Double>>> mdist = new ArrayList<>();
		BitSet action = new BitSet();
		int col, row;
		buildStepGame(csg, null, null, null, s);
		for (row = 0; row < strategies.get(0).size(); row++) {
			mdist.add(row, new ArrayList<>());
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

	/**
	 * Build the matrix game to solve a CSG state s.
	 * This is returned as a list of rows, where each row is a list of values.
	 * Rows correspond to the maximising coalition, columns to the minimising one. 
	 * <br><br>
	 * If argument {@code mmap} is non-null, it is filled with a list of coalition actions
	 * for the coalition for which a strategy is being synthesised.
	 * The list is stored as a map from (ascending integer) indices to
	 * BitSets containing the indices of the actions in the coalition action.
	 * 
	 * @param csg The CSG
	 * @param r The rewards
	 * @param mmap Map in which to store coalition action indices
	 * @param val Array (over states) of values to multiply by when computing matrix values
	 * @param s Index of state to build matrix game for 
	 * @param min Min or max values for the coalition (true=min, false=max)
	 */
	public ArrayList<ArrayList<Double>> buildMatrixGame(CSG<Double> csg, CSGRewards<Double> r, Map<Integer, BitSet> mmap, double[] val, int s, boolean min)
			throws PrismException
	{
		ArrayList<ArrayList<Double>> mgame = new ArrayList<ArrayList<Double>>();
		ArrayList<CSGRewards<Double>> rewards = null;
		Map<BitSet, Integer> imap = new HashMap<BitSet, Integer>();
		Map<Integer, BitSet> rmap;
		BitSet action = new BitSet();
		int col, row;
		// Build matrix game info
		if (r != null) {
			rewards = new ArrayList<>();
			rewards.add(0, r);
		}
		buildStepGame(csg, rewards, imap, val, s);
		// Reverse map for coalition actions (i.e., store mapping from their integer indices,
		// to the coalition actions, represented as BitSets storing the indices of their actions.
		// If requested (for strategy synthesis), store a copy of the coalition actions for one coalition.
		rmap = imap.entrySet().stream().collect(Collectors.toMap(Map.Entry::getValue, Map.Entry::getKey));
		if (mmap != null) {
			if (min) {
				for (col = 0; col < strategies.get(1).size(); col++) {
					mmap.put(col, rmap.get(strategies.get(1).get(col)));
				}
			} else {
				for (row = 0; row < strategies.get(0).size(); row++) {
					mmap.put(row, rmap.get(strategies.get(0).get(row)));
				}
			}
		}
		// For each coalition action of the maximising coalition
		for (row = 0; row < strategies.get(0).size(); row++) {
			mgame.add(row, new ArrayList<Double>());
			action.clear();
			action.set(strategies.get(0).get(row));
			// For each coalition action of the minimising coalition
			for (col = 0; col < strategies.get(1).size(); col++) {
				action.set(strategies.get(1).get(col));
				// Find corresponding matrix value, store 
				if (utilities.containsKey(action)) {
					mgame.get(row).add(col, utilities.get(action).get(0));
				} else
					throw new PrismException("Error in building matrix game");
				action.clear(strategies.get(1).get(col));
			}
		}
		return mgame;
	}

	/**
	 * Build info needed for the matrix game to solve a CSG state s.
	 * A list of all coalition actions (comprising one action, incl. "idle",
	 * for each player in the coalition), used by either coalition in this state,
	 * is computed and stored in {@code imap}, as a mapping from BitSets (of
	 * the indices of actions in the coalition action) to their index in this list. 
	 * Matrix game info is stored in fields:
	 * <ul>
	 * <li> actions (for each coalition, a description of each coalition action in s)
	 * <li> strategies (for each coalition, the indices of all coalition actions in s) 
	 * <li> utilities (matrix game values, stored as mapping from pairs of coalition action
	 *     indices, as a BitSet, to the value, stored in a 1-element ArrayList)
	 * <li> probabilities (CSG distribution used to compute value, again as mapping from BitSet)
	 * <li> varIndex (total number of coalition actions in s)
	 * <li> minEntry (minimum matrix value) 
	 * <li> allEqual (true if all matrix values are equal)
	 * </ul>
	 * @param csg The CSG
	 * @param rewards The rewards (contained in a 1-element list) 
	 * @param imap Map in which to store coalition actions and indices
	 * @param val Array (over states) of values to multiply by when computing matrix values
	 * @param s Index of state to build matrix game for 
	 */
	public void buildStepGame(CSG<Double> csg, List<CSGRewards<Double>> rewards, Map<BitSet, Integer> imap, double[] val, int s) throws PrismException
	{
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
		// For each choice in state s
		for (t = 0; t < csg.getNumChoices(s); t++) {
			jidx = new BitSet();
			joint = csg.getIndexes(s, t);
			// Build bitset of indices of actions (incl. idle) for all players in choice t
			indexes.clear();
			for (p = 0; p < numPlayers; p++) {
				if (joint[p] != -1)
					indexes.set(joint[p]);
				else
					indexes.set(csg.getIdleForPlayer(p));
			}
			// For each coalition: update strategies/action lists.
			// Also build coalition action index pair as bitset to store value 
			for (c = 0; c < 2; c++) {
				v = 0.0;
				// Build bitset of indices of actions (incl. idle) used by coalition c in choice t
				tmp.clear();
				tmp.or(actionIndexes[c]);
				tmp.and(indexes);
				if (tmp.cardinality() != coalitionIndexes[c].cardinality()) {
					// Should be one per player
					throw new PrismException("Error in coalition");
				} else {
					if (!imap.keySet().contains(tmp)) {
						act = "";
						strategies.get(c).add(varIndex);
						for (i = tmp.nextSetBit(0); i >= 0; i = tmp.nextSetBit(i + 1)) {
							act += "[" + csg.getActions().get(i - 1) + "]";
						}
						actions.get(c).add(act);
						jidx.set(varIndex);
						imap.put((BitSet) tmp.clone(), varIndex);
						varIndex++;
					} else {
						jidx.set(imap.get(tmp));
					}
				}
			}
			// Compute matrix value and store
			// (also store distribution and update allEqual/minEntry)
			utilities.put(jidx, new ArrayList<Double>());
			probabilities.put(jidx, new ArrayList<Distribution<Double>>());
			v = 0.0;
			if (val != null) {
				for (Iterator<Map.Entry<Integer, Double>> iter = csg.getTransitionsIterator(s, t); iter.hasNext();) {
					Map.Entry<Integer, Double> e = iter.next();
					v += e.getValue() * val[e.getKey()];
				}
			}
			if (rewards != null)
				v += rewards.get(0).getTransitionReward(s, t);
			if (u != Double.NaN)
				allEqual = allEqual && Double.compare(u, v) == 0;
			utilities.get(jidx).add(0, v);
			probabilities.get(jidx).add(0, csg.getChoice(s, t));
			minEntry = (minEntry > v) ? v : minEntry;
			u = v;
		}
		//System.out.println("## state " + csg.getStatesList().get(s));
		//System.out.println("-- actions " + actions);
		//System.out.println("-- strategies " + strategies);		
		//System.out.println("-- utilities " + utilities);
		//System.out.println("-- imap " + imap);
	}

	/**
	 * Solves a matrix game and return its value.
	 * If requested, store an optimal strategy for the coalition being solved for.   
	 * 
	 * @param lp LpSolve instance to use for solving
	 * @param mgame The matrix
	 * @param strat Storage for strategy (as map from coalition actions to probability of selection)
	 * @param rmap List of coalition actions for the coalition to solve for
	 *             (stored as a map from (ascending integer) indices to
	 *             BitSets containing the indices of the actions in the coalition action)
	 * @param s Index of state matrix game is for 
	 * @param rew Are we solving a reward (true) or probability (false) problem?
	 * @param min Are we minimising or maximising? (dictates which coalition to solve for) 
	 */
	public double val(LpSolve lp, ArrayList<ArrayList<Double>> mgame, List<Map<BitSet, Double>> strat, Map<Integer, BitSet> rmap, int s, boolean rew,
			boolean min) throws PrismException
	{
		long timer = System.currentTimeMillis();
		int nrows = mgame.size(); // Number of rows
		int ncols = mgame.get(0).size(); // Number of columns
		double res = Double.NaN;
		Map<BitSet, Double> d = new HashMap<BitSet, Double>();
		// Special cases
		if (allEqual) {
			if (genStrat) {
				d.put(rmap.get(0), 1.0);
				strat.set(s, d);
			}
			return minEntry;
		} else if (nrows == 1) {
			int srow = 0;
			res = Double.POSITIVE_INFINITY;
			for (int col = 0; col < ncols; col++) {
				if (res > mgame.get(0).get(col)) {
					res = mgame.get(0).get(col);
					srow = (min) ? col : 0;
				}
			}
			if (genStrat) {
				d.put(rmap.get(srow), 1.0); // In case of min, rmap maps columns not rows
				strat.set(s, d);
			}
			return res;
		} else if (ncols == 1) {
			int scol = 0;
			res = Double.NEGATIVE_INFINITY;
			for (int row = 0; row < nrows; row++) {
				if (res < mgame.get(row).get(0)) {
					res = mgame.get(row).get(0);
					scol = (min) ? 0 : row;
				}
			}
			if (genStrat) {
				d.put(rmap.get(scol), 1.0);
				strat.set(s, d);
			}
			return res;
		} else {
			// Should add check for trivial games
			int infty;
			infty = valInfinity(mgame);
			if (infty != -1) {
				res = Double.POSITIVE_INFINITY;
				if (genStrat) {
					d.put(rmap.get(infty), 1.0);
				}
				return res;
			} else if ("SoPlex".equals(lpSolver)) {
				res = valSoPlex(mgame, d, rmap, s, rew, min);
				if (genStrat)
					strat.set(s, d);
			} else {
				try {
					if (min)
						lp.resizeLp(0, maxCols + 1);
					else
						lp.resizeLp(0, maxRows + 1);
					buildLPLpsolve(lp, mgame, rew, min);
				} catch (LpSolveException e1) {
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
						double[] values = (min) ? new double[maxCols + 1] : new double[maxRows + 1];
						lp.getVariables(values);
						if (genStrat) {
							for (int row = 1; row <= ((min) ? ncols : nrows); row++) { // C's actions: columns when min
								if (values[row] > 0)
									d.put(rmap.get(row - 1), values[row]);
							}
							strat.set(s, d);
						}
					} else {
						throw new PrismException("lpSolve could not find an optimal solution for state " + s);
					}
				} catch (Exception e) {
					mainLog.println(
							"Exception raised by lpSolve when computing value for state " + s + ". lpSolve status: " + lp.getStatustext(lp.getStatus()));
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
							double[] values = (min) ? new double[maxCols + 1] : new double[maxRows + 1];
							lp.getVariables(values);
							if (genStrat) {
								for (int row = 1; row <= ((min) ? ncols : nrows); row++) { // C's actions: columns when min
									if (values[row] > 0)
										d.put(rmap.get(row - 1), values[row]);
								}
								strat.set(s, d);
							}
						} else {
							throw new PrismException("Rounding up failed for state " + s + ". Failed to compute solution");
						}
					} catch (LpSolveException e2) {
						throw new PrismException("Rounding up failed for state " + s + ". Failed to compute solution");
					}
					return res;
				}
			}
		}
		timer = System.currentTimeMillis() - timer;
		timerVal += timer;
		return res;
	}

	/**
	 * Solve a matrix game with SoPlex (reusing one SoPlex object across states).
	 * Same LP as {@link #buildLPLpsolve}: the value is bounded to [0,1] for probabilities, free for rewards;
	 * if {@code min} the LP is solved for the column (minimising) coalition, otherwise for the row coalition.
	 * The optimal strategy of that coalition is stored in {@code d} (keys from {@code rmap}).
	 */
	public double valSoPlex(ArrayList<ArrayList<Double>> mgame, Map<BitSet, Double> d, Map<Integer, BitSet> rmap, int s, boolean rew, boolean min) throws PrismException
	{
		int nrows = mgame.size();
		int ncols = mgame.get(0).size();
		double[] A = new double[nrows * ncols];
		for (int i = 0; i < nrows; i++)
			for (int j = 0; j < ncols; j++)
				A[i * ncols + j] = mgame.get(i).get(j);
		if (soplex == null) {
			try {
				soplex = new SoPlex();
				soplex.setScaler(SoPlex.scalerFromName(soplexScaling));
			} catch (UnsatisfiedLinkError | IllegalStateException e) {
				throw new PrismException("Could not load SoPlex (native library soplexj): " + e.getMessage());
			}
		}
		double[] strat = new double[min ? ncols : nrows];
		// Normalise the game to entries in [0,1] (an affine change of the payoffs, which leaves the optimal
		// strategies unchanged): games whose entries differ only slightly (e.g. all close to 1, as in the late
		// iterations of value iteration) are otherwise badly conditioned, and SoPlex can fail on them
		double lo = Double.POSITIVE_INFINITY, hi = Double.NEGATIVE_INFINITY;
		for (double a : A) {
			lo = Math.min(lo, a);
			hi = Math.max(hi, a);
		}
		double range = hi - lo;
		if (range == 0.0) {
			// constant game: any strategy is optimal
			strat[0] = 1.0;
			if (genStrat)
				d.put(rmap.get(0), 1.0);
			return lo;
		}
		for (int i = 0; i < A.length; i++)
			A[i] = (A[i] - lo) / range;
		double vlb = rew ? -SoPlex.INFINITY : (0.0 - lo) / range;
		double vub = rew ? SoPlex.INFINITY : (1.0 - lo) / range;
		double res = soplex.matrixGame(A, nrows, ncols, min, vlb, vub, strat);
		if (Double.isNaN(res))
			throw new PrismException("SoPlex could not find an optimal solution for state " + s + " (status " + soplex.getStatus() + ")");
		res = lo + range * res;
		if (genStrat) {
			for (int k = 0; k < strat.length; k++) {
				if (strat[k] > 0)
					d.put(rmap.get(k), strat[k]);
			}
		}
		return res;
	}

	/**
	 * Deals with infinite cases in solving a matrix game.
	 * If all values in some row are +inf, return the index of that row (it's optimal).
	 * Then remove any column containing a +inf and return -1. 
	 */
	public int valInfinity(ArrayList<ArrayList<Double>> mgame)
	{
		BitSet hasInf = new BitSet();
		int row, col;
		boolean infRow;
		// Check which columns have a +inf value,
		// and whether any row has all values equal to +inf
		for (row = 0; row < mgame.size(); row++) {
			infRow = true;
			for (col = 0; col < mgame.get(0).size(); col++) {
				if (mgame.get(row).get(col) == Double.POSITIVE_INFINITY) {
					hasInf.set(col);
				} else {
					infRow = false;
				}
			}
			if (infRow)
				return row;
		}
		// If any column contains +inf remove it
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

	/**
	 * Build the linear program to solve a matrix game.
	 * 
	 * @param lp LpSolve instance to use for constructing the LP
	 * @param mgame The matrix
	 * @param rew Are we solving a reward (true) or probability (false) problem?
	 * @param min Are we minimising or maximising? (dictates which coalition to solve for) 
	 */
	public void buildLPLpsolve(LpSolve lp, ArrayList<ArrayList<Double>> mgame, boolean rew, boolean min) throws LpSolveException
	{
		int nrows = (min) ? mgame.get(0).size() : mgame.size(); // Number of rows
		int ncols = (min) ? mgame.size() : mgame.get(0).size(); // Number of columns
		int[] vari = new int[nrows + 1]; // Indexes of variables, should be m + 1 for an m x n matrix
		double[] row = new double[nrows + 1];
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
		} else {
			lp.setBounds(1, 0, 1.0);
		}
		lp.setAddRowmode(true);
		int k = 0;
		// Builds each column considering the n + 1 constraints as a matrix
		for (int j = 0; j < ncols; j++) {
			vari[k] = 1;
			row[k] = scaleFactor;
			for (int i = 0; i < nrows; i++) {
				k++;
				vari[k] = k + 1;
				if (min)
					row[k] = -1.0 * scaleFactor * mgame.get(j).get(i);
				else
					row[k] = -1.0 * scaleFactor * mgame.get(i).get(j);
			}
			if (min)
				lp.addConstraintex(nrows + 1, row, vari, LpSolve.GE, 0.0);
			else
				lp.addConstraintex(nrows + 1, row, vari, LpSolve.LE, 0.0);
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
			row[k] = (k > 0) ? 1 : 0;
		}
		lp.addConstraintex(k, row, vari, LpSolve.EQ, 1.0);
		lp.setAddRowmode(false);
		for (k = 0; k < nrows + 1; k++) {
			vari[k] = k + 1;
			row[k] = (k == 0) ? 1 : 0;
		}
		lp.setObjFnex(k, row, vari);
		if (min)
			lp.setMinim();
		else
			lp.setMaxim();
		//lp.printLp();
	}

	/**
	 * Update the strategy
	 *
	 * @param kstrat The strategies for all states computed in iteration k
	 * @param lstrat The overall strategy
	 * @param k Iteration k, if bounded
	 * @param s The index of the state
	 * @param bounded
	 */
	public void updateStrategy(List<Map<BitSet, Double>> kstrat, List<List<List<Map<BitSet, Double>>>> lstrat, int k, int s, boolean bounded)
	{
		// player -> iteration -> state -> indexes -> value
		if (bounded) {
			if (lstrat.get(0).get(k).get(s) == null || !lstrat.get(0).get(k - 1).get(s).equals(kstrat.get(s))) {
				lstrat.get(0).get(k).set(s, kstrat.get(s));
			} else {
				lstrat.get(0).get(k).set(s, lstrat.get(0).get(k - 1).get(s));
			}
		} else {
			if (lstrat.get(0).get(0).get(s) == null) {
				lstrat.get(0).get(0).set(s, kstrat.get(s));
			} else if (!lstrat.get(0).get(0).get(s).equals(kstrat.get(s))) {
				lstrat.get(0).get(0).set(s, kstrat.get(s));
			}
		}
	}

	/**
	 * Constructs the DRA for the CSGxDRA product when computing mixed bounded equilibria.
	 * @param labelA Label for accepting states
	 * @param bounds Bound for either cumulative or instant rewards
	 * @return
	 */
	public static DA<BitSet, AcceptanceRabin> constructDRAForInstant(String labelA, IntegerBound bounds)
	{
		DA<BitSet, AcceptanceRabin> dra;
		List<String> apList = new ArrayList<String>();
		BitSet edge_no, edge_yes;
		BitSet accL, accK;

		int saturation = bounds.getMaximalInterestingValue();

		// [0,saturation] + yes
		int states = saturation + 2;

		apList.add(labelA);

		dra = new DA<BitSet, AcceptanceRabin>(states);
		dra.setAcceptance(new AcceptanceRabin());
		dra.setAPList(apList);
		dra.setStartState(0);

		// edge labeled with the target label
		edge_yes = new BitSet();
		// edge not labeled with the target label
		edge_no = new BitSet();

		edge_yes.set(0); // yes = a, no = !a

		int yes_state = states - 1;
		int next_counter;

		for (int counter = 0; counter <= saturation; counter++) {
			next_counter = counter + 1;
			if (next_counter > saturation)
				next_counter = saturation;

			if (bounds.isInBounds(counter)) {
				dra.addEdge(counter, edge_no, next_counter);
				if (counter <= saturation - 2)
					dra.addEdge(counter, edge_yes, next_counter);
				else {
					dra.addEdge(counter, edge_yes, yes_state);
				}
			} else {
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


}
