//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford)
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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.Stack;
import java.util.Arrays;
import java.math.BigInteger;
import java.lang.Math;
import java.util.Map;
import java.util.HashMap;

import parser.ast.Expression;
import parser.ast.ExpressionPATL;
import parser.ast.ExpressionProb;
import parser.ast.ExpressionTemporal;
import parser.ast.ExpressionUnaryOp;
import parser.ast.ExpressionBinaryOp;
import parser.visitor.ASTTraverse;
import prism.PrismException;
import prism.PrismLangException;
import strat.ExactValueStrategy;
import strat.MultiObjectiveStrategy;
import strat.Strategy;
import explicit.rewards.SMGRewards;

import parser.ast.RewardStruct;
import parser.State;
import explicit.rewards.ConstructRewards;
import explicit.rewards.STPGRewards;
import explicit.PPLSupport;
import explicit.MapMDPSimulator;

import org.apache.commons.math3.fraction.BigFraction;
import parma_polyhedra_library.*;

/**
 * Explicit-state model checker for multi-player stochastic games (SMGs).
 */
public class SMGModelChecker extends STPGModelChecker
{


    protected Map<Integer,Polyhedron> checkMultiObjectiveFormula(Model model, ExpressionPATL exprPATL, boolean min, List<List<Polyhedron>> stochasticStates) throws PrismException
    {
	// dynamically load the Parma Polyhedra Library
	System.loadLibrary("ppl_java");
	// initializa ppl
	Parma_Polyhedra_Library.initialize_library();

	List<STPGRewards> stpgRewards = new ArrayList<STPGRewards>();


	// print model
	int initial_state = model.getFirstInitialState();
	System.out.printf("initial state: %d\n", initial_state);
	System.out.println(((STPG) model));
	// get who is playing against whom
	((SMG) model).setCoalition(exprPATL.getCoalition());
	Expression expr = exprPATL.getExpressionProb().getExpression();


	Expression prob = exprPATL.getExpressionProb().getProb(); // hack to get some bound on maxIter
	double maxIter = 1.0/prob.evaluateDouble(constantValues);

	// Test whether this is a simple path formula (i.e. PCTL)
	// and then pass control to appropriate method.
	if (expr.isSimplePathFormula()) {
	    // #clemens: major hack to get more than two goals
	    //           need to adapt parser anyway
	    List<BitSet> targets = new ArrayList<BitSet>();
	    Expression temp_expr = expr;
	    try {
		while(temp_expr.isSimplePathFormula()){
		    BitSet t = checkExpression(model, ((ExpressionTemporal) temp_expr).getOperand1()).getBitSet();
		    temp_expr = ((ExpressionTemporal) temp_expr).getOperand2();
		    if(temp_expr instanceof ExpressionProb){
			temp_expr = ((ExpressionProb) temp_expr).getExpression();
		    }
		    targets.add(t);
	    }
		BitSet t = checkExpression(model, temp_expr).getBitSet();
		targets.add(t);
	    } catch (Exception e) {
		throw new PrismException("An error occurred trying to extract the goals. Only the following form is supported: P [ goal1 U P [ goal2 U P [ goal3 ... ] ] ]");
	    }
	    
	    // rewards
	    // TODO: properly get reward structure
	    RewardStruct rewStruct;
	    ConstructRewards constructRewards;
	    for(int i = 0; i < modulesFile.getNumRewardStructs(); i++){

		rewStruct = modulesFile.getRewardStruct(i);
		constructRewards = new ConstructRewards(mainLog);
		STPGRewards stpgr = constructRewards.buildSTPGRewardStructure((STPG) model, rewStruct, constantValues);
		stpgRewards.add(stpgr);
	    }

	    /*
	    List<Set<ReachTuple>> result_t = computeReachabilityTuples(min, !min, (STPG) model, targets, maxIter);
	    for(int i = 0; i < result_t.size(); i++){
		System.out.printf("%d: %s\n", i, result_t.get(i));
	    }
	    */

	    // set accuracy


	    long[] accuracy = new long[targets.size()+stpgRewards.size()];
	    long baseline_accuracy = 40;
	    System.err.printf("Accuracy: %d", baseline_accuracy);
	    for(int i = 0; i < targets.size()+stpgRewards.size(); i++) {
		if(i < targets.size()) { // probabilities
		    accuracy[i] = baseline_accuracy;
		} else { // rewards
		    long maxReward = 1;
		    accuracy[i] = baseline_accuracy/maxReward;
		}
	    }

	    //maxIter = stoppingCriterion((STPG) model, 1, stpgRewards, accuracy);
	    maxIter = Double.MAX_VALUE;
	    System.out.printf("maxIter: %e\n", maxIter);

	    long polyTime = System.nanoTime();
	    Map<Integer,Polyhedron> result_p = null;
	    String compute = System.getenv().get("PCOMP");
	    if(compute.equals("compute")) {
		result_p = this.computeReachabilityPolyhedra(min, !min, (STPG) model, stpgRewards, targets, accuracy, maxIter, stochasticStates);
		polyTime = System.nanoTime() - polyTime;
		
		System.out.printf("Polyhedra computation: %.4f ms\n", ((double)polyTime)/1000000.0);

		double[] goal = { 0.6, 0.7, 10.0 };
		
		MultiObjectiveStrategy strategy_mdp = new MultiObjectiveStrategy((STPG) model, initial_state, goal, result_p, stochasticStates, stpgRewards);
		
		MapMDPSimulator mmdps = new MapMDPSimulator((STPG) model, stpgRewards);
		
		mmdps.writeStrategy(strategy_mdp, "mmdps");
	    } else if (compute.equals("recompute")) {
		MapMDPSimulator mmdps = new MapMDPSimulator((STPG) model, stpgRewards);
		mmdps.readStrategy("mmdps");

	        //double[] goal = { 0.019, 0.903, 2.884 };
                //double[] goal = { 0.769, 0.75,  11.69 };
	        double[] goal = { 0.0, 0.576, 17.21 };

		mmdps.recomputeInitial(goal);
		mmdps.writeStrategy("mmdps");

	    } else if (compute.equals("simulate")) {
		MapMDPSimulator mmdps = new MapMDPSimulator((STPG) model, stpgRewards);
		mmdps.readStrategy("mmdps");
		mmdps.evaluateStrategy();
	    } else {
		System.out.println("Invalid command in environment variable PCOMP.");
	    }
	    

	    
	    return result_p;
	    
	}

	throw new PrismException("Explicit engine does not yet handle LTL-style path formulas.");
	
    }




    protected Strategy constructMultiStrategy(Model G, List<Polyhedron> X)
    {
	return null;
    }


	/**
	 * Compute probabilities for the contents of a P operator.
	 */
	protected StateValues checkProbPathFormula(Model model, ExpressionPATL exprPATL, boolean min) throws PrismException
	{

	    //@clemens : don't change this - this works out the player coalition
	    // setting coalition parameter
	    ((SMG) model).setCoalition(exprPATL.getCoalition());
	    
	    Expression expr = exprPATL.getExpressionProb().getExpression();

		// Test whether this is a simple path formula (i.e. PCTL)
		// and then pass control to appropriate method.
		if (expr.isSimplePathFormula()) {
			double p = -1;
			Expression pb = exprPATL.getExpressionProb().getProb();
			if (pb != null) {
			    p = pb.evaluateDouble(constantValues);
			}
			// do the polyhedra computation
			List<List<Polyhedron>> Y = new ArrayList<List<Polyhedron>>(((STPG) model).getStatesList().size());
			Map<Integer,Polyhedron> X = checkMultiObjectiveFormula(model, exprPATL, min, Y);
			

			// here do the standard method that I've basically overridden
			return super.checkProbPathFormulaSimple(model, expr, min, !min, p);
		}

		 /*
		  * TODO implement FG and GF formulae //Test if this is FG if (expr
		  * instanceof ExpressionTemporal) { ExpressionTemporal exprT =
		  * (ExpressionTemporal) expr; if (exprT.getOperator() ==
		  * ExpressionTemporal.P_F) { Expression expr2 = exprT.getOperand2(); if
		  * (expr2 instanceof ExpressionTemporal) { ExpressionTemporal expr2T =
		  * (ExpressionTemporal) expr2; if (expr2T.getOperator() ==
		  * ExpressionTemporal.P_G) { Expression expr3 = expr2T.getOperand2(); if
		  * (!(expr3 instanceof ExpressionTemporal)) { return
		  * super.checkFG(model, expr, min, !min); } } } } }
		  * 
		  * //Test whether this is GF if (expr instanceof ExpressionTemporal) {
		  * ExpressionTemporal exprT = (ExpressionTemporal) expr; if
		  * (exprT.getOperator() == ExpressionTemporal.P_G) { Expression expr2 =
		  * exprT.getOperand2(); if (expr2 instanceof ExpressionTemporal) {
		  * ExpressionTemporal expr2T = (ExpressionTemporal) expr2; if
		  * (expr2T.getOperator() == ExpressionTemporal.P_F) { Expression expr3 =
		  * expr2T.getOperand2(); if (!(expr3 instanceof ExpressionTemporal)) {
		  * return super.checkGF(model, expr, min, !min); } } } } }
		  */

		 // in other case
		 throw new PrismException("Explicit engine does not yet handle LTL-style path formulas except for GF and FG");


	 }

	 /**
	  * Compute rewards for the contents of an R operator.
	  */
	 protected StateValues checkRewardFormula(Model model, SMGRewards modelRewards, ExpressionPATL exprPATL, boolean min)
			 throws PrismException
	 {
		 // setting coalition parameter
		 ((SMG) model).setCoalition(exprPATL.getCoalition());

		 StateValues rewards = null;
		 Expression expr = exprPATL.getExpressionRew().getExpression();

		 if (expr instanceof ExpressionTemporal) {
			 ExpressionTemporal exprTemp = (ExpressionTemporal) expr;
			 switch (exprTemp.getOperator()) {
			 case ExpressionTemporal.R_F:
				 rewards = checkRewardReach(model, modelRewards, exprTemp, min, !min, STPGModelChecker.R_INFINITY);
				 break;
			 case ExpressionTemporal.R_Fc:
				 rewards = checkRewardReach(model, modelRewards, exprTemp, min, !min, STPGModelChecker.R_CUMULATIVE);
				 break;
			 case ExpressionTemporal.R_F0:
				 rewards = checkRewardReach(model, modelRewards, exprTemp, min, !min, STPGModelChecker.R_ZERO);
				 break;
			 default:
				 throw new PrismException("Explicit engine does not yet handle the " + exprTemp.getOperatorSymbol()
						 + " operator in the R operator");
			 }
		 }

		 if (rewards == null)
			 throw new PrismException("Unrecognised operator in R operator");

		 return rewards;
	 }

	 protected StateValues checkExactProbabilityFormula(Model model, ExpressionPATL expr, double p)
			 throws PrismException
	 {
		 if(expr.getExpressionProb().getExpression() instanceof ExpressionTemporal && ((ExpressionTemporal)expr.getExpressionProb().getExpression()).hasBounds())
		 {
			 throw new PrismException(
			 "The exact probability queries are not supported for step-bounded properties");
		 }

		 ((SMG) model).setCoalition(expr.getCoalition());
		 // 1) check whether the game is stopping, if not - terminate
		 // 1.1) find states which have self loops only
		 BitSet terminal = findTerminalStates(model);
		 // 1.2) check whether the minmin prob to reach those states is
		 // 1, if not - terminate, if yes continue to 2)
		 double[] res = ((SMGModelChecker) this).computeUntilProbs((STPG) model, null, terminal, true, true, 1.0).soln;

		 // System.out.println("Terminal states: " + terminal);
		 // System.out.println(Arrays.toString(res));
		 for (int i = 0; i < res.length; i++)
			 if (res[i] < 1.0 - 1e-6)
				 throw new PrismException(
						 "The game is not stopping. The exact probability queries only work for stopping games");

		 // 2) computing minmax and maxmin values for all states
		 double[] minmax = null, maxmin = null; // see the do loop below

		 // 3) removing states from the game which have minmax>maxmin
		 // model.
		 int n = model.getNumStates();
		 boolean repeat;
		 BitSet removed = new BitSet(n), removedNew = new BitSet(n);
		 STPG stpg = ((STPG) model);

		 Strategy minStrat = null;
		 Strategy maxStrat = null;

		 do {
			 // computing minmax and maxmin
			 minmax = this.checkProbPathFormula(model, expr, true).getDoubleArray();
			 if (generateStrategy)
				 minStrat = strategy;
			 maxmin = this.checkProbPathFormula(model, expr, false).getDoubleArray();
			 if (generateStrategy)
				 maxStrat = strategy;

			 repeat = false;
			 // checking which states are marked for removal
			 for (int i = 0; i < n; i++)
				 if (!removed.get(i) && minmax[i] > maxmin[i]) {
					 removed.set(i);
					 removedNew.set(i);
					 repeat = true;
				 }

			 // disabling choices that have transitions to those states
			 removedNew.flip(0, n);
			 for (int i = 0; i < n; i++)
				 for (int j = 0; j < model.getNumChoices(i); j++)
					 if (!stpg.allSuccessorsInSet(i, j, removedNew))
						 stpg.disableChoice(i, j);
			 removedNew.clear();
			 // 4) repeat 2-3 while the set of states from 3 is empty
		 } while (repeat);

		 // 5) if bound is null, return the interval, otherwise check
		 // whether bound is in the interval.
		 BitSet ret = new BitSet(n);
		 for (int i = 0; i < n; i++)
			 ret.set(i, !removed.get(i) && minmax[i] <= p && maxmin[i] >= p);

		 // enabling choices that have been disabled for model checking
		 stpg.enableAllChoices();

		 if (generateStrategy) {
			 strategy = new ExactValueStrategy(minStrat, minmax, maxStrat, maxmin, p, (STPG) model);
		 }

		 return StateValues.createFromBitSet(ret, model);
	 }

	 protected StateValues checkExactRewardFormula(Model model, SMGRewards modelRewards, ExpressionPATL expr, double p)
			 throws PrismException
	 {
		 ((SMG) model).setCoalition(expr.getCoalition());
		 // check if the reward is Fc
		 ExpressionTemporal exprTemp = null;
		 if (expr.getExpressionRew().getExpression() instanceof ExpressionTemporal) {
			 exprTemp = (ExpressionTemporal) expr.getExpressionRew().getExpression();
			 switch (exprTemp.getOperator()) {
			 case ExpressionTemporal.R_Fc:
				 break;
			 case ExpressionTemporal.R_F:
				 throw new PrismException("Only cumulative reward type is supported for exact values.");
			 case ExpressionTemporal.R_F0:
				 throw new PrismException("Only cumulative reward type is supported for exact values.");
			 default:
				 throw new PrismException("Only cumulative reward type is supported for exact values.");
			 }
		 } else {
			 throw new PrismException("Only temporal expression are supported at the moment");
		 }

		 // 1) check whether the game is stopping, if not - terminate
		 // 1.1) find states which have self loops only
		 BitSet terminal = findTerminalStates(model);
		 // 1.2) check whether the minmin prob to reach those states is
		 // 1, if not - terminate, if yes continue to 2)
		 double[] res = ((SMGModelChecker) this).computeUntilProbs((STPG) model, null, terminal, true, true, 1.0).soln;

		 // System.out.println("Terminal states: " + terminal);
		 // System.out.println(Arrays.toString(res));
		 for (int i = 0; i < res.length; i++)
			 if (res[i] < 1.0 - 1e-6)
				 throw new PrismException(
						 "The game is not stopping. The exact probability queries only work for stopping games");

		 // 2) computing minmax and maxmin values for all states
		 double[] minmax = null, maxmin = null; // see the do loop below 

		 // 3) removing states from the game which have minmax>maxmin
		 // model.
		 int n = model.getNumStates();
		 boolean repeat;
		 BitSet removed = new BitSet(n), removedNew = new BitSet(n);
		 STPG stpg = ((STPG) model);

		 Strategy minStrat = null;
		 Strategy maxStrat = null;

		 do {
			 // computing minmax and maxmin
			 minmax = this.checkRewardReach(model, modelRewards, exprTemp, true, false, STPGModelChecker.R_CUMULATIVE).valuesD;
			 if (generateStrategy)
				 minStrat = strategy;
			 maxmin = this.checkRewardReach(model, modelRewards, exprTemp, false, true, STPGModelChecker.R_CUMULATIVE).valuesD;
			 if (generateStrategy)
				 maxStrat = strategy;

			 repeat = false;
			 // checking which states are marked for removal
			 for (int i = 0; i < n; i++)
				 if (!removed.get(i) && minmax[i] > maxmin[i]) {
					 removed.set(i);
					 removedNew.set(i);
					 repeat = true;
				 }

			 // disabling choices that have transitions to those states
			 removedNew.flip(0, n);
			 for (int i = 0; i < n; i++)
				 for (int j = 0; j < model.getNumChoices(i); j++)
					 if (!stpg.allSuccessorsInSet(i, j, removedNew))
						 stpg.disableChoice(i, j);
			 removedNew.clear();
			 // 4) repeat 2-3 while the set of states from 3 is empty
		 } while (repeat);

		 // 5) if bound is null, return the interval, otherwise check
		 // whether bound is in the interval.
		 BitSet ret = new BitSet(n);
		 for (int i = 0; i < n; i++)
			 ret.set(i, !removed.get(i) && minmax[i] <= p && maxmin[i] >= p);

		 // enabling choices that have been disabled for model checking
		 stpg.enableAllChoices();

		 if (generateStrategy) {
			 strategy = new ExactValueStrategy(minStrat, minmax, maxStrat, maxmin, p, (STPG) model);
		 }

		 return StateValues.createFromBitSet(ret, model);
	 }




    private void printReachabilityPolyhedra(Map<Integer,Polyhedron> polyhedra, List<List<Polyhedron>> stochasticStates, int dim)
    {
	for(int s = 0; s < polyhedra.size(); s++){
	    printReachabilityPolyhedron(polyhedra.get(s), dim, s);
	    for(int t = 0; t < stochasticStates.get(s).size(); t++) {
		System.out.printf("    ->");
		printReachabilityPolyhedron(stochasticStates.get(s).get(t), dim, t);
	    }
	}
    }

    public static void printReachabilityPolyhedron(Polyhedron polyhedron, int dim, int s)
    {
	System.out.printf("%d: [", s);
	for(Generator g : polyhedron.minimized_generators()){
	    System.out.printf("[");
	    BigInteger den = g.divisor().getBigInteger();
	    Map<Variable, BigInteger> num = new HashMap<Variable, BigInteger>();
	    PPLSupport.getCoefficientsFromLinearExpression(g.linear_expression(), false, BigInteger.ONE, num);
	    boolean init = true;
	    for(int i = 0; i<dim; i++){
		if(!init){
		    System.out.printf(", ");
		}
		init = false;
		boolean foundvalue = false;
		for(Variable j : num.keySet()){
		    if(j!=null && i==j.id()){
			BigFraction val = new BigFraction(num.get(j), den);
			System.out.printf("%.4f", val.doubleValue());
			foundvalue = true;
			break;
		    }
		}
		if(!foundvalue){
		    System.out.printf("%.4f", 0.0);
		}
	    }
	    System.out.printf("]");
	}
	System.out.printf("]\n");
    }


    private void printMatlab(Map<Integer,Polyhedron> polyhedra, int dim, int iter)
    {
	
	int max_points = 0;
	for(int s = 0; s < polyhedra.size(); s++) {
	    int points_s = polyhedra.get(s).minimized_generators().size();
	    if(points_s > max_points){
		max_points = points_s;
	    }
	    
	}
	System.out.printf("%% maxpoints: %d\n", max_points);
	

	for(int s = 0; s < polyhedra.size(); s++) {
	    //System.out.printf("points{%d, %d} = %d;\n", iter+1, s+1, polyhedra.get(s).minimized_generators().size());
	    System.out.printf("m{%d, %d} = [", iter+1, s+1); // indices must be greater than zero
	     boolean init1 = true;
	     for(Generator g : polyhedra.get(s).minimized_generators()){
		 if(!init1){
		     System.out.printf(" ; ");
		 }
		 init1 = false;
		 BigInteger den = g.divisor().getBigInteger();
		 Map<Variable, BigInteger> num = new HashMap<Variable, BigInteger>();
		 PPLSupport.getCoefficientsFromLinearExpression(g.linear_expression(), false, BigInteger.ONE, num);
		 boolean init2 = true;
		 for(int i = 0; i<dim; i++){
		     if(!init2){
			 System.out.printf(", ");
		     }
		     init2 = false;
		     boolean foundvalue = false;
		     for(Variable j : num.keySet()){
			 if(j!=null && i==j.id()){
			     BigFraction val = new BigFraction(num.get(j), den);
			     System.out.printf("%.4f", val.doubleValue());
			     foundvalue = true;
			     break;
			 }
		     }
		     if(!foundvalue){
			 System.out.printf("%.4f", 0.0);
		     }
		 }
	     }
	     System.out.printf("];\n");
	 }
    }


    private double stoppingCriterion(STPG stpg, int s0, List<STPGRewards> stpgRewards, long[] accuracy)
    {
	int gameSize = stpg.getStatesList().size();

	// first, sort the states topologically
	// ignore cycles all together (justify!)
	BitSet visited = new BitSet(gameSize); // grey
	BitSet marked = new BitSet(gameSize); // black
	Map<Integer,Set<Integer>> backEdges = new HashMap<Integer,Set<Integer>>();
	List<Integer> top = new ArrayList<Integer>(gameSize); // reverse-topologically sorted states
	DFS_visit(stpg, s0, visited, marked, backEdges, top);

	System.out.println(top);
	System.out.println(backEdges);
	
	// now, traverse trough the topological order, and compute stuff
	Map<Integer, Double> distance = new HashMap<Integer, Double>(gameSize); // delta
	Map<Integer, Integer> length = new HashMap<Integer, Integer>(gameSize); // length of path

	for(Integer s : top){
	    // check out all successors, but stop at back edges
	    List<Distribution> dists = ((SMG) stpg).trans.get(s);
	    double max_dis = 0.0;
	    int max_len = 0;
	    for(int d = 0; d < dists.size(); d++){
		for(Integer t : dists.get(d).keySet()){
		    // t is adjacent to s
		    //if(backEdges.containsKey(s) && backEdges.get(s).contains(t)){ // s->t is a back edge
			// do nothing (justify)
		    //} else {

		    // since now i don't disregard back edges, the field for t might not be available yet
		    int next_len = length.containsKey(t) ? length.get(t) + 1 : 1;
		    if(max_len < next_len){
			max_len = next_len;
		    }
		    
		    double next_dis = Math.log(1.0/dists.get(d).get(t)) + (distance.containsKey(t) ? distance.get(t) : 0); 
		    if(max_dis < next_dis){
			max_dis = next_dis;
		    }
		    //}
		}
	    }
	    distance.put(s, max_dis);
	    length.put(s, max_len);
	}

	System.out.println(distance);

	// distances are log(1/dist), so need to take exponent and invert
	for(int i = 0; i < gameSize; i++){
	    distance.put(i, 1.0/Math.exp(distance.get(i)));
	}

	System.out.println(distance);
	
	int L = length.get(s0);
	double delta = distance.get(s0);
	double epsilon = 1.0/((double)accuracy[0]);
	double M = 1.0; // start at one, because terminals incur reward of one - assume that terminals are reached
	for(int s = 0; s < gameSize; s++){
	    for(STPGRewards stpgr : stpgRewards){
		double rew = stpgr.getStateReward(s);
		if(rew > M) M = rew;
	    }
	}
	M = ((double) L) * M / (1.0 - delta);


	System.out.printf("L: %d, delta: %e, epsilon: %f, M: %f ", L, delta, epsilon, M);
	
	return Math.abs(((double) L) * Math.log(epsilon / M) / Math.log(1 - delta));

    }

    private void DFS_visit(STPG stpg, int s, BitSet visited, BitSet marked, Map<Integer, Set<Integer>> backEdges, List<Integer> L)
    {
	visited.set(s);

	List<Distribution> dists = ((SMG) stpg).trans.get(s);
	for(int d = 0; d < dists.size(); d++){
	    for(Integer t : dists.get(d).keySet()){
		// t is adjacent to s
		if(!visited.get(t) && !marked.get(t)){
		    DFS_visit(stpg, t, visited, marked, backEdges, L);
		} else if (visited.get(t) && !marked.get(t)){ // back edge from t to s
		    if(backEdges.containsKey(t)){ // back edge u->w already registered
			backEdges.get(t).add(s);
		    } else {
			Set<Integer> edges = new HashSet<Integer>();
			edges.add(s);
			backEdges.put(t, edges);
		    }
		}
	    }
	}

	marked.set(s);
	L.add(s);
    }


    private double generatorDistance(Generator g1, Generator g2, long[] accuracy)
    {
	Linear_Expression le1 = g1.linear_expression();
	Linear_Expression le2 = g2.linear_expression();
	BigInteger d1 = g1.divisor().getBigInteger();
	BigInteger d2 = g2.divisor().getBigInteger();
	Map<Variable, BigInteger> c1 = new HashMap<Variable, BigInteger>();
	Map<Variable, BigInteger> c2 = new HashMap<Variable, BigInteger>();

	// get coefficients
	PPLSupport.getCoefficientsFromLinearExpression(le1, false, BigInteger.ONE, c1);
	PPLSupport.getCoefficientsFromLinearExpression(le2, false, BigInteger.ONE, c2);

	// index by integers and not variables
	Map<Integer, BigInteger> ic1 = new HashMap<Integer, BigInteger>();
	for(Variable v : c1.keySet()){
	    if(v != null){
		ic1.put(v.id(), c1.get(v));
	    }
	}

	Map<Integer, BigInteger> ic2 = new HashMap<Integer, BigInteger>();
	for(Variable v : c2.keySet()){
	    if(v != null){
		ic2.put(v.id(), c2.get(v));
	    }
	}

	//
	double dist = 0.0;
	Set<Integer> vars = new HashSet<Integer>();
	vars.addAll(ic1.keySet());
	vars.addAll(ic2.keySet());
	for(Integer v : vars){
	    double a = 0.0;
	    if(ic1.containsKey(v)){
		BigFraction temp = new BigFraction(ic1.get(v), d1);
		a = temp.doubleValue();
	    }
	    if(ic2.containsKey(v)){
		BigFraction temp = new BigFraction(ic2.get(v), d2);
		a -= temp.doubleValue();
	    }
	    //dist += (a*a);
	    dist += (a*a*(((double)accuracy[0])/((double)accuracy[v]))*(((double)accuracy[0])/((double)accuracy[v])));
	}

	return Math.sqrt(dist);

    }

    private boolean stop(Map<Integer,Polyhedron> X, Map<Integer,Polyhedron> prevX, long[] accuracy)
    {
	int gameSize = X.size();
	
	for(int s = 0; s < gameSize; s++){
	    double max_dist = 0.0;
	    for(Generator g : X.get(s).minimized_generators()){
		double min_dist = Double.MAX_VALUE;
		for(Generator prevg : prevX.get(s).minimized_generators()){
		    double dist = generatorDistance(g, prevg, accuracy);
		    if(min_dist > dist){
			min_dist = dist;
		    }
		}
		if(max_dist < min_dist){
		    max_dist = min_dist;
		}
	    }
	    if(max_dist >= 1.0/((double)accuracy[0])){ // if any distance is large enough, continue
		System.out.printf("%% Distance of %e >= %e in polyhedron %d found. Continuing...\n", max_dist, 1.0/((double)accuracy[0]), s);
		return false;
	    }
	}
	return true; // if all distances are smaller, stop
    }


    public Map<Integer,Polyhedron> computeReachabilityPolyhedra(boolean min1, boolean min2, STPG stpg, List<STPGRewards> stpgRewards, List<BitSet> targets, long[] accuracy, double maxIter, List<List<Polyhedron>> stochasticStates) throws PrismException
     {

	 int gameSize = stpg.getStatesList().size();
	 Map<Integer,Polyhedron> result = new HashMap<Integer,Polyhedron>(gameSize);
	 Map<Integer,Polyhedron> prev_result = new HashMap<Integer,Polyhedron>(gameSize);

	 // a list of generator systems, one for each state
	 // will be turned into polyhedra later
	 List<Generator_System> gss = new ArrayList<Generator_System>(gameSize);
	 for (int s = 0; s < gameSize; s++)
	     gss.add(new Generator_System());

	 // initialize coordinates
	 List<Variable> dimensions = new ArrayList<Variable>(targets.size()+stpgRewards.size());
	 for (int i = 0; i < targets.size()+stpgRewards.size(); i++){
	     dimensions.add(new Variable(i)); // variable is associated with target i
	 }

	 double[] maxmin;

	 // initialize polyhedra
	 for (int s = 0; s < gameSize; s++){
	     // base point
	     Linear_Expression base = new Linear_Expression_Times(new Coefficient(0), dimensions.get(0));
	     for(int i = 0; i < targets.size()+stpgRewards.size(); i++){
		 base = new Linear_Expression_Sum(base, new Linear_Expression_Times(new Coefficient(0), dimensions.get(i)));
	     }
	     // add base point to generator system
	     gss.get(s).add(Generator.point(base, new Coefficient(1)));

	     // containers for rewards
	     STPGRewards stpgr;
	     BigFraction r;
	     BigInteger num;
	     BigInteger den;

	     // now see if targets are satisfied and if so add corresponding generator
	     for(int i = 0; i < targets.size()+stpgRewards.size(); i++){
		 den = BigInteger.ONE;
		 if(i>=targets.size()){ // rewards
		     stpgr = stpgRewards.get(i-targets.size());
		     r = new BigFraction(stpgr.getStateReward(s), 1.0/1000000000.0, Integer.MAX_VALUE);
		     num = r.getNumerator();
		     den = r.getDenominator();
		 } else { // target
		     num = BigInteger.valueOf(targets.get(i).get(s)?1:0);
		 }
		 // System.out.printf("Add one: state %d, target %d\n", s, i);


		 //System.out.println(num);
		 //System.out.println(es.ascii_dump());

		 if(num.compareTo(BigInteger.ZERO)!=0){
		     List<Generator> to_add = new ArrayList<Generator>();
		     for(Generator g : gss.get(s)){
			 // TODO: don't add but saturate - don't want to get things like 1 + 1
			 Linear_Expression sum;

			 Linear_Expression es = new Linear_Expression_Times(new Coefficient(g.divisor().getBigInteger().multiply(num)), dimensions.get(i));		     
			 if(den.compareTo(BigInteger.ONE)==0){
			     sum = new Linear_Expression_Sum(es, g.linear_expression());
			 } else {
			     sum = new Linear_Expression_Sum(es, g.linear_expression().times(new Coefficient(den)));
			 }
		     
			 to_add.add(Generator.point(sum, new Coefficient(g.divisor().getBigInteger().multiply(den))));
		     }
		     gss.get(s).addAll(to_add);
		 }
	     }

	     // see how poly looks like
	     C_Polyhedron cp = new C_Polyhedron(gss.get(s));
	     //System.out.printf("State: %d\n", s);
	     //System.out.println(cp.ascii_dump());
	     result.put(s,cp);
	 }


	 double step_increase = 4.0;
         double increase_factor = 1.1;
	 int stop_increasing_after = 60;
	 boolean round = true; // round in all iterations

	 maxIter = 20;
	 int last_iter = 0;

	 // ITERATE FUNCTIONAL APPLICATION
	 for(int iter = 0; iter < Math.ceil(maxIter); ) {

	     System.out.printf("%% Starting iteration %d, rounding %s\n", iter, round ? "on" : "off");

	     long itertime = System.nanoTime();
	     stochasticStates.clear();

	     result = ((SMG) stpg).pMultiObjective(min1, min2, result, targets, stpgRewards, accuracy, stochasticStates, prev_result, round);

	     //System.out.printf("time{%d} = %.4f; // ms\n", iter+1, ((double)System.nanoTime()-itertime)/1000000.0);
	     printMatlab(result, targets.size()+stpgRewards.size(), iter);

	     // STOPPING CRITERION
	     if(iter>0 && stop(result, prev_result, accuracy)){
		 break;
	     }
	     for(int s = 0; s < gameSize; s++) {
		 prev_result.put(s,result.get(s));
	     }

	     iter++;
	     if(iter < stop_increasing_after && (double)(iter-last_iter) % step_increase < 1.0 ) { // increase accuracy by incrase_factor every step_increase iterations
                 for(int i = 0; i < targets.size()+stpgRewards.size(); i++) {
                     accuracy[i] *= increase_factor;
                 }
		 step_increase *= increase_factor;
		 last_iter = iter;
                 System.out.printf("%% ACCURACY SET TO %d, increase again after: %f\n", accuracy[0], step_increase);

             }

	 }

	 // formatted output
	 System.out.println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
	 System.out.println("Final Results:");
	 System.out.println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
	 printReachabilityPolyhedra(result, stochasticStates, targets.size()+stpgRewards.size());
	 System.out.println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");

	 // return the set of polyhedra
	 return result;
     }



	 /**
	  * Computes convex hull of all achievable reachability tuples for targets
	  * 
	  * @param min1
	  *            false if player 1 maximising, true otherwise
	  * @param min2
	  *            false if player 2 maximising, true otherwise
	  * @param stpg
	  *            stochastic two player game model
	  * @param targets
	  *            set of reachability targets (only disjoint is supprted atm)
	  * @return
	  * @throws PrismException
	  */

    public List<Set<ReachTuple>> computeReachabilityTuples(boolean min1, boolean min2, STPG stpg, List<BitSet> targets, double maxIter)
			 throws PrismException
    {

		 int gameSize = stpg.getStatesList().size();

		 // checking for disjointness
		 BitSet bs = new BitSet(gameSize);
		 bs.set(0, gameSize - 1, true);
		 for (BitSet tar : targets)
			 bs.and(tar);
		 if (bs.cardinality() != 0)
			 throw new PrismException("Target sets have to be disjoint!");

		 // intialise results storage
		 List<Set<ReachTuple>> result = new ArrayList<Set<ReachTuple>>(gameSize);
		 for (int i = 0; i < gameSize; i++)
			 result.add(new HashSet<ReachTuple>());

		 // compute initial tuples for states
		 // System.out.println("Computing initial values for states..");
		 double[] maxmin;
		 ReachTuple tuple;
		 BitSet ones = new BitSet(gameSize);
		 ones.set(0, gameSize - 1, true);

		 // 0 - from initial intervals
		 // 1 - from 1  // #clemens: this initialization is super weird to me and I think gives wrong results
		 int initValues = 0;

		 for (int t = 0; t < targets.size(); t++) {

		     //@clemens: base case of value iteration
		     // i.e.\ start of value iteration
			 maxmin = computeUntilProbs(stpg, ones, targets.get(t), false, true, 0).soln;
			 System.out.println(Arrays.toString(maxmin));
			 for (int s = gameSize - targets.size(); s < gameSize; s++) {
				 // if(maxmin[s] == 0) continue;
				 tuple = new ReachTuple(targets.size());
				 //tuple.setValue(t, maxmin[s]);
				 tuple.setValue(t, targets.get(t).get(s)?1.0:0.0); // #clemens: removed until precomputation
				 result.get(s).add(tuple);
			 }
			 for (int s = 0; s < gameSize - targets.size(); s++) {
				 // if(maxmin[s] == 0) continue;
				 tuple = new ReachTuple(targets.size());
				 switch (initValues) {
				 case 0:
				         //tuple.setValue(t, maxmin[s]); // #clemens: removed until
				     tuple.setValue(t, targets.get(t).get(s)?1.0:0.0);
					 break;
				 case 1:
					 tuple.setValue(t, 1);
					 break;
				 // case 2 : tuple.setValue(t, 0); break;
				 }
				 result.get(s).add(tuple);
			 }
			 // TODO consider states which are in multiple targets!
		 }

		 // System.out.println(result.toString());
		 // starting value iterationinit
		 BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
		 Set<ReachTuple> sample;
		 for(int i = 0; i < Math.ceil(maxIter); i++){
		     System.out.printf("iteration %d\n", i);
		     /*
			 System.out.println(result.toString());
			 System.out.println("--------------");
			 sample = result.get(0);
			 System.out.println(sample.size());
			 System.out.println(sample);
			 System.out.println(result.get(1));
			 System.out.println("--------------");
		     */
			 // X^{k+1} = f(X^k) for all states (vector with X's is in result)
			 result = ((SMG) stpg).mvMultiObjective(min1, min2, result);
			 // try {
			 // in.readLine();
			 // } catch (IOException e) {
			 // // TODO Auto-generated catch block
			 // e.printStackTrace();
			 // }

			 // perform upward perturbation when half way
			 // if (maxIter == 10)
			 // for (Set<ReachTuple> ts : result)
			 // for (ReachTuple t : ts)
			 // if (initValues == 0)
			 // t.perturbateUp();
			 // else
			 // t.perturbateDown();


			 /*			 System.out.printf("Printing results of iteration %d\n", 5-maxIter);
			 for (int i = 0; i < result.size(); i++) {
			     System.out.println(i + " " + result.get(i));
			 }
			 System.out.println("------------------");*/
				 
		 }

		 return result;
	 }

	 public void computeIntervalSet(boolean min1, boolean min2, STPG stpg, BitSet target, BitSet zero)
			 throws PrismException
	 {
		System.out.println(stpg);
		System.out.println(target);
		System.out.println(zero);

		int it;

		int n = stpg.getNumStates();
		List<List<Interval>> intervals = new ArrayList<List<Interval>>(n);
		List<List<Interval>> intervals_ = new ArrayList<List<Interval>>(n);

		// initialising intervals
		// compute initial tuples for states
		System.out.println("Computing initial values for states..");
		double[] minmin;
		double[] minmax;
		double[] maxmin;
		double[] maxmax;
		BitSet ones = new BitSet(n);
		ones.set(0, n - 1, true);

		minmin = computeUntilProbs(stpg, ones, target, true, true, 0).soln;
		minmax = computeUntilProbs(stpg, ones, target, true, false, 0).soln;
		maxmin = computeUntilProbs(stpg, ones, target, false, true, 0).soln;
		maxmax = computeUntilProbs(stpg, ones, target, false, false, 0).soln;

		List<Interval> init;
		for (int s = 0; s < n; s++) {
			init = new LinkedList<Interval>();
			init.add(new Interval(minmin[s], minmax[s]));
			init.add(new Interval(maxmin[s], maxmax[s]));
			intervals.add(init);
		}
		// TODO consider states which are in multiple targets!

		BufferedReader in = new BufferedReader(new InputStreamReader(System.in));

		// set number of iterations
		it = 1000;
		System.out.println(intervals);
		int i=0;
		while (true) {
			intervals_ = ((SMG) stpg).mvMultIntervals(min1, min2, intervals);

//			try {
//				in.readLine();
//			} catch (IOException e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			}
			
			if(i++ > it) break;
			// do some checks...
			//@clemens: System.out.println(intervals_);

			intervals = intervals_;
		}
	}
}
