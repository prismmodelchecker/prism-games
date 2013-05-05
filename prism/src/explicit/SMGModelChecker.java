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

// TODO: check dependencies
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
import parser.ast.ExpressionFunc;
import parser.ast.ExpressionPATL;
import parser.ast.ExpressionProb;
import parser.ast.ExpressionReward;
import parser.ast.ExpressionTemporal;
import parser.ast.ExpressionUnaryOp;
import parser.ast.ExpressionBinaryOp;
import parser.type.TypeBool;
import parser.type.TypeDouble;
import parser.visitor.ASTTraverse;
import prism.PrismException;
import prism.PrismLangException;
import prism.ModelType;
import strat.ExactValueStrategy;
import strat.MultiObjectiveStrategy;
import strat.Strategy;
import explicit.rewards.SMGRewards;

import parser.ast.RewardStruct;
import parser.State;
import explicit.rewards.ConstructRewards;
import explicit.rewards.STPGRewards;
import explicit.rewards.SMGRewards;
import explicit.rewards.SMGRewardsSimple;
import explicit.PPLSupport;
import explicit.MapMDPSimulator;

import org.apache.commons.math3.fraction.BigFraction;
import parma_polyhedra_library.*;

/**
 * Explicit-state model checker for multi-player stochastic games (SMGs).
 */
public class SMGModelChecker extends STPGModelChecker
{
    protected StateValues checkExpressionMulti(Model model, ExpressionFunc expr) throws PrismException
    {
	// TODO: load from PRISM-properties
	// parameters
	long baseline_accuracy = 1000;
	double maxIter = 100.0;

	System.out.println(model);

	// load the Parma Polyhedra Library (PPL) - used for polyhedra operations
	try {
	    System.loadLibrary("ppl_java");
	    // initialise PPL
	    Parma_Polyhedra_Library.initialize_library();
	} catch (Exception e) {
	    throw new PrismException("Error loading Parma Polyhedra Library. Library properly compiled and linked?");
	}


	// need terminal states to check stopping assumption and reward assumptions
	BitSet terminals = findTerminalStates(model);

	// check stopping assumption and model type
	if(model.getModelType() == ModelType.SMG) {
	    // game is non-stopping if for every strategy pair a terminal state is reached with probability 1
	    // can do single-objective min-min reachability problem for terminals and check if > 0
	    double[] reach_term = computeUntilProbs((STPG)model, null, terminals, true, true, 1.0).soln;
	    for(int i = 0; i < reach_term.length; i++) {
		if(reach_term[i] < 1.0 - 1e-6/*-1e-6*/) {
		    throw new PrismException("The game is not stopping.");
		}
	    }

	} else {
	    throw new PrismException("Only SMGs supported by multi-objective engine.");
	}

	// get the initial state of the model, which is an SMG
	int initial_state;
	if(model.getNumInitialStates() != 1) {
	    throw new PrismException("Multi-objective engine supports only models with a single initial state.");
	} else {
	    initial_state = model.getFirstInitialState();
	}

	// get multi-objective goal
	int n = expr.getNumOperands();
	if(n<=1) {
	    throw new PrismException("Need to have at least two goals for multi-objective engine.");
	}
	List<SMGRewards> rewards = new ArrayList<SMGRewards>(n); // the reward structures
	List<Double> bounds = new ArrayList<Double>(n); // the required bounds

	for(int i = 0; i < n; i++) {
	    Expression expr_i = expr.getOperand(i);
	    if (expr_i instanceof ExpressionProb) {
		String relOp = ((ExpressionProb)expr_i).getRelOp(); // direction and strictness of operator
		boolean minimize = false;
		if(relOp.equals("<") || relOp.equals(">")) {
		    //TODO: properly output log
		    System.out.println("Strict inequalities ignored and turned into nonstrict inequalities.");
		} else if (!relOp.equals(">=") && !relOp.equals("<=")) {
		    throw new PrismException("Only minimization or maximization supported.");
		}
		if(relOp.equals("<") || relOp.equals("<=")) {
		    minimize = true;
		}

		Expression pb = ((ExpressionProb)expr_i).getProb(); // probability bound expression
		double p = -1; // probability bound
		if (pb != null) {
		    p = pb.evaluateDouble(constantValues);
		    if (p < 0 || p > 1)
			throw new PrismException("Invalid probability bound " + p + " in P operator");
		} else {
		    throw new PrismException("Probability bound required");
		}
		bounds.add(minimize ? -p : p); // add probability to vector

		Expression e = ((ExpressionProb)expr_i).getExpression();
		if(e instanceof ExpressionTemporal) {
		    if(((ExpressionTemporal)e).getOperator() == ExpressionTemporal.P_F) {
			BitSet t = checkExpression(model, ((ExpressionTemporal)e).getOperand2()).getBitSet(); // evaluate which states satisfy the property
			// convert target set to reward structure
			SMGRewardsSimple reward = new SMGRewardsSimple(((SMG)model).numStates);
			for(int s = 0; s < ((SMG)model).numStates; s++) {
			    if(t.get(s)) { // attach reward 1.0 to states in target set
				reward.setStateReward(s, minimize ? -1.0 : 1.0);
			    }
			}
			rewards.add(reward);
		    } else {
			// TODO: reduction from LTL to rewards goes here
			throw new PrismException("Invalid property: property " + i + " must be a reachability property.");
		    }
		} else {
		    throw new PrismException("Invalid property: property " + i + " must be an LTL formula.");
		}

	    } else if (expr_i instanceof ExpressionReward) {
		String relOp = ((ExpressionReward)expr_i).getRelOp(); // direction and strictness of operator
		boolean minimize = false;
		if(relOp.equals("<") || relOp.equals(">")) {
		    //TODO: properly output log
		    System.out.println("Strict inequalities ignored and turned into nonstrict inequalities.");
		} else if (!relOp.equals(">=") && !relOp.equals("<=")) {
		    throw new PrismException("Only minimization or maximization supported.");
		}
		if(relOp.equals("<") || relOp.equals("<=")) {
		    minimize = true;
		}

		// evaluate reward bound
		Expression rb = ((ExpressionReward)expr_i).getReward();
		double r = 0.0; // reward bound
		if (rb != null) {
		    r = rb.evaluateDouble(constantValues);
		} else {
		    throw new PrismException("Reward bound required");
		}
		bounds.add(minimize ? -r : r); // add probability to vector

		// check if cumulative reward
		Expression e = ((ExpressionReward)expr_i).getExpression();
		if(e instanceof ExpressionTemporal) {
		    // TODO: parse cumulative reward properly
		    if(((ExpressionTemporal)e).getOperator() == ExpressionTemporal.R_S) { // cumulative reward
			// everything ok here
		    } else {
			throw new PrismException("Only cumulative rewards supported so far.");
		    }
		} else {
		    throw new PrismException("Only temporal expressions supported so far.");
		}

		// get index of reward structure from expression
		Object r_index = ((ExpressionReward)expr_i).getRewardStructIndex();
		// construct state rewards
		RewardStruct rewStruct;
		if(r_index instanceof Integer) {
		    rewStruct = modulesFile.getRewardStruct((Integer)r_index);
		} else if (r_index instanceof String) {
		    rewStruct = modulesFile.getRewardStructByName((String)r_index);
		} else {
		    throw new PrismException("Cannot get reward structure for goal " + i + ".");
		}
		ConstructRewards constructRewards = new ConstructRewards(mainLog);
		SMGRewardsSimple reward = (SMGRewardsSimple)constructRewards.buildSMGRewardStructure((SMG)model, rewStruct, constantValues);
		// check rewards assumptions
		Boolean negative_rewards = null;
		for(int s = 0; s < ((SMG)model).numStates; s++) {
		    if(terminals.get(s) && reward.getStateReward(s)!=0.0) {
			throw new PrismException("Terminal states must have zero rewards. Check state " + s + ".");
		    }
		    if(reward.getStateReward(s)<0) {
			if(negative_rewards==null) {
			    negative_rewards = true;
			} else {
			    if(negative_rewards == false) {
				throw new PrismException("Rewards must be either all negative or all non-positive in each structure.");
			    }
			}
		    } else if(reward.getStateReward(s)>0){
			if(negative_rewards==null) {
			    negative_rewards = false;
			} else {
			    if(negative_rewards == true) {
				throw new PrismException("Rewards must be either all negative or all non-positive in each structure.");
			    }
			}
		    }

		}

		// for minimization need to revert sign:
		if(minimize) {
		    for(int s = 0; s < ((SMG)model).numStates; s++) {
			reward.setStateReward(s, reward.getStateReward(s));
		    }
		}

		// add reward
		rewards.add(reward);

	    } else {
		throw new PrismException("Only the P and R operators are supported so far.");
	    }
	}
	// TODO: set accuracy
	long[] accuracy = new long[n];
	for(int i = 0; i < n; i++) {
	    SMGRewards reward = rewards.get(i);
	    double[] maxmaxreward = (super.computeReachRewardsCumulative((STPG)model, reward, new BitSet(model.getNumStates()), false, false, null, null)).soln;
	    double biggest_reward = 0;
	    for(int s = 0; s < model.getNumStates(); s++) {
		if(biggest_reward < Math.abs(maxmaxreward[s])) {
		    biggest_reward = Math.abs(maxmaxreward[s]);
		}		
	    }
	    if(biggest_reward > 0) {
		accuracy[i] = baseline_accuracy/((long)biggest_reward);
	    } else {
		throw new PrismException("Reward " + i + " is all zeros, remove.");
	    }
	}

	// store polyhedra of stochastic states for strategy construction
	List<List<Polyhedron>> stochasticStates = new ArrayList<List<Polyhedron>>(((SMG)model).numStates);

	// compute polyhedra
	Map<Integer,Polyhedron> result_p = computeParetoSetApproximations((SMG) model, rewards, bounds, terminals, accuracy, maxIter, stochasticStates); // stores Pareto set approximations


	// TODO: Insert actual results here
	//return new StateValues(TypeDouble.getInstance(), new Double(3.141593), model);
	return new StateValues(TypeBool.getInstance(), new Boolean(true), model);
    }


    /*
    // TODO: remove
    protected Map<Integer,Polyhedron> checkMultiObjectiveFormula(Model model, ExpressionPATL exprPATL, boolean min, List<List<Polyhedron>> stochasticStates) throws PrismException
    {

	System.out.println("multiobjcheck");

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


	    // set accuracy


	    long[] accuracy = new long[targets.size()+stpgRewards.size()];
	    long baseline_accuracy = 50;
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

		double[] goal = { 0.75, 0.78 };

		result_p = this.computeReachabilityPolyhedra(min, !min, (STPG) model, stpgRewards, targets, accuracy, maxIter, stochasticStates, goal);
		polyTime = System.nanoTime() - polyTime;
		
		System.out.printf("Polyhedra computation: %.4f ms\n", ((double)polyTime)/1000000.0);

		polyTime = System.nanoTime();
		
		MultiObjectiveStrategy strategy_mdp = new MultiObjectiveStrategy((STPG) model, initial_state, goal, result_p, stochasticStates, stpgRewards);
		
		polyTime = System.nanoTime() - polyTime;

		System.out.printf("Strategy computation: %4f ms\n", ((double)polyTime)/1000000.0);

		MapMDPSimulator mmdps = new MapMDPSimulator((STPG) model, stpgRewards);
		
		mmdps.writeStrategy(strategy_mdp, "mmdps");
	    } else if (compute.equals("recompute")) {
		MapMDPSimulator mmdps = new MapMDPSimulator((STPG) model, stpgRewards);
		mmdps.readStrategy("mmdps");

	        //double[] goal = { 0.019, 0.903, 2.884 };
                //double[] goal = { 0.769, 0.75,  11.69 };
	        //double[] goal = { 0.0, 0.576, 17.21 };
	        double[] goal = { 0.6346, 0.7692, 13.2692 };

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

    */


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
			//TODO get this back
			//Map<Integer,Polyhedron> X = checkMultiObjectiveFormula(model, exprPATL, min, Y);
			

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
	System.out.printf("maxpoints{%d} = %d;\n", iter+1, max_points);
	

	for(int s = 0; s < polyhedra.size(); s++) {
	    //System.out.printf("points{%d, %d} = %d;\n", iter+1, s+1, polyhedra.get(s).minimized_generators().size());
	    System.out.printf("m{%d, %d} = [", iter+1, s+1); // indices must be greater than zero
	     boolean init1 = true;
	     for(Generator g : polyhedra.get(s).minimized_generators()){
		 // ignore rays
		 if(g.type() == Generator_Type.RAY) {
		     continue;
		 }
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


    public Map<Integer,Polyhedron> computeParetoSetApproximations(SMG smg, List<SMGRewards> rewards, List<Double> bounds, BitSet terminals, long[] accuracy, double maxIter, List<List<Polyhedron>> stochasticStates) throws PrismException
    {
	int gameSize = smg.getNumStates();
	int n = rewards.size();
	Map<Integer,Polyhedron> result = new HashMap<Integer,Polyhedron>(gameSize);

	 // a list of generator systems, one for each state - will be turned into polyhedra later
	 List<Generator_System> gss = new ArrayList<Generator_System>(gameSize);

	 // VALUE ITERATION

	 // initialisation: compute polyhedra X_s^0
	 // first precompute MIN(r) for each reward
	 double[][] MIN = new double[n][gameSize];
	 for(int i = 0; i < n; i++) {
	     SMGRewards reward = rewards.get(i);
	     boolean positive_reward = false;
	     for(int s = 0; s < gameSize; s++) {
		 MIN[i][s] = reward.getStateReward(s);
		 if(reward.getStateReward(s) > 0) {
		     positive_reward = true;
		 }
	     }
	     if(positive_reward) {
		 continue;
	     }
	     MIN[i] = (super.computeReachRewardsCumulative((STPG)smg, reward, new BitSet(gameSize), true, true, null, null)).soln;
	     
	 }

	 for(int s = 0; s < gameSize; s++) {
	     Generator_System gs = new Generator_System();

	     // generate corner point from state reward
	     Linear_Expression r_num = null;
	     BigInteger r_den = BigInteger.ONE;
	     // step through rewards and add to r
	     for(int i = 0; i < n; i++) {
		 SMGRewards reward = rewards.get(i);
		 // TODO: What is a good accuracy here?
		 BigFraction ri = new BigFraction(MIN[i][s], 1.0/1000000.0, Integer.MAX_VALUE);

		 BigInteger num = ri.getNumerator();
		 BigInteger den = ri.getDenominator();

		 // add r_num/r_den + num/den:
		 if(r_num==null) {
		     r_num = new Linear_Expression_Times(new Coefficient(num), new Variable(i));
		 } else {
		     Linear_Expression r_num_to_add = new Linear_Expression_Times(new Coefficient(num.multiply(r_den)), new Variable(i));
		     if(den.compareTo(BigInteger.ONE)==0) {
			 // (r_num + num*r_den)/r_den
			 r_num = new Linear_Expression_Sum(r_num, r_num_to_add);
		     } else {
			 // (r_num*den + num*r_den)/(r_den*den)
			 r_num = new Linear_Expression_Sum(r_num.times(new Coefficient(den)), r_num_to_add);
			 r_den = r_den.multiply(den);
		     }
		 }

		 // generate ray for downward closure
		 Linear_Expression ray = new Linear_Expression_Times(new Coefficient((BigInteger.ONE).negate()), new Variable(i));
		 gs.add(Generator.ray(ray));
	     }
	     gs.add(Generator.point(r_num, new Coefficient(r_den)));

	     // generate initial polyhedra: X^0_s
	     C_Polyhedron cp = new C_Polyhedron(gs);
	     result.put(s, cp);
	 }

	 // ITERATE FUNCTIONAL APPLICATION

	 for(int k = 0; k < Math.ceil(maxIter); /* k increased later */ ) {
	     System.out.printf("%% Starting iteration %d with accuracy[0] = %d\n", k, accuracy[0]);

	     // TODO: only store stochastic states at last iteration anyway
	     stochasticStates.clear();
	     
	     result = smg.pMultiObjective(result, rewards, terminals, accuracy, stochasticStates);

	     // TODO: proper logging
	     printMatlab(result, n, k);

	     k++;
	     
	     // increase accuracy
	     // TODO: get increase_factor from properties
	     double increase_factor = 1.1;
	     for(int i = 0; i < n; i++) {
		 accuracy[i] *= increase_factor;
	     }
	 }

	return result;
    }

}

