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

package explicit;

import java.util.Arrays;
import java.util.List;
import java.util.BitSet;
import explicit.rewards.SMGRewards;

import parser.ast.Expression;
import parser.ast.ExpressionFunc;
import prism.PrismException;

/**
 * Gathers parameters used for the multi objective model checker.
 */
public class MultiParameters
{
    // OBJECTIVE TYPES
    public static final int ETCR = 1; // Expected total cumulative reward
    public static final int EAR = 2;  // Expected average reward
    public static final int ERCR = 3;  // Expected ratio cumulative reward
    public static final int PAR = 4;  // Probability (satisfaction) average reward
    public static final int PRCR = 5;  // Probability (satisfaction) ratio cumulative reward

    // GOAL
    List<List<Expression>> expr; // expression this parameter set belongs to
    List<Expression> expressions; // list of all expressions
    List<SMGRewards> rewards; // reward structures
    List<String> reward_names; // reward structure names
    List<SMGRewards> divisors; // reward structures acting as divisors for ratio rewards, null if not available
    List<String> divisor_names; // divisor structure names
    List<Integer> reward_types; // types of rewards
    List<Integer> directions; // direction of inequalities
    List<Double> shifts; // shifts to be applied to rewards
    List<Integer> ji; // indices for ratio MQs
    List<Boolean> nested; // is the R inside a P>=1?
    int CONJUNCTS; // number of conjuncts in CNF
    int[] DISJUNCTS; // number of disjuncts in each conjunct of the CNF
    List<Double> bounds; // bounds of the objective
    int objective_type; // if the objective contains of only one sort of expressions, this is their type

    // PARAMETERS
    boolean checkBounds; // whether to check bounds during iteration and abort once they are met
    boolean generateStrategy; // whether to generate a strategy when finished
    double[][] MIN; // used to start CQ iteration
    long maxCIter; // maximum number of iteration of CQ algorithm
    int maxDIter; // maximum number of iteration of DQ algorithm
    int dIterOffset;
    int maxRIter; // maximum number of iteration of Ratio algorithm
    double varepsilon; // relative termination criterion
    long M; // bounding box for energy objectives
    boolean no_union_with_previous = false; // do not take union with previous Pareto approximation

    // ACCURACY
    boolean rounding; // round
    long baseline_accuracy; // baseline accuracy of CQ algorithm
    double increase_factor;
    double[] baseline_biggest_reward; // biggest (maxmin) reward in each dimension
    double[] baseline_smallest_reward; // smallest (minmin) reward in each dimension
    double[] biggest_reward; // biggest reward in each dimension for each unique reward, up to the sign

    public String toString()
    {
	String s = "";
	s += String.format("expr: %s\n", expr);
	s += String.format("expressions: %s\n", expressions);
	//	s += String.format("rewards: %s\n", rewards);
	s += String.format("reward_names: %s\n", reward_names);
	s += String.format("divisors: %s\n", divisors);
	s += String.format("divisor_names: %s\n", divisor_names);
	s += String.format("reward_types: %s\n", reward_types);
	s += String.format("directions: %s\n", directions);
	s += String.format("shifts: %s\n", shifts);
	s += String.format("ji: %s\n", ji);
	s += String.format("nested: %s\n", nested);
	s += String.format("CONJUNCTS: %d\n", CONJUNCTS);
	s += String.format("DISJUNCTS: %s\n", Arrays.toString(DISJUNCTS));
	s += String.format("bounds: %s\n", bounds);
	s += String.format("objective_type: %d\n", objective_type);
	s += String.format("checkBounds: %s\n", checkBounds);
	s += String.format("generateStrategy: %s\n", generateStrategy);
	//	s += String.format("MIN: %s\n", Arrays.deepToString(MIN));
	s += String.format("maxCIter: %d\n", maxCIter);
	s += String.format("maxDIter: %d\n", maxDIter);
	s += String.format("dIterOffset: %s\n", dIterOffset);
	s += String.format("maxRIter: %d\n", maxRIter);
	s += String.format("varepsilon: %f\n", varepsilon);
	s += String.format("M: %f\n", M);
	s += String.format("rounding: %s\n", rounding);
	s += String.format("baseline_accuracy: %s\n", baseline_accuracy);
	s += String.format("increase_factor: %s\n", increase_factor);
	s += String.format("baseline_biggestReward: %s\n", Arrays.toString(baseline_biggest_reward));
	s += String.format("baseline_smallest_reward: %s\n", Arrays.toString(baseline_smallest_reward));
	s += String.format("biggest_reward: %s\n", Arrays.toString(biggest_reward));

	return s;
    }

    public String getParameterString()
    {
	boolean isConjunction = true;
	for(int i = 0; i < CONJUNCTS; i++)
	    if(DISJUNCTS[i] != 1) {
		isConjunction = false;
		break;
	    }

	boolean containsRatio = false;	
	for(Integer reward_type : reward_types)
	    if(reward_type == ERCR) {
		containsRatio = true;
		break;
	    }
	
	return String.format("\nmaximum C-iterations: %d", maxCIter) +
	    String.format("\n\trelative termination threshold: %f", varepsilon) +
	    String.format("\n\tbounding box: %n", M) +
	    (!isConjunction ? String.format("\nmaximum D-iterations: %d", maxDIter) +
	     String.format("\n\tD-iteration offset: %d", dIterOffset) : "") +
	    (containsRatio ? String.format("\nmaximum R-iterations: %d", maxRIter) : "") +
	    (rounding ? String.format("\nrounding") +
	     String.format("\n\tbaseline accuracy: %d", baseline_accuracy) +
	     String.format("\n\tincrease factor: %f", increase_factor) : "");
    }

    protected void shallow_copy(MultiParameters params) throws PrismException
    {
	expr = params.expr;
	expressions = params.expressions;
	rewards = params.rewards;
	reward_names = params.reward_names;
	divisors = params.divisors;
	divisor_names = params.divisor_names;
	reward_types = params.reward_types;
	directions = params.directions;
	shifts = params.shifts;
	ji = params.ji;
	nested = params.nested;
	CONJUNCTS = params.CONJUNCTS;
	DISJUNCTS = params.DISJUNCTS;
	bounds = params.bounds;
	objective_type = params.objective_type;

	checkBounds = params.checkBounds;
	generateStrategy = params.generateStrategy;
	MIN = params.MIN;
	maxCIter = params.maxCIter;
	maxDIter = params.maxDIter;
	dIterOffset = params.dIterOffset;
	maxRIter = params.maxRIter;
	varepsilon = params.varepsilon;
	M = params.M;
	no_union_with_previous = params.no_union_with_previous;

	rounding = params.rounding;
	baseline_accuracy = params.baseline_accuracy;
	increase_factor = params.increase_factor;
	baseline_biggest_reward = params.baseline_biggest_reward;
	baseline_smallest_reward = params.baseline_smallest_reward;
	biggest_reward = params.biggest_reward;
    }

}
