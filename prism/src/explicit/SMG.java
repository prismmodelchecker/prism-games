// ==============================================================================
//	
// Copyright (c) 2002-
// Authors:
// * Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford)
//	
// ------------------------------------------------------------------------------
//	
// This file is part of PRISM.
//	
// PRISM is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//	
// PRISM is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//	
// You should have received a copy of the GNU General Public License
// along with PRISM; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//	
// ==============================================================================

package explicit;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.math.BigInteger;

import prism.ModelType;
import prism.PrismException;
import explicit.rewards.STPGRewards;
import explicit.rewards.SMGRewards;

import org.apache.commons.math3.fraction.BigFraction;
import parma_polyhedra_library.*;

/**
 * Simple explicit-state representation of a stochastic multi-player game (SMG).
 * States can be labelled arbitrarily with player 1..n, player 0 has a special
 * purpose of scheduling the moves of other players
 */
public class SMG extends STPGExplicit implements STPG
{

	// State labels: states with label i are controlled by player i
	protected List<Integer> stateLabels;

	// Set of players which form a coalition
	protected Set<Integer> coalition;

	// player-integer mapping
	protected Map<String, Integer> players;

	public SMG()
	{
		super();
		stateLabels = new ArrayList<Integer>(numStates);
	}

	public SMG(int n)
	{
		super(n);
		stateLabels = new ArrayList<Integer>(numStates);
	}

	/**
	 * Construct an SMG from an existing one and a state index permutation, i.e.
	 * in which state indexsetPlayer i becomes index permut[i].
	 */
	public SMG(SMG smg, int permut[], Map<String, Integer> players)
	{
		super(smg, permut);
		this.players = players;
		stateLabels = new ArrayList<Integer>(numStates);
		// Create blank array of correct size
		for (int i = 0; i < numStates; i++) {
			stateLabels.add(0);
		}
		// Copy permuted player info
		for (int i = 0; i < numStates; i++) {
			stateLabels.set(permut[i], smg.stateLabels.get(i));
		}
		coalition = new HashSet<Integer>();

	}

	/**
	 * Copy constructor
	 */
	public SMG(SMG smg)
	{
		super(smg);
		this.players = new HashMap<String, Integer>(smg.players);
		stateLabels = new ArrayList<Integer>(smg.stateLabels);
		coalition = new HashSet<Integer>(smg.coalition);
	}

	/**
	 * Returns the list of states that belong to the scheduler
	 * 
	 * @return the list of states that belong to the scheduler
	 */
	public Set<Integer> getSchedulerStates()
	{
		Set<Integer> ret = new HashSet<Integer>();
		return ret;
	}

	/**
	 * Adds one state, assigned to player 0
	 */
	@Override
	public int addState()
	{
		return addState(0);
	}

	/**
	 * Adds specified number of states all assigned to player 0
	 */
	@Override
	public void addStates(int numToAdd)
	{
		for (int i = 0; i < numToAdd; i++)
			stateLabels.add(0);
	}

	/**
	 * Adds state assigned to the specified player
	 * 
	 * @param player
	 *            state owner
	 * @return state id
	 */
	public int addState(int player)
	{
		super.addStates(1);
		stateLabels.add(player);
		return numStates - 1;
	}

	/**
	 * Adds the number of states the same as number of Integer in the list, each
	 * assigned to the corresponding player
	 * 
	 * @param players
	 *            list of players (to which corresponding state belongs)
	 */
	public void addStates(List<Integer> players)
	{
		super.addStates(players.size());
		stateLabels.addAll(players);
	}

	/**
	 * labels the given state with the given player
	 * 
	 * @param s
	 *            state
	 * @param player
	 *            player
	 */
	public void setPlayer(int s, int player)
	{
		if (s < stateLabels.size())
			stateLabels.set(s, player);
	}

	// /**
	// * Sets the coalition (representing Player 1)
	// * @param coalition
	// */
	// public void setCoalition(Set<Integer> coalition) {
	// this.coalition = coalition;
	// }
	//
	/**
	 * Sets the coalition (representing Player 1)
	 * 
	 * @param coalition
	 * @throws PrismException
	 */
	public void setCoalition(Set<String> coalition) throws PrismException
	{

		this.coalition.clear();

		int pl;
		for (String player : coalition) {
			if (players.containsKey(player)) { // get the number of the player
				this.coalition.add(players.get(player));
			} else { // try parsing an integer
				try {
					this.coalition.add(Integer.parseInt(player));
				} catch (NumberFormatException e) {
					throw new PrismException("Player " + player + " is not present in the model");
				}
			}

		}

	}

	/**
	 * Sets the coalition (representing Player 1)
	 * 
	 * @param coalition
	 * @throws PrismException
	 */
	public void setCoalitionInts(Set<Integer> coalition) throws PrismException
	{

		this.coalition = coalition;

	}

	/**
	 * Makes a half-deep (up to one reference level) copy of itself
	 */
	public SMG clone()
	{
		SMG smg = new SMG();
		smg.copyFrom(this);
		smg.actions = new ArrayList<List<Object>>(this.actions);
		smg.allowDupes = this.allowDupes;
		smg.maxNumDistrs = this.maxNumDistrs;
		smg.maxNumDistrsOk = this.maxNumDistrsOk;
		smg.numDistrs = this.numDistrs;
		smg.numTransitions = this.numTransitions;
		smg.stateLabels = new ArrayList<Integer>(this.stateLabels);
		smg.trans = new ArrayList<List<Distribution>>(this.trans);
		return smg;
	}

	// Accessors (for Model)

	@Override
	public ModelType getModelType()
	{
		return ModelType.SMG;
	}

	// Accessors (for STPG)

	/**
	 * Returns the 1 if the state belong to coalition and 2 otherwise
	 */
	@Override
	public int getPlayer(int s)
	{
		return coalition.contains(stateLabels.get(s)) ? PLAYER_1 : PLAYER_2;
	}

	public Set<Integer> getCoalition()
	{
		return this.coalition;
	}

	public Map<String, Integer> getPlayerMapping()
	{
		return this.players;
	}

        public void setPlayerMapping(Map<String, Integer> pl)
	{
		this.players = pl;
	}


    /**
     * take X^k and apply F(X^k)(s) for each state
     */
    public Map<Integer,Polyhedron> pMultiObjective(Map<Integer,Polyhedron> Xk, List<SMGRewards> rewards, BitSet terminals, long[] accuracy, List<List<Polyhedron>> stochasticStates) throws PrismException
    {
	// TODO: load from Prism properties
	boolean gauss_seidel = false;

	Map<Integer,Polyhedron> result;
	if(gauss_seidel) {
	    result = Xk;
	} else {
	    result = new HashMap<Integer, Polyhedron>(Xk.size());
	}

	// iterate for each state separate
	for(int s = 0; s < numStates; s++) {
	    // initialize the polyhedra for the stochastic states of s
	    List<Polyhedron> distPolys = new ArrayList<Polyhedron>(trans.get(s).size());
	    stochasticStates.add(distPolys);

	    // apply F to (X^k)(s)
	    result.put(s, pMultiObjectiveSingle(s, Xk, rewards, terminals, accuracy, distPolys));
	}

	// return X^{k+1}
	return result;
    }

    // distPolys will hold the polyhedra of the stochastic states
    private Polyhedron pMultiObjectiveSingle(int s, Map<Integer,Polyhedron> Xk, List<SMGRewards> rewards, BitSet terminals, long[] accuracy, List<Polyhedron> distPolys) throws PrismException
    {
	double fract_acc = 1.0/10000000.0;
	int fract_iter = Integer.MAX_VALUE;

	int n = rewards.size();

	// ------------------------------------------------------------------------------
	// STOCHASTIC STATE OPERATIONS

	// distributions of successor states
	List<Distribution> dists = trans.get(s);

	// step through all stochastic successors of s
	for(int d = 0; d < dists.size(); d++) {
	    // consider each distribution as a stochastic state
	    // hence get a polyhedron for each distribution
	    Distribution distr = dists.get(d);

	    // the successors of the distribution d
	    ArrayList<Integer> states = new ArrayList<Integer>(distr.keySet());
	    if(states.size() == 0) {
		throw new PrismException("State " + s + " has no successors.");
	    } else if (states.size() == 1) {
		// distribution assigns 1 to first successor
		distPolys.add(Xk.get(states.get(0)));
	    } else {
		// need to compute Minkowski sum

		C_Polyhedron cp = null;
		Linear_Expression lhs, rhs;

		Map<Integer,Polyhedron> m = new HashMap<Integer,Polyhedron>(1);

		// first need to make sure probabilities add to one
		Map<Integer, BigFraction> probs = new HashMap<Integer,BigFraction>(states.size());
		BigFraction residual = BigFraction.ONE;
		for(Integer t : states) {
		    BigFraction prob = new BigFraction(distr.get(t), fract_acc, fract_iter);
		    probs.put(t, prob);
		    residual = residual.subtract(prob);
		}
		//if(s==458) System.out.printf("458 residual: %s vs. %s", residual, probs.get(states.get(0)));
		probs.put(states.get(0), probs.get(states.get(0)).add(residual));

		int supdim = 0;
		for(Integer t : states) {
		    // deep copy!
		    C_Polyhedron p = new C_Polyhedron(Xk.get(t).generators());
		    //if(s==458) System.out.printf("458: %s", p.ascii_dump());
		    p.add_space_dimensions_and_embed(states.size());
		    //if(s==458) System.out.printf("458e: %s", p.ascii_dump());
		    BigFraction prob = probs.get(t);
		    for(int i = 0; i < states.size(); i++) {
			if(i == supdim) {
			    lhs = new Linear_Expression_Times(new Coefficient(prob.getNumerator()), new Variable(n+i));
			    rhs = new Linear_Expression_Coefficient(new Coefficient(prob.getDenominator()));
			} else {
			    lhs = new Linear_Expression_Times(new Coefficient(BigInteger.ONE), new Variable(n+i));
			    rhs = new Linear_Expression_Coefficient(new Coefficient(BigInteger.ZERO));
	
			}
			Constraint c = new Constraint(lhs, Relation_Symbol.EQUAL, rhs);
			p.add_constraint(c);
		    }
		    //if(s==458) System.out.printf("458c: %s", p.ascii_dump());
		    if(cp == null) {
		 	cp = p;
		    } else {
			cp.upper_bound_assign(p);
		    }
		    supdim++;
 		}
		//if(s==458) System.out.printf("458cp: %s", cp.ascii_dump());
		//if(s==458) {m.put(0,cp); SMGModelChecker.printMatlab(m,cp.space_dimension(),0); }
		supdim = 0;
		for(Integer t : states) {
		    lhs = new Linear_Expression_Times(new Coefficient(BigInteger.ONE), new Variable(n+supdim));
		    rhs = new Linear_Expression_Coefficient(new Coefficient(BigInteger.ONE));
		    Constraint c = new Constraint(lhs, Relation_Symbol.EQUAL, rhs);
		    cp.add_constraint(c);
		    supdim++;
		}
		//if(s==458) System.out.printf("458cpc: %s\n", cp.ascii_dump());
		//if(s==458) {m.put(0,cp); SMGModelChecker.printMatlab(m,cp.space_dimension(),0); }
		// project away the unneccessary dimensions
		cp.remove_higher_space_dimensions(n);
		//if(s==458) System.out.printf("458cpp: %s\n", cp.ascii_dump());
		//if(s==458) {m.put(0,cp); SMGModelChecker.printMatlab(m,cp.space_dimension(),0); }
		// now in cp have the minkowski sum for that particular distribution d
		distPolys.add(cp);
	    }
	}

	// ------------------------------------------------------------------------------
	// PLAYER ONE AND PLAYER TWO OPERATIONS
	
	// Xk1s holds the polyhedron of state s
	// need deep copy here because want to retain Minkowski sums
	Polyhedron Xk1s = new C_Polyhedron(distPolys.get(0).generators());
	// restore dimension
	Xk1s.add_space_dimensions_and_project(distPolys.get(0).space_dimension() - Xk1s.space_dimension());
	if (stateLabels.get(s)==PLAYER_1) {
	    // Player 1
	    for (int cp_i = 1; cp_i < distPolys.size(); cp_i++) {
		Xk1s.upper_bound_assign(distPolys.get(cp_i));
	    }
	}
	else { // if (stateLabels.get(s)==PLAYER_2) {
	    // Player 2
	    for (int cp_i = 1; cp_i < distPolys.size(); cp_i++) {
		Xk1s.intersection_assign(distPolys.get(cp_i));
	    }
	} //else {
	//    throw new PrismException("Only two-player games supported yet.");
	//}

	// ------------------------------------------------------------------------------
	// ADD STATE REWARDS


	Generator_System ngs = new Generator_System();
	if(terminals.get(s)){
	    ngs = Xk1s.generators();
	} else {
	    // first set up the reward vector that should be added to each point generator
	    Linear_Expression le = null;
	    Coefficient c = new Coefficient(BigInteger.ONE);
	    for(int i = 0; i < n; i++) {
		SMGRewards reward = rewards.get(i);
		BigFraction r = new BigFraction(reward.getStateReward(s), fract_acc, fract_iter);
		BigInteger num = r.getNumerator();
		BigInteger den = r.getDenominator();
		if(le == null) {
		    le = new Linear_Expression_Times(new Coefficient(num), new Variable(i));
		    c = new Coefficient(den);
		} else {
		    le = new Linear_Expression_Sum(le.times(new Coefficient(den)), new Linear_Expression_Times(new Coefficient(num.multiply(c.getBigInteger())), new Variable(i)));
		    c = new Coefficient(den.multiply(c.getBigInteger()));
		}
	    }
	    // now add reward vector to each point generator
	    
	    for(Generator g : Xk1s.generators()) {
		if(g.type()==Generator_Type.POINT) {
		    Linear_Expression nle = new Linear_Expression_Sum(le.times(g.divisor()), g.linear_expression().times(c));
		    Coefficient nc = new Coefficient(g.divisor().getBigInteger().multiply(c.getBigInteger()));
		    ngs.add(Generator.point(nle, nc));
		} else {
		    ngs.add(g);
		}
	    }
	}

	// ------------------------------------------------------------------------------
	// ROUNDING

	// find baseline accuracy
	long baseline_accuracy = 0;
	for(int i = 0; i < n; i++) {
	    if(accuracy[i] > baseline_accuracy) {
		baseline_accuracy = accuracy[i];
	    }
	}

	Generator_System new_ngs = new Generator_System();
	for(Generator ng : ngs) {
	    if(ng.type()==Generator_Type.POINT) {
		Linear_Expression le = ng.linear_expression();
		Coefficient c = ng.divisor();
		Map<Variable,BigInteger> map = new HashMap<Variable,BigInteger>();
		PPLSupport.getCoefficientsFromLinearExpression(le, false, BigInteger.ONE, map);
		
		// new denominator at baseline accuracy
		Coefficient new_c = new Coefficient(BigInteger.valueOf(baseline_accuracy));
		
		// new linear expression
		Linear_Expression new_le;
		if(map.containsKey(null)) { // there is a coefficient without a variable
		    if(map.get(null).compareTo(BigInteger.ZERO) != 0)
			throw new PrismException("Exception in Polyhedron presentation.");
		    new_le = new Linear_Expression_Coefficient(new Coefficient(map.get(null)));
		} else {
		    new_le = new Linear_Expression_Coefficient(new Coefficient(BigInteger.ZERO));
		}
		for(Variable k : map.keySet()) {
		    if(k != null) {
			if(map.get(k).compareTo(BigInteger.ZERO) < 0) { // negative
			    BigFraction round_test = new BigFraction(map.get(k), c.getBigInteger());
			    long rounded = ((long)(Math.floor(round_test.doubleValue()*accuracy[k.id()])*baseline_accuracy/((double)accuracy[k.id()])));
			    new_le = new_le.subtract(new Linear_Expression_Times(new Coefficient(-rounded), k));
			} else { // positive
			    BigFraction round_test = new BigFraction(map.get(k), c.getBigInteger());
			    long rounded = ((long)(Math.floor(round_test.doubleValue()*accuracy[k.id()])*baseline_accuracy/((double)accuracy[k.id()])));
			    new_le = new_le.sum(new Linear_Expression_Times(new Coefficient(rounded), k));
			}
		    }
		}
		new_ngs.add(Generator.point(new_le, new_c));
	    } else if (ng.type()==Generator_Type.RAY) {
		new_ngs.add(ng);
	    }
	}

	// ------------------------------------------------------------------------------
	// CLEAN UP: DIMENSIONALITY, UNION WITH PREVIOUS RESULT, MINIMIZE REPRESENTATION

	Xk1s = new C_Polyhedron(new_ngs);
	
	// add zero dimensions if rounding deleted them
	if(Xk1s.space_dimension()!=n) {
	    Xk1s.add_space_dimensions_and_project(n - Xk1s.space_dimension());
	}

	// union with previous result
	Xk1s.upper_bound_assign(Xk.get(s));

	// minimize representation
	Xk1s = new C_Polyhedron(Xk1s.minimized_generators());

	// add zero dimensions if minimization deleted them
	if(Xk1s.space_dimension()!=n) {
	    Xk1s.add_space_dimensions_and_project(n - Xk1s.space_dimension());
	}
	
	return Xk1s;
    }

	// Standard methods

	@Override
	public String toString()
	{
		int i, j, n;
		Object o;
		String s = "";
		s = "[ ";
		for (i = 0; i < numStates; i++) {
			if (i > 0)
				s += ", ";
			s += i + "(P-" + stateLabels.get(i) + " " + statesList.get(i) + "): ";
			s += "[";
			n = getNumChoices(i);
			for (j = 0; j < n; j++) {
				if (j > 0)
					s += ",";
				o = getAction(i, j);
				if (o != null)
					s += o + ":";
				s += trans.get(i).get(j);
			}
			s += "]";
		}
		s += " ]\n";
		return s;
	}

	public List<Integer> getStateLabels()
	{
		return stateLabels;
	}

	public void setStateLabels(List<Integer> list)
	{
		stateLabels = list;
	}

}
