// ==============================================================================
//	
// Copyright (c) 2002-
// Authors:
// * Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford)
// * Clemens Wiltsche <clemens.wiltsche@cs.ox.ac.uk> (University of Oxford)
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

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.apache.commons.math3.fraction.BigFraction;

import explicit.rewards.SMGRewards;
import parma_polyhedra_library.C_Polyhedron;
import parma_polyhedra_library.Coefficient;
import parma_polyhedra_library.Constraint;
import parma_polyhedra_library.Generator;
import parma_polyhedra_library.Generator_System;
import parma_polyhedra_library.Generator_Type;
import parma_polyhedra_library.Linear_Expression;
import parma_polyhedra_library.Linear_Expression_Coefficient;
import parma_polyhedra_library.Linear_Expression_Times;
import parma_polyhedra_library.Polyhedron;
import parma_polyhedra_library.Relation_Symbol;
import parma_polyhedra_library.Variable;
import parser.ast.Coalition;
import prism.ModelType;
import prism.PrismException;
import prism.PrismLog;

/**
 * Simple explicit-state representation of a (turn-based) stochastic multi-player game (SMG).
 * States can be labelled arbitrarily with player 1..n, player 0 has a special
 * purpose of scheduling the moves of other players
 */
public class SMG extends STPGExplicit implements STPG
{
	// NB: We re-use the existing stateOwners list in the superclass to assign states to players

	// A definition of the players in the game, i.e., the (integer) index and name of each one.
	// It is stored as a mapping from indices to names.
	// Indices can be arbitrary and do not need to be contiguous.
	// Names are optional and can be null or "" if undefined,
	// but all normally defined names should be unique.
	// The size of this map defines the number of players in the game.
	protected Map<Integer, String> playerNames;
	
	// Optionally, a mapping from player indices to 1 or 2,
	// as induced by a coalition of players (those in the coalition
	// are mapped to 1, those who are not are mapped to 2).
	// The mapped index of the player with original index {@code p}
	// is stored as {@code coalitionPlayerMap[p]}, regardless of the
	// range of player indices actually used.
	protected int[] coalitionPlayerMap;

	// Constructors

	/**
	 * Constructor: empty SMG.
	 */
	public SMG()
	{
		super();
		playerNames = new HashMap<Integer, String>();
		coalitionPlayerMap = null;
	}

	/**
	 * Constructor: new SMG with fixed number of states.
	 */
	public SMG(int numStates)
	{
		super(numStates);
		playerNames = new HashMap<Integer, String>();
		coalitionPlayerMap = null;
	}

	/**
	 * Construct an SMG from an existing one and a state index permutation,
	 * i.e. in which state index i becomes index permut[i].
	 * Player and coalition info is also copied across.
	 */
	public SMG(SMG smg, int permut[])
	{
		super(smg, permut);
		playerNames = new HashMap<Integer, String>(smg.playerNames);
		coalitionPlayerMap = smg.coalitionPlayerMap == null ? null : smg.coalitionPlayerMap.clone(); 
	}

	/**
	 * Copy constructor
	 */
	public SMG(SMG smg)
	{
		super(smg);
		playerNames = new HashMap<Integer, String>(smg.playerNames);
		coalitionPlayerMap = smg.coalitionPlayerMap == null ? null : smg.coalitionPlayerMap.clone(); 
	}

	/**
	 * Remove state by disabling all accessing transitions.
	 * 
	 * @param s State to be removed
	 * @param close_P1 Close game for Player 1 (coalition players)
	 * @param close_P2 Close game for Player 2 (non-coalition players)
	 * 
	 * @return set of newly disabled states
	 */
	public BitSet disableState(int s, boolean close_P1, boolean close_P2)
	{
		throw new UnsupportedOperationException();
		
		// Disabled for now because support disabling transitions has been removed
		// during refactoring of explicit engine model classes
		
		/*
		int numStates = trans.size();
		BitSet disabled = new BitSet(numStates);
		disabled.set(s);
		
		for(int t = 0; t < numStates; t++) {
			List<Distribution> trans_t = trans.get(t);
			int num_trans_t = trans_t.size();
			int disabled_trans_t = 0;
			for(int c = 0; c < num_trans_t; c++) {
				if(trans_t.get(c).getSupport().contains(s)) { // some transition t -> s
					// disable this transition
					disableChoice(t, c);					
				}
				// count how many transitions are disabled;
				if(disabledChoices.get(t) != null && disabledChoices.get(t).get(c))
					disabled_trans_t++;
			}
			
			// test if state became disabled, depending on whether the game should be closed for a type of players
			if(getPlayer(t) == 1 && close_P1 && (num_trans_t == disabled_trans_t) |
					getPlayer(t) == 2 && close_P2 && disabled_trans_t > 0) {
				disabled.or(disableState(t, close_P1, close_P2)); // recursively disable states, and add to return set
			}
		}
		
		return disabled;
		*/
	}

	// Mutators

	/**
	 * Add a new (player 1) state and return its index.
	 */
	@Override
	public int addState()
	{
		return addState(0);
	}

	/**
	 * Add multiple new (player 1) states.
	 */
	@Override
	public void addStates(int numToAdd)
	{
		for (int i = 0; i < numToAdd; i++) {
			addState();
		}
	}

	/**
	 * Set the info about players, i.e., the (integer) index and name of each one.
	 * This is given as a mapping from indices to names.
	 * Indices can be arbitrary and do not need to be contiguous.
	 * Names are optional and can be null or "" if undefined,
	 * but all normally defined names should be unique.
	 */
	public void setPlayerInfo(Map<Integer, String> playerNames)
	{
		this.playerNames = new HashMap<Integer, String>(playerNames);
	}

	/**
	 * Set the info about players, provided as a list of names.
	 * Names are optional and can be null or "" if undefined,
	 * but all normally defined names should be unique.
	 */
	public void setPlayerInfo(List<String> playerNamesList)
	{
		this.playerNames = new HashMap<Integer, String>();
		int numPlayers = playerNamesList.size();
		for (int i = 0; i < numPlayers; i++) {
			this.playerNames.put(i + 1, playerNamesList.get(i));
		}
	}

	/**
	 * Copy player info from another SMG.
	 */
	public void copyPlayerInfo(SMG smg)
	{
		setPlayerInfo(smg.playerNames);
	}

	/**
	 * Set a coalition of players for this SMG
	 * (which effectively makes it an STPG with player 1 representing the coalition and 2 the rest).
	 * Pass null to remove any coalition info from this SMG.
	 * 
	 * @param coalition Coalition info object 
	 */
	public void setCoalition(Coalition coalition) throws PrismException
	{
		// Clear info if coalition is null
		if (coalition == null) {
			coalitionPlayerMap = null;
			return;
		}
		
		// This is first time we really need the {@code playerNames} info,
		// so if it has not been set, we tried to create it based on {@code stateOwners} 
		if (playerNames.isEmpty()) {
			for (int i = 0; i < numStates; i++) {
				int p = stateOwners.get(i);
				if (!playerNames.containsKey(p)) {
					playerNames.put(p, null);
				}
			}
		}
		
		// Find max player index
		int maxIndex = 0;
		for (int index : playerNames.keySet()) {
			maxIndex = Math.max(maxIndex, index);
		}
		
		// Construct mapping
		coalitionPlayerMap = new int[maxIndex + 1];
		for (int i = 0; i < maxIndex + 1; i++) {
			coalitionPlayerMap[i] = -1;
		}
		for (Entry<Integer, String> entry : playerNames.entrySet()) {
			int playerIndex = entry.getKey();
			boolean inCoalition = coalition.isPlayerIndexInCoalition(playerIndex, playerNames);
			// In coalition => player 1; not in coalition (or undefined) => player 2
			coalitionPlayerMap[playerIndex] = inCoalition ? 1 : 2;
		}
	}

	/**
	 * Copy coalition info from another SMG.
	 */
	public void copyCoalitionInfo(SMG smg)
	{
		coalitionPlayerMap = smg.coalitionPlayerMap == null ? null : smg.coalitionPlayerMap.clone();
	}

	/**
	 * Makes a half-deep (up to one reference level) copy of itself
	 */
	public SMG clone()
	{
		SMG smg = new SMG();
		smg.copyFrom(this);
		smg.actions = new ChoiceActionsSimple(this.actions);
		smg.allowDupes = this.allowDupes;
		smg.maxNumDistrs = this.maxNumDistrs;
		smg.maxNumDistrsOk = this.maxNumDistrsOk;
		smg.numDistrs = this.numDistrs;
		smg.numTransitions = this.numTransitions;
		smg.stateOwners = new ArrayList<Integer>(this.stateOwners);
		smg.trans = new ArrayList<List<Distribution>>(this.trans);
		return smg;
	}

	// Accessors (for Model)

	@Override
	public ModelType getModelType()
	{
		return ModelType.SMG;
	}

	// Accessors (for STPG/SMG)

	//@Override
	public int getNumPlayers()
	{
		return playerNames.size();
	}
	
	@Override
	public int getPlayer(int s)
	{
		int playerIndex = stateOwners.get(s); 
		if (coalitionPlayerMap == null) {
			// No coalition: just return index
			return playerIndex;
		} else {
			// Coalition defined: look up if player 1 or 2
			// (note: undefined players are mapped to player 2)
			return playerIndex == -1 ? 2 : coalitionPlayerMap[playerIndex];
		}
	}

	@Override
	public void exportToDotFile(PrismLog out, BitSet mark, boolean states)
	{
		BitSet players = new BitSet(numStates);
		for (int i = 0; i < numStates; i++) {
			if (stateOwners.get(i) == 1) {
				players.set(i);
			}
		}
		super.exportToDotFile(out, mark, players, states);
	}

	/**
	 * take X^k and apply F(X^k)(s) for each state (cf. MFCS'13 and QEST'13)
	 *
	 * arguments:
	 * @param gaussSeidel Gauss Seidel update allowed
	 * @param rounding rounding enabled
	 * @param union_with_previous take union with previous Pareto set
	 * @param cut cut off everything that is strictly above the negative orthant (used for energy objectives)
	 * @param M maximum bound on Pareto sets (quantity is positive)
	 */
	public Pareto[] pMultiObjective(Pareto[] Xk, List<SMGRewards> rewards, boolean gaussSeidel,
					long baseline_accuracy, double[] biggest_reward,
					List<Pareto>[] stochasticStates, boolean rounding,
					boolean union_with_previous, boolean cut, long M)
	    throws PrismException
	{
		Pareto[] result = new Pareto[Xk.length];
		Pareto[] Yk = gaussSeidel ? null : new Pareto[Xk.length]; // if Gauss-Seidel, no memory allocation required
		System.arraycopy(Xk, 0, gaussSeidel ? result : Yk, 0, Xk.length); // if Gauss-Seidel, update result in-place
		// iterate for each state separately
		for (int s = 0; s < numStates; s++) {
			// first, check if cancelled
			// initialize the polyhedra for the stochastic states of s
			List<Pareto> distPolys = new ArrayList<Pareto>(trans.get(s).size());
			// apply F to (X^k)(s)
			//double t0 = (double)System.nanoTime();
			result[s] = pMultiObjectiveSingle(s, gaussSeidel ? result : Yk, rewards, baseline_accuracy, biggest_reward, distPolys, rounding,
							  union_with_previous, cut, M);
			//System.out.printf("total: %f s\n", ((double) (System.nanoTime() - t0)) / 1e9);
			// store stochastic states if requested (by the reference being non-null)
			if (stochasticStates != null)
				stochasticStates[s] = distPolys;
		}

		// return X^{k+1}
		return result;
	}

    private Polyhedron round(Generator_System ngs, long baseline_accuracy, double[] biggest_reward, boolean energy_objective) throws PrismException
	{
		int n = biggest_reward.length;
		// accuracy
		long[] accuracy = new long[n];
		for (int i = 0; i < n; i++) {
		        long tmp_a = energy_objective ? baseline_accuracy : ((long) (((double) baseline_accuracy) / biggest_reward[i]));
			// sanity check to prevent overflow
			accuracy[i] = tmp_a < Long.MAX_VALUE && tmp_a > 0 ? tmp_a : Long.MAX_VALUE;
		}

		Generator_System new_ngs = new Generator_System();
		for (Generator ng : ngs) {
			if (ng.type() == Generator_Type.POINT) {
				Linear_Expression le = ng.linear_expression();
				Coefficient c = ng.divisor();
				Map<Variable, BigInteger> map = new HashMap<Variable, BigInteger>();
				PPLSupport.getCoefficientsFromLinearExpression(le, false, BigInteger.ONE, map);

				// new denominator at baseline accuracy
				Coefficient new_c = new Coefficient(BigInteger.valueOf(baseline_accuracy));

				// new linear expression
				Linear_Expression new_le;
				if (map.containsKey(null)) { // there is a coefficient without a variable
					if (map.get(null).compareTo(BigInteger.ZERO) != 0)
						throw new PrismException("Exception in Polyhedron presentation.");
					new_le = new Linear_Expression_Coefficient(new Coefficient(map.get(null)));
				} else {
					new_le = new Linear_Expression_Coefficient(new Coefficient(BigInteger.ZERO));
				}
				for (Variable k : map.keySet()) {
					if (k != null) {
						BigFraction round_test = new BigFraction(map.get(k), c.getBigInteger());
						long rounded = ((long) (Math.floor(round_test.doubleValue() * accuracy[(int)k.id()]) * baseline_accuracy / ((double) accuracy[(int)k.id()])));
						new_le = new_le.sum(new Linear_Expression_Times(new Coefficient(rounded), k));
					}
				}
				new_ngs.add(Generator.point(new_le, new_c));
			} else if (ng.type() == Generator_Type.RAY) {
				new_ngs.add(ng);
			}
		}

		Polyhedron result = new C_Polyhedron(new_ngs);
		// add zero dimensions if rounding deleted them
		if (result.space_dimension() != n) {
			result.add_space_dimensions_and_project(n - result.space_dimension());
		}

		return result;
	}

    protected Pareto stochasticState(int s, Distribution distr, int d, Pareto[] Xk, List<SMGRewards> rewards, double[] extra_rewards, boolean cut, long M)
			throws PrismException
	{
		int n = rewards.size();

		// the successors of the distribution d
		ArrayList<Integer> states = new ArrayList<Integer>(distr.keySet());
		int b = states.size();

		Pareto cp = null;
		if (b == 0) {
			throw new PrismException("Distribution " + s + ", " + d + " has no successors.");
		} else if (b == 1) {
			// distribution assigns 1 to first successor
			cp = Xk[states.get(0)];
		} else { // need to compute Minkowski sum
			Linear_Expression lhs, rhs;

			// first need to make sure probabilities add to one
			BigFraction[] probs = new BigFraction[b];
			BigFraction residual = BigFraction.ONE;
			int supdim = 0;
			for (Integer t : states) {
				BigFraction prob = new BigFraction(distr.get(t));
				probs[supdim] = prob;
				residual = residual.subtract(prob);
				supdim++;
			}
			probs[0] = probs[0].add(residual); // just add residual to first probability

			supdim = 0;
			for (Integer t : states) {
				C_Polyhedron p = new C_Polyhedron((C_Polyhedron) Xk[t].get()); // deep copy!
				p.add_space_dimensions_and_embed(b);
				BigFraction prob = probs[supdim];
				for (int i = 0; i < b; i++) {
					if (i == supdim) {
						lhs = new Linear_Expression_Times(new Coefficient(prob.getNumerator()), new Variable(n + i));
						rhs = new Linear_Expression_Coefficient(new Coefficient(prob.getDenominator()));
					} else {
						lhs = new Linear_Expression_Times(new Coefficient(BigInteger.ONE), new Variable(n + i));
						rhs = new Linear_Expression_Coefficient(new Coefficient(BigInteger.ZERO));
					}
					Constraint c = new Constraint(lhs, Relation_Symbol.EQUAL, rhs);
					p.add_constraint(c);
				}
				if (cp == null)
					cp = new Pareto(p);
				else
					cp.get().upper_bound_assign(p);

				supdim++;
			}

			supdim = 0;
			for (Integer t : states) {
				lhs = new Linear_Expression_Times(new Coefficient(BigInteger.ONE), new Variable(n + supdim));
				rhs = new Linear_Expression_Coefficient(new Coefficient(BigInteger.ONE));
				Constraint c = new Constraint(lhs, Relation_Symbol.EQUAL, rhs);
				cp.get().add_constraint(c);
				supdim++;
			}

			// project away the unneccessary dimensions
			//double t0 = (double)System.nanoTime();
			cp.get().remove_higher_space_dimensions(n);
			//System.out.printf("rhsd: %f s\n", ((double) (System.nanoTime() - t0)) / 1e9);
			// now in cp have the Minkowski sum for that particular distribution d
		}

		// add rewards
		Polyhedron Yk1 = PPLSupport.add_rewards(cp.get(), s, d, rewards, extra_rewards);

		// cut everything but the negative orthant bounded by -M
		if (cut) PPLSupport.cutBox(Yk1, M);

		// add transition rewards and return polyhedron
		return new Pareto(Yk1);
	}

	// distPolys will hold the polyhedra of the stochastic states
	private Pareto pMultiObjectiveSingle(int s, Pareto[] Xk, List<SMGRewards> rewards, long baseline_accuracy, double[] biggest_reward, List<Pareto> distPolys,
					     boolean rounding, boolean union_with_previous, boolean cut, long M) throws PrismException
	{
		int n = rewards.size();

		// distributions of successor states
		List<Distribution> dists = trans.get(s);

		// ------------------------------------------------------------------------------
		// STOCHASTIC STATE OPERATIONS

		// step through all stochastic successors of s
		int d = 0;
		for (Distribution distr : dists) {
			// add polyhedron to the list of polyhedra in the successors of s
		        distPolys.add(stochasticState(s, distr, d, Xk, rewards, null, cut, M));
			d++;
		}

		// ------------------------------------------------------------------------------
		// PLAYER ONE AND PLAYER TWO OPERATIONS

		// Xk1s holds the polyhedron of state s
		// need deep copy here because want to retain Minkowski sums
		Polyhedron Xk1s;
		if (distPolys.size() > 0) {
			if (getPlayer(s) == 1) {
			        // Player 1
			        Xk1s = new C_Polyhedron(distPolys.get(0).get().generators());
				int cp_start = 0;
			        get_first_cp: for (cp_start = 0; cp_start < distPolys.size(); cp_start++) {
				    if(!distPolys.get(cp_start).get().is_empty()) {
					Xk1s = new C_Polyhedron(distPolys.get(cp_start).get().generators());
					break get_first_cp;
				    }
				}
				if(distPolys.get(0).get().space_dimension() > Xk1s.space_dimension())
				        Xk1s.add_space_dimensions_and_project(distPolys.get(0).get().space_dimension() - Xk1s.space_dimension());

				for (int cp_i = cp_start+1; cp_i < distPolys.size(); cp_i++) {
				    if(!distPolys.get(cp_i).get().is_empty())
					Xk1s.upper_bound_assign(distPolys.get(cp_i).get());
				}
			} else {
				// Player 2
 			        Xk1s = new C_Polyhedron(distPolys.get(0).get().generators());
				Xk1s.add_space_dimensions_and_project(distPolys.get(0).get().space_dimension() - Xk1s.space_dimension());
				for (int cp_i = 1; cp_i < distPolys.size(); cp_i++) {
				    if(!Xk1s.is_empty() && !distPolys.get(cp_i).get().is_empty())
					Xk1s.intersection_assign(distPolys.get(cp_i).get());
				    else if(!Xk1s.is_empty()) // now the other polyhedron must be empty
					Xk1s = new C_Polyhedron(distPolys.get(cp_i).get().generators());
				}
			}
		} else { // deadlock
		        Xk1s = Xk[s].get(); //new C_Polyhedron(new Generator_System()); // empty
		}

		// ------------------------------------------------------------------------------
		// ADD STATE REWARDS
		Xk1s = PPLSupport.add_rewards(Xk1s, s, Integer.MIN_VALUE, rewards, null);

		// ------------------------------------------------------------------------------
		// ROUNDING (if required)
		if (rounding) Xk1s = round(Xk1s.generators(), baseline_accuracy, biggest_reward, cut);

		// ------------------------------------------------------------------------------
		// CLEAN UP: UNION WITH PREVIOUS RESULT OR CUT, MINIMIZE REPRESENTATION, DIMENSIONALITY

		// union with previous result (after rounding)
		if (rounding && union_with_previous) Xk1s.upper_bound_assign(Xk[s].get());
		// cut everything but the negative orthant bounded by -M
		if (cut) PPLSupport.cutBox(Xk1s, M);

		// minimize representation
		Xk1s = new C_Polyhedron(Xk1s.minimized_generators());

		// add zero dimensions if minimization deleted them
		if (Xk1s.space_dimension() != n)
			Xk1s.add_space_dimensions_and_project(n - Xk1s.space_dimension());

		return new Pareto(Xk1s);
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
			if (statesList.size() > i)
				s += i + "(P-" + stateOwners.get(i) + " " + statesList.get(i) + "): ";
			else
				s += i + "(P-" + stateOwners.get(i) + "): ";
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
}
