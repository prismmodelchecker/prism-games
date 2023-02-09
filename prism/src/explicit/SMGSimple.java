//==============================================================================
//	
//	Copyright (c) 2010-
//	Authors:
//	* Dave Parker <d.a.parker@cs.bham.ac.uk> (University of Birmingham)
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

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.fraction.BigFraction;

import explicit.rewards.MDPRewards;
import explicit.rewards.SMGRewards;
import explicit.rewards.STPGRewards;
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
import prism.PlayerInfo;
import prism.PlayerInfoOwner;
import prism.PrismException;

/**
 * Simple explicit-state representation of a (turn-based) stochastic multi-player game (SMG).
 */
public class SMGSimple extends MDPSimple implements SMG
{
	/**
	 * Which player owns each state
	 */
	protected StateOwnersSimple stateOwners;
	
	/**
	 * Player + coalition information
	 */
	protected PlayerInfo playerInfo;
	
	// Constructors

	/**
	 * Constructor: empty SMG.
	 */
	public SMGSimple()
	{
		super();
		stateOwners = new StateOwnersSimple();
		playerInfo = new PlayerInfo();
	}

	/**
	 * Constructor: new SMG with fixed number of states.
	 */
	public SMGSimple(int numStates)
	{
		super(numStates);
		stateOwners = new StateOwnersSimple(numStates);
		playerInfo = new PlayerInfo();
	}

	/**
	 * Copy constructor
	 */
	public SMGSimple(SMGSimple smg)
	{
		super(smg);
		stateOwners = new StateOwnersSimple(smg.stateOwners);
		playerInfo = new PlayerInfo(smg.playerInfo);
	}
	
	/**
	 * Construct an SMG from an existing one and a state index permutation,
	 * i.e. in which state index i becomes index permut[i].
	 * Player and coalition info is also copied across.
	 */
	public SMGSimple(SMGSimple smg, int permut[])
	{
		super(smg, permut);
		stateOwners = new StateOwnersSimple(smg.stateOwners, permut);
		playerInfo = new PlayerInfo(smg.playerInfo);
	}

	// Mutators

	@Override
	public void clearState(int s)
	{
		super.clearState(s);
		stateOwners.clearState(s);
	}

	@Override
	public void addStates(int numToAdd)
	{
		super.addStates(numToAdd);
		// Assume all player 1
		for (int i = 0; i < numToAdd; i++) {
			stateOwners.addState(0);
		}
	}

	/**
	 * Add a new (player {@code p}) state and return its index.
	 * @param p Player who owns the new state (0-indexed)
	 */
	public int addState(int p)
	{
		int s = super.addState();
		stateOwners.setPlayer(s, p);
		return s;
	}

	/**
	 * Set the player that owns state {@code s} to {@code p}.
	 * @param s State to be modified (0-indexed)
	 * @param p Player who owns the state (0-indexed)
	 */
	public void setPlayer(int s, int p)
	{
		stateOwners.setPlayer(s, p);
	}

	/**
	 * Copy the player info from another model
	 */
	public void copyPlayerInfo(PlayerInfoOwner model)
	{
		playerInfo = new PlayerInfo(model.getPlayerInfo());
	}
	
	// Accessors (for Model)
	
	@Override
	public void checkForDeadlocks(BitSet except) throws PrismException
	{
		for (int i = 0; i < numStates; i++) {
			if (trans.get(i).isEmpty() && (except == null || !except.get(i)))
				throw new PrismException("Game has a deadlock in state " + i + (statesList == null ? "" : ": " + statesList.get(i)));
		}
	}
	
	// Accessors (for STPG)
	
	@Override
	public int getPlayer(int s)
	{
		return playerInfo.getPlayer(stateOwners.getPlayer(s));
	}
	
	@Override
	public void prob0step(BitSet subset, BitSet u, boolean forall1, boolean forall2, BitSet result)
	{
		int i;
		boolean b1, b2;
		boolean forall = false;

		for (i = 0; i < numStates; i++) {
			if (subset.get(i)) {
				forall = (getPlayer(i) == 0) ? forall1 : forall2;
				b1 = forall; // there exists or for all
				for (Distribution distr : trans.get(i)) {
					b2 = distr.containsOneOf(u);
					if (forall) {
						if (!b2) {
							b1 = false;
							continue;
						}
					} else {
						if (b2) {
							b1 = true;
							continue;
						}
					}
				}
				result.set(i, b1);
			}
		}
	}
	
	@Override
	public void prob1step(BitSet subset, BitSet u, BitSet v, boolean forall1, boolean forall2, BitSet result)
	{
		int i;
		boolean b1, b2;
		boolean forall = false;

		for (i = 0; i < numStates; i++) {
			if (subset.get(i)) {
				forall = (getPlayer(i) == 0) ? forall1 : forall2;
				b1 = forall; // there exists or for all
				for (Distribution distr : trans.get(i)) {
					b2 = distr.containsOneOf(v) && distr.isSubsetOf(u);
					if (forall) {
						if (!b2) {
							b1 = false;
							continue;
						}
					} else {
						if (b2) {
							b1 = true;
							continue;
						}
					}
				}
				result.set(i, b1);
			}
		}
	}

	@Override
	public void mvMultMinMax(double vect[], boolean min1, boolean min2, double result[], BitSet subset, boolean complement, int adv[])
	{
		int s;
		boolean min = false;
		// Loop depends on subset/complement arguments
		if (subset == null) {
			for (s = 0; s < numStates; s++) {
				min = (getPlayer(s) == 0) ? min1 : min2;
				result[s] = mvMultMinMaxSingle(s, vect, min, adv);
			}
		} else if (complement) {
			for (s = subset.nextClearBit(0); s < numStates; s = subset.nextClearBit(s + 1)) {
				min = (getPlayer(s) == 0) ? min1 : min2;
				result[s] = mvMultMinMaxSingle(s, vect, min, adv);
			}
		} else {
			for (s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s + 1)) {
				min = (getPlayer(s) == 0) ? min1 : min2;
				result[s] = mvMultMinMaxSingle(s, vect, min, adv);
			}
		}
	}

	@Override
	public double mvMultMinMaxSingle(int s, double vect[], boolean min1, boolean min2)
	{
		boolean min = (getPlayer(s) == 0) ? min1 : min2;
		return mvMultMinMaxSingle(s, vect, min, null);
	}

	@Override
	public List<Integer> mvMultMinMaxSingleChoices(int s, double vect[], boolean min1, boolean min2, double val)
	{
		boolean min = (getPlayer(s) == 0) ? min1 : min2;
		return mvMultMinMaxSingleChoices(s, vect, min, val);
	}

	@Override
	public double mvMultGSMinMax(double vect[], boolean min1, boolean min2, BitSet subset, boolean complement, boolean absolute)
	{
		int s;
		double d, diff, maxDiff = 0.0;
		// Loop depends on subset/complement arguments
		if (subset == null) {
			for (s = 0; s < numStates; s++) {
				d = mvMultJacMinMaxSingle(s, vect, min1, min2);
				diff = absolute ? (Math.abs(d - vect[s])) : (Math.abs(d - vect[s]) / d);
				maxDiff = diff > maxDiff ? diff : maxDiff;
				vect[s] = d;
			}
		} else if (complement) {
			for (s = subset.nextClearBit(0); s < numStates; s = subset.nextClearBit(s + 1)) {
				d = mvMultJacMinMaxSingle(s, vect, min1, min2);
				diff = absolute ? (Math.abs(d - vect[s])) : (Math.abs(d - vect[s]) / d);
				maxDiff = diff > maxDiff ? diff : maxDiff;
				vect[s] = d;
			}
		} else {
			for (s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s + 1)) {
				d = mvMultJacMinMaxSingle(s, vect, min1, min2);
				diff = absolute ? (Math.abs(d - vect[s])) : (Math.abs(d - vect[s]) / d);
				maxDiff = diff > maxDiff ? diff : maxDiff;
				vect[s] = d;
			}
		}
		return maxDiff;
	}

	@Override
	public double mvMultJacMinMaxSingle(int s, double vect[], boolean min1, boolean min2)
	{
		boolean min = (getPlayer(s) == 0) ? min1 : min2;
		return mvMultJacMinMaxSingle(s, vect, min, null);
	}

	@Override
	public void mvMultRewMinMax(double vect[], STPGRewards rewards, boolean min1, boolean min2, double result[], BitSet subset, boolean complement, int adv[])
	{
		int s;
		boolean min = false;
		MDPRewards mdpRewards = rewards.buildMDPRewards();
		// Loop depends on subset/complement arguments
		if (subset == null) {
			for (s = 0; s < numStates; s++) {
				min = (getPlayer(s) == 0) ? min1 : min2;
				result[s] = mvMultRewMinMaxSingle(s, vect, mdpRewards, min, adv, 1.0);
			}
		} else if (complement) {
			for (s = subset.nextClearBit(0); s < numStates; s = subset.nextClearBit(s + 1)) {
				min = (getPlayer(s) == 0) ? min1 : min2;
				result[s] = mvMultRewMinMaxSingle(s, vect, mdpRewards, min, adv, 1.0);
			}
		} else {
			for (s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s + 1)) {
				min = (getPlayer(s) == 0) ? min1 : min2;
				result[s] = mvMultRewMinMaxSingle(s, vect, mdpRewards, min, adv, 1.0);
			}
		}
	}

	@Override
	public void mvMultRewMinMax(double vect[], STPGRewards rewards, boolean min1, boolean min2, double result[], BitSet subset, boolean complement, int adv[],
			double disc)
	{
		int s;
		boolean min = false;
		MDPRewards mdpRewards = rewards.buildMDPRewards();
		// Loop depends on subset/complement arguments
		if (subset == null) {
			for (s = 0; s < numStates; s++) {
				min = (getPlayer(s) == 0) ? min1 : min2;
				result[s] = mvMultRewMinMaxSingle(s, vect, mdpRewards, min, adv, disc);
			}
		} else if (complement) {
			for (s = subset.nextClearBit(0); s < numStates; s = subset.nextClearBit(s + 1)) {
				min = (getPlayer(s) == 0) ? min1 : min2;
				result[s] = mvMultRewMinMaxSingle(s, vect, mdpRewards, min, adv, disc);
			}
		} else {
			for (s = subset.nextSetBit(0); s >= 0; s = subset.nextSetBit(s + 1)) {
				min = (getPlayer(s) == 0) ? min1 : min2;
				//System.out.printf("s: %s, min1: %s, min2: %s, min: %s, player: %d\n", s, min1, min2, min, getPlayer(s));
				result[s] = mvMultRewMinMaxSingle(s, vect, mdpRewards, min, adv, disc);
			}
		}
	}

	@Override
	public double mvMultRewMinMaxSingle(int s, double vect[], STPGRewards rewards, boolean min1, boolean min2, int adv[])
	{
		MDPRewards mdpRewards = rewards.buildMDPRewards();
		boolean min = (getPlayer(s) == 0) ? min1 : min2;
		return mvMultRewMinMaxSingle(s, vect, mdpRewards, min, adv);
	}

	@Override
	public List<Integer> mvMultRewMinMaxSingleChoices(int s, double vect[], STPGRewards rewards, boolean min1, boolean min2, double val)
	{
		MDPRewards mdpRewards = rewards.buildMDPRewards();
		boolean min = (getPlayer(s) == 0) ? min1 : min2;
		return mvMultRewMinMaxSingleChoices(s, vect, mdpRewards, min, val);
	}

	/**
	 * Allows discounting
	 * @param s
	 * @param vect
	 * @param mdpRewards
	 * @param min
	 * @param adv
	 * @param disc
	 * @return
	 */
	public double mvMultRewMinMaxSingle(int s, double vect[], MDPRewards mdpRewards, boolean min, int adv[], double disc)
	{
		int j, k, advCh = -1;
		double d, prob, minmax;
		boolean first;
		List<Distribution> step;

		minmax = 0;
		first = true;
		j = -1;
		step = trans.get(s);
		for (Distribution distr : step) {
			j++;
			// Compute sum for this distribution
			d = mdpRewards.getTransitionReward(s, j);

			for (Map.Entry<Integer, Double> e : distr) {
				k = (Integer) e.getKey();
				prob = (Double) e.getValue();
				d += prob * vect[k] * disc;
			}

			// Check whether we have exceeded min/max so far
			if (first || (min && d < minmax) || (!min && d > minmax)) {
				minmax = d;
				// If adversary generation is enabled, remember optimal choice
				if (adv != null) {
					advCh = j;
				}
			}
			first = false;
		}
		// If adversary generation is enabled, store optimal choice
		if (adv != null & !first) {
			// Only remember strictly better choices (required for max)
			if (adv[s] == -1 || (min && minmax < vect[s]) || (!min && minmax > vect[s]) || this instanceof STPG) {
				adv[s] = advCh;
			}
		}

		// Add state reward (doesn't affect min/max)
		minmax += mdpRewards.getStateReward(s);

		return minmax;
	}

	// Accessors (for PlayerInfoOwner)

	@Override
	public PlayerInfo getPlayerInfo()
	{
		return playerInfo;
	}
	
	// Accessors (for SMG)

	@Override
	public void reachpositivestep(BitSet u, boolean forall1, boolean forall2, BitSet result)
	{
		int i;
		boolean forall = false;
		boolean first;
		Set<Integer> u1;
		
		for (i = 0; i < numStates; i++) {
			if (u.get(i)) {
				forall = (getPlayer(i) == 0) ? forall1 : forall2;
				u1 = null; // reach in one step
				first = true;
				for (Distribution distr : trans.get(i)) {
					if(first) { 
						u1 = new HashSet<Integer>(distr.getSupport()); // put all successors in reachable states
						first = false;
					} else if (!first && forall) {
						u1.retainAll(distr.getSupport()); // intersect
					} else if (!first & !forall) {
						u1.addAll(distr.getSupport()); // union
					}
				}
				for(int r : u1)
					result.set(r, true);
			}
		}
	}

    
	@Override
	public void subtreeStep(BitSet u, int closedPlayer, BitSet result)
	{
		for (int i = 0; i < numStates; i++) {
			// go only through states in subtree so far,
			// and only extend subtree if closed for that player,
			// or if the state has only one choice that is enabled
			boolean jump = (getNumChoices(i) == 1) && getPlayer(i) != closedPlayer;
			if (u.get(i) && (getPlayer(i) == closedPlayer || jump)) {
				for (Distribution distr : trans.get(i)) {
					for (int r : distr.getSupport()) { // add all successors (no matter which player)
						result.set(r, true);
					}
				}
			}
		}
	}
	
	@Override
	public Pareto[] pMultiObjective(Pareto[] Xk, List<SMGRewards> rewards, boolean gaussSeidel, long baseline_accuracy, double[] biggest_reward,
			List<Pareto>[] stochasticStates, boolean rounding, boolean union_with_previous, boolean cut, long M) throws PrismException
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
			if (getPlayer(s) == 0) {
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
				s += i + "(P-" + (stateOwners.getPlayer(i)+1) + " " + statesList.get(i) + "): ";
			else
				s += i + "(P-" + (stateOwners.getPlayer(i)+1) + "): ";
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
