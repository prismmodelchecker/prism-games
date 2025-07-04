//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <david.parker@cs.ox.ac.uk> (University of Oxford)
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

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;

import explicit.graphviz.StateOwnerDecorator;
import explicit.rewards.STPGRewards;
import prism.ModelType;
import prism.PlayerInfoOwner;
import prism.PrismException;
import prism.PrismLog;
import prism.PrismUtils;

/**
 * Interface for classes that provide (read) access to an explicit-state stochastic two-player game (STPG).
 * <br><br>
 * These are turn-based STPGs, i.e. at most one player controls each state.
 * Probabilistic states do not need to be stored explicitly; instead, like in an MDP,
 * players have several 'choices', each of which is a probability distribution over successor states.
 */
public interface STPG<Value> extends MDP<Value>, TurnBasedGame
{
	// Accessors (for Model) - default implementations
	
	@Override
	default ModelType getModelType()
	{
		return ModelType.STPG;
	}

	@Override
	default void exportToPrismLanguage(final String filename, int precision) throws PrismException
	{
		throw new UnsupportedOperationException();
	}
	
	// Accessors

	/**
	 * Perform a single step of precomputation algorithm Prob0, i.e., for states i in {@code subset},
	 * set bit i of {@code result} iff, for all/some player 1 choices, for all/some player 2 choices,
	 * there is a transition to a state in {@code u}.
	 * Quantification over player 1/2 choices is determined by {@code forall1}, {@code forall2}.
	 * @param subset Only compute for these states
	 * @param u Set of states {@code u}
	 * @param forall1 For-all or there-exists for player 1 (true=for-all, false=there-exists)
	 * @param forall2 For-all or there-exists for player 2 (true=for-all, false=there-exists)
	 * @param result Store results here
	 */
	public void prob0step(BitSet subset, BitSet u, boolean forall1, boolean forall2, BitSet result);

	/**
	 * Perform a single step of precomputation algorithm Prob1, i.e., for states i in {@code subset},
	 * set bit i of {@code result} iff, for all/some player 1 choices, for all/some player 2 choices,
	 * there is a transition to a state in {@code v} and all transitions go to states in {@code u}.
	 * Quantification over player 1/2 choices is determined by {@code forall1}, {@code forall2}.
	 * @param subset Only compute for these states
	 * @param u Set of states {@code u}
	 * @param v Set of states {@code v}
	 * @param forall1 For-all or there-exists for player 1 (true=for-all, false=there-exists)
	 * @param forall2 For-all or there-exists for player 2 (true=for-all, false=there-exists)
	 * @param result Store results here
	 */
	public void prob1step(BitSet subset, BitSet u, BitSet v, boolean forall1, boolean forall2, BitSet result);

	/**
	 * Do a matrix-vector multiplication followed by two min/max ops, i.e. one step of value iteration,
	 * i.e. for all s: result[s] = min/max_{k1,k2} { sum_j P_{k1,k2}(s,j)*vect[j] }
	 * @param vect Vector to multiply by
	 * @param min1 Min or max for player 1 (true=min, false=max)
	 * @param min2 Min or max for player 2 (true=min, false=max)
	 * @param result Vector to store result in
	 * @param subset Only do multiplication for these rows (ignored if null)
	 * @param complement If true, {@code subset} is taken to be its complement (ignored if {@code subset} is null)
	 * @param adv Storage for adversary choice indices (ignored if null)
	 */
	public void mvMultMinMax(double vect[], boolean min1, boolean min2, double result[], BitSet subset, boolean complement, int adv[]);

	/**
	 * Do a single row of matrix-vector multiplication followed by min/max,
	 * i.e. return min/max_{k1,k2} { sum_j P_{k1,k2}(s,j)*vect[j] }
	 * @param s Row index
	 * @param vect Vector to multiply by
	 * @param min1 Min or max for player 1 (true=min, false=max)
	 * @param min2 Min or max for player 2 (true=min, false=max)
	 */
	public double mvMultMinMaxSingle(int s, double vect[], boolean min1, boolean min2);

	/**
	 * Determine which choices result in min/max after a single row of matrix-vector multiplication.
	 * @param s Row index
	 * @param vect Vector to multiply by
	 * @param min1 Min or max for player 1 (true=min, false=max)
	 * @param min2 Min or max for player 2 (true=min, false=max)
	 * @param val Min or max value to match
	 */
	public List<Integer> mvMultMinMaxSingleChoices(int s, double vect[], boolean min1, boolean min2, double val);

	/**
	 * Do a Gauss-Seidel-style matrix-vector multiplication followed by min/max.
	 * i.e. for all s: vect[s] = min/max_{k1,k2} { (sum_{j!=s} P_{k1,k2}(s,j)*vect[j]) / P_{k1,k2}(s,s) }
	 * and store new values directly in {@code vect} as computed.
	 * The maximum (absolute/relative) difference between old/new
	 * elements of {@code vect} is also returned.
	 * @param vect Vector to multiply by (and store the result in)
	 * @param min1 Min or max for player 1 (true=min, false=max)
	 * @param min2 Min or max for player 2 (true=min, false=max)
	 * @param subset Only do multiplication for these rows (ignored if null)
	 * @param complement If true, {@code subset} is taken to be its complement (ignored if {@code subset} is null)
	 * @param absolute If true, compute absolute, rather than relative, difference
	 * @return The maximum difference between old/new elements of {@code vect}
	 */
	public double mvMultGSMinMax(double vect[], boolean min1, boolean min2, BitSet subset, boolean complement, boolean absolute, int adv[]);

	/**
	 * Do a single row of Jacobi-style matrix-vector multiplication followed by min/max.
	 * i.e. return min/max_{k1,k2} { (sum_{j!=s} P_{k1,k2}(s,j)*vect[j]) / P_{k1,k2}(s,s) }
	 * @param s Row index
	 * @param vect Vector to multiply by
	 * @param min1 Min or max for player 1 (true=min, false=max)
	 * @param min2 Min or max for player 2 (true=min, false=max)
	 */
	public double mvMultJacMinMaxSingle(int s, double vect[], boolean min1, boolean min2, int[] adv);

	/**
	 * Do a matrix-vector multiplication and sum of action reward followed by min/max, i.e. one step of value iteration.
	 * i.e. for all s: result[s] = min/max_{k1,k2} { rew(s) + sum_j P_{k1,k2}(s,j)*vect[j] }
	 * @param vect Vector to multiply by
	 * @param rewards The rewards
	 * @param min1 Min or max for player 1 (true=min, false=max)
	 * @param min2 Min or max for player 2 (true=min, false=max)
	 * @param result Vector to store result in
	 * @param subset Only do multiplication for these rows (ignored if null)
	 * @param complement If true, {@code subset} is taken to be its complement (ignored if {@code subset} is null)
	 * @param adv Storage for adversary choice indices (ignored if null)
	 */
	public void mvMultRewMinMax(double vect[], STPGRewards<Double> rewards, boolean min1, boolean min2, double result[], BitSet subset, boolean complement, int adv[]);

	/**
	 * Do a single row of matrix-vector multiplication and sum of action reward followed by min/max.
	 * i.e. return min/max_{k1,k2} { rew(s) + sum_j P_{k1,k2}(s,j)*vect[j] }
	 * @param s Row index
	 * @param vect Vector to multiply by
	 * @param rewards The rewards
	 * @param min1 Min or max for player 1 (true=min, false=max)
	 * @param min2 Min or max for player 2 (true=min, false=max)
	 * @param adv Storage for adversary choice indices (ignored if null)
	 */
	public double mvMultRewMinMaxSingle(int s, double vect[], STPGRewards<Double> rewards, boolean min1, boolean min2, int adv[]);

	/**
	 * Determine which choices result in min/max after a single row of matrix-vector multiplication and sum of action reward.
	 * @param s Row index
	 * @param vect Vector to multiply by
	 * @param rewards The rewards
	 * @param min1 Min or max for player 1 (true=min, false=max)
	 * @param min2 Min or max for player 2 (true=min, false=max)
	 * @param val Min or max value to match
	 */
	public List<Integer> mvMultRewMinMaxSingleChoices(int s, double vect[], STPGRewards<Double> rewards, boolean min1, boolean min2, double val);

	/**
	 * Do a single row of (discounted) matrix-vector multiplication and sum of action reward followed by min/max.
	 * i.e. return min/max_{k1,k2} { rew(s) + sum_j P_{k1,k2}(s,j)*vect[j] }
	 * @param vect Vector to multiply by
	 * @param rewards The rewards
	 * @param min1 Min or max for player 1 (true=min, false=max)
	 * @param min2 Min or max for player 2 (true=min, false=max)
	 * @param adv Storage for adversary choice indices (ignored if null)
	 * @param disc Discount factor
	 */
	void mvMultRewMinMax(double[] vect, STPGRewards<Double> rewards, boolean min1, boolean min2, double[] result, BitSet subset, boolean complement, int[] adv, double disc);

	/**
	 * Checks  whether all successors of action c in state s are in a given set
	 * @param s state
	 * @param c choice
	 * @param set target set
	 * @return true if all successors are, false otherwise
	 */
	public boolean allSuccessorsInSet(int s, int c, BitSet set);
}
