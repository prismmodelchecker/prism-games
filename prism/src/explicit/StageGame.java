//==============================================================================
//	
//	Copyright (c) 2026-
//	Authors:
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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * A one-shot (normal-form) game solved at a state of a CSG: coalitions, their actions,
 * one payoff per coalition for every joint action, and which coalitions are "free".
 * <p>
 * Joint actions are indexed in row-major order (the last coalition varies fastest);
 * see {@link #jointIndex(int[])} and {@link #jointActions(int)}.
 * A free coalition has no incentive constraints (e.g. it has already reached, or can no
 * longer reach, its objective): its strategy is chosen together with the others'
 * to optimise the criterion.
 *
 * @param <V> payoff type (only Double is currently supported by the solvers)
 */
public class StageGame<V>
{
	private final int numCoalitions;
	private final int[] numActions;
	private final int[] strides;
	private final int numJoint;
	private final List<List<V>> payoffs;
	private final boolean[] free;

	public StageGame(int[] numActions)
	{
		if (numActions.length < 1)
			throw new IllegalArgumentException("A stage game needs at least one coalition");
		this.numCoalitions = numActions.length;
		this.numActions = numActions.clone();
		this.strides = new int[numCoalitions];
		int s = 1;
		for (int c = numCoalitions - 1; c >= 0; c--) {
			if (numActions[c] < 1)
				throw new IllegalArgumentException("Coalition " + c + " has no actions");
			strides[c] = s;
			s *= numActions[c];
		}
		this.numJoint = s;
		this.payoffs = new ArrayList<>(numCoalitions);
		for (int c = 0; c < numCoalitions; c++)
			payoffs.add(new ArrayList<>(java.util.Collections.nCopies(numJoint, (V) null)));
		this.free = new boolean[numCoalitions];
	}

	public int getNumCoalitions()
	{
		return numCoalitions;
	}

	public int getNumActions(int c)
	{
		return numActions[c];
	}

	public int getNumJointActions()
	{
		return numJoint;
	}

	public int jointIndex(int[] actions)
	{
		int j = 0;
		for (int c = 0; c < numCoalitions; c++)
			j += actions[c] * strides[c];
		return j;
	}

	public int[] jointActions(int joint)
	{
		int[] a = new int[numCoalitions];
		for (int c = 0; c < numCoalitions; c++) {
			a[c] = joint / strides[c];
			joint %= strides[c];
		}
		return a;
	}

	/** Joint index obtained from {@code joint} by replacing coalition c's action with {@code action}. */
	public int deviate(int joint, int c, int action)
	{
		int current = (joint / strides[c]) % numActions[c];
		return joint + (action - current) * strides[c];
	}

	public void setPayoff(int c, int joint, V value)
	{
		payoffs.get(c).set(joint, value);
	}

	public void setPayoff(int c, int[] actions, V value)
	{
		setPayoff(c, jointIndex(actions), value);
	}

	public V getPayoff(int c, int joint)
	{
		return payoffs.get(c).get(joint);
	}

	public void setFree(int c, boolean isFree)
	{
		free[c] = isFree;
	}

	public boolean isFree(int c)
	{
		return free[c];
	}

	/** Throws if some payoff has not been set. */
	public void checkComplete()
	{
		for (int c = 0; c < numCoalitions; c++)
			for (int j = 0; j < numJoint; j++)
				if (payoffs.get(c).get(j) == null)
					throw new IllegalStateException("Missing payoff for coalition " + c + " at joint action " + Arrays.toString(jointActions(j)));
	}

	@Override
	public String toString()
	{
		StringBuilder sb = new StringBuilder("StageGame" + Arrays.toString(numActions));
		for (int j = 0; j < numJoint; j++) {
			sb.append("\n  ").append(Arrays.toString(jointActions(j))).append(" -> [");
			for (int c = 0; c < numCoalitions; c++)
				sb.append(c > 0 ? ", " : "").append(payoffs.get(c).get(j));
			sb.append("]");
		}
		return sb.toString();
	}
}
