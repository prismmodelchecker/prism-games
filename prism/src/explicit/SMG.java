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

import java.util.BitSet;
import java.util.List;

import explicit.rewards.SMGRewards;
import prism.ModelType;
import prism.PlayerInfo;
import prism.PlayerInfoOwner;
import prism.PrismException;

/**
 * Interface for classes that provide (read) access to an explicit-state SMG.
 * <br><br>
 * An SMG is a (turn-based) stochastic multi-player game,
 * essentially an MDP with information about which player owns each state.
 * This extends {@link STPG}, which already stores state owners,
 * and adds further info about players, stored as a {@link PlayerInfo} object.
 * <br><br>
 * Although classes implementing this interface only guarantee read-only
 * access to the model info, they can be (temporarily) turned into a 2-player
 * game using {@link #setCoalition(parser.ast.Coalition)}.
 */
public interface SMG extends STPG, PlayerInfoOwner
{
	// Accessors (for Model) - default implementations

	@Override
	default ModelType getModelType()
	{
		return ModelType.SMG;
	}

	// Accessors
	
	/**
	 * @param u
	 * @param forall1
	 * @param forall2
	 * @param result
	 */
	public void reachpositivestep(BitSet u, boolean forall1, boolean forall2, BitSet result);
	
	/**
	 * @param u The subtree so far
	 * @param closedPlayer Player for which subtree is closed
	 * @param result The subtree after extending
	 */
	public void subtreeStep(BitSet u, int closedPlayer, BitSet result);

	/**
	 * Take X^k and apply F(X^k)(s) for each state (cf. MFCS'13 and QEST'13)
	 * @param gaussSeidel Gauss Seidel update allowed
	 * @param rounding rounding enabled
	 * @param union_with_previous take union with previous Pareto set
	 * @param cut cut off everything that is strictly above the negative orthant (used for energy objectives)
	 * @param M maximum bound on Pareto sets (quantity is positive)
	 */
	public Pareto[] pMultiObjective(Pareto[] Xk, List<SMGRewards> rewards, boolean gaussSeidel, long baseline_accuracy, double[] biggest_reward,
			List<Pareto>[] stochasticStates, boolean rounding, boolean union_with_previous, boolean cut, long M) throws PrismException;
}
