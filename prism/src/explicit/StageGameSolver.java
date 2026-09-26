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

import prism.PrismException;

/**
 * Solver for the one-shot games arising at each state during CSG model checking.
 * Returns one selected equilibrium (or the value of a zero-sum game), not a set.
 * <p>
 * Selection (as in PRISM-games):
 * <ul>
 * <li>SOCIAL_WELFARE: maximise the sum of payoffs, then coalition 0's payoff, then coalition 1's, ...</li>
 * <li>SOCIAL_FAIRNESS: minimise (highest payoff - lowest payoff), then maximise coalition 0's payoff, then coalition 1's, ...</li>
 * </ul>
 * With {@code minimise = true} payoffs are costs: equilibria are computed for the negated game
 * (Lemma 1 of Kwiatkowska et al., FMSD 2021), so SOCIAL_WELFARE yields social-cost equilibria and the
 * tie-breaking minimises coalition 0's cost, then coalition 1's, ... Results are reported in the original scale.
 *
 * @param <V> value type
 */
public interface StageGameSolver<V>
{
	enum Concept
	{
		NASH, CORRELATED, ZERO_SUM
	}

	enum Criterion
	{
		SOCIAL_WELFARE, SOCIAL_FAIRNESS
	}

	String getSolverName();

	/**
	 * Solve a stage game.
	 * For ZERO_SUM the game must have two coalitions; coalition 0's payoffs are used, coalition 0 maximises
	 * (minimises if {@code minimise}), and the criterion is ignored.
	 */
	StageGameResult<V> solve(StageGame<V> game, Concept concept, Criterion criterion, boolean minimise) throws PrismException;
}
