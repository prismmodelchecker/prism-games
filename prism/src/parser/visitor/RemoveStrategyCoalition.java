//==============================================================================
//
//	Copyright (c) 2026-
//	Authors:
//	* Dave Parker <d.a.parker@cs.bham.ac.uk> (University of Birmingham)
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

package parser.visitor;

import java.util.BitSet;
import java.util.List;

import parser.ast.Coalition;
import parser.ast.Expression;
import parser.ast.ExpressionProb;
import parser.ast.ExpressionReward;
import parser.ast.ExpressionSS;
import parser.ast.ExpressionStrategy;
import parser.ast.ExpressionTemporal;
import parser.ast.RelOp;
import prism.PrismLangException;

/**
 * Adapt a game property for the model induced by a strategy of some of its players
 * (a coalition C), in which only the other players' choices remain. The property then asks
 * what the strategy guarantees, whatever the other players do.
 * <ul>
 * <li>{@code <<C>>} becomes the empty coalition {@code <<>>} (against all other players)
 * if the induced model is still a game, or is removed if it is an MDP, in which case
 * min/max for C become max/min for the remaining player. Strategy quantifiers for a
 * coalition other than C are not supported.</li>
 * <li>For STPGs, whose properties have no coalitions but give both players' optimisation,
 * e.g. {@code Rminmax=?} (player 1 minimises, player 2 maximises),
 * only the remaining player's part is kept, for the induced MDP.</li>
 * </ul>
 * Bounds (e.g. {@code P>=0.5}) are unchanged, since for MDPs they hold for all strategies.
 */
public class RemoveStrategyCoalition extends ASTTraverseModify
{
	/** Players of the strategy (i.e., the coalition C) */
	private BitSet stratPlayers;
	/** Names of the players of the original model */
	private List<String> playerNames;
	/** Is the induced model an MDP (rather than a game)? */
	private boolean inducedMDP;
	/** Are we inside a strategy quantifier {@code <<C>>} (that is being removed)? */
	private boolean inCoalition = false;

	/**
	 * Constructor.
	 * @param stratPlayers Players of the strategy (indices, from 0)
	 * @param playerNames Names of the players of the original model
	 * @param inducedMDP Is the induced model an MDP (rather than a game)?
	 */
	public RemoveStrategyCoalition(BitSet stratPlayers, List<String> playerNames, boolean inducedMDP)
	{
		this.stratPlayers = stratPlayers;
		this.playerNames = playerNames;
		this.inducedMDP = inducedMDP;
	}

	@Override
	public Object visit(ExpressionStrategy e) throws PrismLangException
	{
		if (!e.isThereExists() || e.getNumCoalitions() != 1 || e.getNumOperands() != 1) {
			throw new PrismLangException("Only strategy quantifiers <<C>> with a single operand are supported here", e);
		}
		if (!isStrategyCoalition(e.getCoalition())) {
			throw new PrismLangException("Coalition " + e.getCoalition() + " does not match the players of the strategy", e);
		}
		if (!inducedMDP) {
			e.setCoalition(new Coalition());
			return e;
		}
		inCoalition = true;
		try {
			return e.getOperand(0).accept(this);
		} finally {
			inCoalition = false;
		}
	}

	@Override
	public void visitPost(ExpressionProb e) throws PrismLangException
	{
		e.setRelOp(adaptRelOp(e.getRelOp()));
	}

	@Override
	public void visitPost(ExpressionReward e) throws PrismLangException
	{
		if (inducedMDP && containsFc(e.getExpression())) {
			throw new PrismLangException("The Fc reward operator is not supported for MDPs", e);
		}
		e.setRelOp(adaptRelOp(e.getRelOp()));
	}

	@Override
	public void visitPost(ExpressionSS e) throws PrismLangException
	{
		e.setRelOp(adaptRelOp(e.getRelOp()));
	}

	/**
	 * Does an expression use the Fc reward operator (which is only supported for games)?
	 */
	public static boolean containsFc(Expression expr) throws PrismLangException
	{
		try {
			expr.accept(new ASTTraverse()
			{
				@Override
				public void visitPost(ExpressionTemporal e) throws PrismLangException
				{
					if (e.getOperator() == ExpressionTemporal.R_Fc) {
						throw new PrismLangException("Fc");
					}
				}
			});
		} catch (PrismLangException e) {
			return true;
		}
		return false;
	}

	/**
	 * Is this coalition the same set of players as the strategy's?
	 */
	private boolean isStrategyCoalition(Coalition coalition)
	{
		int numPlayers = playerNames.size();
		for (int p = 0; p < numPlayers; p++) {
			boolean inCoalition = coalition.isAllPlayers() || coalition.isPlayerIndexInCoalition(p, playerNames);
			if (inCoalition != stratPlayers.get(p)) {
				return false;
			}
		}
		return true;
	}

	/**
	 * Get the optimisation (if any) for the remaining player, in an induced MDP.
	 */
	private RelOp adaptRelOp(RelOp relOp) throws PrismLangException
	{
		if (!inducedMDP || relOp == null) {
			return relOp;
		}
		// Inside <<C>>: the opposite of the coalition's
		if (inCoalition) {
			if (relOp == RelOp.MIN) {
				return RelOp.MAX;
			} else if (relOp == RelOp.MAX) {
				return RelOp.MIN;
			}
			return relOp;
		}
		// Otherwise, for two-player (STPG) operators, the remaining player's part
		boolean remainingIsFirst = !stratPlayers.get(0);
		switch (relOp) {
		case MINMIN:
			return RelOp.MIN;
		case MAXMAX:
			return RelOp.MAX;
		case MINMAX:
			return remainingIsFirst ? RelOp.MIN : RelOp.MAX;
		case MAXMIN:
			return remainingIsFirst ? RelOp.MAX : RelOp.MIN;
		default:
			return relOp;
		}
	}
}
