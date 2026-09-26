//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Gabriel Santos
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

package parser.ast;

import parser.EvaluateContext;
import parser.visitor.ASTVisitor;
import parser.visitor.DeepCopy;
import prism.PrismLangException;

/**
 * Qualitative (sure/almost/limit) coalition-strategy expression, e.g. {@code sure [ G phi ]}.
 * Appears as the child expression of an {@link ExpressionStrategy} ({@code <<C>>}/{@code [[C]]})
 * node, as an alternative to {@link ExpressionProb}/{@link ExpressionReward}/
 * {@link ExpressionMultiNash}. Strictly zero-sum (one coalition vs. its complement) — no
 * min/max, no equilibrium options; see ExpressionStrategy.EquilibriumType/EquilibriumCriterion,
 * which do not apply here and should be rejected (post-parse) when this is the child operand.
 */

public class ExpressionStrategyQual extends Expression
{
	/** Qualitative winning mode */
	protected QualMode mode;

	public enum QualMode { SURE, ALMOST, LIMIT };

	/** The wrapped path formula, e.g. G phi / F phi / G F phi / F G phi.
	 *  Legality of the shape (safety/reachability/Buchi/co-Buchi only, for now)
	 *  is checked at model-check time, not here. */
	protected Expression expression = null;

	// Constructors

	public ExpressionStrategyQual()
	{
	}

	public ExpressionStrategyQual(QualMode mode, Expression e)
	{
		this.mode = mode;
		expression = e;
	}

	// Set methods

	public void setMode(QualMode mode)
	{
		this.mode = mode;
	}

	public void setExpression(Expression e)
	{
		expression = e;
	}

	// Get methods

	public QualMode getMode()
	{
		return mode;
	}

	public Expression getExpression()
	{
		return expression;
	}

	// Methods required for Expression:

	@Override
	public boolean isConstant()
	{
		return false;
	}

	@Override
	public boolean isProposition()
	{
		return false;
	}

	@Override
	public Object evaluate(EvaluateContext ec) throws PrismLangException
	{
		throw new PrismLangException("Cannot evaluate a qualitative strategy operator without a model");
	}

	@Override
	public boolean returnsSingleValue()
	{
		return false;
	}

	// Methods required for ASTElement:

	@Override
	public Object accept(ASTVisitor v) throws PrismLangException
	{
		return v.visit(this);
	}

	@Override
	public ExpressionStrategyQual deepCopy(DeepCopy copier) throws PrismLangException
	{
		expression = copier.copy(expression);

		return this;
	}

	@Override
	public ExpressionStrategyQual clone()
	{
		return (ExpressionStrategyQual) super.clone();
	}

	// Standard methods

	@Override
	public String toString()
	{
		String s = "";

		switch (mode) {
		case SURE:
			s += "sure";
			break;
		case ALMOST:
			s += "almost";
			break;
		case LIMIT:
			s += "limit";
			break;
		}
		s += " [ " + expression + " ]";

		return s;
	}

	@Override
	public int hashCode()
	{
		final int prime = 31;
		int result = 1;
		result = prime * result + ((mode == null) ? 0 : mode.hashCode());
		result = prime * result + ((expression == null) ? 0 : expression.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj)
	{
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		ExpressionStrategyQual other = (ExpressionStrategyQual) obj;
		if (mode != other.mode)
			return false;
		if (expression == null) {
			if (other.expression != null)
				return false;
		} else if (!expression.equals(other.expression))
			return false;
		return true;
	}
}
