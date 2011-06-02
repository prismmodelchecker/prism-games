//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford, formerly University of Birmingham)
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

import java.util.HashSet;
import java.util.Set;

import parser.*;
import parser.visitor.*;
import prism.PrismLangException;

public class ExpressionCoalnOp extends Expression
{
	Set coalition = null;
	Expression expression = null;
	// Note: this "old-style" filter is just for display purposes
	// The parser creates an (invisible) new-style filter around this expression
	Filter filter = null;

	// Constructors

	public ExpressionCoalnOp()
	{
	}

	public ExpressionCoalnOp(Set c, Expression e)
	{

		coalition = c;
		expression = e;
	}

	// Set methods

	public void setCoalition(Set c)
	{
		this.coalition = c;
	}
	
	public void setExpression(Expression e)
	{
		expression = e;
	}

	public void setFilter(Filter f)
	{
		filter = f;
	}

	// Get methods

	public Expression getExpression()
	{
		return expression;
	}
	
	public Set getCoalition()
	{
		return coalition;
	}

	public Filter getFilter()
	{
		return filter;
	}

	// Methods required for Expression:

	/**
	 * Is this expression constant?
	 */
	public boolean isConstant()
	{
		return false;
	}

	/**
	 * Evaluate this expression, return result.
	 * Note: assumes that type checking has been done already.
	 */
	public Object evaluate(EvaluateContext ec) throws PrismLangException
	{
		throw new PrismLangException("Cannot evaluate a P operator without a model");
	}

	/**
	  * Get "name" of the result of this expression (used for y-axis of any graphs plotted)
	  */
	public String getResultName()
	{
		return "Result for coalition " + coalition.toString();
	}

	@Override
	public boolean returnsSingleValue()
	{
		return false;
	}

	// Methods required for ASTElement:

	/**
	 * Visitor method.
	 */
	public Object accept(ASTVisitor v) throws PrismLangException
	{
		return v.visit(this);
	}

	/**
	 * Convert to string.
	 */
	public String toString()
	{
		return "<<" + coalition + ">> " + expression.toString();
	}

	/**
	 * Perform a deep copy.
	 */
	public Expression deepCopy()
	{
		ExpressionCoalnOp expr = new ExpressionCoalnOp();
		expr.setCoalition(new HashSet(coalition));
		expr.setExpression(expression == null ? null : expression.deepCopy());
		expr.setFilter(filter == null ? null : (Filter) filter.deepCopy());
		expr.setType(type);
		expr.setPosition(this);
		return expr;
	}
}

//------------------------------------------------------------------------------
