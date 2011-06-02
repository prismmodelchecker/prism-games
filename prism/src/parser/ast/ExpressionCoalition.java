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

public class ExpressionCoalition extends ExpressionProb
{
	Set coalition = null;

	// Constructors

	public ExpressionCoalition()
	{
		super();
	}

	public ExpressionCoalition(Set c, Expression e)
	{
		super(((ExpressionProb)e).getExpression(), ((ExpressionProb)e).getRelOp(), ((ExpressionProb)e).getProb());
		coalition = c;
	}

	// Set methods

	public void setCoalition(Set c)
	{
		this.coalition = c;
	}
	
	public void setExpressionProb(ExpressionProb expr)
	{
		setRelOp(expr.relOp);
		setExpression(expr.expression);
		setProb(expr.getProb());
	}


	// Get methods
	
	public Set getCoalition()
	{
		return coalition;
	}

	// Methods required for Expression:

	/**
	  * Get "name" of the result of this expression (used for y-axis of any graphs plotted)
	  */
	public String getResultName()
	{
		return "Result for coalition " + coalition.toString() + " " + super.toString();
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
		return "<<" + coalition + ">> " + super.toString();
	}

	/**
	 * Perform a deep copy.
	 */
	public Expression deepCopy()
	{
		ExpressionCoalition expr = new ExpressionCoalition();
		expr.setCoalition(new HashSet(coalition));
		expr.setProb(super.getProb().deepCopy());
		expr.setRelOp(super.relOp);
		expr.setExpression(expression == null ? null : expression.deepCopy());
		expr.setFilter(filter == null ? null : (Filter) filter.deepCopy());
		expr.setType(type);
		expr.setPosition(this);
		return expr;
	}
}

//------------------------------------------------------------------------------
