//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford, formerly University of Birmingham)
//	* Aistis Simaitis <aistis.simaitis@cs.ox.ac.uk> (University of Oxford)
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

import parser.EvaluateContext;
import parser.visitor.ASTVisitor;
import prism.PrismLangException;

public class ExpressionPATL extends Expression {

	public static final int PRB = 0; 
	public static final int REW = 1;
	
	private int exprType;
	
	private ExpressionProb exprProb; // expression for probabilistic PATL formula
	private ExpressionReward exprRew; // expression for reward-based PATL formula
	
	private Set<String> coalition;
	
	public ExpressionPATL()
	{}
	
	public ExpressionPATL(Set<String> coal, Expression expr, int type)
	{
		coalition = coal;
		switch(type)
		{
		case PRB:
			this.exprType = type;
			exprProb = (ExpressionProb) expr;
			break;
		case REW:
			this.exprType = type;
			exprRew = (ExpressionReward) expr;
			break;
		}
	}
	
	// getters
	public int getExpressionType()
	{return exprType;}
	
	public ExpressionProb getExpressionProb()
	{return exprProb;}
	
	public ExpressionReward getExpressionRew()
	{return exprRew;}
	
	public Set<String> getCoalition()
	{return coalition;}
	
	// setters
	public void setExpressionType(int type)
	{this.exprType=type;}
	
	public void setExpressionProb(ExpressionProb expr)
	{exprProb = expr;}
	
	public void setExpressionRew(ExpressionReward expr)
	{exprRew = expr;}
	
	public void setCoalition(Set<String> coal)
	{coalition=coal;}
	

	// Overriden methods
	@Override
	public Object evaluate(EvaluateContext ec) throws PrismLangException {
		throw new PrismLangException("Cannot evaluate a <<...>> operator without a model");
	}

	@Override
	public boolean isConstant() {
		return false;
	}

	@Override
	public boolean isProposition()
	{
		return false;
	}
	
	@Override
	public boolean returnsSingleValue() {
		return false;
	}

	@Override
	public Object accept(ASTVisitor v) throws PrismLangException
	{
		return v.visit(this);
	}

	@Override
	public String toString() {
		String coalitionText = coalition.toString();
		return "<<"+coalitionText.substring(1, coalitionText.length() - 1)+">> " + (exprType==PRB?exprProb.toString():exprType==REW?exprRew.toString():"");
	}
	
	@Override
	public String getResultName()
	{
		return "Result for coalition " + coalition.toString();
	}
	
	@Override
	public Expression deepCopy() {
		
		ExpressionPATL expr = new ExpressionPATL();
		
		expr.setCoalition(new HashSet<String>(coalition));
		expr.setExpressionType(exprType);
		expr.setType(type);
		expr.setPosition(this);
		expr.setExpressionProb(exprProb==null?null:(ExpressionProb)exprProb.deepCopy());
		expr.setExpressionRew(exprRew==null?null:(ExpressionReward)exprRew.deepCopy());
		
		return expr;
	}

}
