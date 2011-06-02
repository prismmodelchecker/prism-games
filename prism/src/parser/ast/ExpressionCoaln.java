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

import parser.EvaluateContext;
import parser.visitor.ASTVisitor;
import prism.PrismLangException;

public class ExpressionCoaln extends Expression
{
	Set members = new HashSet();

	/**
	 * Adds a member to the coalition
	 * @param m member
	 */
	public void addMember(Object m)
	{
		members.add(m);
	}

	/**
	 * 
	 * @return members of the coalition
	 */
	public Set getMembers()
	{
		return members;
	}

	@Override
	public Expression deepCopy()
	{
		ExpressionCoaln ret = new ExpressionCoaln();
		for (Object m : members)
			ret.addMember(m);
		return ret;
	}

	@Override
	public Object evaluate(EvaluateContext ec) throws PrismLangException
	{
		throw new PrismLangException("Cannot evaluate a <<C>> operator without a model");
	}

	@Override
	public boolean isConstant()
	{
		return false;
	}

	@Override
	public boolean returnsSingleValue()
	{
		return false;
	}

	@Override
	public Object accept(ASTVisitor v) throws PrismLangException
	{
		return v.visit(this);
	}

	@Override
	public String toString()
	{
		String out = "";
		for(Object o: members)
			out += o+",";
		return out.substring(0,out.length()-1);
	}

}
