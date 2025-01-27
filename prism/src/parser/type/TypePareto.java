//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Clemens Wiltschhe <clemens.wiltsche@cs.ox.ac.uk> (University of Oxford)
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

package parser.type;

import explicit.Pareto;
import parser.EvaluateContext;
import prism.PrismLangException;

public class TypePareto extends Type 
{
	private static TypePareto singleton;
	
	static
	{
		singleton = new TypePareto();
	}	
	
	private TypePareto()
	{		
	}

	public static TypePareto getInstance()
	{
		return singleton;
	}

	// Methods required for Type:

	@Override
	public String getTypeString()
	{
		return "Pareto";
	}
	
	@Override
	public Object defaultValue()
	{
	    return new Pareto();
	}
	
	@Override
	public boolean canCastTypeTo(Type type)
	{
		return (type instanceof TypePareto);
	}
	
	@Override
	public Pareto castValueTo(Object value) throws PrismLangException
	{
		// only identity casting allowed
		if (value instanceof Pareto) {
			return (Pareto) value;
		} else {
		    throw new PrismLangException("Can't convert " + value.getClass() + " to type " + getTypeString());
		}
	}

	@Override
	public Pareto castValueTo(Object value, EvaluateContext.EvalMode evalMode) throws PrismLangException
	{
		// only identity casting allowed
		if (value instanceof Pareto) {
			return (Pareto) value;
		} else {
			throw new PrismLangException("Can't convert " + value.getClass() + " to type " + getTypeString());
		}
	}

	// Standard methods:

	public boolean equals(Object o)
	{
		return (o instanceof TypePareto);
	}
}
