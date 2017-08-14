//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford)
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

import java.util.ArrayList;
import java.util.List;

import parser.visitor.ASTVisitor;
import prism.PrismLangException;

/**
 * Class representing a player definition in a model file.
 */
public class Player extends ASTElement
{
	// Name
	private String name;
	// Modules (names)
	private List<String> modules;
	// Actions (names)
	private List<String> actions;

	// Constructor

	public Player(String name)
	{
		this.name = name;
		modules = new ArrayList<String>();
		actions = new ArrayList<String>();
	}

	// Set methods

	public void addModule(String moduleName)
	{
		modules.add(moduleName);
	}

	public void addAction(String actionName)
	{
		actions.add(actionName);
	}

	// Get methods

	public String getName()
	{
		return name;
	}
	
	public List<String> getModules()
	{
		return modules;
	}

	public List<String> getActions()
	{
		return actions;
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
		String s = "";
		boolean first1, first2;

		s += "player " + name + "\n\t";
		first1 = true;
		for (String m : modules) {
			if (first1)
				first1 = false;
			else
				s += ", ";
			s += m;
		}
		first2 = true;
		for (String a : actions) {
			if (first2) {
				first2 = false;
				if (!first1)
					s += ",\n\t";
			} else {
				s += ", ";
			}
			s += "[" + a + "]";
		}
		s += "\nendplayer";

		return s;
	}

	/**
	 * Perform a deep copy.
	 */
	public Player deepCopy()
	{
		Player ret = new Player(name);
		for (String m : modules)
			ret.addModule(m);
		for (String a : actions)
			ret.addAction(a);
		ret.setPosition(this);
		return ret;
	}
}

// ------------------------------------------------------------------------------
