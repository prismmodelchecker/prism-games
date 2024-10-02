//==============================================================================
//
//	Copyright (c) 2024-
//	Authors:
//	* Dave Parker <david.parker@cs.ox.ac.uk> (University of Oxford)
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

package prism;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Class representing a joint action in a concurrent game model
 */
public class JointAction extends ArrayList<Object>
{
	// Specification of idle actions
	public static class IdleAction { public String toString() { return IDLE_ACTION_STRING; } }
	public static final Object IDLE_ACTION = new IdleAction();
	public static final String IDLE_ACTION_STRING = "-";

	/**
	 * Construct a new empty joint action.
	 */
	public JointAction()
	{
		super();
	}

	/**
	 * Construct a new joint action with {@code size} undefined (null) actions.
	 */
	public JointAction(int size)
	{
		super(size);
		addAll(Collections.nCopies(size, null));
	}

	/**
	 * Construct a new joint action based on (1-indexed) indices into a list of actions.
	 * An index of -1 denotes the idle action.
	 */
	public JointAction(int[] indexes, List<Object> actions)
	{
		super(indexes.length);
		for (int i : indexes) {
			add(i >= 0 ? actions.get(i - 1) : JointAction.IDLE_ACTION);
		}
	}

	/**
	 * Copy constructor.
	 */
	public JointAction(JointAction joint)
	{
		// Shallow copy
		super(joint);
	}

	@Override
	public String toString()
	{
		// Comma-separated list of action strings, with "" for undefined, e.g. "[a,b,,c,-,d]"
		return "[" + stream().map(a -> a == null ? "" : a.toString()).collect(Collectors.joining(",")) + "]";
	}
}
