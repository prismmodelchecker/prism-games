//==============================================================================
//	
//	Copyright (c) 2014-
//	Authors:
//	* Clemens Wiltsche <clemens.wiltsche@cs.ox.ac.uk> (University of Oxford)
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


package userinterface.graph;

import java.util.Map;
import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;
import java.util.Set;


public class Slice
{
    private Map<Integer, Double> slice; // dimensions to be sliced, together with the value at which they are sliced
    private List<Integer> project; // dimensions to be projected out

    public boolean swapdimensions = false; // the first element is X if false, and Y if true

    public Slice()
    {
	slice = new HashMap<Integer, Double>();
	project = new ArrayList<Integer>();
	swapdimensions = false;
    }
    
    public void project(int index)
    {
	project.add(index);
    }

    public List<Integer> getProject()
    {
	return project;
    }

    public Double put(int index, double value)
    {
	return slice.put(index, value);
    }

    public int fullSize()
    {
	return slice.size() + project.size();
    }

    /**
     * Returns whether the index i is removed due to slicing or projection
     **/
    public boolean removed(int i)
    {
	return (slice.containsKey(i) || project.contains(i));
    }

    public void clear()
    {
	project.clear();
	slice.clear();
	swapdimensions = false;
    }

    public Set<Integer> keySet()
    {
	return slice.keySet();
    }

    public Double get(Integer key)
    {
	return slice.get(key);
    }

    public boolean containsKey(Integer key)
    {
	return slice.containsKey(key);
    }

    public String toString()
    {
	return String.format("slice: %s, project %s, swap: %s", slice, project, swapdimensions);
    }
}
