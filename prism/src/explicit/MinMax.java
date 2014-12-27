//==============================================================================
//	
//	Copyright (c) 2014-
//	Authors:
//	* Dave Parker <d.a.parker@cs.bham.ac.uk> (University of Birmingham/Oxford)
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

package explicit;

import parser.ast.Coalition;

/**
 * Class to store info about types of probabilities that are to be computed
 * (typically how to quantify over strategies, e.g. "min" or "max"). 
 */
public class MinMax
{
	// Info about quantification over a single class of strategies/adversaries
	
	protected boolean min;
	
	public void setMin(boolean min)
	{
		this.min = min;
	}
	
	public boolean isMin()
	{
		return min;
	}
	
	public boolean isMax()
	{
		return !min;
	}
	
	// Info about quantification over a two classes of strategies (e.g. for 2-player games)
	
	protected boolean min1;
	protected boolean min2;
	protected double bound;
	protected Coalition coalition;
	
	public void setMinMin(boolean min1, boolean min2)
	{
		this.min1 = min1;
		this.min2 = min2;
	}
	
	public void setBound(double bound)
	{
		this.bound = bound;
	}
	
	public void setCoalition(Coalition coalition)
	{
		this.coalition = coalition;
	}

	public boolean isMin1()
	{
		return min1;
	}
	
	public boolean isMin2()
	{
		return min2;
	}
	
	public double getBound()
	{
		return bound;
	}
	
	public Coalition getCoalition()
	{
		return coalition;
	}

	// Create a new instance by applying some operation
	
	public MinMax negate()
	{
		MinMax neg = new MinMax();
		neg.setMin(!isMin());
		neg.setMinMin(!isMin1(), !isMin2());
		neg.setBound(bound == -1 ? -1 : 1.0 - bound);
		return neg;
	}
	
	// Utility methods to create instances of this class
	
	public static MinMax blank()
	{
		return new MinMax();
	}
	
	public static MinMax min()
	{
		MinMax minMax = new MinMax();
		minMax.setMin(true);
		return minMax;
	}
	
	public static MinMax max()
	{
		MinMax minMax = new MinMax();
		minMax.setMin(false);
		return minMax;
	}
	
	public static MinMax minMin(boolean min1, boolean min2)
	{
		MinMax minMax = new MinMax();
		minMax.setMinMin(min1, min2);
		return minMax;
	}
}
