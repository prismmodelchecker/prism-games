//==============================================================================
//	
//	Copyright (c) 20014-
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

package prism;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;


public class Arrow
{	
    
    private Point tip;
    private Point tail;

    public Arrow(Point tip, Point tail)
    {
	this.tip = tip;
	this.tail = tail;
    }

    public Arrow toRealProperties(OpsAndBoundsList obl)
    {
	return new Arrow(tip.toRealProperties(obl), tail.toRealProperties(obl));
    }

    public double getTipX()
    {
	return tip.getCoord(0);
    }

    public double getTipY()
    {
	return tip.getCoord(1);
    }

    // angle in radians
    public double getAngle()
    {
	double angle = -Math.atan2(tip.getCoord(1)-tail.getCoord(1), tip.getCoord(0)-tail.getCoord(0));
	return angle;
    }

}
