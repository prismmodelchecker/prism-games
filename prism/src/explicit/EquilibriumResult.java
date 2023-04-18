//==============================================================================
//	
//	Copyright (c) 2020-
//	Authors:
//	* Gabriel Santos <gabriel.santos@cs.ox.ac.uk> (University of Oxford)
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

import explicit.CSGModelCheckerEquilibria.CSGResultStatus;

import java.util.ArrayList;

public class EquilibriumResult
{
	private CSGResultStatus status;
	private ArrayList<Double> payoffVector;
	private ArrayList<Distribution<Double>> strategy;

	public CSGResultStatus getStatus()
	{
		return status;
	}

	public void setStatus(CSGResultStatus status)
	{
		this.status = status;
	}

	public ArrayList<Double> getPayoffVector()
	{
		return payoffVector;
	}

	public void setPayoffVector(ArrayList<Double> payoffVector)
	{
		this.payoffVector = payoffVector;
	}

	public ArrayList<Distribution<Double>> getStrategy()
	{
		return strategy;
	}

	public void setStrategy(ArrayList<Distribution<Double>> strategy)
	{
		this.strategy = strategy;
	}
}