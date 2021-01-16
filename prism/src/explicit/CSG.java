//==============================================================================
//	
//	Copyright (c) 2020-
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford)
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

import java.util.BitSet;
import java.util.Vector;

import prism.ModelType;
import prism.PlayerInfoOwner;

public interface CSG extends MDP, PlayerInfoOwner
{
	// Accessors (for Model) - default implementations
	
	@Override
	public default ModelType getModelType()
	{
		return ModelType.CSG;
	}
	
	// Accessors (for CSG)
	
	public BitSet[] getIndexes();
	
	public Vector<String> getActions();
	
	public int getIdleForPlayer(int p);
	
	public String[] getActions(int s, int i);
	
	public BitSet getIndexesForPlayer(int s, int p);
	
	public int[] getIndexes(int s, int i);
	
	public int[] getIdles();
	
	public BitSet getConcurrentPlayers(int s);
	
	// Temp:

	public Distribution getChoice(int s, int i);
}
