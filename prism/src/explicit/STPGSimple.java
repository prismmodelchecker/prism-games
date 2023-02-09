//==============================================================================
//	
//	Copyright (c) 2020-
//	Authors:
//	* Dave Parker <d.a.parker@cs.bham.ac.uk> (University of Birmingham)
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

import prism.ModelType;

/**
 * Simple explicit-state representation of a (turn-based) stochastic two-player game (STPG).
 * 
 * This is nothing more than a specific case of an SMG (stochastic multi-player game)
 * where the number of players can be assumed to be two. An SMG already implements the
 * {@link STPG} interface, via coalitions, so this class also implements {@link STPG}
 * but without the need to specify a coalition.
 * 
 * This class should rarely be needed - it is just to support the corner case of a
 * PRISM language model being defined directly as type "stpg", rather than "smg".
 */
public class STPGSimple extends SMGSimple
{
	@Override
	public ModelType getModelType()
	{
		return ModelType.STPG;
	}

	/**
	 * Constructor: empty STPG.
	 */
	public STPGSimple()
	{
		super();
	}

	/**
	 * Construct an STPG from an existing one and a state index permutation,
	 * i.e. in which state index i becomes index permut[i].
	 * Player and coalition info is also copied across.
	 */
	public STPGSimple(STPGSimple stpg, int permut[])
	{
		super(stpg, permut);
	}
}
