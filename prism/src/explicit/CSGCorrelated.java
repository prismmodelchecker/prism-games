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

import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;

public interface CSGCorrelated {
	
	/**
	 * Computes correlated equilibria for a normal form game
	 * 
	 * @param utilities Utility table as a map from action indexes to utility arrays
	 * @param ce_constraints An auxiliary mapping of actions to utilities (player i -> action a -> actions of -i -> utility) used to set the constraints for correlated equilibria 
	 * @param strategies Array with individual actions for each player
	 * @param ce_var_map Map for joint actions
	 * @param crit Criterion (social-welfare = 3, fair = 4)
	 */
	public EquilibriumResult computeEquilibrium(HashMap<BitSet, ArrayList<Double>> utilities, 
												ArrayList<ArrayList<HashMap<BitSet, Double>>> ce_constraints,
												ArrayList<ArrayList<Integer>> strategies,
												HashMap<BitSet, Integer> ce_var_map, int crit);
	
	/**
	 * Gets the solver's name
	 */
	public String getSolverName();
	
	public void clear();
		
	public void printModel();
	
}
