// ==============================================================================
//	
// Copyright (c) 2013-
// Authors:
// * Clemens Wiltsche <clemens.wiltsche@stx.ox.ac.uk> (University of Oxford)
//	
// ------------------------------------------------------------------------------
//	
// This file is part of PRISM.
//	
// PRISM is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//	
// PRISM is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//	
// You should have received a copy of the GNU General Public License
// along with PRISM; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//	
// ==============================================================================

package explicit;

import java.util.Vector;
import java.util.Map;
import java.math.BigInteger;

import parma_polyhedra_library.*;

public class PPLSupport
{
    /*
    static protected Vector<BigInteger> getCoefficientsFromGenerator(Generator g)
    {
	Vector<BigInteger> result = new Vector<BigInteger>();

	Linear_Expression le = g.linear_expression();
	
	Map<Variable, BigInteger> map = new Map<Variable, BigInteger>();
	getCoefficientsFromLinearExpression(le, false, BigInteger.ONE, map);

	
	for(int i = 
	
    }
    */

    static void getCoefficientsFromLinearExpression(Linear_Expression le, boolean minus, BigInteger coefficient, Map<Variable, BigInteger> result)
    {
	
	if(le instanceof Linear_Expression_Coefficient){
	    if(result.containsKey(null)){
		result.put(null,
			   minus
			   ?
			   result.get(null).multiply(((Linear_Expression_Coefficient)le).argument().getBigInteger().multiply(coefficient).negate())
			   :
			   result.get(null).multiply(((Linear_Expression_Coefficient)le).argument().getBigInteger().multiply(coefficient)));
	    }
	    else {
		result.put(null,
			   minus
			   ?
			   ((Linear_Expression_Coefficient)le).argument().getBigInteger().multiply(coefficient).negate()
			   :
			   ((Linear_Expression_Coefficient)le).argument().getBigInteger().multiply(coefficient));
	    }
	} else if (le instanceof Linear_Expression_Difference){
	    getCoefficientsFromLinearExpression(((Linear_Expression_Difference)le).left_hand_side(), minus, coefficient, result);
	    getCoefficientsFromLinearExpression(((Linear_Expression_Difference)le).right_hand_side(), !minus, coefficient, result);
	} else if (le instanceof Linear_Expression_Sum){
	    getCoefficientsFromLinearExpression(((Linear_Expression_Sum)le).left_hand_side(), minus, coefficient, result);
	    getCoefficientsFromLinearExpression(((Linear_Expression_Sum)le).right_hand_side(), minus, coefficient, result);
	} else if (le instanceof Linear_Expression_Times){
	    getCoefficientsFromLinearExpression(((Linear_Expression_Times)le).linear_expression(), minus, coefficient.multiply(((Linear_Expression_Times)le).coefficient().getBigInteger()), result);
	} else if (le instanceof Linear_Expression_Unary_Minus) {
	    getCoefficientsFromLinearExpression(((Linear_Expression_Unary_Minus)le).argument(),!minus, coefficient, result);
	} else if (le instanceof Linear_Expression_Variable) {
	    if(result.containsKey(((Linear_Expression_Variable)le).argument())){
		result.put(((Linear_Expression_Variable)le).argument(),
			   minus
			   ?
			   result.get(((Linear_Expression_Variable)le).argument()).multiply(coefficient).negate()
			   :
			   result.get(((Linear_Expression_Variable)le).argument()).multiply(coefficient));
	    } else {
		result.put(((Linear_Expression_Variable)le).argument(),
			   minus
			   ?
			   coefficient.negate()
			   :
			   coefficient);
	    }
	}

    }

}
