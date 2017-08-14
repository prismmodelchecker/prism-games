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

import prism.PrismException;
import java.util.Arrays;

public class HyperplaneSelector
{

    private double[][] y;
    private int[] exclude;
    private double start;
    private double step;
    private int CONJUNCTS;
    private int[] DISJUNCTS;

    public HyperplaneSelector(int CONJUNCTS, int[] DISJUNCTS)
    {
	this.CONJUNCTS = CONJUNCTS;
	this.DISJUNCTS = DISJUNCTS;
	this.start = 0.0;
	this.step = 1.0;
	init_next_hyperplane();
    }

    public double getStart()
    {
	return start;
    }
    
    public int[] getExclude()
    {
	int[] result = new int[exclude.length];
	System.arraycopy(exclude, 0, result, 0, exclude.length);
	return result;
    }
    
    private double[][] init_next_hyperplane()
    {
	// initialise exclusion
	exclude = new int[CONJUNCTS]; // note: initialised to zero

	// initialise vector to be returned
	double[][] x = new double[CONJUNCTS][];
	for(int i = 0; i < CONJUNCTS; i++) {
	    x[i] = new double[DISJUNCTS[i]];
	    Arrays.fill(x[i], 1.0-start);
	}
	y = x;
	return x;
    }

    public double[][] next_hyperplane() throws PrismException
    {
	try {
	    return next_hyperplane(y);
	} catch (PrismException e) { // no more iterations needed at this accuracy
	    // increase accuracy
	    if(start>0.0) {
		step = start;
		start = start/2.0;
	    } else {
		step = 1;
		start = 0.5;
	    }
	    // reinitialise hyperplane
	    return init_next_hyperplane();
	}
    }

    private double[][] next_hyperplane(double[][] x) throws PrismException
    { 
	// first evaluate overall length of x
	int overall_length = 0;
	for(int i = 0; i < x.length; i++) {
	    overall_length += x[i].length;
	}
	
	// now create next x
	boolean carry = true;
	int i = 0;
	int j = 0;
	int count = 0;
	boolean ex = true; // check exclusion
	while(ex) {
	    ex = false;
	    if(exclude[i]==j) { // current index is excluded
		j++;
		count++;
		if(j>=x[i].length) { // moved j too far, need to move i
		    j=0;
		    i++;
		    ex = true; // need to re-check exclusion after moving i
		    if(i >= x.length) { // moved i too far, need to move exclusion
			exclude[0]++;
			boolean exx = true;
			int ii = 0;
			while(exx) {
			    exx = false;
			    if(exclude[ii]>=x[ii].length) { // moved exclusion too far, need to move the next one
				exclude[ii] = 0;
				ii++;
				if(ii>=x.length) { // moved ii too far, done iterating
				    throw new PrismException("Done iterating");
				}
				exclude[ii]++; // move next one
				exx = true; // re-check moving exclusion
			    }
			}
			// after moving exclusion, return new hyperplane
			return next_hyperplane(x);
		    }
		}
	    }
	}
	while(carry) {
	    carry = false;
	    x[i][j] -= step;
	    // TODO: replace constant
	    if(x[i][j]<0.0-1.0e-8) { // if too big need to reset
		carry = true;
		x[i][j] = 1.0-start; // reset to start
		
		// 1. move j to next one
		j++;
		count++;
		if(j>=x[i].length) { // moved j too far, need to move i
		    j=0;
		    i++;
		    if(i>=x.length) { // moved i too far, need to move exclusion
			exclude[0]++;
			boolean exx = true;
			int ii = 0;
			while(exx) {
			    exx = false;
			    if(exclude[ii]>=x[ii].length) { // moved exclusion too far, need to move the next one
				exclude[ii] = 0;
				ii++;
				if(ii>=x.length) { // moved ii too far, done iterating
				    throw new PrismException("Done iterating");
				}
				exclude[ii]++; // move next one
				exx = true; // re-check moving exclusion
			    }
			}
			// after moving exclusion, return new hyperplane
			return next_hyperplane(x);
		    }
		}
		// at this point have moved j, and potentially i
		// at this point have not moved i too far, and have not moved exclusion too far

		// 2. check exclusion
		ex = true; // check exclusion
		while(ex) {
		    ex = false;
		    if(exclude[i]==j) { // current index is excluded
			j++;
			count++;
			if(j>=x[i].length) { // moved j too far, need to move i
			    j=0;
			    i++;
			    ex = true; // need to re-check exclusion after moving i
			    if(i >= x.length) { // moved i too far, need to move exclusion
				exclude[0]++;
				boolean exx = true;
				int ii = 0;
				while(exx) {
				    exx = false;
				    if(exclude[ii]>=x[ii].length) { // moved exclusion too far, need to move the next one
					exclude[ii] = 0;
					ii++;
					if(ii>=x.length) { // moved ii too far, done iterating
					    throw new PrismException("Done iterating");
					}
					exclude[ii]++; // move next one
					exx = true; // re-check moving exclusion
				    }
				}
				// after moving exclusion, return new hyperplane
				return next_hyperplane(x);
			    }
			}
		    }
		}
	    }
	}
	return x;
    }


}
