//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Clemens Wiltschhe <clemens.wiltsche@cs.ox.ac.uk> (University of Oxford)
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
import java.util.List;
import java.math.BigInteger;

import org.apache.commons.math3.fraction.BigFraction;

import prism.PrismException;

import parma_polyhedra_library.Polyhedron;
import parma_polyhedra_library.Coefficient;
import parma_polyhedra_library.Linear_Expression_Variable;
import parma_polyhedra_library.Variable;
import parma_polyhedra_library.C_Polyhedron;
import parma_polyhedra_library.Generator_Type;
import parma_polyhedra_library.Constraint;
import parma_polyhedra_library.Linear_Expression;
import parma_polyhedra_library.Linear_Expression_Sum;
import parma_polyhedra_library.Linear_Expression_Times;
import parma_polyhedra_library.Generator;
import parma_polyhedra_library.Generator_System;

public class Pareto
{
    private List<Polyhedron> sets;

    // default constructor - empty set
    public Pareto()
    {
	sets = new ArrayList<Polyhedron>();
    }

    public Pareto(int capacity)
    {
	sets = new ArrayList<Polyhedron>(capacity);
    }

    /**
     * Deep copy constructor for Pareto set. Empty Polyhedra are ignored and not added.
     **/
    public Pareto(Pareto pareto)
    {
	sets = new ArrayList<Polyhedron>(pareto.size());
	for(Polyhedron p : pareto.getSets())
	    if(!p.is_empty()) // only nonempty polyhedra
		sets.add(new C_Polyhedron((C_Polyhedron)p)); // deep copy
    }

    public Pareto(Polyhedron set) {
	this(1);
	add(set);
    }

    /**
     * Deep copy constructor for list of Polyhedra.
     **/
    public Pareto(List<Polyhedron> sets)
    {
	this.sets = new ArrayList<Polyhedron>(sets.size());
	for(Polyhedron p : sets) {
	    this.sets.add(new C_Polyhedron((C_Polyhedron)p)); // deep copy
	}
    }

    public void replace(int i, Polyhedron set)
    {
	sets.set(i, set);
    }

    public void add(Pareto set) {
	for(Polyhedron p : set.getSets()) {
	    add(p);
	}
    }

    /**
     * Assuming all sets have the same dimension, returns it.
     **/
    public long getDimension()
    {
	return sets.get(0).space_dimension();
    }

    /**
     * add a new Polyhedron to the list
     * ensures that no Polyhedra are completely contained in others
     **/
    public void add(Polyhedron set) {
	List<Polyhedron> to_remove = new ArrayList<Polyhedron>();
	to_remove.clear();
	for(Polyhedron p : this.sets) {
	    if(p.contains(set)) { // already covered - doesn't need to be added
		return;
	    } else if (set.contains(p)) { // p contains a set already in the list - can be removed
		to_remove.add(p);
	    }
	}
	this.sets.removeAll(to_remove);
	this.sets.add(set);
    }

    public void interior(double varepsilon)
    {
	if(varepsilon < 0) return;

	Pareto result = new Pareto(this.sets.size());
	for(Polyhedron p : this.sets) {
	    long n = p.space_dimension();

	    Generator_System ngs = new Generator_System();
	    // first set up the reward vector that should be added to each point generator
	    
	    BigFraction r = new BigFraction(-varepsilon);
	    BigInteger num = r.getNumerator();
	    BigInteger den = r.getDenominator();

	    // prepare vector pointing in direction (-varepsilon, -varepsilon, ...)
	    Linear_Expression le = new Linear_Expression_Times(new Coefficient(num), new Variable(0));
	    Coefficient c = new Coefficient(den);
	    for (int i = 1; i < n; i++) {
		le = new Linear_Expression_Sum(le, new Linear_Expression_Times(new Coefficient(num), new Variable(i)));
	    }
	    
	    // now add reward vector to each point generator 
	    for (Generator g : p.generators()) {
		if (g.type() == Generator_Type.POINT) {
		    Linear_Expression nle = new Linear_Expression_Sum(le.times(g.divisor()), g.linear_expression().times(c));
		    Coefficient nc = new Coefficient(g.divisor().getBigInteger().multiply(c.getBigInteger()));
		    ngs.add(Generator.point(nle, nc));
		} else {
		    ngs.add(g);
		}
	    }
	    result.sets.add(new C_Polyhedron(ngs));
	}
	this.sets = result.getSets();
    }

    /**
     * turn the Pareto set stored in this with the directions turned in the right way
     **/
    public void turn(List<Variable> negative_dimensions)
    {
	this.sets = turned(negative_dimensions).getSets();
    }
    /**
     * Return the Pareto set with directions turned the right way.
     * Makes a deep copy.
     **/
    public Pareto turned(List<Variable> negative_dimensions)
    {
	Pareto result = new Pareto(this.sets.size());
	for(Polyhedron p : this.sets) {
	    Polyhedron q = new C_Polyhedron(p.generators()); // deep copy
	    for(Variable v : negative_dimensions) {
		q.affine_image(v, new Linear_Expression_Variable(v), new Coefficient(-1));
	    }
	    result.sets.add(q);
	}
	return result;
    }

    public Polyhedron get() throws PrismException
    {
	if(sets.size() == 1) {
	    return sets.get(0);
	} else {
	    throw new PrismException("Pareto set does not consist of a single Polyhedron");
	}
    }

    public Polyhedron get(int i)
    {
	return sets.get(i);
    }

    public List<Polyhedron> getSets()
    {
	return sets;
    }
    
    public int size()
    {
	return sets.size();
    }

    public boolean isConvex()
    {
	if(sets.size()<=1)
	    return true;
	else {
	    // TODO: test if union of individual sets remains convex
	}

	return false;
    }
    
}
