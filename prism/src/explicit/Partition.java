//==============================================================================
//	
//	Copyright (c) 2002-
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

package explicit;

import java.util.BitSet;
import java.util.Map;
import java.util.List;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.Iterator;

import prism.PrismException;

public class Partition
{
	private List<List<Integer>> partition;
	private Map<Integer, Integer> f; // maps concrete state to abstract state
        private List<Double> weightsP1; // (optionally) assigns weights to states in P1 measure
        private List<Double> weightsP2; // (optionally) assigns weights to states in P2 measure
        public BitSet marked; // marks abstract states

	public Partition()
	{
	    partition = new ArrayList<List<Integer>>();
	    f = new HashMap<Integer, Integer>();
	    marked = new BitSet();
	}

	public Partition(int capacity)
	{
	    partition = new ArrayList<List<Integer>>(capacity);
	    f = new HashMap<Integer, Integer>();
	    marked = new BitSet(capacity);
	}

       /**
         * deep copy constructors
         **/
        public Partition(Partition partition)
        {
	    this(partition.partition);
	    this.marked = new BitSet(partition.marked.size());
	    this.marked.or(partition.marked);
	}
	public Partition(List<List<Integer>> partition)
	{
	    this();
	    int i = 0;
	    for(List<Integer> as : partition) {
		List<Integer> pi = new ArrayList<Integer>();
		for(Integer s : as) {
		    f.put(s, i);
		    pi.add(s);
		}
		this.partition.add(pi);
		i++;
	    }
	}

	public void clear()
	{
	    partition.clear();
	    f.clear();
	}

	public int size()
	{
	    return partition.size();
	}

	public int size(int i)
	{
	    return partition.get(i).size();
	}

	public Iterator<Integer> getComponent(int i) throws PrismException
	{
	    if(partition.get(i) == null || partition.get(i).size()==0) throw new PrismException("Empty component");
	    return partition.get(i).iterator();
	}

	public Iterator<List<Integer>> iterator()
	{
	    return partition.iterator();
	}
    
        public List<Double> getWeightsP1()
        {
	    return weightsP1;
	}
        public List<Double> getWeightsP2()
        {
	    return weightsP2;
	}

        public void setWeightsP1(List<Double> weightsP1)
        {
	    this.weightsP1 = weightsP1;
	}
        public void setWeightsP2(List<Double> weightsP2)
        {
	    this.weightsP2 = weightsP2;
	}

	public int get(int i, int j)
	{
	    return partition.get(i).get(j);
	}
	public List<Integer> get(int i)
	{
	    return partition.get(i);
	}

    public void merge()
    {
	Partition H = new Partition();
	int j = -1;
	int i = 0;
	for(List<Integer> as : partition) { // first must be unmarked (if it exists)
	    if(!marked.get(i)) {
		H.add(as);
		j = i;
		break;
	    }
	    i++;
	}

	if(j < 0) { // all marked
	    for(List<Integer> as : partition) {
		if(H.size()==0)
		    H.add(as);
		else
		    H.add(0, as);
	    }
	} else { // some unmarked
	    i = 0;
	    for(List<Integer> as : partition) {
		if(marked.get(i) && H.size()==1) {
		    H.add(as);
		} else if(marked.get(i)) {
		    H.add(1, as);
		} else if (i != j) {
		    H.add(0, as);
		}
		i++;
	    }
	}
	    
	// now set
	this.weightsP1.clear();
	this.weightsP2.clear();
	this.partition = H.partition;
	this.f = H.f;
    }

        // all abstract states that are not marked are merged into one state
        public void mergeUnmarked()
        {
	    Partition H = new Partition();
	    int j = -1;
	    int i = 0;
	    for(List<Integer> as : partition) { // first must be unmarked (if exists)
		if(!marked.get(i)) {
		    H.add(as);
		    j = i;
		    break;
		}
		i++;
	    }

	    i = 0;
	    for(List<Integer> as : partition) { // first must be unmarked (if exists)
		if(marked.get(i)) {
		    H.add(as);
		} else if (i != j) {
		    H.add(0, as);
		}
		i++;
	    }

	    // now set
	    if(weightsP1!=null) this.weightsP1.clear();
	    if(weightsP2!=null) this.weightsP2.clear();
	    this.partition = H.partition;
	    this.f = H.f;
	}

	public Map<Integer, Integer> getf()
	{
	    return f;
	}

	/**
	 *  assume S is a refinement of Pi' is a refinement of Pi, i.e. S <= Pi' <= Pi
	 *  f maps S to Pi (this partition)
	 *  g maps S to Pi'
	 *  return: map from Pi' to Pi
	 **/
	public Map<Integer, Integer> getf(Map<Integer, Integer> g) throws PrismException
	{
	    Map<Integer, Integer> h = new HashMap<Integer, Integer>();
	    for(Map.Entry<Integer, Integer> s_pi : f.entrySet()) {
		int pi = s_pi.getValue();
		int pi_prime = g.get(s_pi.getKey());

		if(h.get(pi_prime) == null)
		    h.put(pi_prime, pi);
		else if (h.get(pi_prime) != pi)
		    throw new PrismException("Partition is not a valid refinement");
	    }

	    return h;
	}

	public boolean add(List<Integer> pi)
	{
	    for(Integer s : pi) {
		f.put(s, partition.size());
	    }
	    return partition.add(pi);
	}

    public boolean add(int i, List<Integer> pi)
    {
	for(Integer s : pi) {
	    f.put(s, i);
	}
	return partition.get(i).addAll(pi);
    }

	public boolean add(int i, int s)
	{
	    f.put(s, i);
	    if(partition.size() <= i) {
		for(int j = partition.size(); j <= i; j++) {
		    partition.add(new ArrayList<Integer>());
		}
	    }
	    return partition.get(i).add(s);
	}

	/**
	 * split state i by partition into H
	 * appends all but one new partition to the end
	 * optionally, add the weights of H (default false)
	 **/
        public void split(int i, Partition H)
        {
	    split(i, H, false);
	}
        public void split(int i, Partition H, boolean addWeights)
        {
	    Iterator<Double> wP1_iterator = null, wP2_iterator = null;
	    if(addWeights) {
		if(weightsP1==null) {
		    weightsP1 = new ArrayList<Double>();
		    for(int j = 0; j < partition.size() + H.size()-1; j++) {
			weightsP1.add(0.0); // fill with zeros
		    }
		}
		if(weightsP2==null) {
		    weightsP2 = new ArrayList<Double>();
		    for(int j = 0; j < partition.size() + H.size()-1; j++) {
			weightsP2.add(0.0); // fill with zeros
		    }
		}
		weightsP1.remove(i);
		weightsP2.remove(i);

		wP1_iterator = H.weightsP1.iterator();
		wP2_iterator = H.weightsP2.iterator();
	    }

	    partition.remove(i);
	    Iterator<List<Integer>> H_iterator = H.iterator();
	    boolean first = true;
	    while(H_iterator.hasNext()) {
		List<Integer> as = H_iterator.next();
		if(first==true) {
		    partition.add(i, as);
		    if(addWeights) {
			weightsP1.add(i, wP1_iterator.next());
			weightsP2.add(i, wP2_iterator.next());
		    }
		    first = false;
		    for(Integer s : as) { // not sure if this is necessary
			f.put(s, i);
		    }
		} else {
		    // append to end
		    if(addWeights) {
			weightsP1.add(wP1_iterator.next());
			weightsP2.add(wP2_iterator.next());
		    }
		    partition.add(as);
		    for(Integer s : as) { // not sure if this is necessary
			f.put(s, partition.size()-1);
		    }
		}
	    }
	}

	public String toString()
	{
	    String result = "[";
	    boolean first = true; 
	    for(int as = 0; as < partition.size(); as++) {
		List<Integer> abstract_state = partition.get(as);
		if(!first)
		    result += ", ";
		else
		    first = false;
		String mark = marked.get(as)?"*":"";
		if(weightsP1!=null && as < weightsP1.size() && weightsP2!=null && as < weightsP2.size())
		    result += String.format("%d(%f|%f)%s:%s", as, weightsP1.get(as), weightsP2.get(as), mark, abstract_state.toString());
		else if(weightsP1!=null && as < weightsP1.size())
		    result += String.format("%d(%f|_)%s:%s", as, weightsP1.get(as), mark, abstract_state.toString());
		else if(weightsP2!=null && as < weightsP2.size())
		    result += String.format("%d(_|%f)%s:%s", as, weightsP2.get(as), mark, abstract_state.toString());
		else
		    result += String.format("%d%s:%s", as, mark, abstract_state.toString());
	    }
	    return result + "]";
	}
}
