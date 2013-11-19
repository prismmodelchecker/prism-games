//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Mateusz Ujma <mateusz.ujma@cs.ox.ac.uk> (University of Oxford)
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

import java.util.BitSet;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;

import explicit.SCCComputerTarjan;
import explicit.SCCComputer;
import explicit.ECComputer;
import explicit.SMG;

// Class for generating detailed statistics about the model

public class ModelStatistics
{

	private SMG model;
	
	public ModelStatistics(SMG model) {
		this.model = model;
	}
	
	public void generateStatistics(PrismLog log) throws PrismException{
		//players, actions, outcomes
		int player1 = 0;
		int player1ActionMax = 0;
		int player1ActionMin = Integer.MAX_VALUE;
		int player1ActionSum = 0;
		int player1OutcomesMax = 0;
		int player1OutcomesMin = Integer.MAX_VALUE;
		int player1OutcomesSum = 0;
		double player1OutcomesAvg = 0;
		double player1ActionAvg = 0;
		int player2 = 0;
		int player2ActionMax = 0;
		int player2ActionMin = Integer.MAX_VALUE;
		int player2ActionSum = 0;
		double player2ActionAvg = 0;
		int player2OutcomesMax = 0;
		int player2OutcomesMin = Integer.MAX_VALUE;
		int player2OutcomesSum = 0;
		double player2OutcomesAvg = 0;
		
		for(int i=0;i<model.getNumStates();i++) {
			if(model.getPlayer(i) == SMG.PLAYER_1) {
				player1++;
				player1ActionMax = Math.max(player1ActionMax, model.getNumChoices(i));
				player1ActionMin = Math.min(player1ActionMin, model.getNumChoices(i));
				player1ActionSum += model.getNumChoices(i);
				
				for(int j=0;j<model.getNumChoices(i);j++) {
					player1OutcomesMax = Math.max(player1OutcomesMax, model.getChoice(i, j).getSupport().size());
					player1OutcomesMin = Math.min(player1OutcomesMin, model.getChoice(i, j).getSupport().size());
					player1OutcomesSum += model.getChoice(i, j).getSupport().size();
				}
			} else {
				player2ActionMax = Math.max(player2ActionMax, model.getNumChoices(i));
				player2ActionMin = Math.min(player2ActionMin, model.getNumChoices(i));
				player2ActionSum += model.getNumChoices(i);
				
				for(int j=0;j<model.getNumChoices(i);j++) {
					player2OutcomesMax = Math.max(player2OutcomesMax, model.getChoice(i, j).getSupport().size());
					player2OutcomesMin = Math.min(player2OutcomesMin, model.getChoice(i, j).getSupport().size());
					player2OutcomesSum += model.getChoice(i, j).getSupport().size();
				}
			}
		}
		player2 = model.getNumStates() - player1;
		
		player1ActionAvg = (double)player1ActionSum/(double)player1;
		player2ActionAvg = (double)player2ActionSum/(double)player2;
		
		player1OutcomesAvg = (double)player1OutcomesSum/(double)player1ActionSum;
		player2OutcomesAvg = (double)player2OutcomesSum/(double)player2ActionSum;
		
		log.println("Player1 states: " + player1);
		log.println("Player2 states: " + player2);
		log.println("Player1 actions max: " + player1ActionMax);
		log.println("Player1 actions min: " + player1ActionMin);
		log.println("Player1 actions avg: " + player1ActionAvg);
		log.println("Player2 actions max: " + player2ActionMax);
		log.println("Player2 actions min: " + player2ActionMin);
		log.println("Player2 actions avg: " + player2ActionAvg);
		log.println("Player1 outcomes max: " + player1OutcomesMax);
		log.println("Player1 outcomes min: " + player1OutcomesMin);
		log.println("Player1 outcomes avg: " + player1OutcomesAvg);
		log.println("Player2 outcomes max: " + player2OutcomesMax);
		log.println("Player2 outcomes min: " + player2OutcomesMin);
		log.println("Player2 outcomes avg: " + player2OutcomesAvg);
		log.println("Diameter: " + computeModelDiameter(model));
		log.println("SCCs " + SCCStatistics(model));
		log.println("End components" + MECStatistics(model));
	}
	
	private String SCCStatistics(SMG model) throws PrismException{
		SCCComputer sccComp = new SCCComputerTarjan(null, model);
		sccComp.computeSCCs();
		List<BitSet> sccs = sccComp.getSCCs();
		int countNonTrivial = 0;
		int nonTrivialSize = 0;
		for(int i=0;i<sccs.size();i++) {
			if(sccs.get(i).cardinality() > 1) {
				countNonTrivial++;
				nonTrivialSize += sccs.get(i).cardinality();
			}
		}
		return "SCCs: " + sccs.size()  + "\nNon-trivial SCCs: " + countNonTrivial + " \nAverage size: " + (double)nonTrivialSize/(double)countNonTrivial;
	}
	
	private String MECStatistics(SMG model) throws PrismException{
		ECComputer eccComp = ECComputer.createECComputer(null, model);
		eccComp.computeMECStates();
		List<BitSet> eccs = eccComp.getMECStates();
		int countNonTrivial = 0;
		int nonTrivialSize = 0;
		for(int i=0;i<eccs.size();i++) {
			if(eccs.get(i).cardinality() > 1) {
				countNonTrivial++;
				nonTrivialSize += eccs.get(i).cardinality();
			}
		}
		return "MECs: " + eccs.size()  + "\nNon-trivial MECs: " + countNonTrivial + " \nAverage size: " + (double)nonTrivialSize/(double)countNonTrivial;
	}
	
	private int computeModelDiameter(SMG model) {
		computeShortestPathDijkstra(model);
		Iterator<Vertex> it = state2Vertex.values().iterator();
		int maxDist = 0;
		while(it.hasNext()) {
			maxDist = Math.max(maxDist, it.next().distance);
		}
		
		return maxDist;
	}
	
	private Map<Integer,Vertex> computeShortestPathDijkstra(SMG model) {
		
		Vertex initialState = getVertex(model.getFirstInitialState());
		initialState.setDistance(0);
		
        PriorityQueue<Vertex> vertexQueue = new PriorityQueue<Vertex>();
      	vertexQueue.add(initialState);

		while (!vertexQueue.isEmpty()) {
		    Vertex u = vertexQueue.poll();
		    Iterator<Integer> it = model.getSuccessorsIterator(u.getState());
		    
		    while(it.hasNext()) {
		    	Vertex v = getVertex(it.next());
		    	int dist = u.getDistance() + 1;
		    	if (dist < v.getDistance()) {
	    		    vertexQueue.remove(v);
	    		    v.setDistance(dist) ;
	    		    vertexQueue.add(v);
	    		}
		    	
		    }
		}	            
		return state2Vertex;  
	}
	
	private Vertex getVertex(int state) {
		Vertex v = state2Vertex.get(state);
		if(v == null) {
			v = new Vertex(state);
			state2Vertex.put(state, v);
		} 
		return v;
	}
	
	private Map<Integer,Vertex> state2Vertex = new HashMap<Integer,Vertex>();
	
	class Vertex implements Comparable<Vertex>{
		private int state;
		private int distance = Integer.MAX_VALUE;
		
		public Vertex(int s) {
			state = s;
		}
		
		public int getDistance() {
			return distance;
		}
		
		public void setDistance(int d) {
			distance = d;
		}
		
		public int getState() {
			return state;
		}
		
		public int compareTo(Vertex v)
	    {
	        return Double.compare(getDistance(), v.getDistance());
	    }
	}
	
	

}

//------------------------------------------------------------------------------
