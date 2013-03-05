package explicit;

import java.io.IOException;
//import java.io.ClassNotFoundException;
import java.io.FileOutputStream;
import java.io.FileInputStream;
import java.io.ObjectOutputStream;
import java.io.ObjectInputStream;

import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;

import prism.PrismException;
import parser.State;
import strat.MultiObjectiveStrategy;
import explicit.rewards.STPGRewards;
import explicit.STPG;

public class MapMDPSimulator
{

    FileInputStream fi;
    FileOutputStream fo;
    ObjectOutputStream oo;
    ObjectInputStream oi;

    MultiObjectiveStrategy strat = null;
    List<STPGRewards> stpgRewards;
    STPG model;

    public MapMDPSimulator(STPG model, List<STPGRewards> stpgRewards)
    {
	this.stpgRewards = stpgRewards;
	this.model = model;
    }

    @Override
    protected void finalize() throws Throwable
    {
	try {
	    fi.close();
	    oi.close();
	    fo.close();
	    oo.close();
	} finally {
	    super.finalize();
	}
    }

    public void writeStrategy(MultiObjectiveStrategy strategy, String filename)
    {
	try {
	    fo = new FileOutputStream(filename);
	    oo = new ObjectOutputStream(fo);
	    oo.writeObject(strategy);
	    fo.close();
	    oo.close();
	} catch (IOException ie) {
	    ie.printStackTrace();
	}
    }

    public boolean readStrategy(String filename)
    {
	try {
	    fi = new FileInputStream(filename);
	    oi = new ObjectInputStream(fi);
	    strat = (MultiObjectiveStrategy) oi.readObject();
	    fi.close();
	    oi.close();
	    return true;
	} catch (IOException ie) {
	    ie.printStackTrace();
	    return false;
	} catch (ClassNotFoundException c) {
	    System.out.println("Multi Objective Strategy class not found");
	    c.printStackTrace();
	    return false;
	}
    }

    public void evaluateStrategy() throws PrismException
    {
	
	System.out.println("---- SAMPLES: ----");
	List<List<State>> samples = strat.simulateMDP(10000);
	//for(int sample = 0; sample < samples.size(); sample++) {
	//	System.out.printf("%d: %s\n", sample, samples.get(sample).toString());
	//}
	
	System.out.println("COLLATED PATHS:");
	Map<Integer,Double> expected_reward = new HashMap<Integer,Double>(stpgRewards.size());
	Map<List<State>,Double> c_samples = collateSamples(samples);
	List<State> max_sample = samples.get(0);
	double max_sample_prob = 0.0;
	for(List<State> c_sample : c_samples.keySet()) {
	    if(c_samples.get(c_sample) > max_sample_prob) {
		max_sample_prob = c_samples.get(c_sample);
		max_sample = c_sample;
	    }
	    //System.out.printf("%.4f: %s\n", c_samples.get(c_sample), c_sample.toString());
	    for(int r = 0; r < stpgRewards.size(); r++) {
		if(expected_reward.containsKey(r)) {
		    expected_reward.put(r, expected_reward.get(r) + getPathReward((STPG) model, c_sample, stpgRewards.get(r))*c_samples.get(c_sample));
		} else {
		    expected_reward.put(r, getPathReward((STPG) model, c_sample, stpgRewards.get(r))*c_samples.get(c_sample));	    
		}
	    }
	}
	System.out.println("MOST LIKELY PATH:");
	System.out.printf("%.4f: %s\n", c_samples.get(max_sample), max_sample.toString());

	System.out.println("REWARDS:");
	for(int r = 0; r < stpgRewards.size(); r++) {
	    System.out.printf("E_%d = %f\n", r, expected_reward.get(r));
	}
	
	System.out.println("COLLATED TERMINALS:");
	double cumulative = 0.0;
	Map<State,Double> c_terminals = collateTerminals(samples);
	for(State c_sample : c_terminals.keySet()) {
	    cumulative += c_terminals.get(c_sample);
	    System.out.printf("%.4f: %s\n", c_terminals.get(c_sample), c_sample.toString());
	}
	System.out.printf("%.4f cumulative\n", cumulative);
	
	System.out.println("COLLATED CAR POSITIONS:");
	int count = 1;
	Map<Integer,Double> c_positions = collateCarPositions(samples);
	for(Integer car_position : c_positions.keySet()) {
	    System.out.printf("%d: %.4f\t", car_position, c_positions.get(car_position));
	    if(count % 5 == 0) {
		System.out.printf("\n");
	    }
	    count++;
	}

    }


    private double getPathReward(STPG G, List<State> path, STPGRewards stpgReward)
    {
	double reward = 0.0;
	List<State> S = G.getStatesList();
	for(State s : path) {
	    if(s!=null) {
		reward += stpgReward.getStateReward(S.indexOf(s));
	    }
	}
	return reward;
    }


    private Map<Integer,Double> collateCarPositions(List<List<State>> samples) {
	Map<Integer,Double> result = new HashMap<Integer,Double>();
	for(List<State> sample : samples) {
	    for(State s : sample) {
	        int car_position = (Integer) s.varValues[1];
		if(result.containsKey(car_position)) {
		    result.put(car_position, result.get(car_position) + 1.0);
		} else {
		    result.put(car_position, 1.0);
		}
	    }
	}
	double total = 0.0;
	for(Integer car_position : result.keySet()) {
	    total += result.get(car_position);
	}
	for(Integer car_position : result.keySet()) {
	    result.put(car_position, result.get(car_position) / total);
	}

	return result;
    }


    private Map<List<State>,Double> collateSamples(List<List<State>> samples) {
	Map<List<State>,Double> result = new HashMap<List<State>,Double>();
	go_through_samples:
	for(List<State> sample : samples) {
	    // first, check if sample is alrady assigned a probability in the results
	    check_if_already_assigned:
	    for(List<State> c_sample : result.keySet()) {
		if(c_sample.size() != sample.size()) {
		    // definitely not equal
		    continue check_if_already_assigned;
		} else { // go through all states and check
		    for(int s = 0; s < c_sample.size() ; s++) {
			if(c_sample.get(s) != sample.get(s)) {
			    // definitely not equal
			    continue check_if_already_assigned;
			}
		    }
		    // fall through only if equal
		    result.put(c_sample, result.get(c_sample) + 1.0);
		    continue go_through_samples;
		}
	    }
	    // if not already assigned
	    result.put(sample, 1.0);
	}
	// divide by total number of samples
	double num_samples = samples.size();
	for(List<State> c_sample : result.keySet()) {
	    result.put(c_sample, result.get(c_sample) / num_samples);
	}
	return result;
    }


    private Map<State,Double> collateTerminals(List<List<State>> samples) {
	Map<State,Double> result = new HashMap<State,Double>();
	go_through_samples:
	for(List<State> sample : samples) {
	    check_if_already_assigned:
	    for(State c_sample : result.keySet()) {
		if(c_sample != sample.get(sample.size()-1)) {
		    continue check_if_already_assigned;
		} else {
		    result.put(c_sample, result.get(c_sample) + 1.0);
		    continue go_through_samples;
		}
	    }
	    // if not already assigned
	    result.put(sample.get(sample.size()-1), 1.0);
	}
	// divide by total number of samples
	double num_samples = samples.size();
	for(State c_sample : result.keySet()) {
	    result.put(c_sample, result.get(c_sample) / num_samples);
	}
	return result;
    }

    
}
