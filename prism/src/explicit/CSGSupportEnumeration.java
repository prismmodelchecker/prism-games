package explicit;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Map;

import explicit.CSGModelCheckerEquilibria.CSGResultStatus;
import prism.Pair;

public interface CSGSupportEnumeration {
	
	public void init();
	
	public void setNumPlayers(int n);
	
	public void setIndexes(ArrayList<ArrayList<Integer>> a);
	
	public Pair<CSGResultStatus, ArrayList<Double>> computeEquilibria(BitSet supp, HashMap<Integer, int[]> map, int s);
	
	public void computeConstraints(BitSet supp);
	
	public void computeSupport(BitSet supp, HashMap<Integer, int[]> map);

	public void translateAssertions(HashMap<Integer, HashMap<Integer, ArrayList<Pair<BitSet, Double>>>> assertionsIdx, HashMap<Integer, int[]> map);	

	public void setGradient(HashMap<Integer, HashMap<Integer, ArrayList<Pair<BitSet, Double>>>> gradient);
	
	public void setAssertions(HashMap<Integer, HashMap<Integer, ArrayList<Pair<BitSet, Double>>>> assertions);
	
	public ArrayList<Distribution> getStrat();

}
