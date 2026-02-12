package explicit;

import java.util.BitSet;
import java.util.Map;

public class ModelCheckerResultMulti extends ModelCheckerResult
{
	// Constructor
	public ModelCheckerResultMulti() 
	{
		super();
	}
	
	public Map<BitSet, ModelCheckerResultMulti> subGames = null;
	public double[][] solnMulti = null;
	public BitSet D = null;
	public BitSet E = null;
	
}
