package explicit;

import java.util.ArrayList;

public interface CSGLabeledPolytopes {
	
	public void update(int nrows, int ncols, double[][] a, double[][] b);
	
	public void compEq();
	
	public void compPayoffs();
	
	public ArrayList<ArrayList<Distribution>> getStrat();
	
	public double[] getP1p(); 
	
	public double[] getP2p(); 
	
	public int getNeq();

	public void clear();
	
}
