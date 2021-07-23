package explicit;

import java.util.ArrayList;

import prism.PrismException;

public interface CSGLabeledPolytopes {
	
    public String getSolverName();
    
	public void update(int nrows, int ncols, double[][] a, double[][] b);
	
	public void compEq() throws PrismException;
	
	public void compPayoffs();
	
	public ArrayList<ArrayList<Distribution>> getStrat();
	
	public double[] getP1p(); 
	
	public double[] getP2p(); 
	
	public int getNeq();
	
}
