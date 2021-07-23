package explicit;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import com.sri.yices.Config;
import com.sri.yices.Context;
import com.sri.yices.Status;
import com.sri.yices.Terms;
import com.sri.yices.Version;
import com.sri.yices.Yices;

import prism.PrismException;

public class CSGLabeledPolytopesYicesStack implements CSGLabeledPolytopes
{
	private String solverName = "Yices";

	private Context ctx;
	
	private int realType;
	private int boolType;
	
	private int zero;
	private int one;
	private int ytrue;
	private int yfalse;
	
	private int eq;
	
	private int curr;
	
	private int[] xlabels;
	private int[] ylabels;
		
	private int[] vars;
	private int[] xexps;
	private int[] yexps;
	
	private int[] tmpc;
	private int[] tmpr;
	
	private com.sri.yices.Model model;
	
	private int xctr;
	private int yctr;
	private int c1;
	private int c2; 
	
	private String[] lvp1;
	private String[] lvp2;
	
	private double[][] a;
	private double[][] b;

	private double[] p1p;
	private double[] p2p;
		
	private int neq;
	private int nrows;
	private int ncols;
	
	private HashMap<String,ArrayList<Double>> eqs;
    private ArrayList<ArrayList<Distribution>> strat;
	
	protected String eqPolicy;
    
    public CSGLabeledPolytopesYicesStack() {
    	
    }
    
    public CSGLabeledPolytopesYicesStack(int nrows, int ncols) throws PrismException
    {
    	initSolver();
        eqs = new HashMap<String,ArrayList<Double>>();
        realType = Yices.realType();
        boolType = Yices.boolType();
		zero = Terms.intConst(0);
		one = Terms.intConst(1);
		vars = new int[nrows+ncols];
		lvp1 = new String[nrows];
		lvp2 = new String[ncols];
		tmpc = new int[ncols-1];
		tmpr = new int[nrows-1];
		xlabels = new int[nrows+ncols];
		ylabels = new int[nrows+ncols];
		yexps = new int[nrows];
		xexps = new int[ncols];
		xctr = zero;
		yctr = zero;
		ytrue = Yices.mkTrue();
		yfalse = Yices.mkFalse();
		int i = 0, j = 0;
		for(; i < nrows; i++) {
			vars[i] = Terms.newUninterpretedTerm(realType);
			lvp1[i] = "x_" + i;
		}
		for(; j < ncols; j++) {
			vars[i] = Terms.newUninterpretedTerm(realType);
			lvp2[j] = "y_" + j;
			i++;
		}
		for(int v : vars) {
			ctx.assertFormula(Terms.arithLeq(v, one));
			ctx.assertFormula(Terms.arithGeq(v, zero));
		}
    }
    
    public void clear() throws PrismException
    {
    	initSolver();
        eqs = new HashMap<String,ArrayList<Double>>();
        realType = Yices.realType();
        boolType = Yices.boolType();
		zero = Terms.intConst(0);
		one = Terms.intConst(1);
		vars = new int[nrows+ncols];
		lvp1 = new String[nrows];
		lvp2 = new String[ncols];
		tmpc = new int[ncols-1];
		tmpr = new int[nrows-1];
		xlabels = new int[nrows+ncols];
		ylabels = new int[nrows+ncols];
		yexps = new int[nrows];
		xexps = new int[ncols];
		xctr = zero;
		yctr = zero;
		ytrue = Yices.mkTrue();
		yfalse = Yices.mkFalse();
		int i = 0, j = 0;
		for(; i < nrows; i++) {
			vars[i] = Terms.newUninterpretedTerm(realType);
			lvp1[i] = "x_" + i;
		}
		for(; j < ncols; j++) {
			vars[i] = Terms.newUninterpretedTerm(realType);
			lvp2[j] = "y_" + j;
			i++;
		}
		for(int v : vars) {
			ctx.assertFormula(Terms.arithLeq(v, one));
			ctx.assertFormula(Terms.arithGeq(v, zero));
		}
    }
    
    /**
     * Initialise the solver
     */
    private void initSolver() throws PrismException
    {
        try {
        	Config cfg = new Config("QF_LRA");
    		cfg.set("mode", "push-pop");
    		ctx = new Context(cfg);
    		System.out.println("NEW");
    		cfg.close();
        	solverName = "Yices " + Version.versionString;
        } catch (UnsatisfiedLinkError e) {
        	throw new PrismException("Could not initialise Yices: " + e.getMessage());
        }
    }
    
    @Override
    public String getSolverName()
    {
    	return solverName;
    }
    
	public void update(int nrows, int ncols, double[][] a, double[][] b) {
		this.nrows = nrows;
		this.ncols = ncols;
		this.a = a;
		this.b = b;
	}
	
	private void xLabels() {
		int l = 0;
		int j = 0;
		for(int i = 0; i < nrows+ncols; i++) {
			if(i < nrows) {
				xlabels[i] = Terms.arithEq(vars[i], zero);
			}
			else {
				for(int k = 0; k < ncols; k++) {
					if(j != k) {
						tmpc[l] = Terms.arithGeq(xexps[j], xexps[k]);
						l++;
					}
				}
				xlabels[i] = tmpc[0];
				if(ncols-1 > 1) {
					for(int m = 1; m < ncols-1; m++) {
						xlabels[i] = Terms.and(xlabels[i], tmpc[m]);
					}
				}
				j++;
				l = 0;
			}
		}
	}
	
	private void yLabels() {
		int l = 0;
		int j = 0;
		for(int i = 0; i < nrows+ncols; i++) {
			if(i < nrows) {
				for(int k = 0; k < nrows; k++) {
					if(j != k) {
						tmpr[l] = Terms.arithGeq(yexps[j], yexps[k]);
						l++;
					}
				}
				ylabels[i] = tmpr[0];
				if(nrows-1 > 1) {
					for(int m = 1; m < nrows-1; m++) {
						ylabels[i] = Terms.and(ylabels[i], tmpr[m]);
					}
				}
				j++;
				l = 0;
			}
			else {
				ylabels[i] = Terms.arithEq(vars[i], zero);
			}
		}
	}
	
	private void vMult() {
		for(int i = 0; i < nrows; i++) {
			curr = zero;
			for(int j = 0; j < ncols; j++) {
				try {
					curr = Terms.add(curr, Terms.mul(vars[nrows+j], 
													Terms.parseFloat(String.valueOf(a[i][j]))));
				}
				catch(Exception e) {
					e.printStackTrace();
				}
			}
			yexps[i] = curr;
		}
		for(int j = 0; j < ncols; j++) {
			curr = zero;
			for(int i = 0; i < nrows; i++) {
				try {		
					curr = Terms.add(curr, Terms.mul(vars[i], 
													Terms.parseFloat(String.valueOf(b[i][j]))));
				}
				catch(Exception e) {
					e.printStackTrace();
				}
			}
			xexps[j] = curr;
		}
	}

	public void compEq() throws PrismException
	{
		ArrayList<Distribution> dists;
		Distribution dist1;
		Distribution dist2;
		double p;
		int i, j;	
		clear();
		vMult();
		xLabels();
		yLabels();
		eq = ytrue;
		xctr = zero;
		yctr = zero;
		for(i = 0; i < nrows+ncols; i++) {
			eq = Terms.and(eq, Terms.or(xlabels[i],ylabels[i]));
		}
		for(i = 0; i < nrows; i++) {
			xctr = Terms.add(xctr, vars[i]);
		}
		xctr = Terms.arithEq(xctr, one);
		for(j = i; j < nrows+ncols; j++) {
			yctr = Terms.add(yctr, vars[j]);
		}
		yctr = Terms.arithEq(yctr, one);			
		ctx.assertFormula(xctr);
		ctx.assertFormula(yctr);
		ctx.assertFormula(eq);
        strat = new ArrayList<ArrayList<Distribution>>();
		eqs.clear();
        j = 0; 
        while (ctx.check() == Status.SAT) {
        	model = ctx.getModel();        
            c1 = ytrue;
            c2 = ytrue;
            dists = new ArrayList<Distribution>();
            dist1 = new Distribution();
            dist2 = new Distribution();
    		//System.out.println("---");
            //String prt = "(";
            for (i = 0; i < nrows+ncols; i++) {            	
        		//p = model.bigRationalValue(vars[i]).doubleValue();
            	p = model.doubleValue(vars[i]);
        		//prt += p;
        		//if (i < vars.length - 1)
        			//prt += ",";
        		//System.out.println(p);
        		if(p > 0) {
        			if(i < nrows)
        				dist1.add(i, p);
        			else 
        				dist2.add(i - nrows, p);
        		}
                if (j == 0) {
        			if (i < nrows) {
                		eqs.put(lvp1[i], new ArrayList<Double>());
                		eqs.get(lvp1[i]).add(p);
                	} 
        			else {
        	       		eqs.put(lvp2[i-nrows], new ArrayList<Double>());
        	       		eqs.get(lvp2[i-nrows]).add(p);
                	}
                }
                else {
        			if (i < nrows) {
                		eqs.get(lvp1[i]).add(p);
                	} 
        			else {
        				eqs.get(lvp2[i-nrows]).add(p);		
               		}
       			}
            }
            //prt += ")";
            //System.out.println(prt);
            if (eqPolicy == "first_found") {
                dists.add(0, dist1);
                dists.add(1, dist2);
                strat.add(j, dists);
        		j++;
        		break;
            }
            for(i = 0; i < nrows+ncols; i++) {
				if (i < nrows) {
        			if (Double.compare(eqs.get(lvp1[i]).get(j), 0.0) != 0) {
        				c1 = Terms.and(c1, Terms.not(Terms.arithEq(vars[i], zero)));
        			}
            		else {
            			c1 = Terms.and(c1, Terms.arithEq(vars[i], zero));
            		}
				} 
				else {
    				if (Double.compare(eqs.get(lvp2[i-nrows]).get(j), 0.0) != 0) {
    					c2 = Terms.and(c2, Terms.not(Terms.arithEq(vars[i], zero)));
            		}
            		else {
            			c2 = Terms.and(c2, Terms.arithEq(vars[i], zero));
            		}
				}
            }
            dists.add(0, dist1);
            dists.add(1, dist2);
            strat.add(j, dists);
    		j++;
            ctx.assertFormula(Terms.or(Terms.not(c1), Terms.not(c2)));
        }
        model.close();
        neq = j;
        Yices.reset();
        compPayoffs();
	}

	public void compPayoffs() {
		p1p = new double[neq];
		p2p = new double[neq];
		Arrays.fill(p1p, 0.0);
		Arrays.fill(p2p, 0.0);
		for(int e = 0; e < neq; e++) {
        	for(int i = 0; i < nrows; i++) {
        		for(int j = 0; j < ncols; j++) {
        			p1p[e] += eqs.get(lvp1[i]).get(e) * eqs.get(lvp2[j]).get(e) * a[i][j];
        			p2p[e] += eqs.get(lvp1[i]).get(e) * eqs.get(lvp2[j]).get(e) * b[i][j];
        		}
        	}
        }			
	}

	public ArrayList<ArrayList<Distribution>> getStrat() {
		return strat;
	}

	public double[] getP1p() {
		return p1p;
	}

	public double[] getP2p() {
		return p2p;
	}
	
	public int getNeq() {
		return neq;
	}

	public void setPolicy(String s) {
		this.eqPolicy = s;
	}

}
