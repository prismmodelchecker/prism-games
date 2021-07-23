package explicit;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import com.microsoft.z3.AlgebraicNum;
import com.microsoft.z3.ArithExpr;
import com.microsoft.z3.BoolExpr;
import com.microsoft.z3.Context;
import com.microsoft.z3.Expr;
import com.microsoft.z3.IntExpr;
import com.microsoft.z3.Model;
import com.microsoft.z3.RatNum;
import com.microsoft.z3.RealExpr;
import com.microsoft.z3.Solver;
import com.microsoft.z3.Status;
import com.microsoft.z3.Version;

import prism.PrismException;

public class CSGLabeledPolytopesZ3Stack implements CSGLabeledPolytopes
{
	private String solverName = "Z3";
	
	private RealExpr[] payvars;
	private ArithExpr[] payoffs;
	private RealExpr[] vars;
	private RealExpr[] p1vars;
	private RealExpr[] p2vars;
	private IntExpr zero;
	private IntExpr one;
	private BoolExpr ytrue;
	private BoolExpr yfalse;
	private BoolExpr[] xlabels;
	private BoolExpr[] ylabels;
	private ArithExpr[] xexps;
	private ArithExpr[] yexps;
	private ArithExpr curr;

	private BoolExpr eq;
	private int nrows = 0;
	private int ncols = 0;
	
	private BoolExpr[] tmpc;
	private BoolExpr[] tmpr;
	
	private Model model;
	
	private Expr xctr;
	private Expr yctr;
	private BoolExpr c1;
	private BoolExpr c2; 
	
	private int neq = 0;
	private double[] p1p;
	private double[] p2p;
	private double[][] a;
	private double[][] b;
	
	private String[] lvp1;
    private String[] lvp2;
    
    private Context ctx;
    private Solver s;
	
	private HashMap<String,ArrayList<Double>> eqs;
    private ArrayList<ArrayList<Distribution>> strat;
    
    public CSGLabeledPolytopesZ3Stack(int nrows, int ncols) throws PrismException
    {
    	initSolver();
        s = ctx.mkSolver(); 
        eqs = new HashMap<String,ArrayList<Double>>();
		zero = ctx.mkInt(0);
	    one = ctx.mkInt(1);
		vars = new RealExpr[nrows+ncols];
		lvp1 = new String[nrows];
		lvp2 = new String[ncols];
		tmpc = new BoolExpr[ncols-1];
		tmpr = new BoolExpr[nrows-1];
		xlabels = new BoolExpr[nrows+ncols];
		ylabels = new BoolExpr[nrows+ncols];
		yexps = new ArithExpr[nrows];
		xexps = new ArithExpr[ncols];
		xctr = zero;
		yctr = zero;
		ytrue = ctx.mkTrue();
		yfalse = ctx.mkFalse();
		int i = 0, j = 0;
		for(; i < nrows; i++) {
			vars[i] = ctx.mkRealConst("x" + i);
			lvp1[i] = "x_" + i;
		}
		for(; j < ncols; j++) {
			vars[i] = ctx.mkRealConst("y" + j);
			lvp2[j] = "y_" + j;
			i++;
		}
		for(RealExpr v : vars) {
			s.add(ctx.mkLe(v, one));
			s.add(ctx.mkGe(v, zero));
		}
    }
    
    /**
     * Initialise the solver
     */
    private void initSolver() throws PrismException
    {
    	HashMap<String, String> cfg = new HashMap<String, String>();
        cfg.put("model", "true");
        cfg.put("auto_config", "true");
        try {
        	ctx = new Context(cfg);
        	solverName = Version.getFullVersion();
        } catch (UnsatisfiedLinkError e) {
        	throw new PrismException("Could not initialise Z3: " + e.getMessage());
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
				xlabels[i] = ctx.mkEq(vars[i], zero);
			}
			else {
				for(int k = 0; k < ncols; k++) {
					if(j != k) {
						tmpc[l] = ctx.mkGe(xexps[j], xexps[k]);
						l++;
					}
				}
				xlabels[i] = tmpc[0];
				if(ncols-1 > 1) {
					for(int m = 1; m < ncols-1; m++) {
						xlabels[i] = ctx.mkAnd(xlabels[i], tmpc[m]);
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
						tmpr[l] = ctx.mkGe(yexps[j], yexps[k]);
						l++;
					}
				}
				ylabels[i] = tmpr[0];
				if(nrows-1 > 1) {
					for(int m = 1; m < nrows-1; m++) {
						ylabels[i] = ctx.mkAnd(ylabels[i], tmpr[m]);
					}
				}
				j++;
				l = 0;
			}
			else {
				ylabels[i] = ctx.mkEq(vars[i], zero);
			}
		}
	}
    
	private void vMult() {
		for(int i = 0; i < nrows; i++) {
			curr = zero;
			for(int j = 0; j < ncols; j++) {
				try {
					curr = ctx.mkAdd(curr, ctx.mkMul(vars[nrows+j], ctx.mkReal(String.valueOf(a[i][j]))));
				}
				catch(Exception e) {
					System.out.println(a[i][j]);
					System.out.println(String.valueOf(a[i][j]));
					e.printStackTrace();
				}
			}
			yexps[i] = curr;
		}
		for(int j = 0; j < ncols; j++) {
			curr = zero;
			for(int i = 0; i < nrows; i++) {
				try {
					curr = ctx.mkAdd(curr, ctx.mkMul(vars[i], ctx.mkReal(String.valueOf(b[i][j]))));
				}
				catch(Exception e) {
					System.out.println(a[i][j]);
					System.out.println(String.valueOf(a[i][j]));
					e.printStackTrace();
				}
			}
			xexps[j] = curr;
		}
	}
	
	public void compEq() {
		ArrayList<Distribution> dists;
		Distribution dist1;
		Distribution dist2;
		double p;
		int i, j;		
		vMult();
		xLabels();
		yLabels();
		s.push();
		eq = ytrue;
		xctr = zero;
		yctr = zero;
		for(i = 0; i < nrows+ncols; i++) {
			/*
			if (st == 1)
				System.out.println(xlabels[i]);
			*/
			eq = ctx.mkAnd(eq, ctx.mkOr(xlabels[i], ylabels[i]));
		}
		for(i = 0; i < nrows; i++) {
			xctr = ctx.mkAdd((ArithExpr) xctr, vars[i]);
		}
		xctr = ctx.mkEq(xctr, one);
		for(j = i; j < nrows+ncols; j++) {
			yctr = ctx.mkAdd((ArithExpr) yctr, vars[j]);
		}
		yctr = ctx.mkEq(yctr, one);			
		s.add((BoolExpr) xctr);
		s.add((BoolExpr) yctr);
		s.add(eq);
        strat = new ArrayList<ArrayList<Distribution>>();
		eqs.clear();
        j = 0; 
        while (Status.SATISFIABLE == s.check()) {
            model = s.getModel();
            c1 = ytrue;
            c2 = ytrue;
            dists = new ArrayList<Distribution>();
            dist1 = new Distribution();
            dist2 = new Distribution();
    		//System.out.println("---");
            for (i = 0; i < nrows+ncols; i++) {            	
        		//p = model.bigRationalValue(vars[i]).doubleValue();
            	p = getDoubleValue(model, vars[i]);
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
            for(i = 0; i < nrows+ncols; i++) {
				if (i < nrows) {
        			//if (Double.compare(model.bigRationalValue(vars[i]).doubleValue(), 0.0) != 0) {
        			if (Double.compare(eqs.get(lvp1[i]).get(j), 0.0) != 0) {
        				c1 = ctx.mkAnd(c1, ctx.mkNot(ctx.mkEq(vars[i], zero)));
        			}
            		else {
            			c1 = ctx.mkAnd(c1, ctx.mkEq(vars[i], zero));
            		}
				} 
				else {
    				//if (Double.compare(model.bigRationalValue(vars[i]).doubleValue(), 0.0) != 0) {
    				if (Double.compare(eqs.get(lvp2[i-nrows]).get(j), 0.0) != 0) {
    					c2 = ctx.mkAnd(c2, ctx.mkNot(ctx.mkEq(vars[i], zero)));
            		}
            		else {
            			c2 = ctx.mkAnd(c2, ctx.mkEq(vars[i], zero));
            		}
				}
            }
            dists.add(0, dist1);
            dists.add(1, dist2);
            strat.add(j, dists);
    		j++;
            s.add(ctx.mkOr(ctx.mkNot(c1), ctx.mkNot(c2)));
        }
		//System.out.println(eqs);
        s.pop();
        //YicesWrapper.garbage_collect();
        //YicesWrapper.yices_exit();
        neq = j;
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

	public double getDoubleValue(Model model, Expr expr) {
		RatNum v1;
		AlgebraicNum v2;
		if(model.getConstInterp(expr) instanceof RatNum) {
			v1 = (RatNum) model.getConstInterp(expr);
			return (Double) (v1.getBigIntNumerator().doubleValue() / v1.getBigIntDenominator().doubleValue());
		}
		else if (model.getConstInterp(expr) instanceof AlgebraicNum) {
			v2 = (AlgebraicNum) model.getConstInterp(expr);
			v1 = v2.toUpper(12);
			return (Double) (v1.getBigIntNumerator().doubleValue() / v1.getBigIntDenominator().doubleValue());
		}
		else
			return Double.NaN;
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

	public void clear() {
		
	}
	
}
