//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford)
//  * Gabriel Santos <gabriel.santos@cs.ox.ac.uk> (University of Oxford)
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
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.List;

import org.apache.commons.math3.fraction.BigFraction;

import com.microsoft.z3.*;
import com.microsoft.z3.Model;

import explicit.Distribution;
import org.apache.commons.math3.util.Precision;
import java.math.BigDecimal;

public class CSGLabeledPolytopesZ3 implements CSGLabeledPolytopes {
	
	private RealExpr[] payvars;
	private ArithExpr[] payoffs;

	private RealExpr[] vars;
	private RealExpr[] p1vars;
	private RealExpr[] p2vars;
	private IntExpr zero;
	private IntExpr one;
	private BoolExpr[] xlabels;
	private BoolExpr[] ylabels;
	private ArithExpr[] xexps;
	private ArithExpr[] yexps;
	private BoolExpr eq;
	private int nrows = 0;
	private int ncols = 0;
	private HashMap<String,ArrayList<Double>> eqs;
	private int neq = 0;
	private double[] p1p;
	private double[] p2p;
	private double[][] a;
	private double[][] b;
    
	private String[] lvp1;
    private String[] lvp2;
    
    private HashMap<String, String> cfg;
    private Context ctx;
    private Solver s;
    
    private ArithExpr cp1;
    private ArithExpr cp2;
    private double time;
        
    private int ndigits = 9;
    
    private ArrayList<ArrayList<Distribution>> strat;
    
    public CSGLabeledPolytopesZ3() {
    	cfg = new HashMap<String, String>();
        cfg.put("model", "true");
        ctx = new Context(cfg);
        s = ctx.mkSolver();        
    }
    
    public CSGLabeledPolytopesZ3(int nrows, int ncols, double[][] a, double[][] b) {
    	cfg = new HashMap<String, String>();
        cfg.put("model", "true");
        ctx = new Context(cfg);
        s = ctx.mkSolver();
        this.nrows = nrows;
        this.ncols = ncols;
        this.a = a;
        this.b = b;
        initialize();
    }
        
    public void update(int nrows, int ncols, double[][] a, double[][] b) {
        s.reset();
        this.nrows = nrows;
        this.ncols = ncols;
        this.a = a;
        this.b = b;
        /*
        System.out.println("-- matrix a");
        for (int row = 0; row < nrows; row++) {
            System.out.println("-- row " + row + " " + Arrays.toString(a[row]));
        }
        System.out.println("-- matrix b");
        for (int row = 0; row < nrows; row++) {
            System.out.println("-- row " + row + " " + Arrays.toString(b[row]));
        }
        */
        initialize();
    }
    
	private void initialize() {
		zero = ctx.mkInt(0);
	    one = ctx.mkInt(1);
		vars = new RealExpr[nrows+ncols];
		p1vars = new RealExpr[nrows];
		p2vars = new RealExpr[ncols];
		lvp1 = new String[nrows];
		lvp2 = new String[ncols];
		strat = new ArrayList<ArrayList<Distribution>>();
		int i = 0, j = 0;
		for(; i < nrows; i++) {
			vars[i] = ctx.mkRealConst("x" + i);
			p1vars[i] = vars[i];
			lvp1[i] = "x" + i;
		}
		for(; j < ncols; j++) {
			vars[i] = ctx.mkRealConst("y" + j);
			p2vars[j] = vars[i];
			lvp2[j] = "y" + j;
			i++;
		}
		for(RealExpr v : vars) {
			s.add(ctx.mkLe(v, one));
			s.add(ctx.mkGe(v, zero));
		}
	}
	
	@Override
	public void clear() {
		// TODO Auto-generated method stub
		
	}
	
	private void xLabels() {
		/*
		System.out.println("-- xexps");
		for (int l = 0; l < xexps.length; l++) {
			System.out.println(l + ": " + xexps[l]);
		}
		*/
		BoolExpr[] tmp = new BoolExpr[ncols-1];
		int l = 0;
		int j = 0;
		xlabels = new BoolExpr[nrows+ncols];
		for(int i = 0; i < nrows+ncols; i++) {
			if(i < nrows) {
				xlabels[i] = ctx.mkEq(vars[i], zero);
			}
			else {
				for(int k = 0; k < ncols; k++) {
					if(j != k) {
						tmp[l] = ctx.mkGe(xexps[j], xexps[k]);
						l++;
					}
				}
				xlabels[i] = tmp[0];
				if(tmp.length > 1) {
					for(int m = 1; m < tmp.length; m++) {
						xlabels[i] = ctx.mkAnd(xlabels[i], tmp[m]);
					}
				}
				j++;
				l = 0;
			}
		}
	}
	
	private void yLabels() {
		/*
		System.out.println("-- yexps");
		for (int l = 0; l < yexps.length; l++) {
			System.out.println(l + ": " + yexps[l]);
		}
		*/
		BoolExpr[] tmp = new BoolExpr[nrows-1];
		int l = 0;
		int j = 0;
		ylabels = new BoolExpr[nrows+ncols];
		for(int i = 0; i < nrows+ncols; i++) {
			if(i < nrows) {
				for(int k = 0; k < nrows; k++) {
					if(j != k) {
						tmp[l] = ctx.mkGe(yexps[j], yexps[k]);
						l++;
					}
				}
				ylabels[i] = tmp[0];
				if(tmp.length > 1) {
					for(int m = 1; m < tmp.length; m++) {
						ylabels[i] = ctx.mkAnd(ylabels[i], tmp[m]);
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
		yexps = new ArithExpr[nrows];
		xexps = new ArithExpr[ncols];
		ArithExpr curr;
		ArithExpr prod;
		for(int i = 0; i < nrows; i++) {
			curr = zero;
			for(int j = 0; j < ncols; j++) {
				try {
					curr = ctx.mkAdd(curr, ctx.mkMul(vars[nrows+j], getRealValue(a[i][j])));
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
					curr = ctx.mkAdd(curr, ctx.mkMul(vars[i], getRealValue(b[i][j])));
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
        Model model = null;                
        BoolExpr rstr;
        BoolExpr prstr;
        long start, end;
		int i,j;
        
        BoolExpr c1;
        BoolExpr c2;
        
        BoolExpr p1;
        BoolExpr p2; 
		
        //in order to check for indifference
        BoolExpr xeqrst = ctx.mkTrue();
        BoolExpr yeqrst = ctx.mkTrue();
        
        start = new Date().getTime();
        
		vMult();
		xLabels();
		yLabels();
		
		eq = ctx.mkTrue();
		for(i = 0; i < nrows+ncols; i++) {
			eq = ctx.mkAnd(eq, ctx.mkOr(xlabels[i], ylabels[i]));
		}
		
		Expr xctr = zero;
		Expr yctr = zero;
		
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
	
        j = 0; 

        strat = new ArrayList<ArrayList<Distribution>>();
        eqs = new HashMap<String,ArrayList<Double>> ();
        
        /*
    	System.out.println("-- xlabels");
        for (int l = 0; l < xlabels.length; l++) {
        	System.out.println(l + ": " + xlabels[l]);
        }
        
    	System.out.println("-- ylabels");
        for (int l = 0; l < ylabels.length; l++) {
        	System.out.println(l + ": " + ylabels[l]);        	
        }
        */
        
        /*
        for(int nv = 0; nv < p1vars.length - 1; nv++) {
        	xeqrst = ctx.mkAnd(xeqrst, ctx.mkEq(p1vars[nv], p1vars[nv+1]));
        }
        
        for(int nv = 0; nv < p2vars.length - 1; nv++) {
        	yeqrst = ctx.mkAnd(yeqrst, ctx.mkEq(p2vars[nv], p2vars[nv+1]));
        }
        */
       
        /*** Too strong? ***/
       
       /* 
        for(int nv = 0; nv < p1vars.length; nv++) {
            for(int nnv = 0; nnv < p1vars.length; nnv++) {
            	xeqrst = ctx.mkAnd(xeqrst, 
            					   ctx.mkOr(
            							   	ctx.mkNot(ctx.mkEq(p1vars[nv], p1vars[nnv])),
            							   	ctx.mkAnd(ctx.mkEq(p1vars[nv], zero),
            							   			  ctx.mkEq(p1vars[nnv], zero))	
            							   )
            						);
            }
        }
        
        for(int nv = 0; nv < p2vars.length; nv++) {
            for(int nnv = 0; nnv < p2vars.length; nnv++) {
            	xeqrst = ctx.mkAnd(xeqrst, 
            					   ctx.mkOr(
            							    ctx.mkNot(ctx.mkEq(p2vars[nv], p2vars[nnv])),
            							    ctx.mkAnd(ctx.mkEq(p2vars[nv], zero),
            							    		  ctx.mkEq(p2vars[nnv], zero))	
            							   )
            						);
            }
        }
        */
        
        /*** Too strong? ***/
        
        xeqrst = ctx.mkNot(xeqrst);
        yeqrst = ctx.mkNot(yeqrst);
        
        boolean indif = false;
        
        addPayoffs();
        
        /*
        Optimize opt = ctx.mkOptimize();
        opt.Add(ctx.mkEq(payoffs[0],  payvars[0]));
        opt.Add(ctx.mkEq(payoffs[1],  payvars[1]));
        opt.MkMaximize(payvars[0]);
        Optimize.Handle mx = opt.MkMaximize(payvars[0]);
        Optimize.Handle my = opt.MkMaximize(payvars[1]);
        */
        
        while (Status.SATISFIABLE == s.check()) {
            model = s.getModel();
            rstr = ctx.mkFalse();

            c1 = ctx.mkTrue();
            c2 = ctx.mkTrue();
                     
            /*
            if(j > 0) {          	
            	if(((BoolExpr) model.eval(xeqrst, true)).isFalse()) {
            		System.out.println("## Player 1 is indifferent. ##");
            		indif = true;
            		//break;
            	}
            	if(((BoolExpr) model.eval(yeqrst, true)).isFalse()) {
            		System.out.println("## Player 2 is indifferent. ##");
            		indif = true;
            		//break;	
            	}      	            	
            }
            */
            
            //System.out.println(getDoubleValue(model, payvars[0]));
            //System.out.println(getDoubleValue(model, payvars[1]));
            //compPayoffs(model);
                        
            dists = new ArrayList<Distribution>();
            dist1 = new Distribution();
            dist2 = new Distribution();
            
            for(i = 0; i < vars.length; i++) {
            		p = getDoubleValue(model, vars[i]);
            		if(p > 0) {
            			if(i < nrows)
            				dist1.add(i, p);
            			else 
            				dist2.add(i - nrows, p);
            		}
            		if(j == 0) {
            			eqs.put(vars[i].getSExpr(), new ArrayList<Double> ());
            			eqs.get(vars[i].getSExpr()).add(p);
            		}
            		else {
            			//System.out.println((RatNum) model.getConstInterp(vars[i]));
            			eqs.get(vars[i].getSExpr().toString()).add(p);	
            		}
            	//rstr = ctx.mkOr(rstr, ctx.mkNot(ctx.mkEq(vars[i], (RealExpr) model.getConstInterp(vars[i]))));
            }
            
            for(i = 0; i < p1vars.length; i++) {
            		if(((BoolExpr) model.eval(ctx.mkEq(p1vars[i], zero), true)).isFalse()) {
            			c1 = ctx.mkAnd(c1, ctx.mkNot(ctx.mkEq(p1vars[i], zero)));
            		}
            		else {
            			c1 = ctx.mkAnd(c1, ctx.mkEq(p1vars[i], zero));
            		}
            }
            
            for(i = 0; i < p2vars.length; i++) {
            		if(((BoolExpr) model.eval(ctx.mkEq(p2vars[i], zero), true)).isFalse()) {
            			c2 = ctx.mkAnd(c2, ctx.mkNot(ctx.mkEq(p2vars[i], zero)));
            		}
            		else {
            			c2 = ctx.mkAnd(c2, ctx.mkEq(p2vars[i], zero));
            		}
            }            
            
            /*
            p1 = ctx.mkAnd(c2, ctx.mkEq(payvars[0], model.getConstInterp(payvars[0])));
            p2 = ctx.mkAnd(c1, ctx.mkEq(payvars[1], model.getConstInterp(payvars[1])));
            */

            /*
        	prstr = ctx.mkAnd(ctx.mkEq(payvars[0], model.getConstInterp(payvars[0])), ctx.mkEq(payvars[1], model.getConstInterp(payvars[1])));
            */
            
            //System.out.println(opt.Check()); 
            //System.out.println(mx.getValue());
            //System.out.println(my.getValue());
            
            /*
			System.out.println("--");
			for(String v : eqs.keySet()) {
				System.out.println(v + " " + eqs.get(v).get(j));
			}
			System.out.println("--");
            */
            
            //System.out.println("$$ dist1 " + dist1);
            //System.out.println("$$ dist2 " + dist2);
            
            
            dists.add(0,dist1);
            dists.add(1,dist2);
            strat.add(j,dists);
    		
            j++;
    		s.add(ctx.mkOr(ctx.mkNot(c1), ctx.mkNot(c2)));
    		
    		//s.add(rstr);
            //s.add(ctx.mkNot(p2));
            //s.add(ctx.mkNot(p1));
        }
        end = new Date().getTime();
        time = (end - start)/1000.000;
        //System.out.println("Time: " + time); 
        if(!indif) {
        	//System.out.println(j + " equilibrium(a)");
        }
        else {
        	//System.out.println(j-1 + " equilibrium(a)");
        	//System.out.println("## One of the players is indifferent. ##");
        }
        neq = j;
	}

	public void addPayoffs() {
		payvars = new RealExpr[2];
		payoffs = new ArithExpr[2];
		payvars[0] = ctx.mkRealConst("pp1");
		payvars[1] = ctx.mkRealConst("pp2");
		payoffs[0] = ctx.mkInt(0);
		payoffs[1] = ctx.mkInt(0);
		for(int i = 0; i < nrows; i++) {
			for(int j = 0; j < ncols; j++) {
				payoffs[0] = ctx.mkAdd(payoffs[0], ctx.mkMul(p1vars[i], p2vars[j], getRealValue(a[i][j])));
				payoffs[1] = ctx.mkAdd(payoffs[1], ctx.mkMul(p1vars[i], p2vars[j], getRealValue(b[i][j])));
			}
		}
		s.add(ctx.mkEq(payvars[0], payoffs[0]));
		s.add(ctx.mkEq(payvars[1], payoffs[1]));
	}	

	public void compPayoffs(Model model) {	
		double p1p = 0.0;
		double p2p = 0.0;
        for(int i = 0; i < nrows; i++) {
        	for(int j = 0; j < ncols; j++) {
        		p1p += getDoubleValue(model, p1vars[i])*
        			   getDoubleValue(model, p2vars[j])*a[i][j];
        		p2p += getDoubleValue(model, p1vars[i])*
        			   getDoubleValue(model, p2vars[j])*b[i][j];
        	}
        }	
        cp1 = getRealValue(p1p);
        cp2 = getRealValue(p2p);
	}
	
	public void compPayoffs() {	
		p1p = new double[neq];
		p2p = new double[neq];
		Arrays.fill(p1p, 0.0);
		Arrays.fill(p2p, 0.0);
		for(int e = 0; e < neq; e++) {
        	for(int i = 0; i < nrows; i++) {
        		for(int j = 0; j < ncols; j++) {
           			//p1p[e] += eqs.get(lvp1[i].toString()).get(e) * eqs.get(lvp2[j].toString()).get(e) * Precision.round(a[i][j], 9, BigDecimal.ROUND_HALF_EVEN);
        			//p2p[e] += eqs.get(lvp1[i].toString()).get(e) * eqs.get(lvp2[j].toString()).get(e) * Precision.round(b[i][j], 9, BigDecimal.ROUND_HALF_EVEN);
        			p1p[e] += eqs.get(p1vars[i].toString()).get(e) * eqs.get(p2vars[j].toString()).get(e) * a[i][j];
        			p2p[e] += eqs.get(p1vars[i].toString()).get(e) * eqs.get(p2vars[j].toString()).get(e) * b[i][j];
        		}
        	}
        }		
	}
	
	public void print() {
			for(int i = 0; i < neq; i++) {
				for(String v : eqs.keySet()) {
					System.out.println(v + " " + eqs.get(v).get(i));
				}
				System.out.println("p1 " + p1p[i]);
				System.out.println("p2 " + p2p[i]);
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
			v1 = v2.toUpper(9);
			return (Double) (v1.getBigIntNumerator().doubleValue() / v1.getBigIntDenominator().doubleValue());
		}
		else
			return Double.NaN;
	}
	
	public ArithExpr getRealValue(double v) {
		BigFraction fr;
		ArithExpr result = null;
		result = ctx.mkReal((long) (Math.pow(10, ndigits) * v));
		result = ctx.mkDiv(result, ctx.mkReal((long) (Math.pow(10, ndigits))));
		/*
		try {
			fr = new BigFraction(v, 0.0000000000001D, 10000);
			result = ctx.mkReal(fr.getNumeratorAsLong());
			result = ctx.mkDiv(result, ctx.mkReal(fr.getDenominatorAsLong()));	
		}
		catch (Exception e) {
			e.printStackTrace();
		}
		*/
		return result;
	}
	
	public RealExpr[] getP1vars() {
		return p1vars;
	}

	public RealExpr[] getP2vars() {
		return p2vars;
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

	public ArrayList<ArrayList<Distribution>> getStrat() {
		return strat;
	}
	
	public HashMap<String, ArrayList<Double>> getEqs() {
		return eqs;
	}

	public String[] getLvp1() {
		return lvp1;
	}

	public String[] getLvp2() {
		return lvp2;
	}
	
}
