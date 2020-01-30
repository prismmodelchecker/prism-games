package explicit;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;

import com.microsoft.z3.AlgebraicNum;
import com.microsoft.z3.ArithExpr;
import com.microsoft.z3.BoolExpr;
import com.microsoft.z3.Context;
import com.microsoft.z3.Expr;
import com.microsoft.z3.IntNum;
import com.microsoft.z3.Model;
import com.microsoft.z3.Params;
import com.microsoft.z3.RatNum;
import com.microsoft.z3.Solver;
import com.microsoft.z3.Status;

import explicit.CSGModelCheckerEquilibria.CSGResultStatus;
import prism.Pair;

public class CSGSupportEnumerationZ3 implements CSGSupportEnumeration {

	private HashMap<String, String> cfg;
	private Context ctx;
    private Solver s;
    
	private ArrayList<ArithExpr> payoffs;	
	private ArrayList<ArrayList<ArithExpr>> strategies;
	private HashMap<Integer, HashMap<Integer, ArithExpr>> assertions;
	private ArrayList<ArrayList<Integer>> indexes;

	private IntNum zero;
	private IntNum one;
	
	private BitSet players;
	private BoolExpr support;
	
	private int numPlayers;
	private Params params;
	
    private ArrayList<Distribution> strat;
    
    private int ndigits = 9;
	
	public CSGSupportEnumerationZ3() {
    	cfg = new HashMap<String, String>();
        cfg.put("model", "true");
        ctx = new Context(cfg);
        s = ctx.mkSolver();        
        
        params = ctx.mkParams();
		params.add("timeout", 20);
		s.setParameters(params);
		
		assertions = new HashMap<Integer, HashMap<Integer, ArithExpr>>();
		strategies = new ArrayList<ArrayList<ArithExpr>>();
		payoffs = new ArrayList<ArithExpr>();
		players = new BitSet();
		strat = new ArrayList<Distribution>();
	}
	
	public void init() {
		s.reset();
		s.setParameters(params);
		
		players.clear();
		players.set(0, numPlayers);
		strategies.clear();
		payoffs.clear();
		zero = ctx.mkInt(0);
		one = ctx.mkInt(1);
		ArithExpr cnstr1;
		for(int p = 0; p < numPlayers; p++) {
			cnstr1 = zero;			
			strategies.add(p, new ArrayList<ArithExpr>());
			payoffs.add(p, ctx.mkRealConst("v_" + p));
			for(int a = 0; a < indexes.get(p).size(); a++) {
				strategies.get(p).add(a, ctx.mkRealConst("p_" + p + "_" + a));
				s.add(ctx.mkGe(strategies.get(p).get(a), zero));
				s.add(ctx.mkLe(strategies.get(p).get(a), one));
				cnstr1 = ctx.mkAdd(cnstr1, strategies.get(p).get(a));
			}
			s.add(ctx.mkEq(cnstr1, one));
		}
		//System.out.println(numPlayers);
		//System.out.println(indexes);
		//System.out.println(strategies);
	}
	
	public void computeConstraints(BitSet supp) {
		for(int p = 0; p < indexes.size(); p++) {
			for(int a = 0; a < indexes.get(p).size(); a++) {
				if(supp.get(indexes.get(p).get(a))) {
					s.add(ctx.mkEq(assertions.get(p).get(a), payoffs.get(p)));
				}
				else {
					s.add(ctx.mkLe(assertions.get(p).get(a), payoffs.get(p)));
				}
			}
		}
	}

	public void computeSupport(BitSet supp, HashMap<Integer, int[]> map) {
		support = ctx.mkTrue();
		for(int p = 0; p < indexes.size(); p++) {
			for(int a = 0; a < indexes.get(p).size(); a++) {
				if(supp.get(indexes.get(p).get(a))) {
					support = ctx.mkAnd(support, ctx.mkGt(strategies.get(p).get(a), zero));
				}
				else {
					support = ctx.mkAnd(support, ctx.mkEq(strategies.get(p).get(a), zero));
				}
			}
		}
		s.add(support);
	}
	
	public Pair<CSGResultStatus, ArrayList<Double>> computeEquilibria(BitSet supp, HashMap<Integer, int[]> map, int st) {
		s.push();
		computeConstraints(supp);
		computeSupport(supp, map);
		/*
		System.out.println("-- to prove");
		for(int i = 0; i < s.getNumAssertions(); i++) {
			System.out.println(s.getAssertions()[i]);
		}
		*/
		Model model;
		ArrayList<Double> eq = new ArrayList<>();
		Distribution d;
		double v;
		
		Status result = s.check();
		if(result == Status.SATISFIABLE) {
			model = s.getModel();
			strat.clear();
			for(int p = 0; p < numPlayers; p++) {
				d = new Distribution();
				eq.add(p, getDoubleValue(model, payoffs.get(p)));
				for (int a = 0; a < strategies.get(p).size(); a++) {
					v = getDoubleValue(model, strategies.get(p).get(a));
					if (v > 0)
						d.add(a, v);
				}
				strat.add(d);
			}
		}
		s.pop();
		return new Pair<CSGResultStatus, ArrayList<Double>>((result == Status.SATISFIABLE)? CSGResultStatus.SAT : 
															(result == Status.UNKNOWN)? CSGResultStatus.UNKNOWN :
																						CSGResultStatus.UNSAT, new ArrayList<Double>(eq));	
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
  	
	public void translateAssertions(HashMap<Integer, HashMap<Integer, ArrayList<Pair<BitSet, Double>>>> assertionsIdx, HashMap<Integer, int[]> map) {
		assertions.clear();
		ArithExpr sum, prod;
		int a, p;
		for(p = 0; p < assertionsIdx.keySet().size(); p++) {
			assertions.put(p, new HashMap<Integer, ArithExpr>());
			for(a = 0; a < assertionsIdx.get(p).keySet().size(); a++) {
				sum = zero;
				for(Pair<BitSet, Double> product : assertionsIdx.get(p).get(a)) {
					prod = getRealValue(product.second);
					for(int id = product.first.nextSetBit(0); id >= 0; id = product.first.nextSetBit(id + 1)) {
						prod = ctx.mkMul(prod, strategies.get(map.get(id)[0]).get(map.get(id)[1]));
					}
					sum = ctx.mkAdd(sum, prod);
				}
				assertions.get(p).put(a, sum);
			}
		}
	}
	
	public ArithExpr getRealValue(double v) {
		ArithExpr result;
		result = ctx.mkReal((long) (Math.pow(10, ndigits) * v));
		result = ctx.mkDiv(result, ctx.mkReal((long) (Math.pow(10, ndigits))));
		return result;
	}
	
	public void setNumPlayers(int n) {
		this.numPlayers = n;
	}
	
	public void setIndexes(ArrayList<ArrayList<Integer>> a) {
		this.indexes = a;
	}

	public ArrayList<Distribution> getStrat() {
		return new ArrayList<Distribution>(this.strat);
	}
	
	@Override
	public void setGradient(HashMap<Integer, HashMap<Integer, ArrayList<Pair<BitSet, Double>>>> gradient) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setAssertions(HashMap<Integer, HashMap<Integer, ArrayList<Pair<BitSet, Double>>>> assertions) {
		// TODO Auto-generated method stub
		
	}

}
