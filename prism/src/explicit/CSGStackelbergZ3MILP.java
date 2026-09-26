package explicit;

import java.util.ArrayList;
import java.util.HashMap;

import com.microsoft.z3.AlgebraicNum;
import com.microsoft.z3.ArithExpr;
import com.microsoft.z3.Context;
import com.microsoft.z3.Expr;
import com.microsoft.z3.Model;
import com.microsoft.z3.Optimize;
import com.microsoft.z3.RatNum;
import com.microsoft.z3.RealExpr;
import com.microsoft.z3.Status;
import com.microsoft.z3.Version;

public class CSGStackelbergZ3MILP {
	
	private HashMap<String, String> cfg;
	private Context ctx;
	private Optimize s;
	private String name;
	
	public CSGStackelbergZ3MILP() {
		cfg = new HashMap<String, String>();
		cfg.put("model", "true");
		ctx = new Context(cfg);
		s = ctx.mkOptimize();
		name = Version.getFullVersion();
	}
	
	public EquilibriumResult computeEquilibrium(ArrayList<ArrayList<ArrayList<Double>>> bmgame) {
		EquilibriumResult result = new EquilibriumResult();
		int nrows = bmgame.get(0).size();
		int ncols = bmgame.get(0).get(0).size();
		RealExpr[] payoff_vars = new RealExpr[2];
		RealExpr[] p1_strategy_vars = new RealExpr[nrows];
		RealExpr[] p2_strategy_vars = new RealExpr[ncols];
		ArithExpr expr1, expr2;
		ArrayList<Double> payoffs_result = new ArrayList<Double>();
		int c, p, r;
		
		for (p = 0; p < 2; p++) {
			payoff_vars[p] = ctx.mkRealConst("p" + p);
		}
		
		expr1 = ctx.mkInt(0);
		for (r = 0; r < nrows; r++) {
			p1_strategy_vars[r] = ctx.mkRealConst("x" + r);
			s.Add(ctx.mkGe(p1_strategy_vars[r], ctx.mkInt(0)));
			expr1 = ctx.mkAdd(expr1, p1_strategy_vars[r]);
		}
		s.Add(ctx.mkEq(expr1, ctx.mkInt(1)));
		
		expr2 = ctx.mkInt(0);
		for (c = 0; c < ncols; c++) {
			p2_strategy_vars[c] = ctx.mkRealConst("y" + c);
			s.Add(ctx.mkOr(ctx.mkEq(p2_strategy_vars[c], ctx.mkInt(0)), ctx.mkEq(p2_strategy_vars[c], ctx.mkInt(1))));
			expr2 = ctx.mkAdd(expr2, p2_strategy_vars[c]);
		}
		s.Add(ctx.mkEq(expr2, ctx.mkInt(1)));
	
		for (c = 0; c < ncols; c++) {
			expr1 = ctx.mkInt(0);
			expr2 = ctx.mkInt(0);
			for (r = 0; r < nrows; r++) {
				expr1 = ctx.mkAdd(expr1, ctx.mkMul(p1_strategy_vars[r], ctx.mkReal(String.valueOf(bmgame.get(0).get(r).get(c)))));
				expr2 = ctx.mkAdd(expr2, ctx.mkMul(p1_strategy_vars[r], ctx.mkReal(String.valueOf(bmgame.get(1).get(r).get(c)))));
			}
			s.Add(ctx.mkLe(expr2, payoff_vars[1]));
			s.Add(ctx.mkImplies(ctx.mkGt(p2_strategy_vars[c], ctx.mkInt(0)), ctx.mkEq(payoff_vars[0], expr1)));
			s.Add(ctx.mkImplies(ctx.mkGt(p2_strategy_vars[c], ctx.mkInt(0)), ctx.mkEq(payoff_vars[1], expr2)));
		}
		
		s.MkMaximize(payoff_vars[0]);
		s.MkMaximize(payoff_vars[1]);
		
		if (s.Check() == Status.SATISFIABLE) {
			System.out.println(s.getModel());
			for (p = 0; p < 2; p++) {
				payoffs_result.add(getDoubleValue(s.getModel(), payoff_vars[p]));
			}
			result.setStatus(CSGModelCheckerEquilibria.CSGResultStatus.SAT);
			result.setPayoffVector(payoffs_result);
		}
		else {
			result.setStatus(CSGModelCheckerEquilibria.CSGResultStatus.UNSAT);
		}
		return result;
	}

	/**
	 * Return a double value for a given expression (usually variable), converting from BigInt fractions
	 * @param model The SMT model
	 * @param expr The SMT expression
	 * @return
	 */
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
	
	public String getSolverName() {
		return name;
	}
	
}

/*
 * Stage-level use of this solver, formerly CSGModelCheckerEquilibria.stepStackelbergEquilibria (never called;
 * removed there because the semantics for the whole game/logic were never worked out). Kept for reference:

	public double[] stepStackelbergEquilibria(CSG<Double> csg, List<CSGRewards<Double>> rewards, List<Map<Integer, BitSet>> mmap,
											  List<List<Map<BitSet, Double>>> strats, double[][] val, int s, boolean min) throws PrismException {
		EquilibriumResult result;
		ArrayList<ArrayList<Integer>> nmap;
		ArrayList<ArrayList<ArrayList<Double>>> bmgame;
		double[] eqs = new double[3];
		int i, p;
		
		mmap = new ArrayList<Map<Integer, BitSet>>();
		nmap = new ArrayList<ArrayList<Integer>>();
		for (p = 0; p < 2; p++) {
			mmap.add(p, new HashMap<Integer, BitSet>());
			nmap.add(p, new ArrayList<Integer>());
		}
		bmgame = buildBimatrixGame(csg, rewards.get(0), rewards.get(1), mmap, nmap, val, s, min);	
		
		CSGStackelbergZ3MILP stackelbergSolver = new CSGStackelbergZ3MILP();
		result = stackelbergSolver.computeEquilibrium(bmgame);
		
		if (result.getStatus() == CSGResultStatus.SAT) {
			eqs[0] = 0.0;
			for (Double d : result.getPayoffVector()) {
				eqs[0] += d;
			}
			for (p = 0; p < numCoalitions; p++) {
				eqs[p+1] = result.getPayoffVector().get(p);
			}
		}
		else {
			throw new PrismException(stackelbergSolver.getSolverName() + " could not find an optimal solution for state " + s);
		}
		
		if (rewards != null) {
			addStateRewards(eqs, rewards, s, min);
			// keep eqs[0] = sum of the payoffs
			eqs[0] = 0.0;
			for (i = 1; i < eqs.length; i++)
				eqs[0] += eqs[i];
		}
		if (min) {
			for (i = 0; i < eqs.length; i++)
				eqs[i] = -1.0 * eqs[i];
		}
		return eqs;
	}

*/
