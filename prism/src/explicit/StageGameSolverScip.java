//==============================================================================
//	
//	Copyright (c) 2026-
//	Authors:
//	* Gabriel Santos <gabriel.santos@cs.ox.ac.uk> (University of Oxford)
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
import java.util.Collections;
import java.util.List;

import jscip.Constraint;
import jscip.Expression;
import jscip.SCIP_ParamSetting;
import jscip.SCIP_Status;
import jscip.SCIP_Vartype;
import jscip.Scip;
import jscip.Solution;
import jscip.Variable;
import prism.PrismException;

/**
 * {@link StageGameSolver} using SCIP (via JSCIPOpt), in floating point.
 * <ul>
 * <li>ZERO_SUM: two LPs (one per player).</li>
 * <li>CORRELATED: LP over the joint distribution; social fairness via two auxiliary variables
 *     t >= u_c >= s, minimising t - s (no binaries needed).</li>
 * <li>NASH: values u_c and regrets r_{c,a} = u_c - U_c(a, p_{-c}) >= 0 with complementarity
 *     p_{c,a} * r_{c,a} = 0 (SOS1 or big-M). Linear (MILP) for two coalitions, polynomial
 *     (nonconvex MINLP, solved globally by spatial branch-and-bound) for more.</li>
 * </ul>
 * Tie-breaking objectives are handled by sequential solves: each stage fixes the previous
 * objectives at their optimum (up to a relative tolerance) and optimises the next.
 * Free coalitions (see {@link StageGame#setFree}: in PRISM, those that have reached their target or failed)
 * keep their actions but have no incentive constraints and are not part of any objective: they cooperate
 * with the others. Their payoffs are still reported (computed from the solution).
 */
public class StageGameSolverScip implements StageGameSolver<Double>
{
	public enum Complementarity
	{
		SOS1, BIG_M
	}

	private static boolean libraryLoaded = false;

	/** Time limit per solve() call in seconds (<= 0: none) */
	private double timeLimit = 0.0;
	/** Relative tolerance used when fixing an optimal objective value in later stages */
	private double fixTolerance = 1e-9;
	/** Probabilities below this are set to zero in results */
	private double zeroTolerance = 1e-9;
	private Complementarity complementarity = Complementarity.SOS1;
	/** Social fairness tie-breaking: false = gap, then each coalition's payoff (as for correlated equilibria);
	 *  true = gap, then the sum of payoffs, then each coalition's payoff (as for Nash equilibria) */
	private boolean fairnessWelfareTieBreak = false;
	/** Whether SCIP's primal heuristics are switched off (default: true) */
	private boolean heuristicsOff = true;

	public StageGameSolverScip() throws PrismException
	{
		loadLibrary();
	}

	private static synchronized void loadLibrary() throws PrismException
	{
		if (libraryLoaded)
			return;
		try {
			System.loadLibrary("jscip");
			libraryLoaded = true;
		} catch (UnsatisfiedLinkError e) {
			throw new PrismException("Could not load SCIP (JSCIPOpt native library jscip): " + e.getMessage());
		}
	}

	public void setTimeLimit(double seconds)
	{
		timeLimit = seconds;
	}

	public void setFixTolerance(double tol)
	{
		fixTolerance = tol;
	}

	public void setZeroTolerance(double tol)
	{
		zeroTolerance = tol;
	}

	public void setComplementarity(Complementarity c)
	{
		complementarity = c;
	}

	public void setHeuristicsOff(boolean b)
	{
		heuristicsOff = b;
	}

	public void setFairnessWelfareTieBreak(boolean b)
	{
		fairnessWelfareTieBreak = b;
	}

	@Override
	public String getSolverName()
	{
		return "SCIP (JSCIPOpt)";
	}

	// ------------------------------------------------------------------------------------------

	@Override
	public StageGameResult<Double> solve(StageGame<Double> game, Concept concept, Criterion criterion, boolean minimise) throws PrismException
	{
		game.checkComplete();
		int n = game.getNumCoalitions();
		int nj = game.getNumJointActions();
		double sign = minimise ? -1.0 : 1.0;
		double[][] P = new double[n][nj];
		for (int c = 0; c < n; c++)
			for (int j = 0; j < nj; j++)
				P[c][j] = sign * game.getPayoff(c, j);

		StageGameResult<Double> res;
		if (concept == Concept.ZERO_SUM) {
			if (n != 2)
				throw new PrismException("Zero-sum stage games must have two coalitions");
			res = solveZeroSum(game, P);
		} else {
			res = solveStaged(game, P, concept, criterion);
		}
		if (minimise && res.hasSolution()) {
			List<Double> neg = new ArrayList<>();
			for (Double v : res.getPayoffs())
				neg.add(-v);
			res = new StageGameResult<>(res.getStatus(), neg, res.getStrategies(), res.getJointDistribution(), res.getGap(), res.getMessage());
		}
		return res;
	}

	// ---------------------------------------------------------------- staged (lexicographic) solving

	/** One objective (always maximised): sum_i coef_i * var_i over the model's variables. */
	private interface Objective
	{
		void apply(Model m, double[] coefsOut);
	}

	private static final int SUM = 0, GAP = 1, PLAYER = 2;

	private StageGameResult<Double> solveStaged(StageGame<Double> game, double[][] P, Concept concept, Criterion criterion) throws PrismException
	{
		// Free coalitions (no incentives: done or failed) are not part of any objective
		List<Integer> active = new ArrayList<>();
		for (int c = 0; c < game.getNumCoalitions(); c++)
			if (!game.isFree(c))
				active.add(c);
		int na = active.size();
		List<int[]> stages = new ArrayList<>(); // {kind, coalition}
		if (criterion == Criterion.SOCIAL_WELFARE) {
			stages.add(new int[] { SUM, -1 });
			for (int i = 0; i < na - 1; i++) // with the sum fixed, the last payoff is determined
				stages.add(new int[] { PLAYER, active.get(i) });
		} else if (fairnessWelfareTieBreak) {
			stages.add(new int[] { GAP, -1 });
			stages.add(new int[] { SUM, -1 });
			for (int i = 0; i < na - 1; i++)
				stages.add(new int[] { PLAYER, active.get(i) });
		} else {
			stages.add(new int[] { GAP, -1 });
			for (int i = 0; i < na; i++)
				stages.add(new int[] { PLAYER, active.get(i) });
		}
		boolean needGap = criterion == Criterion.SOCIAL_FAIRNESS;
		if (FREE_TRANSFORM != null)
			return solveStagedReuse(game, P, concept, stages, needGap);
		long start = System.nanoTime();
		double[] optima = new double[stages.size()];
		StageGameResult<Double> last = null;
		for (int k = 0; k < stages.size(); k++) {
			Model m = concept == Concept.NASH ? buildNash(game, P, needGap) : buildCorrelated(game, P, needGap);
			try {
				for (int i = 0; i < k; i++) {
					double tol = fixTolerance * Math.max(1.0, Math.abs(optima[i]));
					m.addObjectiveBound(stages.get(i), optima[i] - tol);
				}
				m.setObjective(stages.get(k));
				StageOutcome o = solveStage(m, k, start);
				if (o.done != null)
					return o.done.apply(last);
				optima[k] = o.optimum;
				last = o.result;
			} finally {
				m.free();
			}
		}
		return last;
	}

	/**
	 * Same stages as solveStaged, on one SCIP model: after each stage the transformed problem is freed
	 * (SCIPfreeTransform, which keeps the best solutions as starting solutions), the objective just optimised is
	 * bounded by its optimum (within the tolerance) and the next objective is set.
	 */
	private StageGameResult<Double> solveStagedReuse(StageGame<Double> game, double[][] P, Concept concept, List<int[]> stages, boolean needGap)
			throws PrismException
	{
		long start = System.nanoTime();
		StageGameResult<Double> last = null;
		double lastOptimum = 0.0;
		Model m = concept == Concept.NASH ? buildNash(game, P, needGap) : buildCorrelated(game, P, needGap);
		try {
			for (int k = 0; k < stages.size(); k++) {
				if (k > 0) {
					freeTransform(m.scip);
					int[] prev = stages.get(k - 1);
					m.clearObjective(prev);
					double tol = fixTolerance * Math.max(1.0, Math.abs(lastOptimum));
					m.addObjectiveBound(prev, lastOptimum - tol);
				}
				m.setObjective(stages.get(k));
				StageOutcome o = solveStage(m, k, start);
				if (o.done != null)
					return o.done.apply(last);
				lastOptimum = o.optimum;
				last = o.result;
			}
		} finally {
			m.free();
		}
		return last;
	}

	/** Result of solving one stage: either an optimum (and its solution), or a final result given the previous stage's */
	private static final class StageOutcome
	{
		double optimum;
		StageGameResult<Double> result;
		java.util.function.Function<StageGameResult<Double>, StageGameResult<Double>> done;
	}

	private StageOutcome solveStage(Model m, int k, long start)
	{
		StageOutcome o = new StageOutcome();
		if (timeLimit > 0) {
			double left = timeLimit - (System.nanoTime() - start) / 1e9;
			if (left <= 0) {
				o.done = last -> timeLimited(last, "time limit reached before stage " + k);
				return o;
			}
			m.scip.setRealParam("limits/time", left);
		}
		m.scip.solve();
		SCIP_Status st = m.scip.getStatus();
		boolean hasSol = m.scip.getNSols() > 0;
		if (st == SCIP_Status.SCIP_STATUS_OPTIMAL && hasSol) {
			Solution sol = m.scip.getBestSol();
			o.optimum = m.scip.getSolOrigObj(sol);
			o.result = m.extract(sol, StageGameResult.Status.OPTIMAL, 0.0, null);
		} else if (hasSol) {
			StageGameResult<Double> r = m.extract(m.scip.getBestSol(), StageGameResult.Status.TIME_LIMIT, m.scip.getGap(), "stopped at stage " + k + " (" + st + ")");
			o.done = last -> r;
		} else if (st == SCIP_Status.SCIP_STATUS_INFEASIBLE) {
			if (k == 0)
				o.done = last -> StageGameResult.failed(StageGameResult.Status.INFEASIBLE, "no equilibrium found (infeasible)");
			else
				// fixing previous optima made the model infeasible: numerical issue; keep the previous stage's solution
				o.done = last -> new StageGameResult<>(StageGameResult.Status.OPTIMAL, last.getPayoffs(), last.getStrategies(), last.getJointDistribution(), 0.0,
						"tie-breaking stopped at stage " + k + " (infeasible after fixing previous optima)");
		} else {
			o.done = last -> k == 0 ? StageGameResult.failed(StageGameResult.Status.FAILED, "SCIP status " + st)
					: timeLimited(last, "stopped at stage " + k + " (" + st + ")");
		}
		return o;
	}

	/** Scip.freeTransform() if the JSCIPOpt build provides it (added to JSCIPOpt for PRISM), else null */
	private static final java.lang.invoke.MethodHandle FREE_TRANSFORM = findFreeTransform();

	private static java.lang.invoke.MethodHandle findFreeTransform()
	{
		try {
			return java.lang.invoke.MethodHandles.publicLookup().findVirtual(Scip.class, "freeTransform", java.lang.invoke.MethodType.methodType(void.class));
		} catch (ReflectiveOperationException | RuntimeException e) {
			return null;
		}
	}

	private static void freeTransform(Scip scip) throws PrismException
	{
		try {
			FREE_TRANSFORM.invokeExact(scip);
		} catch (Throwable e) {
			throw new PrismException("SCIP: could not free the transformed problem: " + e.getMessage());
		}
	}

	private static StageGameResult<Double> timeLimited(StageGameResult<Double> last, String msg)
	{
		if (last == null)
			return StageGameResult.failed(StageGameResult.Status.TIME_LIMIT, msg);
		// criterion optimal, tie-breaking incomplete
		return new StageGameResult<>(StageGameResult.Status.TIME_LIMIT, last.getPayoffs(), last.getStrategies(), last.getJointDistribution(), 0.0, msg);
	}

	// ---------------------------------------------------------------- model building

	/** A SCIP model plus the variables needed to set objectives and read solutions. */
	private final class Model
	{
		final Scip scip = new Scip();
		final StageGame<Double> game;
		final double[][] P;
		final List<Variable> vars = new ArrayList<>();
		Variable[] u;          // coalition values
		Variable t, s;         // highest / lowest value (social fairness)
		Variable[][] p;        // Nash: mixed strategies
		Variable[] q;          // correlated: joint distribution

		Model(String name, StageGame<Double> game, double[][] P)
		{
			scip.create(name);
			scip.hideOutput(true);
			scip.setMaximize();
			if (heuristicsOff)
				scip.setHeuristics(SCIP_ParamSetting.SCIP_PARAMSETTING_OFF, true);
			this.game = game;
			this.P = P;
		}

		Variable var(String name, double lb, double ub, SCIP_Vartype type)
		{
			Variable v = scip.createVar(name, lb, ub, 0.0, type);
			vars.add(v);
			return v;
		}

		Variable cont(String name, double lb, double ub)
		{
			return var(name, lb, ub, SCIP_Vartype.SCIP_VARTYPE_CONTINUOUS);
		}

		void linear(String name, List<Variable> vs, List<Double> cs, double lhs, double rhs)
		{
			Constraint cons = scip.createConsLinear(name, vs.toArray(new Variable[0]), toArray(cs), lhs, rhs);
			scip.addCons(cons);
			scip.releaseCons(cons);
		}

		void addFairnessVars()
		{
			double inf = scip.infinity();
			t = cont("t", -inf, inf);
			s = cont("s", -inf, inf);
			for (int c = 0; c < u.length; c++) {
				if (game.isFree(c))
					continue;
				linear("t" + c, List.of(t, u[c]), List.of(1.0, -1.0), 0.0, inf);   // t >= u_c
				linear("s" + c, List.of(s, u[c]), List.of(1.0, -1.0), -inf, 0.0);  // s <= u_c
			}
		}

		/** Objective of a stage as (variables, coefficients), to be maximised. */
		void objectiveTerms(int[] stage, List<Variable> vs, List<Double> cs)
		{
			switch (stage[0]) {
			case SUM:
				for (int c = 0; c < u.length; c++)
					if (!game.isFree(c)) { vs.add(u[c]); cs.add(1.0); }
				break;
			case GAP:
				vs.add(t); cs.add(-1.0);
				vs.add(s); cs.add(1.0);
				break;
			default:
				vs.add(u[stage[1]]); cs.add(1.0);
			}
		}

		void setObjective(int[] stage)
		{
			List<Variable> vs = new ArrayList<>();
			List<Double> cs = new ArrayList<>();
			objectiveTerms(stage, vs, cs);
			for (int i = 0; i < vs.size(); i++)
				scip.changeVarObj(vs.get(i), cs.get(i));
		}

		void clearObjective(int[] stage)
		{
			List<Variable> vs = new ArrayList<>();
			List<Double> cs = new ArrayList<>();
			objectiveTerms(stage, vs, cs);
			for (Variable v : vs)
				scip.changeVarObj(v, 0.0);
		}

		void addObjectiveBound(int[] stage, double lb)
		{
			List<Variable> vs = new ArrayList<>();
			List<Double> cs = new ArrayList<>();
			objectiveTerms(stage, vs, cs);
			linear("fix" + stage[0] + "_" + stage[1], vs, cs, lb, scip.infinity());
		}

		StageGameResult<Double> extract(Solution sol, StageGameResult.Status status, double gap, String msg)
		{
			int n = game.getNumCoalitions();
			List<Double> payoffs = new ArrayList<>(Collections.nCopies(n, 0.0));
			if (p != null) {
				List<List<Double>> strategies = new ArrayList<>();
				double[][] x = new double[n][];
				for (int c = 0; c < n; c++) {
					x[c] = new double[p[c].length];
					for (int a = 0; a < p[c].length; a++)
						x[c][a] = scip.getSolVal(sol, p[c][a]);
					x[c] = clean(x[c]);
					List<Double> l = new ArrayList<>();
					for (double v : x[c])
						l.add(v);
					strategies.add(l);
				}
				for (int j = 0; j < game.getNumJointActions(); j++) {
					int[] a = game.jointActions(j);
					double pr = 1.0;
					for (int c = 0; c < n && pr != 0.0; c++)
						pr *= x[c][a[c]];
					for (int c = 0; c < n; c++)
						payoffs.set(c, payoffs.get(c) + pr * P[c][j]);
				}
				return new StageGameResult<>(status, payoffs, strategies, null, gap, msg);
			} else {
				double[] y = new double[q.length];
				for (int j = 0; j < q.length; j++)
					y[j] = scip.getSolVal(sol, q[j]);
				y = clean(y);
				List<Double> joint = new ArrayList<>();
				for (int j = 0; j < y.length; j++) {
					joint.add(y[j]);
					for (int c = 0; c < n; c++)
						payoffs.set(c, payoffs.get(c) + y[j] * P[c][j]);
				}
				return new StageGameResult<>(status, payoffs, null, joint, gap, msg);
			}
		}

		void free()
		{
			for (Variable v : vars)
				scip.releaseVar(v);
			scip.free();
		}
	}

	/** Clamp to [0,1], zero tiny entries and renormalise. */
	private double[] clean(double[] x)
	{
		double sum = 0.0;
		for (int i = 0; i < x.length; i++) {
			x[i] = Math.min(1.0, Math.max(0.0, x[i]));
			if (x[i] < zeroTolerance)
				x[i] = 0.0;
			sum += x[i];
		}
		if (sum > 0)
			for (int i = 0; i < x.length; i++)
				x[i] /= sum;
		return x;
	}

	private Model buildCorrelated(StageGame<Double> game, double[][] P, boolean fairness)
	{
		int n = game.getNumCoalitions(), nj = game.getNumJointActions();
		Model m = new Model("ce", game, P);
		double inf = m.scip.infinity();
		m.q = new Variable[nj];
		List<Variable> all = new ArrayList<>();
		List<Double> ones = new ArrayList<>();
		for (int j = 0; j < nj; j++) {
			m.q[j] = m.cont("q" + j, 0.0, 1.0);
			all.add(m.q[j]);
			ones.add(1.0);
		}
		m.linear("dist", all, ones, 1.0, 1.0);
		m.u = new Variable[n];
		for (int c = 0; c < n; c++) {
			m.u[c] = m.cont("u" + c, -inf, inf);
			List<Variable> vs = new ArrayList<>();
			List<Double> cs = new ArrayList<>();
			for (int j = 0; j < nj; j++)
				if (P[c][j] != 0.0) { vs.add(m.q[j]); cs.add(P[c][j]); }
			vs.add(m.u[c]); cs.add(-1.0);
			m.linear("u" + c, vs, cs, 0.0, 0.0);  // u_c = sum_j q_j P_c(j)
		}
		// incentive constraints: following recommendation a is at least as good as deviating to b
		for (int c = 0; c < n; c++) {
			if (game.isFree(c))
				continue;
			for (int a = 0; a < game.getNumActions(c); a++) {
				for (int b = 0; b < game.getNumActions(c); b++) {
					if (a == b)
						continue;
					List<Variable> vs = new ArrayList<>();
					List<Double> cs = new ArrayList<>();
					for (int j = 0; j < nj; j++) {
						if (game.jointActions(j)[c] != a)
							continue;
						double d = P[c][j] - P[c][game.deviate(j, c, b)];
						if (d != 0.0) { vs.add(m.q[j]); cs.add(d); }
					}
					if (!vs.isEmpty())
						m.linear("ic" + c + "_" + a + "_" + b, vs, cs, 0.0, inf);
				}
			}
		}
		if (fairness)
			m.addFairnessVars();
		return m;
	}

	private Model buildNash(StageGame<Double> game, double[][] P, boolean fairness)
	{
		int n = game.getNumCoalitions(), nj = game.getNumJointActions();
		Model m = new Model("ne", game, P);
		Scip scip = m.scip;
		double inf = scip.infinity();
		m.p = new Variable[n][];
		for (int c = 0; c < n; c++) {
			int k = game.getNumActions(c);
			m.p[c] = new Variable[k];
			List<Variable> vs = new ArrayList<>();
			List<Double> cs = new ArrayList<>();
			for (int a = 0; a < k; a++) {
				m.p[c][a] = m.cont("p" + c + "_" + a, 0.0, 1.0);
				vs.add(m.p[c][a]);
				cs.add(1.0);
			}
			m.linear("simplex" + c, vs, cs, 1.0, 1.0);
		}
		m.u = new Variable[n];
		for (int c = 0; c < n; c++)
			m.u[c] = game.isFree(c) ? m.cont("u" + c, 0.0, 0.0) : m.cont("u" + c, -inf, inf); // free: unused

		// one expression per strategy variable, shared by all nonlinear terms
		Expression[][] px = null;
		List<Expression> toRelease = new ArrayList<>();
		if (n > 2) {
			px = new Expression[n][];
			for (int c = 0; c < n; c++) {
				px[c] = new Expression[m.p[c].length];
				for (int a = 0; a < m.p[c].length; a++) {
					px[c][a] = scip.createExprVar(m.p[c][a]);
					toRelease.add(px[c][a]);
				}
			}
		}

		for (int c = 0; c < n; c++) {
			if (game.isFree(c))
				continue; // no incentives, and not part of the objectives
			double bigM = max(P[c]) - min(P[c]);
			for (int a = 0; a < game.getNumActions(c); a++) {
				Variable r = m.cont("r" + c + "_" + a, 0.0, inf);
				// u_c - r_{c,a} - U_c(a, p_{-c}) = 0
				if (n <= 2) {
					List<Variable> vs = new ArrayList<>();
					List<Double> cs = new ArrayList<>();
					vs.add(m.u[c]); cs.add(1.0);
					vs.add(r); cs.add(-1.0);
					if (n == 2) {
						int o = 1 - c;
						for (int b = 0; b < game.getNumActions(o); b++) {
							int[] act = new int[2];
							act[c] = a;
							act[o] = b;
							double v = P[c][game.jointIndex(act)];
							if (v != 0.0) { vs.add(m.p[o][b]); cs.add(-v); }
						}
						m.linear("reg" + c + "_" + a, vs, cs, 0.0, 0.0);
					} else {
						m.linear("reg" + c + "_" + a, vs, cs, P[c][a], P[c][a]);
					}
				} else {
					List<Expression> ch = new ArrayList<>();
					List<Double> cf = new ArrayList<>();
					for (int j = 0; j < nj; j++) {
						int[] act = game.jointActions(j);
						if (act[c] != a || P[c][j] == 0.0)
							continue;
						Expression[] fs = new Expression[n - 1];
						int f = 0;
						for (int d = 0; d < n; d++)
							if (d != c)
								fs[f++] = px[d][act[d]];
						Expression prod = scip.createExprProduct(fs, P[c][j]);
						toRelease.add(prod);
						ch.add(prod);
						cf.add(-1.0);
					}
					Expression eu = scip.createExprVar(m.u[c]);
					Expression er = scip.createExprVar(r);
					toRelease.add(eu);
					toRelease.add(er);
					ch.add(eu); cf.add(1.0);
					ch.add(er); cf.add(-1.0);
					Expression sum = scip.createExprSum(ch.toArray(new Expression[0]), toArray(cf), 0.0);
					toRelease.add(sum);
					Constraint cons = scip.createConsNonlinear("reg" + c + "_" + a, sum, 0.0, 0.0);
					scip.addCons(cons);
					scip.releaseCons(cons);
				}
				// complementarity p_{c,a} * r_{c,a} = 0
				if (complementarity == Complementarity.SOS1) {
					Constraint cons = scip.createConsSOS1("cmp" + c + "_" + a, new Variable[] { m.p[c][a], r }, new double[] { 1.0, 2.0 });
					scip.addCons(cons);
					scip.releaseCons(cons);
				} else {
					Variable b = m.var("b" + c + "_" + a, 0.0, 1.0, SCIP_Vartype.SCIP_VARTYPE_BINARY);
					m.linear("cmpa" + c + "_" + a, List.of(m.p[c][a], b), List.of(1.0, -1.0), -inf, 0.0);   // p <= b
					m.linear("cmpb" + c + "_" + a, List.of(r, b), List.of(1.0, bigM), -inf, bigM);           // r <= M (1 - b)
				}
			}
		}
		for (Expression e : toRelease)
			scip.releaseExpr(e);
		if (fairness)
			m.addFairnessVars();
		return m;
	}



	// ---------------------------------------------------------------- zero-sum

	private StageGameResult<Double> solveZeroSum(StageGame<Double> game, double[][] P) throws PrismException
	{
		int m0 = game.getNumActions(0), m1 = game.getNumActions(1);
		double[][] A = new double[m0][m1];
		for (int a = 0; a < m0; a++)
			for (int b = 0; b < m1; b++)
				A[a][b] = P[0][game.jointIndex(new int[] { a, b })];
		double[] x = new double[m0], y = new double[m1];
		double v = matrixGame(A, x, false);
		matrixGame(A, y, true);
		List<List<Double>> strategies = new ArrayList<>();
		strategies.add(toList(clean(x)));
		strategies.add(toList(clean(y)));
		List<Double> payoffs = new ArrayList<>();
		payoffs.add(v);
		payoffs.add(-v);
		return new StageGameResult<>(StageGameResult.Status.OPTIMAL, payoffs, strategies, null, 0.0, null);
	}

	/**
	 * Row player (column=false): max v s.t. sum_i A[i][j] x_i >= v for all j.
	 * Column player (column=true): min w s.t. sum_j A[i][j] y_j <= w for all i.
	 */
	private double matrixGame(double[][] A, double[] strat, boolean column) throws PrismException
	{
		int rows = A.length, cols = A[0].length;
		int k = column ? cols : rows, h = column ? rows : cols;
		Model m = new Model("matrix", null, null);
		try {
			Scip scip = m.scip;
			double inf = scip.infinity();
			Variable[] z = new Variable[k];
			List<Variable> all = new ArrayList<>();
			List<Double> ones = new ArrayList<>();
			for (int i = 0; i < k; i++) {
				z[i] = m.cont("z" + i, 0.0, 1.0);
				all.add(z[i]);
				ones.add(1.0);
			}
			m.linear("simplex", all, ones, 1.0, 1.0);
			Variable v = m.cont("v", -inf, inf);
			if (column)
				scip.setMinimize();
			scip.changeVarObj(v, 1.0);
			for (int o = 0; o < h; o++) {
				List<Variable> vs = new ArrayList<>();
				List<Double> cs = new ArrayList<>();
				for (int i = 0; i < k; i++) {
					double a = column ? A[o][i] : A[i][o];
					if (a != 0.0) { vs.add(z[i]); cs.add(a); }
				}
				vs.add(v); cs.add(-1.0);
				if (column)
					m.linear("c" + o, vs, cs, -inf, 0.0);
				else
					m.linear("c" + o, vs, cs, 0.0, inf);
			}
			if (timeLimit > 0)
				scip.setRealParam("limits/time", timeLimit);
			scip.solve();
			if (scip.getStatus() != SCIP_Status.SCIP_STATUS_OPTIMAL || scip.getNSols() == 0)
				throw new PrismException("SCIP could not solve matrix game (status " + scip.getStatus() + ")");
			Solution sol = scip.getBestSol();
			for (int i = 0; i < k; i++)
				strat[i] = scip.getSolVal(sol, z[i]);
			return scip.getSolVal(sol, v);
		} finally {
			m.free();
		}
	}

	// ---------------------------------------------------------------- helpers

	private static double[] toArray(List<Double> l)
	{
		double[] a = new double[l.size()];
		for (int i = 0; i < a.length; i++)
			a[i] = l.get(i);
		return a;
	}

	private static List<Double> toList(double[] a)
	{
		List<Double> l = new ArrayList<>(a.length);
		for (double v : a)
			l.add(v);
		return l;
	}

	private static double max(double[] a)
	{
		double m = Double.NEGATIVE_INFINITY;
		for (double v : a)
			m = Math.max(m, v);
		return m;
	}

	private static double min(double[] a)
	{
		double m = Double.POSITIVE_INFINITY;
		for (double v : a)
			m = Math.min(m, v);
		return m;
	}
}
