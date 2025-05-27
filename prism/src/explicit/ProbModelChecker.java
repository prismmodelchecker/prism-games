//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford)
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

import static prism.PrismSettings.DEFAULT_EXPORT_MODEL_PRECISION;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.List;

import explicit.rewards.CSGRewards;
import explicit.rewards.ConstructRewards;
import explicit.rewards.MCRewards;
import explicit.rewards.MDPRewards;
import explicit.rewards.Rewards;
import explicit.rewards.RewardsExplicit;
import explicit.rewards.STPGRewards;
import parser.BooleanUtils;
import parser.ast.Coalition;
import parser.ast.Expression;
import parser.ast.ExpressionFunc;
import parser.ast.ExpressionLiteral;
import parser.ast.ExpressionMultiNash;
import parser.ast.ExpressionMultiNashProb;
import parser.ast.ExpressionMultiNashReward;
import parser.ast.ExpressionProb;
import parser.ast.ExpressionQuant;
import parser.ast.ExpressionReward;
import parser.ast.ExpressionSS;
import parser.ast.ExpressionStrategy;
import parser.ast.ExpressionTemporal;
import parser.ast.ExpressionUnaryOp;
import parser.type.TypeBool;
import parser.type.TypeDouble;
import parser.type.TypePathBool;
import parser.type.TypePathDouble;
import prism.AccuracyFactory;
import prism.Evaluator;
import prism.IntegerBound;
import prism.ModelType;
import prism.OpRelOpBound;
import prism.Prism;
import prism.PrismComponent;
import prism.PrismException;
import prism.PrismLog;
import prism.PrismNotSupportedException;
import prism.PrismSettings;
import prism.RewardGenerator;

/**
 * Super class for explicit-state probabilistic model checkers.
 */
public class ProbModelChecker extends NonProbModelChecker
{
	// Flags/settings
	// (NB: defaults do not necessarily coincide with PRISM)

	// Method used to solve linear equation systems
	protected LinEqMethod linEqMethod = LinEqMethod.GAUSS_SEIDEL;
	// Method used to solve MDPs
	protected MDPSolnMethod mdpSolnMethod = MDPSolnMethod.GAUSS_SEIDEL;
	// Method used to solve STPGs
	protected STPGSolnMethod stpgSolnMethod = STPGSolnMethod.GAUSS_SEIDEL;
	// Method used to solve IMDPs (and IDTMCs)
	protected IMDPSolnMethod imdpSolnMethod = IMDPSolnMethod.GAUSS_SEIDEL;
	// Iterative numerical method termination criteria
	protected TermCrit termCrit = TermCrit.RELATIVE;
	// Parameter for iterative numerical method termination criteria
	protected double termCritParam = 1e-8;
	// Max iterations for numerical solution
	protected int maxIters = 100000;
	// Resolution for POMDP fixed grid approximation algorithm
	protected int gridResolution = 10;
	// Use precomputation algorithms in model checking?
	protected boolean precomp = true;
	protected boolean prob0 = true;
	protected boolean prob1 = true;
	// should we suppress log output during precomputations?
	protected boolean silentPrecomputations = false;
	// Use predecessor relation? (e.g. for precomputation)
	protected boolean preRel = true;
	// Direction of convergence for value iteration (lfp/gfp)
	protected ValIterDir valIterDir = ValIterDir.BELOW;
	// Method used for numerical solution
	protected SolnMethod solnMethod = SolnMethod.VALUE_ITERATION;
	// Is non-convergence of an iterative method an error?
	protected boolean errorOnNonConverge = true;

	protected boolean useDiscounting = false;
	protected double discountFactor = 1.0;

	// Delay between occasional updates for slow processes, e.g. numerical solution (milliseconds)
	public static final int UPDATE_DELAY = 5000;

	// Enums for flags/settings

	// Method used for numerical solution
	public enum LinEqMethod {
		POWER, JACOBI, GAUSS_SEIDEL, BACKWARDS_GAUSS_SEIDEL, JOR, SOR, BACKWARDS_SOR;
		public String fullName()
		{
			switch (this) {
			case POWER:
				return "Power method";
			case JACOBI:
				return "Jacobi";
			case GAUSS_SEIDEL:
				return "Gauss-Seidel";
			case BACKWARDS_GAUSS_SEIDEL:
				return "Backwards Gauss-Seidel";
			case JOR:
				return "JOR";
			case SOR:
				return "SOR";
			case BACKWARDS_SOR:
				return "Backwards SOR";
			default:
				return this.toString();
			}
		}
	};

	// Method used for solving MDPs
	public enum MDPSolnMethod {
		VALUE_ITERATION, GAUSS_SEIDEL, POLICY_ITERATION, MODIFIED_POLICY_ITERATION, LINEAR_PROGRAMMING;
		public String fullName()
		{
			switch (this) {
			case VALUE_ITERATION:
				return "Value iteration";
			case GAUSS_SEIDEL:
				return "Gauss-Seidel";
			case POLICY_ITERATION:
				return "Policy iteration";
			case MODIFIED_POLICY_ITERATION:
				return "Modified policy iteration";
			case LINEAR_PROGRAMMING:
				return "Linear programming";
			default:
				return this.toString();
			}
		}
	};

	// Method used for solving STPGs
	public enum STPGSolnMethod {
		VALUE_ITERATION, GAUSS_SEIDEL;
		public String fullName()
		{
			switch (this) {
				case VALUE_ITERATION:
					return "Value iteration";
				case GAUSS_SEIDEL:
					return "Gauss-Seidel";
				default:
					return this.toString();
			}
		}
	};

	// Method used for solving IMDPs (and IDTMCs)
	public enum IMDPSolnMethod {
		VALUE_ITERATION, GAUSS_SEIDEL;
		public String fullName()
		{
			switch (this) {
			case VALUE_ITERATION:
				return "Value iteration";
			case GAUSS_SEIDEL:
				return "Gauss-Seidel";
			default:
				return this.toString();
			}
		}
	};

	// Iterative numerical method termination criteria
	public enum TermCrit {
		ABSOLUTE, RELATIVE
	};

	// Direction of convergence for value iteration (lfp/gfp)
	public enum ValIterDir {
		BELOW, ABOVE
	};

	// Method used for numerical solution
	public enum SolnMethod {
		VALUE_ITERATION, GAUSS_SEIDEL, POLICY_ITERATION, MODIFIED_POLICY_ITERATION, LINEAR_PROGRAMMING
	};

	/**
	 * Create a new ProbModelChecker, inherit basic state from parent (unless null).
	 */
	public ProbModelChecker(PrismComponent parent) throws PrismException
	{
		super(parent);

		// If present, initialise settings from PrismSettings
		if (settings != null) {
			String s;
			// PRISM_LIN_EQ_METHOD
			s = settings.getString(PrismSettings.PRISM_LIN_EQ_METHOD);
			if (s.equals("Power")) {
				setLinEqMethod(LinEqMethod.POWER);
			} else if (s.equals("Jacobi")) {
				setLinEqMethod(LinEqMethod.JACOBI);
			} else if (s.equals("Gauss-Seidel")) {
				setLinEqMethod(LinEqMethod.GAUSS_SEIDEL);
			} else if (s.equals("Backwards Gauss-Seidel")) {
				setLinEqMethod(LinEqMethod.BACKWARDS_GAUSS_SEIDEL);
			} else if (s.equals("JOR")) {
				setLinEqMethod(LinEqMethod.JOR);
			} else if (s.equals("SOR")) {
				setLinEqMethod(LinEqMethod.SOR);
			} else if (s.equals("Backwards SOR")) {
				setLinEqMethod(LinEqMethod.BACKWARDS_SOR);
			} else {
				throw new PrismNotSupportedException("Explicit engine does not support linear equation solution method \"" + s + "\"");
			}
			// PRISM_MDP_SOLN_METHOD
			s = settings.getString(PrismSettings.PRISM_MDP_SOLN_METHOD);
			if (s.equals("Value iteration")) {
				setMDPSolnMethod(MDPSolnMethod.VALUE_ITERATION);
			} else if (s.equals("Gauss-Seidel")) {
				setMDPSolnMethod(MDPSolnMethod.GAUSS_SEIDEL);
			} else if (s.equals("Policy iteration")) {
				setMDPSolnMethod(MDPSolnMethod.POLICY_ITERATION);
			} else if (s.equals("Modified policy iteration")) {
				setMDPSolnMethod(MDPSolnMethod.MODIFIED_POLICY_ITERATION);
			} else if (s.equals("Linear programming")) {
				setMDPSolnMethod(MDPSolnMethod.LINEAR_PROGRAMMING);
			} else {
				throw new PrismNotSupportedException("Explicit engine does not support MDP solution method \"" + s + "\"");
			}
			// PRISM_STPG_SOLN_METHOD
			s = settings.getString(PrismSettings.PRISM_STPG_SOLN_METHOD);
			if (s.equals("Value iteration")) {
				setSTPGSolnMethod(STPGSolnMethod.VALUE_ITERATION);
			} else if (s.equals("Gauss-Seidel")) {
				setSTPGSolnMethod(STPGSolnMethod.GAUSS_SEIDEL);
			} else {
				throw new PrismNotSupportedException("Explicit engine does not support STPG solution method \"" + s + "\"");
			}
			// PRISM_IMDP_SOLN_METHOD
			s = settings.getString(PrismSettings.PRISM_IMDP_SOLN_METHOD);
			if (s.equals("Value iteration")) {
				setIMDPSolnMethod(IMDPSolnMethod.VALUE_ITERATION);
			} else if (s.equals("Gauss-Seidel")) {
				setIMDPSolnMethod(IMDPSolnMethod.GAUSS_SEIDEL);
			} else {
				throw new PrismNotSupportedException("Explicit engine does not support IMDP solution method \"" + s + "\"");
			}
			// PRISM_TERM_CRIT
			s = settings.getString(PrismSettings.PRISM_TERM_CRIT);
			if (s.equals("Absolute")) {
				setTermCrit(TermCrit.ABSOLUTE);
			} else if (s.equals("Relative")) {
				setTermCrit(TermCrit.RELATIVE);
			} else {
				throw new PrismNotSupportedException("Unknown termination criterion \"" + s + "\"");
			}
			// PRISM_TERM_CRIT_PARAM
			setTermCritParam(settings.getDouble(PrismSettings.PRISM_TERM_CRIT_PARAM));
			// PRISM_MAX_ITERS
			setMaxIters(settings.getInteger(PrismSettings.PRISM_MAX_ITERS));
			// PRISM_GRID_RESOLUTION
			setGridResolution(settings.getInteger(PrismSettings.PRISM_GRID_RESOLUTION));
			// PRISM_PRECOMPUTATION
			setPrecomp(settings.getBoolean(PrismSettings.PRISM_PRECOMPUTATION));
			// PRISM_PROB0
			setProb0(settings.getBoolean(PrismSettings.PRISM_PROB0));
			// PRISM_PROB1
			setProb1(settings.getBoolean(PrismSettings.PRISM_PROB1));
			// PRISM_USE_PRE
			setPreRel(settings.getBoolean(PrismSettings.PRISM_PRE_REL));
			// PRISM_FAIRNESS
			if (settings.getBoolean(PrismSettings.PRISM_FAIRNESS)) {
				throw new PrismNotSupportedException("The explicit engine does not support model checking MDPs under fairness");
			}
		}
	}

	// Settings methods

	/**
	 * Inherit settings (and the log) from another ProbModelChecker object.
	 * For model checker objects that inherit a PrismSettings object, this is superfluous
	 * since this has been done already.
	 */
	public void inheritSettings(ProbModelChecker other)
	{
		super.inheritSettings(other);
		setLinEqMethod(other.getLinEqMethod());
		setMDPSolnMethod(other.getMDPSolnMethod());
		setSTPGSolnMethod(other.getSTPGSolnMethod());
		setIMDPSolnMethod(other.getIMDPSolnMethod());
		setTermCrit(other.getTermCrit());
		setTermCritParam(other.getTermCritParam());
		setMaxIters(other.getMaxIters());
		setGridResolution(other.getGridResolution());
		setPrecomp(other.getPrecomp());
		setProb0(other.getProb0());
		setProb1(other.getProb1());
		setValIterDir(other.getValIterDir());
		setSolnMethod(other.getSolnMethod());
		setErrorOnNonConverge(other.geterrorOnNonConverge());
	}

	/**
	 * Print summary of current settings.
	 */
	public void printSettings()
	{
		super.printSettings();
		mainLog.print("linEqMethod = " + linEqMethod + " ");
		mainLog.print("mdpSolnMethod = " + mdpSolnMethod + " ");
		mainLog.print("stpgSolnMethod = " + stpgSolnMethod + " ");
		mainLog.print("imdpSolnMethod = " + imdpSolnMethod + " ");
		mainLog.print("termCrit = " + termCrit + " ");
		mainLog.print("termCritParam = " + termCritParam + " ");
		mainLog.print("maxIters = " + maxIters + " ");
		mainLog.print("gridResolution = " + gridResolution + " ");
		mainLog.print("precomp = " + precomp + " ");
		mainLog.print("prob0 = " + prob0 + " ");
		mainLog.print("prob1 = " + prob1 + " ");
		mainLog.print("valIterDir = " + valIterDir + " ");
		mainLog.print("solnMethod = " + solnMethod + " ");
		mainLog.print("errorOnNonConverge = " + errorOnNonConverge + " ");
	}

	// Set methods for flags/settings

	/**
	 * Set verbosity level, i.e. amount of output produced.
	 */
	public void setVerbosity(int verbosity)
	{
		this.verbosity = verbosity;
	}

	/**
	 * Set flag for suppressing log output during precomputations (prob0, prob1, ...)
	 * @param value silent?
	 * @return the previous value of this flag
	 */
	public boolean setSilentPrecomputations(boolean value)
	{
		boolean old = silentPrecomputations;
		silentPrecomputations = value;
		return old;
	}

	/**
	 * Set method used to solve linear equation systems.
	 */
	public void setLinEqMethod(LinEqMethod linEqMethod)
	{
		this.linEqMethod = linEqMethod;
	}

	/**
	 * Set method used to solve MDPs.
	 */
	public void setMDPSolnMethod(MDPSolnMethod mdpSolnMethod)
	{
		this.mdpSolnMethod = mdpSolnMethod;
	}

	/**
	 * Set method used to solve STPGs.
	 */
	public void setSTPGSolnMethod(STPGSolnMethod stpgSolnMethod)
	{
		this.stpgSolnMethod = stpgSolnMethod;
	}

	/**
	 * Set method used to solve IMDPs (and IDTMCs).
	 */
	public void setIMDPSolnMethod(IMDPSolnMethod imdpSolnMethod)
	{
		this.imdpSolnMethod = imdpSolnMethod;
	}

	/**
	 * Set termination criteria type for numerical iterative methods.
	 */
	public void setTermCrit(TermCrit termCrit)
	{
		this.termCrit = termCrit;
	}

	/**
	 * Set termination criteria parameter (epsilon) for numerical iterative methods.
	 */
	public void setTermCritParam(double termCritParam)
	{
		this.termCritParam = termCritParam;
	}

	/**
	 * Set maximum number of iterations for numerical iterative methods.
	 */
	public void setMaxIters(int maxIters)
	{
		this.maxIters = maxIters;
	}

	/**
	 * Set resolution for POMDP fixed grid approximation algorithm.
	 */
	public void setGridResolution(int gridResolution)
	{
		this.gridResolution = gridResolution;
	}

	/**
	 * Set whether or not to use precomputation (Prob0, Prob1, etc.).
	 */
	public void setPrecomp(boolean precomp)
	{
		this.precomp = precomp;
	}

	/**
	 * Set whether or not to use Prob0 precomputation
	 */
	public void setProb0(boolean prob0)
	{
		this.prob0 = prob0;
	}

	/**
	 * Set whether or not to use Prob1 precomputation
	 */
	public void setProb1(boolean prob1)
	{
		this.prob1 = prob1;
	}

	/**
	 * Set whether or not to use pre-computed predecessor relation
	 */
	public void setPreRel(boolean preRel)
	{
		this.preRel = preRel;
	}

	/**
	 * Set direction of convergence for value iteration (lfp/gfp).
	 */
	public void setValIterDir(ValIterDir valIterDir)
	{
		this.valIterDir = valIterDir;
	}

	/**
	 * Set method used for numerical solution.
	 */
	public void setSolnMethod(SolnMethod solnMethod)
	{
		this.solnMethod = solnMethod;
	}

	/**
	 * Set whether non-convergence of an iterative method an error
	 */
	public void setErrorOnNonConverge(boolean errorOnNonConverge)
	{
		this.errorOnNonConverge = errorOnNonConverge;
	}

	// Get methods for flags/settings

	public int getVerbosity()
	{
		return verbosity;
	}

	public LinEqMethod getLinEqMethod()
	{
		return linEqMethod;
	}

	public MDPSolnMethod getMDPSolnMethod()
	{
		return mdpSolnMethod;
	}

	public STPGSolnMethod getSTPGSolnMethod()
	{
		return stpgSolnMethod;
	}

	public IMDPSolnMethod getIMDPSolnMethod()
	{
		return imdpSolnMethod;
	}

	public TermCrit getTermCrit()
	{
		return termCrit;
	}

	public double getTermCritParam()
	{
		return termCritParam;
	}

	public int getMaxIters()
	{
		return maxIters;
	}

	public int getGridResolution()
	{
		return gridResolution;
	}

	public boolean getPrecomp()
	{
		return precomp;
	}

	public boolean getProb0()
	{
		return prob0;
	}

	public boolean getProb1()
	{
		return prob1;
	}

	public boolean getPreRel()
	{
		return preRel;
	}

	public ValIterDir getValIterDir()
	{
		return valIterDir;
	}

	public SolnMethod getSolnMethod()
	{
		return solnMethod;
	}

	/**
	 * Is non-convergence of an iterative method an error?
	 */
	public boolean geterrorOnNonConverge()
	{
		return errorOnNonConverge;
	}

	// Model checking functions

	@Override
	public StateValues checkExpression(Model<?> model, Expression expr, BitSet statesOfInterest) throws PrismException
	{
		StateValues res;
		// <<>> or [[]] operator
		if (expr instanceof ExpressionStrategy) {
			res = checkExpressionStrategy(model, (ExpressionStrategy) expr, statesOfInterest);
		}
		// P operator
		else if (expr instanceof ExpressionProb) {
			res = checkExpressionProb(model, (ExpressionProb) expr, statesOfInterest);
		}
		// R operator
		else if (expr instanceof ExpressionReward) {
			res = checkExpressionReward(model, (ExpressionReward) expr, statesOfInterest);
		}
		// S operator
		else if (expr instanceof ExpressionSS) {
			res = checkExpressionSteadyState(model, (ExpressionSS) expr);
		}
		// Functions (for multi-objective)
		else if (expr instanceof ExpressionFunc) {
			res = checkExpressionFunc(model, (ExpressionFunc) expr, statesOfInterest);
		}
		// Otherwise, use the superclass
		else {
			res = super.checkExpression(model, expr, statesOfInterest);
		}

		return res;
	}

	/**
	 * Model check a <<>> or [[]] operator expression and return the values for the statesOfInterest.
	 * * @param statesOfInterest the states of interest, see checkExpression()
	 */
	protected StateValues checkExpressionStrategy(Model<?> model, ExpressionStrategy expr, BitSet statesOfInterest) throws PrismException
	{
		// Only support <<>> right now, not [[]]
		if (!expr.isThereExists())
			throw new PrismNotSupportedException("The " + expr.getOperatorString() + " operator is not yet supported");

		// Only support <<>> for MDPs/SMGs right now
		if (!(this instanceof MDPModelChecker || this instanceof SMGModelChecker || this instanceof CSGModelChecker))
			throw new PrismNotSupportedException("The " + expr.getOperatorString() + " operator is only supported for MDPs and SMGs currently");

		// Will we be quantifying universally or existentially over strategies/adversaries?
		boolean forAll = !expr.isThereExists();
		
		// Multiple coalitions only supported for CSGs
		if (expr.getNumCoalitions() > 1 && model.getModelType() != ModelType.CSG) {
			throw new PrismNotSupportedException("The " + expr.getOperatorString() + " operator can only contain one coalition");
		}
		// Extract coalition info
		Coalition coalition = expr.getCoalition();

		// For non-games (i.e., models with a single player), deal with the coalition operator here and then remove it
		if (coalition != null && !model.getModelType().multiplePlayers()) {
			if (coalition.isEmpty()) {
				// An empty coalition negates the quantification ("*" has no effect)
				forAll = !forAll;
			}
			coalition = null;
		}

		// For now, just support a single expression (which may encode a Boolean combination of objectives)
		List<Expression> exprs = expr.getOperands();
		if (exprs.size() > 1) {
			throw new PrismException("Cannot currently check strategy operators wth lists of expressions");
		}

		Expression exprSub = exprs.get(0);
		// Pass onto relevant method:
		// P operator
		if (exprSub instanceof ExpressionProb) {
			return checkExpressionProb(model, (ExpressionProb) exprSub, forAll, coalition, statesOfInterest);
		}
		// R operator
		else if (exprSub instanceof ExpressionReward) {
			if (((ExpressionReward) exprSub).getDiscount() != null) {
				useDiscounting = true;
				discountFactor = ((Expression) ((ExpressionReward) exprSub).getDiscount()).evaluateDouble();
			}
			return checkExpressionReward(model, (ExpressionReward) exprSub, forAll, coalition, statesOfInterest);
		}
		// Equilibria
		else if (exprSub instanceof ExpressionMultiNash) {
			return checkExpressionMultiNash(model, (ExpressionMultiNash) exprSub, expr.getCoalitions(), expr.getEquilibriumType(), expr.getEquilibriumCriterion());
		}
		// Anything else is treated as multi-objective
		else {
			return checkExpressionMultiObjective(model, expr, forAll, coalition);
		}
	}

	protected StateValues checkExpressionMultiNash(Model<?> model, ExpressionMultiNash expr, List<Coalition> coalitions, ExpressionStrategy.EquilibriumType equilibriumType, ExpressionStrategy.EquilibriumCriterion equilibriumCriterion) throws PrismException
	{
		ModelCheckerResult res = new ModelCheckerResult();
		List<ExpressionQuant> formulae = expr.getOperands();
		List<CSGRewards<Double>> rewards = new ArrayList<>();
		BitSet[] remain = new BitSet[coalitions.size()];
		BitSet[] targets = new BitSet[coalitions.size()];
		BitSet bounded  = new BitSet();
		BitSet unbounded = new BitSet();
		int[] bounds = new int[coalitions.size()];
		Expression expr1;
		IntegerBound bound;
		int e, p, r;
		boolean type = true;
		boolean rew = false;
		boolean min = expr.getRelOp().isMin();

		int eqType = CSGModelCheckerEquilibria.NASH;
		if (equilibriumType == ExpressionStrategy.EquilibriumType.NASH) {
			eqType = CSGModelCheckerEquilibria.NASH;
		} else if (equilibriumType == ExpressionStrategy.EquilibriumType.CORRELATED) {
			eqType = CSGModelCheckerEquilibria.CORR;
		}
		int crit = CSGModelCheckerEquilibria.SWEQ;
		if (equilibriumCriterion == ExpressionStrategy.EquilibriumCriterion.SOCIAL) {
			crit = CSGModelCheckerEquilibria.SWEQ;
		} else if (equilibriumCriterion == ExpressionStrategy.EquilibriumCriterion.FAIR) {
			crit = CSGModelCheckerEquilibria.FAIR;
		}

		Arrays.fill(remain, null);
		Arrays.fill(targets, null);
		Arrays.fill(bounds, -1);

		/*
		System.out.println("-- RelOp");
		System.out.println("isMax " + expr.getRelOp().isMax());
		System.out.println("isMin " + expr.getRelOp().isMin());
		System.out.println("isLowerBound " + expr.getRelOp().isLowerBound());
		System.out.println("isUpperBound " + expr.getRelOp().isUpperBound());
		System.out.println("isStrict " + expr.getRelOp().isStrict());
		*/

		OpRelOpBound opInfo = expr.getRelopBoundInfo(constantValues);
		MinMax minMax = opInfo.getMinMax(model.getModelType());

		/*
		System.out.println("-- minMax");
		System.out.println("isMax " + minMax.isMax());
		System.out.println("isMin " + minMax.isMin());
		System.out.println("isMin1 " + minMax.isMin1());
		System.out.println("isMin2 " + minMax.isMin2());
		*/

		if (coalitions.size() != formulae.size())
			throw new PrismException("The number of coalitions and objectives must be equal");

		for (e = 0; e < formulae.size() -1 ; e++) {
			type = type && formulae.get(e).getClass() == formulae.get(e+1).getClass();
		}

		if (!type)
			throw new PrismException("Mixing P and R operators is not yet supported");

		List<ExpressionTemporal> exprs = new ArrayList<ExpressionTemporal>();

		for (p = 0; p < formulae.size() ; p++) {
			expr1 = formulae.get(p);
			if (expr1 instanceof ExpressionMultiNashProb) {
				expr1 = ((ExpressionMultiNashProb) (expr1)).getExpression();
				expr1 = Expression.convertSimplePathFormulaToCanonicalForm(expr1);
				if (expr1 instanceof ExpressionTemporal) {
		 			ExpressionTemporal exprTemp = (ExpressionTemporal) expr1;
		 			exprs.add(p, exprTemp);
					switch (exprTemp.getOperator()) {
						case ExpressionTemporal.P_X:
							bounded.set(p);
							bounds[p] = 1;
							break;
						case ExpressionTemporal.P_U:
							if (exprTemp.hasBounds()) {
								bounded.set(p);
								bound = IntegerBound.fromExpressionTemporal((ExpressionTemporal) expr1, constantValues, true);
								if (bound.hasUpperBound())
									bounds[p] = bound.getHighestInteger();
								else
									throw new PrismException("Only upper bounds are allowed");
							} else {
								unbounded.set(p);
							}
							if (!((ExpressionTemporal) expr1).getOperand1().equals(ExpressionLiteral.True()))
								remain[p] = checkExpression(model, ((ExpressionTemporal) expr1).getOperand1(), null).getBitSet();
							break;
						default:
							throw new PrismNotSupportedException("The reward operator " + exprTemp.getOperatorSymbol() + " is not yet supported for equilibria-based properties");
					}
				}
				targets[p] = checkExpression(model, ((ExpressionTemporal) expr1).getOperand2(), null).getBitSet();
			}
			else if (expr1 instanceof ExpressionMultiNashReward) {
				r = ExpressionReward.getRewardStructIndexByIndexObject(((ExpressionMultiNashReward) expr1).getRewardStructIndex(), rewardGen, constantValues);
				expr1 = ((ExpressionMultiNashReward) (expr1)).getExpression();
				rewards.add(p, (CSGRewards<Double>) constructRewards(model, r));
				rew = true;
				ExpressionTemporal exprTemp = (ExpressionTemporal) expr1;
	 			exprs.add(p, exprTemp);
				switch (exprTemp.getOperator()) {
					case ExpressionTemporal.P_F:
						targets[p] = checkExpression(model, exprTemp.getOperand2(), null).getBitSet();
						unbounded.set(p);
						break;
					case ExpressionTemporal.R_I:
						bounded.set(p);
						bounds[p] = exprTemp.getUpperBound().evaluateInt(constantValues);
						break;
					case ExpressionTemporal.R_C:
						if (exprTemp.hasBounds()) {
							bounded.set(p);
							bound = IntegerBound.fromExpressionTemporal((ExpressionTemporal) expr1, constantValues, true);
							if (bound.hasUpperBound())
								bounds[p] = bound.getHighestInteger();
							else
								throw new PrismException("Only upper bounds are allowed");
						} else {
							throw new PrismNotSupportedException("Total rewards " + exprTemp.getOperatorSymbol() + " is not yet supported for equilibria-based properties");
						}
						break;
					default:
						throw new PrismNotSupportedException("The reward operator " + exprTemp.getOperatorSymbol() + " is not yet supported for equilibria-based properties");
				}
			}
		}
		if (!rew)
			rewards = null;
		if (coalitions.size() == 1) {
			throw new PrismNotSupportedException("Equilibria-based properties require at least two coalitions");
		}
		else if (coalitions.size() == 2) {
			if (unbounded.cardinality() == formulae.size()) {
				if (rew) {
					res = ((CSGModelChecker) this).computeRewReachEquilibria((CSG<Double>) model, coalitions, rewards, targets, eqType, crit, min);
				}
				else {
					res = ((CSGModelChecker) this).computeProbReachEquilibria((CSG<Double>) model, coalitions, targets, remain, eqType, crit, min);
				}
			}
			else if (bounded.cardinality() == formulae.size()) {
				if (rew) {
					res = ((CSGModelChecker) this).computeRewBoundedEquilibria((CSG<Double>) model, coalitions, rewards, exprs, bounds, eqType, crit, min);
				}
				else {
					res = ((CSGModelChecker) this).computeProbBoundedEquilibria((CSG<Double>) model, coalitions, exprs, targets, remain, bounds, eqType, crit, min);
				}
			}
			else {
				res = ((CSGModelChecker) this).computeMixedEquilibria((CSG<Double>) model, coalitions, rewards, exprs, bounded, targets, remain, bounds, eqType, crit, min);
			}
		}
		else if (coalitions.size() > 2) {
			// throw new PrismNotSupportedException("Equilibria-based properties with more than two coalitions are not yet supported");
			if (unbounded.cardinality() == formulae.size()) {
				if (rew) {
					res = ((CSGModelChecker) this).computeMultiRewReachEquilibria((CSG) model, coalitions, rewards, targets, eqType, crit, min);
				}
				else {
					res = ((CSGModelChecker) this).computeMultiProbReachEquilibria((CSG) model, coalitions, targets, remain, eqType, crit, min);
				}
			}
			else if (bounded.cardinality() == formulae.size()) {
				int first_bound = bounds[0];
				for (int i = 1; i < bounds.length; i++) {
					if (first_bound != bounds[i])
						throw new PrismNotSupportedException("Multi-coalitional properties with different bounds are not yet supported");
				}
				if (rew) {
					res = ((CSGModelChecker) this).computeMultiRewBoundedEquilibria((CSG) model, coalitions, rewards, exprs, bounds, eqType, crit, min);
				}
				else {
					res = ((CSGModelChecker) this).computeMultiProbBoundedEquilibria((CSG) model, coalitions, targets, remain, bounds, eqType, crit, min);
				}
			}
			else {
				throw new PrismNotSupportedException("Multi-coalitional mixed properties are not yet supported");
			}
		}

		result.setStrategy(res.strat);
		StateValues sv = StateValues.createFromDoubleArray(res.soln, model);

		// Print out probabilities
		if (getVerbosity() > 5) {
			mainLog.print("\nProbabilities (non-zero only) for all states:\n");
			sv.print(mainLog);
		}

		// For =? properties, just return values; otherwise compare against bound
		if (!opInfo.isNumeric()) {
			sv.applyPredicate(v -> opInfo.apply((double) v, sv.getAccuracy()));
		}
		return sv;
	}

	/**
	 * Model check a P operator expression and return the values for the statesOfInterest.
 	 * @param statesOfInterest the states of interest, see checkExpression()
	 */
	protected StateValues checkExpressionProb(Model<?> model, ExpressionProb expr, BitSet statesOfInterest) throws PrismException
	{
		// Use the default semantics for a standalone P operator
		// (i.e. quantification over all strategies, and no game-coalition info)
		return checkExpressionProb(model, expr, true, null, statesOfInterest);
	}
	
	/**
	 * Model check a P operator expression and return the values for the states of interest.
	 * @param model The model
	 * @param expr The P operator expression
	 * @param forAll Are we checking "for all strategies" (true) or "there exists a strategy" (false)? [irrelevant for numerical (=?) queries] 
	 * @param coalition If relevant, info about which set of players this P operator refers to (null if irrelevant)
	 * @param statesOfInterest the states of interest, see checkExpression()
	 */
	protected StateValues checkExpressionProb(Model<?> model, ExpressionProb expr, boolean forAll, Coalition coalition, BitSet statesOfInterest) throws PrismException
	{
		// For now, need separate handling of S and C operator for SMGs
		if (expr.getExpression() instanceof ExpressionReward) {
			Expression e = ((ExpressionReward) expr.getExpression()).getExpression();
			if (e.getType() instanceof TypePathDouble) {
				ExpressionTemporal eTemp = (ExpressionTemporal) e;
				if (model.getModelType() == ModelType.SMG) {
					switch (eTemp.getOperator()) {
					case ExpressionTemporal.R_S: // average rewards
						return ((SMGModelChecker) this).checkExpressionMultiObjective(model,
								BooleanUtils.convertToCNFLists(expr), coalition);
					case ExpressionTemporal.R_C: // total or ratio rewards
						if (!eTemp.hasBounds()) {
							return ((SMGModelChecker) this).checkExpressionMultiObjective(model,
									BooleanUtils.convertToCNFLists(expr), coalition);
						}
					}
				}
			}
		}

		// Get info from P operator
		OpRelOpBound opInfo = expr.getRelopBoundInfo(constantValues);
		MinMax minMax = opInfo.getMinMax(model.getModelType(), forAll, coalition);

		// Compute probabilities
		StateValues probs = checkProbPathFormula(model, expr.getExpression(), minMax, statesOfInterest);

		// Print out probabilities
		if (getVerbosity() > 5) {
			mainLog.print("\nProbabilities (non-zero only) for all states:\n");
			probs.print(mainLog);
		}

		// For =? properties, just return values; otherwise compare against bound
		if (!opInfo.isNumeric()) {
			probs.applyPredicate(v -> opInfo.apply((double) v, probs.getAccuracy()));
		}
		return probs;
	}

	/**
	 * Compute probabilities for the contents of a P operator.
	 * @param statesOfInterest the states of interest, see checkExpression()
	 */
	protected StateValues checkProbPathFormula(Model<?> model, Expression expr, MinMax minMax, BitSet statesOfInterest) throws PrismException
	{
		// If the path formula is of the form R(path)~r, deal with it
		if (expr instanceof ExpressionReward && "path".equals(((ExpressionReward) expr).getModifier())) {
			return checkProbBoundedRewardFormula(model, (ExpressionReward) expr, minMax, statesOfInterest);
		}

		// Test whether this is a simple path formula (i.e. PCTL)
		// and whether we want to use the corresponding algorithms
		boolean useSimplePathAlgo = expr.isSimplePathFormula();

		if (useSimplePathAlgo &&
		    settings.getBoolean(PrismSettings.PRISM_PATH_VIA_AUTOMATA) &&
		    LTLModelChecker.isSupportedLTLFormula(model.getModelType(), expr)) {
			// If PRISM_PATH_VIA_AUTOMATA is true, we want to use the LTL engine
			// whenever possible
			useSimplePathAlgo = false;
		}

		if (useSimplePathAlgo) {
			return checkProbPathFormulaSimple(model, expr, minMax, statesOfInterest);
		} else {
			// Some model checkers will behave differently for cosafe vs full LTL
			if (Expression.isCoSafeLTLSyntactic(expr, true)) {
				return checkProbPathFormulaCosafeLTL(model, expr, false, minMax, statesOfInterest);
			} else {
				return checkProbPathFormulaLTL(model, expr, false, minMax, statesOfInterest);
			}
		}
	}

	/**
	 * Compute probabilities for a simple, non-LTL path operator.
	 */
	protected StateValues checkProbPathFormulaSimple(Model<?> model, Expression expr, MinMax minMax, BitSet statesOfInterest) throws PrismException
	{
		boolean negated = false;
		StateValues probs = null;

		expr = Expression.convertSimplePathFormulaToCanonicalForm(expr);

		// Negation
		if (expr instanceof ExpressionUnaryOp &&
		    ((ExpressionUnaryOp)expr).getOperator() == ExpressionUnaryOp.NOT) {
			negated = true;
			minMax = minMax.negate();
			expr = ((ExpressionUnaryOp)expr).getOperand();
		}

		if (expr instanceof ExpressionTemporal) {
 			ExpressionTemporal exprTemp = (ExpressionTemporal) expr;
			// Next
			if (exprTemp.getOperator() == ExpressionTemporal.P_X) {
				probs = checkProbNext(model, exprTemp, minMax, statesOfInterest);
			}
			// Until
			else if (exprTemp.getOperator() == ExpressionTemporal.P_U) {
				if (exprTemp.hasBounds()) {
					probs = checkProbBoundedUntil(model, exprTemp, minMax, statesOfInterest);
				} else {
					probs = checkProbUntil(model, exprTemp, minMax, statesOfInterest);
				}
			}
		}

		if (probs == null)
			throw new PrismException("Unrecognised path operator in P operator");

		if (negated) {
			// Subtract from 1 for negation
			probs.applyFunction(TypeDouble.getInstance(), v -> 1.0 - (double) v);
		}

		return probs;
	}

	/**
	 * Compute probabilities for a next operator.
	 */
	@SuppressWarnings("unchecked")
	protected StateValues checkProbNext(Model<?> model, ExpressionTemporal expr, MinMax minMax, BitSet statesOfInterest) throws PrismException
	{
		// Model check the operand for all states
		BitSet target = checkExpression(model, expr.getOperand2(), null).getBitSet();

		// Compute/return the probabilities
		ModelCheckerResult res = null;
		switch (model.getModelType()) {
		case CTMC:
			res = ((CTMCModelChecker) this).computeNextProbs((CTMC<Double>) model, target);
			break;
		case CSG:
			res = ((CSGModelChecker) this).computeNextProbs((CSG<Double>) model, target, minMax.isMin1(), minMax.isMin2(), minMax.getCoalition());
			break;
		case DTMC:
			res = ((DTMCModelChecker) this).computeNextProbs((DTMC<Double>) model, target);
			break;
		case MDP:
			res = ((MDPModelChecker) this).computeNextProbs((MDP<Double>) model, target, minMax.isMin());
			break;
		case STPG:
			res = ((STPGModelChecker) this).computeNextProbs((STPG<Double>) model, target, minMax.isMin1(), minMax.isMin2());
			break;
		case SMG:
			res = ((SMGModelChecker) this).computeNextProbs((SMG<Double>) model, target, minMax.isMin1(), minMax.isMin2(), minMax.getCoalition());
			break;
		case IDTMC:
			res = ((IDTMCModelChecker) this).computeNextProbs((IDTMC<Double>) model, target, minMax);
			break;
		case IMDP:
			res = ((IMDPModelChecker) this).computeNextProbs((IMDP<Double>) model, target, minMax);
			break;
		default:
			throw new PrismNotSupportedException("Cannot model check " + expr + " for " + model.getModelType() + "s");
		}
		result.setStrategy(res.strat);
		return StateValues.createFromArrayResult(res, model);
	}

	/**
	 * Compute probabilities for a bounded until operator.
	 */
	@SuppressWarnings("unchecked")
	protected StateValues checkProbBoundedUntil(Model<?> model, ExpressionTemporal expr, MinMax minMax, BitSet statesOfInterest) throws PrismException
	{
		// This method just handles discrete time
		// Continuous-time model checkers will override this method

		// Get info from bounded until
		Integer lowerBound;
		IntegerBound bounds;
		int i;

		// get and check bounds information
		bounds = IntegerBound.fromExpressionTemporal(expr, constantValues, true);

		// Model check operands for all states
		BitSet remain = checkExpression(model, expr.getOperand1(), null).getBitSet();
		BitSet target = checkExpression(model, expr.getOperand2(), null).getBitSet();

		if (bounds.hasLowerBound()) {
			lowerBound = bounds.getLowestInteger();
		} else {
			lowerBound = 0;
		}

		Integer windowSize = null;  // unbounded

		if (bounds.hasUpperBound()) {
			windowSize = bounds.getHighestInteger() - lowerBound;
		}

		// compute probabilities for Until<=windowSize
		StateValues sv = null;

		if (windowSize == null) {
			// unbounded
			ModelCheckerResult res = null;
			switch (model.getModelType()) {
			case DTMC:
				res = ((DTMCModelChecker) this).computeUntilProbs((DTMC<Double>) model, remain, target);
				break;
			case MDP:
				res = ((MDPModelChecker) this).computeUntilProbs((MDP<Double>) model, remain, target, minMax.isMin());
				break;
			case STPG:
				res = ((STPGModelChecker) this).computeUntilProbs((STPG<Double>) model, remain, target, minMax.isMin1(), minMax.isMin2(), minMax.getBound());
				break;
			case SMG:
				res = ((SMGModelChecker) this).computeUntilProbs((SMG<Double>) model, remain, target, minMax.isMin1(), minMax.isMin2(), minMax.getCoalition());
				break;
			case CSG:
				res = ((CSGModelChecker) this).computeUntilProbs((CSG<Double>) model, remain, target, minMax.isMin1(), minMax.isMin2(), minMax.getCoalition());
				break;
			case IDTMC:
				res = ((IDTMCModelChecker) this).computeUntilProbs((IDTMC<Double>) model, remain, target, minMax);
				break;
			case IMDP:
				res = ((IMDPModelChecker) this).computeUntilProbs((IMDP<Double>) model, remain, target, minMax);
				break;
			default:
				throw new PrismException("Cannot model check " + expr + " for " + model.getModelType() + "s");
			}
			result.setStrategy(res.strat);
			sv = StateValues.createFromArrayResult(res, model);
		} else if (windowSize == 0) {
			// A trivial case: windowSize=0 (prob is 1 in target states, 0 otherwise)
			sv = StateValues.createFromBitSetAsDoubles(target, model);
		} else {
			// Otherwise: numerical solution
			ModelCheckerResult res = null;
			switch (model.getModelType()) {
			case DTMC:
				res = ((DTMCModelChecker) this).computeBoundedUntilProbs((DTMC<Double>) model, remain, target, windowSize);
				break;
			case MDP:
				res = ((MDPModelChecker) this).computeBoundedUntilProbs((MDP<Double>) model, remain, target, windowSize, minMax.isMin());
				break;
			case STPG:
				res = ((STPGModelChecker) this).computeBoundedUntilProbs((STPG<Double>) model, remain, target, windowSize, minMax.isMin1(), minMax.isMin2());
				break;
			case SMG:
				res = ((SMGModelChecker) this).computeBoundedUntilProbs((SMG<Double>) model, remain, target, windowSize, minMax.isMin1(), minMax.isMin2(), minMax.getCoalition());
				break;
			case CSG:
				res = ((CSGModelChecker) this).computeBoundedUntilProbs((CSG<Double>) model, remain, target, windowSize, minMax.isMin1(), minMax.isMin2(), minMax.getCoalition());
				break;
			case IDTMC:
				res = ((IDTMCModelChecker) this).computeBoundedUntilProbs((IDTMC<Double>) model, remain, target, windowSize, minMax);
				break;
			case IMDP:
				res = ((IMDPModelChecker) this).computeBoundedUntilProbs((IMDP<Double>) model, remain, target, windowSize, minMax);
				break;
			default:
				throw new PrismNotSupportedException("Cannot model check " + expr + " for " + model.getModelType() + "s");
			}
			result.setStrategy(res.strat);
			sv = StateValues.createFromArrayResult(res, model);
		}

		// perform lowerBound restricted next-step computations to
		// deal with lower bound.
		if (lowerBound > 0) {
			double[] probs = sv.getDoubleArray();
			for (i = 0; i < lowerBound; i++) {
				switch (model.getModelType()) {
				case DTMC:
					probs = ((DTMCModelChecker) this).computeRestrictedNext((DTMC<Double>) model, remain, probs);
					break;
				case MDP:
					probs = ((MDPModelChecker) this).computeRestrictedNext((MDP<Double>) model, remain, probs, minMax.isMin());
					break;
				case STPG:
					// TODO (JK): Figure out if we can handle lower bounds for STPG in the same way
					throw new PrismNotSupportedException("Lower bounds not yet supported for STPGModelChecker");
				case SMG:
					throw new PrismException("Lower bounds not yet supported for SMGs");
				default:
					throw new PrismNotSupportedException("Cannot model check " + expr + " for " + model.getModelType() + "s");
				}
			}

			sv = StateValues.createFromDoubleArray(probs, model);
			sv.setAccuracy(AccuracyFactory.boundedNumericalIterations());
		}

		return sv;
	}

	/**
	 * Compute probabilities for an (unbounded) until operator.
	 */
	@SuppressWarnings("unchecked")
	protected StateValues checkProbUntil(Model<?> model, ExpressionTemporal expr, MinMax minMax, BitSet statesOfInterest) throws PrismException
	{
		// Model check operands for all states
		BitSet remain = checkExpression(model, expr.getOperand1(), null).getBitSet();
		BitSet target = checkExpression(model, expr.getOperand2(), null).getBitSet();

		// Compute/return the probabilities
		ModelCheckerResult res = null;
		switch (model.getModelType()) {
		case CTMC:
			res = ((CTMCModelChecker) this).computeUntilProbs((CTMC<Double>) model, remain, target);
			break;
		case DTMC:
			res = ((DTMCModelChecker) this).computeUntilProbs((DTMC<Double>) model, remain, target);
			break;
		case MDP:
			res = ((MDPModelChecker) this).computeUntilProbs((MDP<Double>) model, remain, target, minMax.isMin());
			break;
		case POMDP:
			res = ((POMDPModelChecker) this).computeReachProbs((POMDP<Double>) model, remain, target, minMax.isMin(), statesOfInterest);
			break;
		case STPG:
			res = ((STPGModelChecker) this).computeUntilProbs((STPG<Double>) model, remain, target, minMax.isMin1(), minMax.isMin2(), minMax.getBound());
			break;
		case SMG:
			res = ((SMGModelChecker) this).computeUntilProbs((SMG<Double>) model, remain, target, minMax.isMin1(), minMax.isMin2(), minMax.getCoalition());
			break;
		case CSG:
			res = ((CSGModelChecker) this).computeUntilProbs((CSG<Double>) model, remain, target, minMax.isMin1(), minMax.isMin2(), minMax.getCoalition());
			break;
		case IDTMC:
			res = ((IDTMCModelChecker) this).computeUntilProbs((IDTMC<Double>) model, remain, target, minMax);
			break;
		case IMDP:
			res = ((IMDPModelChecker) this).computeUntilProbs((IMDP<Double>) model, remain, target, minMax);
			break;
		default:
			throw new PrismNotSupportedException("Cannot model check " + expr + " for " + model.getModelType() + "s");
		}
		result.setStrategy(res.strat);
		return StateValues.createFromArrayResult(res, model);
	}

	/**
	 * Compute probabilities for an LTL path formula
	 */
	protected StateValues checkProbPathFormulaLTL(Model<?> model, Expression expr, boolean qual, MinMax minMax, BitSet statesOfInterest) throws PrismException
	{
		// To be overridden by subclasses
		throw new PrismNotSupportedException("Computation not implemented yet");
	}

	/**
	 * Compute probabilities for a co-safe LTL path formula
	 */
	protected StateValues checkProbPathFormulaCosafeLTL(Model<?> model, Expression expr, boolean qual, MinMax minMax, BitSet statesOfInterest) throws PrismException
	{
		// Just treat as an arbitrary LTL formula by default
		return checkProbPathFormulaLTL(model, expr, qual, minMax, statesOfInterest);
	}

	/**
	 * Compute probabilities for an LTL path formula
	 */
	protected StateValues checkProbBoundedRewardFormula(Model model, ExpressionReward expr, MinMax minMax, BitSet statesOfInterest) throws PrismException
	{
		// No support for this by default
		throw new PrismNotSupportedException("Reward-bounded path formulas not yet supported");
	}

	/**
	 * Model check an R operator expression and return the values for all states.
	 */
	protected StateValues checkExpressionReward(Model<?> model, ExpressionReward expr, BitSet statesOfInterest) throws PrismException
	{
		// Use the default semantics for a standalone R operator
		// (i.e. quantification over all strategies, and no game-coalition info)
		return checkExpressionReward(model, expr, true, null, statesOfInterest);
	}
	
	/**
	 * Model check an R operator expression and return the values for all states.
	 */
	protected StateValues checkExpressionReward(Model<?> model, ExpressionReward expr, boolean forAll, Coalition coalition, BitSet statesOfInterest) throws PrismException
	{

		// For now, need separate handling of S and C operator for SMGs
		Expression e = expr.getExpression();
		if (e.getType() instanceof TypePathDouble) {
			ExpressionTemporal eTemp = (ExpressionTemporal) e;
			if (model.getModelType() == ModelType.SMG) {
				switch (eTemp.getOperator()) {
				case ExpressionTemporal.R_S: // average rewards
					return ((SMGModelChecker) this).checkExpressionMultiObjective(model,
							BooleanUtils.convertToCNFLists(expr), coalition);
				case ExpressionTemporal.R_C: // total or ratio rewards
					if (!eTemp.hasBounds()) {
						return ((SMGModelChecker) this).checkExpressionMultiObjective(model,
								BooleanUtils.convertToCNFLists(expr), coalition);
					}
				}
			}
		}

		// Check if ratio reward
		if (expr.getRewardStructIndexDiv() != null) {
			throw new PrismException("Ratio rewards not supported with the selected engine and module type.");
		}
		OpRelOpBound opInfo = expr.getRelopBoundInfo(constantValues);
		MinMax minMax = opInfo.getMinMax(model.getModelType(), forAll, coalition);

		// Build rewards
		int r = expr.getRewardStructIndexByIndexObject(rewardGen, constantValues);
		mainLog.println("Building reward structure...");
		Rewards<?> rewards = Expression.usesInstantaneousReward(expr.getExpression()) ? constructRewards(model, r) : constructExpectedRewards(model, r);

		// Compute rewards
		StateValues rews = checkRewardFormula(model, rewards, expr.getExpression(), minMax, statesOfInterest);

		// Print out rewards
		if (getVerbosity() > 5) {
			mainLog.print("\nRewards (non-zero only) for all states:\n");
			rews.print(mainLog);
		}

		// For =? properties, just return values; otherwise compare against bound
		if (!opInfo.isNumeric()) {
			rews.applyPredicate(v -> opInfo.apply((double) v, rews.getAccuracy()));
		}
		return rews;
	}

	/**
	 * Compute rewards for the contents of an R operator.
	 */
	protected StateValues checkRewardFormula(Model<?> model, Rewards<?> modelRewards, Expression expr, MinMax minMax, BitSet statesOfInterest) throws PrismException
	{
		StateValues rewards = null;

		if (expr.getType() instanceof TypePathDouble) {
			ExpressionTemporal exprTemp = (ExpressionTemporal) expr;
			switch (exprTemp.getOperator()) {
			case ExpressionTemporal.R_F0:
			case ExpressionTemporal.R_Fc:
				rewards = checkRewardReach(model, modelRewards, exprTemp, minMax, statesOfInterest);
				break;
			case ExpressionTemporal.R_I:
				rewards = checkRewardInstantaneous(model, modelRewards, exprTemp, minMax, statesOfInterest);
				break;
			case ExpressionTemporal.R_C:
				if (exprTemp.hasBounds()) {
					rewards = checkRewardCumulative(model, modelRewards, exprTemp, minMax);
				} else {
					rewards = checkRewardTotal(model, modelRewards, exprTemp, minMax);
				}
				break;
			case ExpressionTemporal.R_S:
				rewards = checkRewardSteady(model, modelRewards);
				break;
			default:
				throw new PrismNotSupportedException("Explicit engine does not yet handle the " + exprTemp.getOperatorSymbol() + " reward operator");
			}
		} else if (expr.getType() instanceof TypePathBool || expr.getType() instanceof TypeBool) {
			rewards = checkRewardPathFormula(model, modelRewards, expr, minMax, statesOfInterest);
		}

		if (rewards == null)
			throw new PrismException("Unrecognised operator in R operator");

		return rewards;
	}

	/**
	 * Compute rewards for an instantaneous reward operator.
	 */
	@SuppressWarnings("unchecked")
	protected StateValues checkRewardInstantaneous(Model<?> model, Rewards<?> modelRewards, ExpressionTemporal expr, MinMax minMax, BitSet statesOfInterest) throws PrismException
	{
		// Compute/return the rewards
		ModelCheckerResult res = null;
		int k = -1;
		double t = -1;
		if (model.getModelType().continuousTime()) {
			t = expr.getUpperBound().evaluateDouble(constantValues);
		} else {
			k = expr.getUpperBound().evaluateInt(constantValues);
		}
		switch (model.getModelType()) {
		case DTMC:
			res = ((DTMCModelChecker) this).computeInstantaneousRewards((DTMC<Double>) model, (MCRewards<Double>) modelRewards, k, statesOfInterest);
			break;
		case CTMC:
			res = ((CTMCModelChecker) this).computeInstantaneousRewards((CTMC<Double>) model, (MCRewards<Double>) modelRewards, t);
			break;
		case MDP:
			res = ((MDPModelChecker) this).computeInstantaneousRewards((MDP<Double>) model, (MDPRewards<Double>) modelRewards, k, minMax.isMin());
			break;
		case STPG:
			res = ((STPGModelChecker) this).computeInstantaneousRewards((STPG<Double>) model, (STPGRewards<Double>) modelRewards, k, minMax.isMin1(), minMax.isMin2());
			break;
		case SMG:
			res = ((SMGModelChecker) this).computeInstantaneousRewards((SMG<Double>) model, (Rewards<Double>) modelRewards, k, minMax.isMin1(), minMax.isMin2(), minMax.getCoalition());
			break;
		case CSG:
			res = ((CSGModelChecker) this).computeInstantaneousRewards((CSG<Double>) model, (CSGRewards<Double>) modelRewards, minMax.coalition, k, minMax.isMin1(), minMax.isMin2());
			break;
		default:
			throw new PrismNotSupportedException("Explicit engine does not yet handle the " + expr.getOperatorSymbol() + " reward operator for " + model.getModelType()
					+ "s");
		}
		result.setStrategy(res.strat);
		return StateValues.createFromArrayResult(res, model);
	}

	/**
	 * Compute rewards for a cumulative reward operator.
	 */
	@SuppressWarnings("unchecked")
	protected StateValues checkRewardCumulative(Model<?> model, Rewards<?> modelRewards, ExpressionTemporal expr, MinMax minMax) throws PrismException
	{
		int timeInt = -1;
		double timeDouble = -1;

		// Check that there is an upper time bound
		if (expr.getUpperBound() == null) {
			throw new PrismNotSupportedException("This is not a cumulative reward operator");
		}

		// Get time bound
		if (model.getModelType().continuousTime()) {
			timeDouble = expr.getUpperBound().evaluateDouble(constantValues);
			if (timeDouble < 0) {
				throw new PrismException("Invalid time bound " + timeDouble + " in cumulative reward formula");
			}
		} else {
			timeInt = expr.getUpperBound().evaluateInt(constantValues);
			if (timeInt < 0) {
				throw new PrismException("Invalid time bound " + timeInt + " in cumulative reward formula");
			}
		}

		// Compute/return the rewards
		// A trivial case: "C<=0" (prob is 1 in target states, 0 otherwise)
		if (timeInt == 0 || timeDouble == 0) {
			StateValues res = StateValues.createFromSingleValue(TypeDouble.getInstance(), 0.0, model);
			res.setAccuracy(AccuracyFactory.doublesFromQualitative());
			return res;
		}
		// Otherwise: numerical solution
		ModelCheckerResult res = null;
		switch (model.getModelType()) {
		case DTMC:
			res = ((DTMCModelChecker) this).computeCumulativeRewards((DTMC<Double>) model, (MCRewards<Double>) modelRewards, timeInt);
			break;
		case CTMC:
			res = ((CTMCModelChecker) this).computeCumulativeRewards((CTMC<Double>) model, (MCRewards<Double>) modelRewards, timeDouble);
			break;
		case MDP:
			res = ((MDPModelChecker) this).computeCumulativeRewards((MDP<Double>) model, (MDPRewards<Double>) modelRewards, timeInt, minMax.isMin());
			break;
		case STPG:
			res = ((STPGModelChecker) this).computeCumulativeRewards((STPG<Double>) model, (STPGRewards<Double>) modelRewards, timeInt, minMax.isMin1(), minMax.isMin2());
			break;
		case SMG:
			res = ((SMGModelChecker) this).computeCumulativeRewards((SMG<Double>) model, (Rewards<Double>) modelRewards, timeInt, minMax.isMin1(), minMax.isMin2(), minMax.getCoalition());
			break;
		case CSG:
			res = ((CSGModelChecker) this).computeCumulativeRewards((CSG<Double>) model, (CSGRewards<Double>) modelRewards, minMax.getCoalition(), timeInt, minMax.isMin1(), minMax.isMin2(), false);
			break;
		default:
			throw new PrismNotSupportedException("Explicit engine does not yet handle the " + expr.getOperatorSymbol() + " reward operator for " + model.getModelType()
					+ "s");
		}
		result.setStrategy(res.strat);
		return StateValues.createFromArrayResult(res, model);
	}

	/**
	 * Compute expected rewards for a total reward operator.
	 */
	@SuppressWarnings("unchecked")
	protected StateValues checkRewardTotal(Model<?> model, Rewards<?> modelRewards, ExpressionTemporal expr, MinMax minMax) throws PrismException
	{
		// Check that there is no upper time bound
		if (expr.getUpperBound() != null) {
			throw new PrismException("This is not a total reward operator");
		}

		// Compute/return the rewards
		ModelCheckerResult res = null;
		switch (model.getModelType()) {
		case DTMC:
			res = ((DTMCModelChecker) this).computeTotalRewards((DTMC<Double>) model, (MCRewards<Double>) modelRewards);
			break;
		case CTMC:
			res = ((CTMCModelChecker) this).computeTotalRewards((CTMC<Double>) model, (MCRewards<Double>) modelRewards);
			break;
		case MDP:
			res = ((MDPModelChecker) this).computeTotalRewards((MDP<Double>) model, (MDPRewards<Double>) modelRewards, minMax.isMin());
			break;
		case CSG:
			res = ((CSGModelChecker) this).computeTotalRewards((CSG<Double>) model, (CSGRewards<Double>) modelRewards, minMax.isMin1(), minMax.isMin2(), minMax.getCoalition());
			break;
		default:
			throw new PrismNotSupportedException("Explicit engine does not yet handle the " + expr.getOperatorSymbol() + " reward operator for " + model.getModelType()
					+ "s");
		}
		result.setStrategy(res.strat);
		return StateValues.createFromArrayResult(res, model);
	}

	/**
	 * Compute expected rewards for a steady-state reward operator.
	 */
	@SuppressWarnings("unchecked")
	protected StateValues checkRewardSteady(Model<?> model, Rewards<?> modelRewards) throws PrismException
	{
		// Compute/return the rewards
		ModelCheckerResult res = null;
		switch (model.getModelType()) {
		case DTMC:
			res = ((DTMCModelChecker) this).computeSteadyStateRewards((DTMC<Double>) model, (MCRewards<Double>) modelRewards);
			break;
		case CTMC:
			res = ((CTMCModelChecker) this).computeSteadyStateRewards((CTMC<Double>) model, (MCRewards<Double>) modelRewards);
			break;
		default:
			throw new PrismNotSupportedException("Explicit engine does not yet handle the steady-state reward operator for " + model.getModelType() + "s");
		}
		result.setStrategy(res.strat);
		return StateValues.createFromArrayResult(res, model);
	}

	/**
	 * Compute rewards for a path formula in a reward operator.
	 */
	protected StateValues checkRewardPathFormula(Model<?> model, Rewards<?> modelRewards, Expression expr, MinMax minMax, BitSet statesOfInterest) throws PrismException
	{
		if (Expression.isReach(expr)) {
			return checkRewardReach(model, modelRewards, (ExpressionTemporal) expr, minMax, statesOfInterest);
		}
		else if (Expression.isCoSafeLTLSyntactic(expr, true)) {
			return checkRewardCoSafeLTL(model, modelRewards, expr, minMax, statesOfInterest);
		}
		throw new PrismException("R operator contains a path formula that is not syntactically co-safe: " + expr);
	}
	
	/**
	 * Compute rewards for a reachability reward operator.
	 */
	@SuppressWarnings("unchecked")
	protected StateValues checkRewardReach(Model<?> model, Rewards<?> modelRewards, ExpressionTemporal expr, MinMax minMax, BitSet statesOfInterest) throws PrismException
	{
		// Non-game models don't yet support other variants of R[F]
		if (expr.getOperator() != ExpressionTemporal.P_F) {
			if (!(model.getModelType() == ModelType.STPG || model.getModelType() == ModelType.SMG || model.getModelType() == ModelType.CSG)) {
				throw new PrismException("The " + expr.getOperatorSymbol() + " reward operator only works for game models");
			}
		}

		// No time bounds allowed
		if (expr.hasBounds()) {
			throw new PrismNotSupportedException("R operator cannot contain a bounded F operator: " + expr);
		}
		
		// Model check the operand for all states
		BitSet target = checkExpression(model, expr.getOperand2(), null).getBitSet();

		// Compute/return the rewards
		ModelCheckerResult res = null;
		switch (model.getModelType()) {
		case DTMC:
			res = ((DTMCModelChecker) this).computeReachRewards((DTMC<Double>) model, (MCRewards<Double>) modelRewards, target);
			break;
		case CTMC:
			res = ((CTMCModelChecker) this).computeReachRewards((CTMC<Double>) model, (MCRewards<Double>) modelRewards, target);
			break;
		case MDP:
			res = ((MDPModelChecker) this).computeReachRewards((MDP<Double>) model, (MDPRewards<Double>) modelRewards, target, minMax.isMin());
			break;
		case POMDP:
			res = ((POMDPModelChecker) this).computeReachRewards((POMDP<Double>) model, (MDPRewards<Double>) modelRewards, target, minMax.isMin(), statesOfInterest);
			break;
		case STPG:
			res = ((STPGModelChecker) this).computeReachRewards((STPG<Double>) model, (STPGRewards<Double>) modelRewards, target, minMax.isMin1(), minMax.isMin2());
			break;
		case SMG:
			switch (expr.getOperator()) {
			case ExpressionTemporal.P_F:
				res = ((SMGModelChecker) this).computeReachRewards((SMG<Double>) model, (Rewards<Double>) modelRewards, target, STPGModelChecker.R_INFINITY, minMax.isMin1(), minMax.isMin2(), minMax.getCoalition());
				break;
			case ExpressionTemporal.R_Fc:
				res = ((SMGModelChecker) this).computeReachRewards((SMG<Double>) model, (Rewards<Double>) modelRewards, target, STPGModelChecker.R_CUMULATIVE, minMax.isMin1(), minMax.isMin2(), minMax.getCoalition());
				break;
			case ExpressionTemporal.R_F0:
				res = ((SMGModelChecker) this).computeReachRewards((SMG<Double>) model, (Rewards<Double>) modelRewards, target, STPGModelChecker.R_ZERO, minMax.isMin1(), minMax.isMin2(), minMax.getCoalition());
				break;
			}
			break;
		case CSG:
			switch (expr.getOperator()) {
			case ExpressionTemporal.P_F:
				res = ((CSGModelChecker) this).computeReachRewards((CSG<Double>) model, (CSGRewards<Double>) modelRewards, target, CSGModelChecker.R_INFINITY, minMax.isMin1(), minMax.isMin2(), minMax.getCoalition());
				break;
			case ExpressionTemporal.R_Fc:
				res = ((CSGModelChecker) this).computeReachRewards((CSG<Double>) model, (CSGRewards<Double>) modelRewards, target, CSGModelChecker.R_CUMULATIVE, minMax.isMin1(), minMax.isMin2(), minMax.getCoalition());
				break;
			}
			break;
		case IDTMC:
			res = ((IDTMCModelChecker) this).computeReachRewards((IDTMC<Double>) model, (MCRewards<Double>) modelRewards, target, minMax);
			break;
		case IMDP:
			res = ((IMDPModelChecker) this).computeReachRewards((IMDP<Double>) model, (MDPRewards<Double>) modelRewards, target, minMax);
			break;
		default:
			throw new PrismNotSupportedException("Explicit engine does not yet handle the " + expr.getOperatorSymbol() + " reward operator for " + model.getModelType()
					+ "s");
		}
		result.setStrategy(res.strat);
		return StateValues.createFromArrayResult(res, model);
	}

	/**
	 * Compute rewards for a co-safe LTL reward operator.
	 */
	protected StateValues checkRewardCoSafeLTL(Model<?> model, Rewards<?> modelRewards, Expression expr, MinMax minMax, BitSet statesOfInterest) throws PrismException
	{
		// To be overridden by subclasses
		throw new PrismException("Computation not implemented yet");
	}

	/**
	 * Model check an S operator expression and return the values for all states.
	 */
	protected StateValues checkExpressionSteadyState(Model<?> model, ExpressionSS expr) throws PrismException
	{
		// Get info from S operator
		OpRelOpBound opInfo = expr.getRelopBoundInfo(constantValues);
		MinMax minMax = opInfo.getMinMax(model.getModelType(), true, null);

		// Compute probabilities
		StateValues probs = checkSteadyStateFormula(model, expr.getExpression(), minMax);

		// Print out probabilities
		if (getVerbosity() > 5) {
			mainLog.print("\nProbabilities (non-zero only) for all states:\n");
			probs.print(mainLog);
		}

		// For =? properties, just return values; otherwise compare against bound
		if (!opInfo.isNumeric()) {
			probs.applyPredicate(v -> opInfo.apply((double) v, probs.getAccuracy()));
		}
		return probs;
	}

	/**
	 * Compute steady-state probabilities for an S operator.
	 */
	@SuppressWarnings("unchecked")
	protected StateValues checkSteadyStateFormula(Model<?> model, Expression expr, MinMax minMax) throws PrismException
	{
		// Model check operand for all states
		BitSet b = checkExpression(model, expr, null).getBitSet();

		// Compute/return the probabilities
		switch (model.getModelType()) {
		case DTMC:
			return ((DTMCModelChecker) this).computeSteadyStateFormula((DTMC<Double>) model, b);
		case CTMC:
			return ((CTMCModelChecker) this).computeSteadyStateFormula((CTMC<Double>) model, b);
		default:
			throw new PrismNotSupportedException("Explicit engine does not yet handle the S operator for " + model.getModelType() + "s");
		}
	}

	protected StateValues checkExpressionMultiObjective(Model<?> model, ExpressionStrategy expr, boolean forAll, Coalition coalition) throws PrismException
	{
		// Copy expression because we will modify it
		expr = (ExpressionStrategy) expr.deepCopy();

		// Strip any outer parentheses in operand
		Expression exprSub = expr.getOperand(0);
		while (Expression.isParenth(exprSub)) {
			exprSub = ((ExpressionUnaryOp) exprSub).getOperand();
		}

		// Boolean
		if (exprSub.getType() instanceof TypeBool) {
			// We will solve an existential query, so negate if universal
			if (forAll) {
				exprSub = Expression.Not(exprSub);
			}
			// Convert to CNF
			List<List<Expression>> cnf = BooleanUtils.convertToCNFLists(exprSub);
			// Check all "propositions" of CNF are valid
			for (List<Expression> disjunction : cnf) {
				for (Expression prop : disjunction) {
					if (Expression.isNot(prop)) {
						prop = ((ExpressionUnaryOp) prop).getOperand();
					}
					if (!(prop instanceof ExpressionQuant)) {
						throw new PrismException("Expression " + prop.getClass() + " is not allowed in a multi-objective query");
					}
				}
			}
			// Push negation inside objectives
			for (List<Expression> disjunction : cnf) {
				for (int j = 0; j < disjunction.size(); j++) {
					Expression prop = disjunction.get(j);
					if (Expression.isNot(prop)) {
						ExpressionQuant exprQuant = (ExpressionQuant) ((ExpressionUnaryOp) prop).getOperand();
						exprQuant.setRelOp(exprQuant.getRelOp().negate());
						disjunction.set(j, exprQuant);
					}
				}
			}
			// Print reduced query
			mainLog.println("\nReducing multi-objective query to CNF: " + BooleanUtils.convertCNFListsToExpression(cnf));

			// TODO: handle negation for 'forall'
			return checkExpressionMultiObjective(model, cnf, coalition);
		}
		else if (exprSub.getType() instanceof TypeDouble) {
			throw new PrismException("Multi-objective model checking not supported for: " + exprSub);
		} else {
			throw new PrismException("Multi-objective model checking not supported for: " + exprSub);
		}
	}

	public StateValues checkExpressionMultiObjective(Model<?> model, List<List<Expression>> cnf, Coalition coalition) throws PrismException
	{
		throw new PrismException("Not implemented");
	}

	/**
	 * Model check a multi-objective expression and return the values for all states.
	 */
	protected StateValues checkExpressionMulti(Model<?> model, Expression expr, boolean forAll, Coalition coalition) throws PrismException
	{
		// Assume "there exists" for now
		if (forAll)
			throw new PrismException("Nor support for forall in multi-objective queries yet");

		// For now, assume this is a function
		if (!(expr instanceof ExpressionFunc))
			throw new PrismException("Unsupported format for multi-objective query");

		switch (((ExpressionFunc) expr).getNameCode()) {
		case ExpressionFunc.MULTI:
			throw new PrismException("Properties with \"multi\" no longer supported.");
		default:
			throw new PrismException("Unsupported format for multi-objective query");
		}
	}

	/**
	 * Model check a function.
	 */
	protected StateValues checkExpressionFunc(Model<?> model, ExpressionFunc expr, BitSet statesOfInterest) throws PrismException
	{
		switch (expr.getNameCode()) {
		case ExpressionFunc.MULTI:
			throw new PrismNotSupportedException("Properties with \"multi\" not supported.");
		default:
			return super.checkExpressionFunc(model, expr, statesOfInterest);
		}
	}

	/**
	 * Finds states of the model which only have self-loops
	 *
	 * @param model
	 *            model
	 * @return the bitset indicating which states are terminal
	 */
	protected BitSet findTerminalStates(Model<?> model)
	{
		int n = model.getStatesList().size();
		BitSet ret = new BitSet(n), bs = new BitSet(n);
		for (int i = 0; i < n; i++) {
			bs.clear();
			bs.set(i);
			if (model.allSuccessorsInSet(i, bs))
				ret.set(i);
		}
		return ret;
	}

	// Utility methods for probability distributions

	/**
	 * Generate a probability distribution, stored as a StateValues object, from a file.
	 * If {@code distFile} is null, so is the return value.
	 */
	public StateValues readDistributionFromFile(File distFile, Model<?> model) throws PrismException
	{
		StateValues dist = null;

		if (distFile != null) {
			mainLog.println("\nImporting probability distribution from file \"" + distFile + "\"...");
			dist = StateValues.createFromFile(TypeDouble.getInstance(), distFile, model);
		}

		return dist;
	}

	/**
	 * Build a probability distribution, stored as a StateValues object,
	 * from the initial states info of the current model: either probability 1 for
	 * the (single) initial state or equiprobable over multiple initial states.
	 */
	public StateValues buildInitialDistribution(Model<?> model) throws PrismException
	{
		int numInitStates = model.getNumInitialStates();
		if (numInitStates == 1) {
			int sInit = model.getFirstInitialState();
			return StateValues.create(TypeDouble.getInstance(), s -> s == sInit ? 1.0 : 0.0, model);
		} else {
			double pInit = 1.0 / numInitStates;
			return StateValues.create(TypeDouble.getInstance(), s -> model.isInitialState(s) ? pInit : 0.0, model);
		}
	}

	/**
	 * Wrap a reward structure, replacing zeros with another small value {@code epsilon}.
	 */
	protected RewardsExplicit<Double> replaceZeroRewards(Rewards<Double> rewards, double epsilon)
	{
		return new RewardsExplicit<>() {

			@Override
			public Evaluator<Double> getEvaluator()
			{
				return rewards.getEvaluator();
			}

			@Override
			public boolean hasStateRewards()
			{
				return rewards.hasStateRewards();
			}

			@Override
			public boolean hasTransitionRewards()
			{
				return rewards.hasTransitionRewards();
			}

			@Override
			public Double getStateReward(int s)
			{
				double r = rewards.getStateReward(s);
				return r == 0 ? epsilon : r;
			}

			@Override
			public Double getTransitionReward(int s, int i)
			{
				double r = rewards.getTransitionReward(s, i);
				return r == 0 ? epsilon : r;
			}

			@Override
			public RewardsExplicit<Double> liftFromModel(Product<?> product)
			{
				return (RewardsExplicit<Double>) rewards.liftFromModel(product);
			}
		};
	}
}
