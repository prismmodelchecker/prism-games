package prism;

import parser.ast.RelOp;
import explicit.MinMax;

/**
 * Class to represent info (operator, relational operator, bound, etc.) found in a P/R/S operator.
 */
public class OpRelOpBound
{
	protected String op;
	protected RelOp relOp;
	protected boolean numeric;
	protected double bound;

	public OpRelOpBound(String op, RelOp relOp, Double boundObject)
	{
		this.op = op;
		this.relOp = relOp;
		numeric = (boundObject == null);
		if (boundObject != null)
			bound = boundObject.doubleValue();
	}

	public boolean isProbabilistic()
	{
		return "P".equals(op);
	}

	public boolean isReward()
	{
		return "R".equals(op);
	}

	public RelOp getRelOp()
	{
		return relOp;
	}

	public boolean isNumeric()
	{
		return numeric;
	}

	public double getBound()
	{
		return bound;
	}

	public boolean isQualitative()
	{
		return !isNumeric() && op.equals("P") && (bound == 0 || bound == 1);
	}

	public boolean isTriviallyTrue()
	{
		if (!isNumeric() && op.equals("P")) {
			// >=0
			if (bound == 0 && relOp == RelOp.GEQ)
				return true;
			// <=1
			if (bound == 1 && relOp == RelOp.LEQ)
				return true;
		}
		return false;
	}

	public boolean isTriviallyFalse()
	{
		if (!isNumeric() && op.equals("P")) {
			// <0
			if (bound == 0 && relOp == RelOp.LT)
				return true;
			// >1
			if (bound == 1 && relOp == RelOp.GT)
				return true;
		}
		return false;
	}

	public boolean isExact()
	{
		return relOp == RelOp.EQ && !numeric; 
	}
	
	public MinMax getMinMax(ModelType modelType) throws PrismLangException
	{
		return getMinMax(modelType, true);
	}

	public MinMax getMinMax(ModelType modelType, boolean forAll) throws PrismLangException
	{
		MinMax minMax = MinMax.blank();
		if (modelType.nondeterministic()) {
			if (!modelType.multiplePlayers()) {
				if (!(modelType == ModelType.MDP || modelType == ModelType.CTMDP)) {
					throw new PrismLangException("Don't know how to model check " + getTypeOfOperator() + " properties for " + modelType + "s");
				}
				if (isNumeric()) {
					if (relOp == RelOp.EQ) {
						throw new PrismLangException("Can't use \"" + op + "=?\" for nondeterministic models; use e.g. \"" + op + "min=?\" or \"" + op + "max=?\"");
					}
					minMax = relOp.isMin() ? MinMax.min() : MinMax.max();
				} else {
					if (forAll) {
						minMax = (relOp.isLowerBound()) ? MinMax.min() : MinMax.max();
					} else {
						minMax = (relOp.isLowerBound()) ? MinMax.max() : MinMax.min();
					}
				}
			} else {
				if (modelType == ModelType.SMG) {
					if (relOp == RelOp.EQ && isNumeric()) {
						throw new PrismLangException("Can't use \"" + op + "=?\" for SMGs; use e.g. \"" + op + "min=?\" or \"" + op + "max=?\"");
					}
					if (relOp.isMin() || (forAll && relOp.isLowerBound()) || (!forAll && relOp.isUpperBound())) {
						minMax = MinMax.minMin(true, false);
					} else {
						minMax = MinMax.minMin(false, true);
					}
				} else if (modelType == ModelType.STPG) {
					if (relOp == RelOp.EQ && isNumeric()) {
						throw new PrismLangException("Can't use \"" + op + "=?\" for STPGs; use e.g. \"" + op + "minmax=?\"");
					}
					if (relOp == RelOp.MINMIN) {
						minMax = MinMax.minMin(true, true);
					} else if (relOp == RelOp.MINMAX) {
						minMax = MinMax.minMin(true, false);
					} else if (relOp == RelOp.MAXMIN) {
						minMax = MinMax.minMin(false, true);
					} else if (relOp == RelOp.MAXMAX) {
						minMax = MinMax.minMin(false, false);
					} else {
						throw new PrismLangException("Use e.g. \"Rminmax=?\" for stochastic games");
					}
				}
			}
		}

		return minMax;
	}

	public String getTypeOfOperator()
	{
		String s = "";
		s += op + relOp;
		s += isNumeric() ? "?" : "p"; // TODO: always "p"?
		return s;
	}

	public String relOpBoundString()
	{
		return relOp.toString() + bound;
	}

	@Override
	public String toString()
	{
		return op + relOp.toString() + bound;
	}
}
