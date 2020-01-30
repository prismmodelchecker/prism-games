package parser.ast;

import java.util.ArrayList;
import java.util.List;

import param.BigRational;
import parser.EvaluateContext;
import parser.Values;
import parser.visitor.ASTVisitor;
import prism.OpRelOpBound;
import prism.PrismLangException;

public class ExpressionMultiNash extends Expression
{

	/** The attached relational operator (e.g. "&lt;" in "P&lt;0.1"). */
	protected RelOp relOp = null;
	/** The attached (probability/reward) bound, as an expression (e.g. "p" in "P&lt;p"). Null if absent (e.g. "P=?"). */
	protected Expression bound = null;
	// Operands
	protected ArrayList<ExpressionQuant> operands;

	// Constructors
	
	public ExpressionMultiNash()
	{
		operands = new ArrayList<>();
	}

	// Set methods
	
	/**
	 * Set the attached relational operator.
	 * Uses the enum {@link RelOp}. For example: {@code setRelOp(RelOp.GT);}
	 */
	public void setRelOp(RelOp relOp)
	{
		this.relOp = relOp;
	}

	/**
	 * Set the attached relational operator.
	 * The operator is passed as a string, e.g. "&lt;" or "&gt;=".
	 */
	public void setRelOp(String relOpString)
	{
		relOp = RelOp.parseSymbol(relOpString);
	}

	/**
	 * Set the attached bound, as an expression. Should be null if absent (e.g. "=?").
	 */
	public void setBound(Expression bound)
	{
		this.bound = bound;
	}

	public void addOperand(ExpressionQuant e)
	{
		operands.add(e);
	}

	public void setOperand(int i, ExpressionQuant e)
	{
		operands.set(i, e);
	}

	// Get methods
	
	/**
	 * Get the attached relational operator, as a {@link RelOp}.
	 */
	public RelOp getRelOp()
	{
		return relOp;
	}

	/**
	 * Get the attached bound, as an expression. Should be null if absent (e.g. "=?").
	 */
	public Expression getBound()
	{
		return bound;
	}

	public int getNumOperands()
	{
		return operands.size();
	}

	public ExpressionQuant getOperand(int i)
	{
		return operands.get(i);
	}

	public List<ExpressionQuant> getOperands()
	{
		return operands;
	}
	
	public OpRelOpBound getRelopBoundInfo(Values constantValues) throws PrismLangException
	{
		if (getBound() != null) {
			double boundVal = getBound().evaluateDouble(constantValues);
			return new OpRelOpBound("", getRelOp(), boundVal);
		} else {
			return new OpRelOpBound("", getRelOp(), null);
		}
	}

	// Methods required for Expression:

	@Override
	public boolean isConstant()
	{
		return false;
	}

	@Override
	public boolean isProposition()
	{
		return false;
	}

	@Override
	public Object evaluate(EvaluateContext ec) throws PrismLangException
	{
		throw new PrismLangException("Cannot evaluate a Nash operator without a model");
	}

	@Override
	public BigRational evaluateExact(EvaluateContext ec) throws PrismLangException
	{
		throw new PrismLangException("Cannot evaluate a Nash operator without a model");
	}

	@Override
	public boolean returnsSingleValue()
	{
		return false;
	}

	// Methods required for ASTElement:

	@Override
	public Object accept(ASTVisitor v) throws PrismLangException
	{
		return v.visit(this);
	}

	@Override
	public Expression deepCopy()
	{
		ExpressionMultiNash expr = new ExpressionMultiNash();
		expr.setRelOp(getRelOp());
		expr.setBound(getBound() == null ? null : getBound().deepCopy());
		int n = getNumOperands();
		for (int i = 0; i < n; i++) {
			expr.addOperand((ExpressionQuant) getOperand(i).deepCopy());
		}
		expr.setType(type);
		expr.setPosition(this);
		return expr;
	}
	
	// Standard methods
	
	@Override
	public String toString()
	{
		String s = "";
		s += getRelOp();
		s += (getBound() == null) ? "?" : getBound().toString();
		s += " (";
		int n = operands.size();
		boolean first = true;
		for (int i = 0; i < n; i++) {
			if (!first)
				s += " + ";
			else
				first = false;
			s = s + getOperand(i);
		}
		s += ")";
		return s;
	}

	@Override
	public int hashCode()
	{
		final int prime = 31;
		int result = 1;
		result = prime * result + ((bound == null) ? 0 : bound.hashCode());
		result = prime * result + ((operands == null) ? 0 : operands.hashCode());
		result = prime * result + ((relOp == null) ? 0 : relOp.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj)
	{
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		ExpressionMultiNash other = (ExpressionMultiNash) obj;
		if (bound == null) {
			if (other.bound != null)
				return false;
		} else if (!bound.equals(other.bound))
			return false;
		if (operands == null) {
			if (other.operands != null)
				return false;
		} else if (!operands.equals(other.operands))
			return false;
		if (relOp != other.relOp)
			return false;
		return true;
	}
}
