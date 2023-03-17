package parser.ast;

import param.BigRational;
import parser.EvaluateContext;
import parser.Values;
import parser.visitor.ASTVisitor;
import parser.visitor.DeepCopy;
import prism.OpRelOpBound;
import prism.PrismException;
import prism.PrismLangException;

public class ExpressionMultiNashProb extends ExpressionQuant
{

	protected Expression expr;

	public ExpressionMultiNashProb()
	{

	}

	public void setExpression(Expression e)
	{
		this.expr = e;
	}

	public Expression getExpression()
	{
		return this.expr;
	}

	@Override
	public OpRelOpBound getRelopBoundInfo(Values constantValues) throws PrismException
	{
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean isConstant()
	{
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean isProposition()
	{
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public Object evaluate(EvaluateContext ec) throws PrismLangException
	{
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public BigRational evaluateExact(EvaluateContext ec) throws PrismLangException
	{
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean returnsSingleValue()
	{
		// TODO Auto-generated method stub
		return false;
	}

	// Methods required for ASTElement:

	@Override
	public Object accept(ASTVisitor v) throws PrismLangException
	{
		return v.visit(this);
	}

	@Override
	public ExpressionMultiNashProb deepCopy(DeepCopy copier) throws PrismLangException
	{
		super.deepCopy(copier);
		expr = copier.copy(expr);
		return this;
	}

	@Override
	public ExpressionMultiNashProb clone()
	{
		return (ExpressionMultiNashProb) super.clone();
	}

	// Standard methods

	@Override
	public String toString()
	{
		// TODO Auto-generated method stub
		String s = "P[" + expr.toString() + "]";
		return s;
	}

}
