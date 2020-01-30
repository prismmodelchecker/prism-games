package parser.ast;

import param.BigRational;
import parser.EvaluateContext;
import parser.Values;
import parser.visitor.ASTVisitor;
import prism.OpRelOpBound;
import prism.PrismException;
import prism.PrismLangException;

public class ExpressionMultiNashProb extends ExpressionQuant {
	
	protected Expression expr;

	public ExpressionMultiNashProb() {
		
	}
	
	public void setExpression(Expression e) {
		this.expr = e;
	}
	
	public Expression getExpression() {
		return this.expr;
	}

	@Override
	public OpRelOpBound getRelopBoundInfo(Values constantValues) throws PrismException {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean isConstant() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean isProposition() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public Object evaluate(EvaluateContext ec) throws PrismLangException {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public BigRational evaluateExact(EvaluateContext ec) throws PrismLangException {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean returnsSingleValue() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public Expression deepCopy() {
		ExpressionMultiNashProb expr = new ExpressionMultiNashProb();
		expr.setExpression(getExpression() == null ? null : getExpression().deepCopy());
		expr.setType(type);
		expr.setPosition(this);
		return expr;
	}

	@Override
	public Object accept(ASTVisitor v) throws PrismLangException {
		return v.visit(this);
	}
	
	@Override
	public String toString() {
		// TODO Auto-generated method stub
		String s = "P[" + expr.toString() + "]";
		return s;
	}
	
}
