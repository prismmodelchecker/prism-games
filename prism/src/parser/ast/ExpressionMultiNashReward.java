package parser.ast;

import param.BigRational;
import parser.EvaluateContext;
import parser.Values;
import parser.visitor.ASTVisitor;
import prism.OpRelOpBound;
import prism.PrismException;
import prism.PrismLangException;

public class ExpressionMultiNashReward extends ExpressionReward {

	protected Expression expr;
	
	public ExpressionMultiNashReward() {

	}
	
	public void setExpression(Expression e) {
		expr = e;
	}
	
	public Expression getExpression() {
		return expr;
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
		ExpressionMultiNashReward expr = new ExpressionMultiNashReward();
		expr.setExpression(getExpression() == null ? null : getExpression().deepCopy());
		expr.setRewardStructIndex(getRewardStructIndex());
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
		String s = "R";
		s += "{";
		if (rewardStructIndex != null) {
			if (rewardStructIndex instanceof Expression) s += rewardStructIndex;
			else if (rewardStructIndex instanceof String) s += "\"" + rewardStructIndex + "\"";
		}
		s += "}";
		s += "[" + expr.toString() + "]";
		return s;
	}

}
