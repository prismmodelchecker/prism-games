package parser.visitor;

import java.util.List;

import parser.ast.ExpressionVar;
import prism.PrismLangException;

public class GetPrimedVars extends ASTTraverse {

	private List<String> v;
	
	public GetPrimedVars(List<String> v)
	{
		this.v = v;
	}
	
	public void visitPost(ExpressionVar e) throws PrismLangException
	{
		if (!v.contains(e.getName()) && e.getPrime()) {
			v.add(e.getName());
		}
	}
}
