package parser.visitor;

import java.util.Vector;
import parser.ast.*;
import prism.PrismLangException;

public class GetPrimedVars extends ASTTraverse {

	private Vector<String> v;
	
	public GetPrimedVars(Vector<String> v)
	{
		this.v = v;
	}
	
	public void visitPost(ExpressionVar e) throws PrismLangException
	{
		if (!v.contains(e.getName()) && e.getPrime()) {
			v.addElement(e.getName());
		}
	}
	
}
