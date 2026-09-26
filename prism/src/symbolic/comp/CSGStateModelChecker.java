package symbolic.comp;

import jdd.JDD;
import jdd.JDDNode;
import jdd.JDDVars;
import parser.VarList;
import parser.Values;
import parser.ast.ExpressionVar;
import prism.Prism;
import prism.PrismException;
import prism.PrismLangException;
import symbolic.states.StateValues;
import symbolic.states.StateValuesMTBDD;

/**
 * A {@link StateModelChecker} specialised for translating CSG guard/update expressions to
 * ADDs, where (unlike every other PRISM model type) an expression may reference the
 * <i>next-state</i> ("primed") value of a variable, not just its current-state value -- e.g.
 * the RHS of an update such as {@code (s1'=c'?1:0)}. The primed-variable dependency
 * acyclicity check that this relies on being well-defined is not implemented by this class --
 * this class only provides the translation itself.
 *
 * <p>{@link StateModelChecker#checkExpressionVar} is the existing base case for an ordinary
 * (non-primed) variable reference: it decodes the variable's value from its <i>row</i>
 * (current-state) ADD variables via {@code SetVectorElement}. This subclass overrides just
 * that one method: when {@link ExpressionVar#getPrime()} is true, it performs the identical
 * construction over the variable's <i>column</i> (next-state) ADD variables instead: return
 * an ADD {@code f} where {@code f(s,s')} is the value of the primed variable in {@code s'}.
 * Every other expression construct (arithmetic, boolean operators, ITE, function calls, ...)
 * is untouched and simply falls through to the superclass's existing, tested recursion.
 */
public class CSGStateModelChecker extends StateModelChecker
{
	/** Column (next-state) ADD variables per model variable -- not held by the superclass. */
	private final JDDVars[] varDDColVars;

	public CSGStateModelChecker(Prism prism, VarList varList, JDDVars allDDRowVars, JDDVars[] varDDRowVars, JDDVars[] varDDColVars, Values constantValues)
			throws PrismException
	{
		super(prism, varList, allDDRowVars, varDDRowVars, constantValues);
		this.varDDColVars = varDDColVars;
	}

	@Override
	protected StateValues checkExpressionVar(ExpressionVar expr, JDDNode statesOfInterest) throws PrismException
	{
		if (!expr.getPrime()) {
			return super.checkExpressionVar(expr, statesOfInterest);
		}

		// Primed reference: decode the variable's value from its COLUMN (next-state) ADD
		// variables instead of its row variables -- otherwise identical to the superclass's
		// own checkExpressionVar. statesOfInterest is deliberately ignored, matching the
		// superclass's own documented behaviour ("more efficient not to restrict... here").
		JDD.Deref(statesOfInterest);

		String s = expr.getName();
		int v = varList.getIndex(s);
		if (v == -1) {
			throw new PrismLangException("Unknown variable \"" + expr.getName() + "\"", expr);
		}
		int l = varList.getLow(v);
		int h = varList.getHigh(v);

		JDDNode dd = JDD.Constant(0);
		for (int i = l; i <= h; i++) {
			dd = JDD.SetVectorElement(dd, varDDColVars[v], i - l, i);
		}

		return new StateValuesMTBDD(dd, model);
	}
}
