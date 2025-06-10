package explicit;

import prism.PrismComponent;
import prism.PrismException;

public class UCSGModelChecker extends ProbModelChecker
{
	/**
	 * Create a new UCSGModelChecker, inherit basic state from parent (unless null).
	 */
	public UCSGModelChecker(PrismComponent parent) throws PrismException
	{
		super(parent);
//		mcMDP = new MDPModelChecker(this);
//		mcMDP.inheritSettings(this);
	}
}
