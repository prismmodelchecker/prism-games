package strat;

/**
 * Exception to be thrown when the strategy is not used in a proper way, i.e.,
 * next move function is undefined for a given pair or memory and state.
 * 
 * @author aistis
 * 
 */
public class InvalidStrategyStateException extends Exception
{

	public InvalidStrategyStateException(String string)
	{
		super(string);
	}

	/**
	 * 
	 */
	private static final long serialVersionUID = 6297385672024098036L;

}
