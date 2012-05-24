package strat;

/**
 * Class contains methods to work with strategies
 * 
 * @author aistis
 *
 */
public class Strategies
{
	private Strategies()
	{
		throw new AssertionError("This class should not be initialised.");
	}

	/**
	 * Loads the strategy from the given file
	 * @param filename name/path of the file
	 * @return the generated strategy
	 * @throws IllegalArgumentException if the file format is not recognised 
	 */
	public static Strategy loadStrategyFromFile(String filename) throws IllegalArgumentException
	{
		return null;
	}
}
