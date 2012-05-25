package strat;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

/**
 * Class contains methods to work with strategies
 * 
 * @author aistis
 *
 */
public class Strategies
{
	public static final String FORMAT_STRING_MD_STRAT = "$MD.strat4355t4fre#";
	public static final String FORMAT_STRING_STEP_BOUNDED_STRAT = "$SB.strat12sdd3";
	public static final String FORMAT_STRING_BOUNDED_REW_STRAT = "$RB.strat453dfs";
	public static final String FORMAT_STRING_EXACT_VALUE_MD_STRAT = "$EVMD.strat679csxc";

	private Strategies()
	{
		throw new AssertionError("This class should not be initialised.");
	}

	/**
	 * Loads the strategy from the given file
	 * @param filename name/path of the file
	 * @return the generated strategy
	 * @throws IllegalArgumentException if the file format is not recognised or file is not found
	 */
	public static Strategy loadStrategyFromFile(String filename) throws IllegalArgumentException
	{
		try {
			Scanner scan = new Scanner(new File(filename));
			try {
				String type = scan.nextLine();
				if (type.equals(FORMAT_STRING_MD_STRAT)) {
					return new MemorylessDeterministicStrategy(scan);
				} else if (type.equals(FORMAT_STRING_STEP_BOUNDED_STRAT)) {
					return new StepBoundedDeterministicStrategy(scan);
				} else if (type.equals(FORMAT_STRING_BOUNDED_REW_STRAT)) {
					return new BoundedRewardDeterministicStrategy(scan);
				} else if (type.equals(FORMAT_STRING_EXACT_VALUE_MD_STRAT)) {
					return new ExactValueStrategy(scan);
				}
				throw new IllegalArgumentException("Format not supported");
			} finally {
				scan.close();
			}
		} catch (FileNotFoundException error) {
			throw new IllegalArgumentException("File not found.");
		}
	}

	public static void main(String[] args)
	{
		String fn = "/users/aissim/prism-games/prism-games/prism/md.adv";
		String fn2 = "/users/aissim/prism-games/prism-games/prism/md2.adv";
		Strategy mdstrat = new MemorylessDeterministicStrategy(new int[] { 1, 2, 4, 6, 2 });
		mdstrat.exportToFile(fn);
		Strategy mdstrat2 = Strategies.loadStrategyFromFile(fn);
		mdstrat2.exportToFile(fn2);
		mdstrat = mdstrat2 = null;

		int[][] choices = { { 30, 1, 28, 2 }, { 25, 1, 24, 2 } };
		int bound = 25;

		String sbfn = "/users/aissim/prism-games/prism-games/prism/sb.adv";
		String sbfn2 = "/users/aissim/prism-games/prism-games/prism/sb2.adv";
		StepBoundedDeterministicStrategy sbstrat = new StepBoundedDeterministicStrategy(choices, bound);
		sbstrat.exportToFile(sbfn);
		Strategy sbstrat2 = Strategies.loadStrategyFromFile(sbfn);
		sbstrat2.exportToFile(sbfn2);
	}

}
