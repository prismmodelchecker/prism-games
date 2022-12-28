package strat;

import explicit.NondetModel;
import prism.PrismLog;

public class StepBoundedDeterministicStrategy extends StrategyExplicit
{
	// memory: the number of steps currently made
	protected int memory;

	// bound: maximum number of steps to be made
	protected int bound;

	// choices of a strategy
	protected int[][] choices;
	protected int chSize;

	// information
	protected String info = "No information available";

	/**
	 * Initialises the strategy
	 * 
	 * @param model
	 *            Model associated with the strategy
	 * @param choices
	 *            strategy choice table: for every state contains an array of
	 *            integers structured as follows {bk, ck, bk-1, ck-1,..., b0 c0}
	 *            where ck represents the choice to be made by the strategy when
	 *            B-k steps have elapsed.
	 * @param bound
	 *            maximum number of steps to be made
	 */
	public StepBoundedDeterministicStrategy(NondetModel model, int[][] choices, int bound)
	{
		super(model);
		this.choices = choices;

		if (bound < 0)
			throw new IllegalArgumentException("The bound should be positive.");

		this.bound = bound;

		// computing the size of the choice function and validating the format
		chSize = 0;
		int prev;
		for (int i = 0; i < choices.length; i++) {
			prev = bound;
			for (int j = 0; j < choices[i].length; j++) {
				chSize++;

				// performing validation
				if (choices[i][j] < 0)
					throw new IllegalArgumentException(
							"The format of choices is invalid: array cannot contain negative numbers.");

				// adjusting the choices to be at most the bound
				if (j % 2 == 0 && choices[i][j] > bound) {
					choices[i][j] = bound;
					prev = bound;
				} else if (j == 0 && choices[i][j] < bound) {
					throw new IllegalArgumentException(
							"The format of choices is invalid: the first pivot has to be >= than the bound.");
				}

				// checking if ordering is correct
				if (j % 2 == 0)
					if (choices[i][j] > prev)
						throw new IllegalArgumentException(
								"The format of choices is invalid: pivots have to be in decreasing order.");
					else
						prev = choices[i][j];
			}
		}
	}

	/**
	 * Creates a StepBoundedDeterministicStrategy.
	 *
	 * @param scan
	 */
	/*public StepBoundedDeterministicStrategy(Scanner scan)
	{
		for (int i = 0; i < 6; i++)
			scan.nextLine();
		chSize = 0;
		Scanner local;
		List<int[]> choices = new LinkedList<int[]>();
		List<Integer> single;
		int[] singleA;
		String nextLine;
		while (scan.hasNext()) {
			nextLine = scan.nextLine();
			if (nextLine.startsWith("Rewards:"))
				break;
			local = new Scanner(nextLine);
			local.nextInt();
			bound = local.nextInt();
			single = new LinkedList<Integer>();
			single.add(bound);
			while (local.hasNext()) {
				chSize++;
				single.add(local.nextInt());
			}
			singleA = new int[single.size()];
			for (int i = 0; i < single.size(); i++)
				singleA[i] = single.get(i);
			choices.add(singleA);
		}
		this.choices = choices.toArray(new int[][] {});
	}*/

	@Override
	public Memory memory()
	{
		return Memory.FINITE;
	}
	
	@Override
	public Object getChoiceAction(int s, int m)
	{
		int c = getChoiceIndex(s, m);
		return c >= 0 ? model.getAction(s, c) : Strategy.UNDEFINED;
	}
	
	@Override
	public int getChoiceIndex(int s, int m)
	{
		int[] actions = choices[s];
		int c = 0;
		for (int i = 0; i < actions.length; i += 2)
			if (actions[i] >= m)
				c = actions[i + 1];
			else
				break;
		return c;
	}
	
	@Override
	public int getMemorySize()
	{
		return bound;
	}
	
	@Override
	public int getInitialMemory(int sInit)
	{
		return bound;
	}
	
	@Override
	public int getUpdatedMemory(int m, Object action, int sNext)
	{
		return m > 0 ? m - 1 : m;
	}
	
	//@Override
	public String getDescription()
	{
		String desc = "";
		desc += "Finite memory deterministic strategy\n";
		desc += "Size of memory: " + bound + "\n";
		desc += "Size of next move function: " + chSize + " \n";
		//desc += "Memory state: " + memory;
		return desc;
	}

	/*@Override
	public void exportToFile(String file)
	{
		// Print adversary
		FileWriter out = null;
		try {
			out = new FileWriter(new File(file));

			out.write(Strategies.FORMAT_STRING_STEP_BOUNDED_STRAT + "\n");
			out.write("// Strategy for step-bounded properties\n");
			out.write("// format: stateId, b1, c1, b2, c2,..., bn, cn\n");
			out.write("// (b1>b2>...>bn)\n");
			out
					.write("// where: ci  (1<=i<n )is the choice taken when the number of steps remaining before the bound is exceeded is >=bi and <bi+1\n");
			out.write("// cn is the choice taken after bn or less steps remain until bound is exceeded.\n");
			out.write("Strategy:\n");
			for (int i = 0; i < choices.length; i++) {
				out.write("" + i);
				for (int j = 0; j < choices[i].length; j++) {
					out.write(" " + choices[i][j]);
				}
				out.write("\n");
			}
			out.flush();
		} catch (IOException error) {
			// TODO Auto-generated catch block
			error.printStackTrace();
		} finally {
			if (out != null)
				try {
					out.close();
				} catch (IOException error) {
					// nothing we can do
				}

		}

	}*/

	/*public static void main(String[] args) throws InvalidStrategyStateException
	{
		int[][] choices = { { 30, 1, 28, 2 }, { 25, 1, 24, 2 } };
		int bound = 25;

		StepBoundedDeterministicStrategy strat = new StepBoundedDeterministicStrategy(null, choices, bound);
		strat.init(0);

		for (int i = 0; i < 25; i++) {
			System.out.println("i = " + i);
			System.out.println(strat.getNextMove(0) + ", " + strat.getNextMove(1));
			strat.updateMemory(0, 0);
		}
	}*/


	public void export(PrismLog out) {}

	@Override
	public void exportActions(PrismLog out)
	{
		int n = (int) model.getNumStates();
		for (int s = 0; s < n; s++) {
			for (int m = 0; m <= bound; m++) {
				Object action = getChoiceAction(s, m);
				if (action!=null)
					out.println(s + "," + m + ":" + action.toString() + "," + (m<bound ? m+1 : m));
			}
		}
	}

	@Override
	public void clear()
	{
		// TODO Auto-generated method stub
	}

	@Override
	public void exportIndices(PrismLog out)
	{
		int n = (int) model.getNumStates();
		for (int s = 0; s < n; s++) {
			for (int m = 0; m <= bound; m++) {
				int c = getChoiceIndex(s, m);
				if (c>=0)
					out.println(s + "," + m + ":" + c + "," + (m<bound ? m+1 : m));
			}
		}
	}

	@Override
	public void exportInducedModel(PrismLog out, int precision)
	{
		// TODO Auto-generated method stub
		
	}

	@Override
	public void exportDotFile(PrismLog out, int precision)
	{
		// TODO Auto-generated method stub
		
	}

	@Override
	public int getNumStates()
	{
		// TODO Auto-generated method stub
		return 0;
	};
}
