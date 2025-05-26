package strat;

import java.util.ArrayList;
import java.util.List;

import explicit.NondetModel;
import explicit.rewards.STPGRewards;

/**
 * 
 * @author aissim
 * @version
 */
public class BoundedRewardDeterministicStrategy<Value> extends StepBoundedDeterministicStrategy<Value>
{

	// rewards
	private List<Value> rewards;

	// number of states in the game for which the strategy is defined
	private int nStates;

	public BoundedRewardDeterministicStrategy(NondetModel<Value> model, int[][] choices, int bound, STPGRewards<Value> rewards)
	{
		super(model, choices, bound);
		nStates = choices.length;
		this.rewards = new ArrayList<>(nStates);
		for (int i = 0; i < nStates; i++) {
			this.rewards.add(rewards.getStateReward(i));
		}
	}

	/**
	 * Creates a BoundedRewardDeterministicStrategy.
	 *
	 * @param scan
	 */
	/*public BoundedRewardDeterministicStrategy(Scanner scan)
	{
		super(scan);
		nStates = choices.length;
		// parse state rewards
		rewards = new double[nStates];
		int i = 0;
		while (scan.hasNext())
			rewards[i++] = scan.nextDouble();
	}*/

	@Override
	public int getInitialMemory(int sInit)
	{
		return bound - (int) (double) rewards.get(sInit);
	}
	
	@Override
	public int getUpdatedMemory(int m, Object action, int sNext)
	{
		int memory = m - (int) (double) rewards.get(sNext);
		return memory < 0 ? 0 : memory;
	}
	
	/*@Override
	public void exportToFile(String file)
	{
		// Print adversary
		FileWriter out = null;
		try {
			out = new FileWriter(new File(file));

			out.write(Strategies.FORMAT_STRING_BOUNDED_REW_STRAT + "\n");
			out.write("// Strategy for F0 reward properties\n");
			out.write("// format: stateId, b1, c1, b2, c2,..., bn, cn\n");
			out.write("// (b1>b2>...>bn)\n");
			out.write("// where: ci  (1<=i<n )is the choice taken when the reward left to accumulate before the bound is reached is >=bi and <bi+1\n");
			out.write("// cn is the choice taken after bn or less remain to accummulate until bound is reached.\n");
			out.write("Strategy:\n");
			for (int i = 0; i < choices.length; i++) {
				out.write("" + i);
				for (int j = 0; j < choices[i].length; j++) {
					out.write(" " + choices[i][j]);
				}
				out.write("\n");
			}
			out.write("Rewards:\n");
			for (int i = 0; i < nStates; i++)
				out.write(" " + rewards[i]);
			out.flush();
		} catch (IOException error) {
			// TODO Auto-generated catch block
			error.printStackTrace();
		} finally {
			if (out != null)
				try {
					out.close();
				} catch (IOException error) {
					// do nothings
				}
		}
	}*/
}
