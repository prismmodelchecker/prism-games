package strat;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

import parser.State;
import prism.PrismException;
import explicit.Distribution;
import explicit.Model;
import explicit.NondetModel;
import explicit.SMG;
import explicit.rewards.STPGRewards;

/**
 * 
 * @author aissim
 * @version
 */
public class BoundedRewardDeterministicStrategy extends StepBoundedDeterministicStrategy
{

	// Model associated with the strategy
	private NondetModel model;
		
	// rewards
	private double[] rewards;

	// number of states in the game for which the strategy is defined
	private int nStates;

	public BoundedRewardDeterministicStrategy(NondetModel model, int[][] choices, int bound, STPGRewards rewards)
	{
		super(model, choices, bound);
		nStates = choices.length;
		this.rewards = new double[nStates];
		for (int i = 0; i < nStates; i++)
			this.rewards[i] = rewards.getStateReward(i);
	}

	/**
	 * Creates a BoundedRewardDeterministicStrategy.
	 *
	 * @param scan
	 */
	public BoundedRewardDeterministicStrategy(Scanner scan)
	{
		super(scan);
		nStates = choices.length;
		// parse state rewards
		rewards = new double[nStates];
		int i = 0;
		while (scan.hasNext())
			rewards[i++] = scan.nextDouble();
	}

	@Override
	public void init(int state) throws InvalidStrategyStateException
	{
		memory = bound - (int) rewards[state];
	}

	@Override
	public void updateMemory(int action, int state) throws InvalidStrategyStateException
	{
		memory -= rewards[state];
		if (memory < 0)
			memory = 0;
	}

	@Override
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
	}

	/**
	 * 
	 * @param model
	 * @return
	 * @throws PrismException
	 */
	@Override
	public Model buildProduct(Model model) throws PrismException
	{

		if (model.getClass().equals(SMG.class)) {
			return this.buildProductSMG((SMG) model);
		}

		throw new UnsupportedOperationException("The product building is not supported for this class of models");
	}

	private Model buildProductSMG(SMG model) throws PrismException
	{
		// construct a new SMG of size ModelSize * MemorySize
		SMG smg = new SMG(model.getStatesList().size() * bound);
		smg.copyPlayerInfo(model);
		smg.copyCoalitionInfo(model);
		int n = smg.getNumStates();

		List<State> oldStates = model.getStatesList();

		// creating helper states for constructing the product
		State[] mem = new State[bound];
		for (int i = bound; i >= 1; i--) {
			mem[bound - i] = new State(1);
			mem[bound - i].setValue(0, i);
		}

		// creating product state list
		List<State> newStates = new ArrayList<State>(n);
		int count = 0;
		for (int j = 0; j < bound; j++)
			for (int i = 0; i < oldStates.size(); i++) {
				newStates.add(new State(oldStates.get(i), mem[j]));
				smg.setPlayer(count++, model.getPlayer(i));
			}

		// setting the states list to SMG
		smg.setStatesList(newStates);

		// adding choices for the product SMG
		// adding transitions to the state with the next memory element
		Distribution distr, newDistr;
		int initial = -1;
		for (int j = bound; j >= 1; j--) {
			// setting memory
			this.memory = j;
			for (int i = 0; i < oldStates.size(); i++) {
				if (i == 0 && j == bound - rewards[i])
					initial = (bound - j) * oldStates.size() + i;
				// if the state belongs to player 1 retrieving choice chosen by
				// the optimal strategy
				if (model.getPlayer(i) == 1) {
					try {
						distr = model.getChoice(i, this.getNextMove(i).keySet().iterator().next());

						// create a new distribution for the product
						newDistr = new Distribution();
						for (Integer succ : distr.keySet())
							// adding transition to the state with the memory
							// element one larger smaller than the current one
							// (j)
							// except for the case where j==1, when we add
							// transition to the same
							newDistr.add(oldStates.size() * (bound - j + (j == 1 ? 0 : (int) rewards[succ])) + succ, distr.get(succ));

						// adding the choice
						smg.addChoice(oldStates.size() * (bound - j) + i, newDistr);

					} catch (InvalidStrategyStateException error) {
						// TODO Auto-generated catch block
						error.printStackTrace();
					}
				} else // otherwise copying all distributions
				{
					for (int k = 0; k < model.getNumChoices(i); k++) {
						distr = model.getChoice(i, k);

						// create a new distribution for the product
						newDistr = new Distribution();
						for (Integer succ : distr.keySet())
							// adding transition to the state with the memory
							// element one larger smaller than the current one
							// (j)
							// except for the case where j==1, when we add
							// transition to the same
							newDistr.add(oldStates.size() * (bound - j + j == 1 ? 0 : (int) rewards[succ]) + succ, distr.get(succ));

						// adding the choice
						smg.addChoice(oldStates.size() * (bound - j) + i, newDistr);
					}
				}
			}
		}

		// setting initial state for the game
		smg.addInitialState(initial);

		return smg;

	}

	@Override
	public int getInitialStateOfTheProduct(int s)
	{
		return bound - (int) rewards[s % nStates];
	}
}
