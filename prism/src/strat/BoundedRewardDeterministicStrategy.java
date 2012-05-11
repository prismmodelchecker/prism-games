package strat;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import parser.State;
import prism.PrismException;
import prism.PrismFileLog;
import prism.PrismLog;
import explicit.Distribution;
import explicit.MDPSimple;
import explicit.MDPSparse;
import explicit.Model;
import explicit.SMG;
import explicit.STPGExplicit;
import explicit.rewards.STPGRewards;

/**
 * 
 * @author aissim
 * @version
 */
public class BoundedRewardDeterministicStrategy extends
		StepBoundedDeterministicStrategy {

	private STPGRewards rewards;

	public BoundedRewardDeterministicStrategy(int[][] choices, int bound,
			STPGRewards rewards) {
		super(choices, bound);
		this.rewards = rewards;
	}

	@Override
	public void init(int state) throws InvalidStrategyStateException {
		memory = bound - (int) rewards.getStateReward(state);
	}

	@Override
	public void updateMemory(int action, int state)
			throws InvalidStrategyStateException {
		memory -= rewards.getStateReward(state);
		if (memory < 0)
			memory = 0;
	}

	@Override
	public void exportToFile(String file) {
		// Print adversary
		PrismLog out = new PrismFileLog(file);
		out.print("// Strategy for F0 reward properties\n");
		out.print("// format: stateId, b1, c1, b2, c2,..., bn, cn\n");
		out.print("// (b1>b2>...>bn)\n");
		out.print("// where: ci  (1<=i<n )is the choice taken when the reward left to accumulate before the bound is reached is >=bi and <bi+1\n");
		out.print("// cn is the choice taken after bn or less remain to accummulate until bound is reached.\n");
		out.print("Strategy:\n");
		for (int i = 0; i < choices.length; i++) {
			out.print(i);
			for (int j = 0; j < choices[i].length; j++) {
				out.print(", " + choices[i][j]);
			}
			out.println();
		}
		out.flush();
	}

	/**
	 * 
	 * @param model
	 * @return
	 * @throws PrismException
	 */
	@Override
	public Model buildProduct(Model model) throws PrismException {

		if (model.getClass().equals(SMG.class)) {
			return this.buildProductSMG((SMG) model);
		}

		throw new UnsupportedOperationException(
				"The product building is not supported for this class of models");
	}

	private Model buildProductSMG(SMG model) throws PrismException {
		// construct a new SMG of size ModelSize * MemorySize
		SMG smg = new SMG(model.getStatesList().size() * bound);
		smg.setPlayerMapping(new HashMap<String, Integer>(model
				.getPlayerMapping()));
		smg.setCoalitionInts(new HashSet<Integer>(model.getCoalition()));
		int n = smg.getNumStates();

		List<Integer> stateLabels = new ArrayList<Integer>(n);

		List<State> oldStates = model.getStatesList();

		// creating helper states for constructing the product
		State[] mem = new State[bound];
		for (int i = bound; i >= 1; i--) {
			mem[bound - i] = new State(1);
			mem[bound - i].setValue(0, i);
		}

		// creating product state list
		List<State> newStates = new ArrayList<State>(n);
		for (int j = 0; j < bound; j++)
			for (int i = 0; i < oldStates.size(); i++) {
				newStates.add(new State(oldStates.get(i), mem[j]));
				stateLabels.add(model.getPlayer(i));
			}

		// setting the states list to SMG
		smg.setStatesList(newStates);
		smg.setStateLabels(stateLabels);

		// adding choices for the product SMG
		// adding transitions to the state with the next memory element
		Distribution distr, newDistr;
		for (int j = bound; j >= 1; j--) {
			// setting memory
			this.memory = j;
			for (int i = 0; i < oldStates.size(); i++) {
				// if the state belongs to player 1 retrieving choice chosen by
				// the optimal strategy
				if (model.getPlayer(i) == SMG.PLAYER_1) {
					try {
						distr = model.getChoice(i, this.getNextMove(i).keySet()
								.iterator().next());

						// create a new distribution for the product
						newDistr = new Distribution();
						for (Integer succ : distr.keySet())
							// adding transition to the state with the memory
							// element one larger smaller than the current one
							// (j)
							// except for the case where j==1, when we add
							// transition to the same
							newDistr.add(
									oldStates.size()
											* (bound - j + j == 1 ? 0
													: (int) rewards
															.getStateReward(succ))
											+ succ, distr.get(succ));

						// adding the choice
						smg.addChoice(oldStates.size() * (bound - j) + i,
								newDistr);

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
							newDistr.add(
									oldStates.size()
											* (bound - j + j == 1 ? 0
													: (int) rewards
															.getStateReward(succ))
											+ succ, distr.get(succ));

						// adding the choice
						smg.addChoice(oldStates.size() * (bound - j) + i,
								newDistr);
					}
				}
			}
		}

		// setting initial state for the game
		smg.addInitialState(0);

		return smg;

	}

}
