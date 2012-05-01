package strat;

import cern.colt.Arrays;
import prism.PrismFileLog;
import prism.PrismLog;
import explicit.Distribution;
import explicit.MDPSimple;
import explicit.MDPSparse;
import explicit.Model;

/**
 * Implementation of the memoryless deterministic strategy
 * 
 * @author aistis
 * 
 */
public class MemorylessDeterministicStrategy implements Strategy {

	private Distribution[] choices;

	public MemorylessDeterministicStrategy(int[] choices) {
		this.choices = new Distribution[choices.length];
		Distribution dist;
		for (int i = 0; i < choices.length; i++) {
			dist = new Distribution();
			dist.add(choices[i], 1);
			this.choices[i] = dist;
		}
	}

	@Override
	public void init(int state) throws InvalidStrategyStateException {
		// do nothing
	}

	@Override
	public void updateMemory(int action, int state)
			throws InvalidStrategyStateException {
		// do nothing
	}

	@Override
	public Distribution getNextMove(int state)
			throws InvalidStrategyStateException {

		if (choices == null || state >= choices.length || state < 0)
			throw new InvalidStrategyStateException(
					"Strategy not defined for state " + state + ".");

		return choices[state];
	}

	@Override
	public void reset() {
		// do nothing
	}

	@Override
	public Model buildProduct(Model model) {
		
		// checking for supported model types
		if (model instanceof MDPSimple) {
			return this.buildProductMDPSimple((MDPSimple) model);
		}
		if (model instanceof MDPSparse) {
			return this.buildProductMDPSparse((MDPSparse) model);
		}

		throw new UnsupportedOperationException(
				"The product building is not supported for this class of models");

	}

	@Override
	public void exportToFile(String file) {
		// Print adversary
		PrismLog out = new PrismFileLog(file);
		out.print("Adv:");
		for (int i = 0; i < choices.length; i++) {
			out.print(" " + i + ":" + choices[i]);
		}
		out.println();
	}

	@Override
	public String toString() {
		return Arrays.toString(choices);
	}

	/**
	 * Builds product of MDPSimple and a given strategy
	 * 
	 * @param model
	 *            model
	 */
	public MDPSimple buildProductMDPSimple(MDPSimple model) {
		MDPSimple mdp = new MDPSimple(model);
		int n = mdp.getNumStates();
		int c;
		Distribution distr;

		for (int s = 0; s < n; s++) {

			// getting the choice of a strategy for this state
			try {
				c = this.getNextMove(s).keySet().iterator().next();

				// if for adversary choice is undefined, taking the first one as
				// default
				if (c < 0)
					c = 0;
			} catch (InvalidStrategyStateException e) {
				// strategy undefined for this state -- keep calm and carry
				// on
				continue;
			}

			// replacing the choices with the one prescribed by the strategy
			distr = mdp.getChoice(s, c);
			mdp.clearState(s);
			mdp.addChoice(s, distr);
		}
		return mdp;
	}

	/**
	 * Builds product of MDPSparse and a given strategy
	 * 
	 * @param model
	 *            model
	 */
	public MDPSparse buildProductMDPSparse(MDPSparse model) {
		return new MDPSparse(buildProductMDPSimple(new MDPSimple(model)));
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
