package strat;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Scanner;

import explicit.Distribution;
import explicit.MDPSimple;
import explicit.MDPSparse;
import explicit.Model;
import explicit.NondetModel;
import explicit.SMG;
import explicit.STPGExplicit;
import parser.State;
import prism.Prism.StrategyExportType;
import prism.PrismException;
import prism.PrismLog;

public class StepBoundedDeterministicStrategy implements Strategy
{
	// Model associated with the strategy
	private NondetModel model;
	
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
		this.model = model;
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
	public StepBoundedDeterministicStrategy(Scanner scan)
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
	}

	@Override
	public void init(int state) throws InvalidStrategyStateException
	{
		memory = bound;
	}

	@Override
	public void updateMemory(int action, int state) throws InvalidStrategyStateException
	{
		if (memory > 0)
			memory--;
	}

	@Override
	public Distribution getNextMove(int state) throws InvalidStrategyStateException
	{

		if (state > choices.length)
			throw new InvalidStrategyStateException("The strategy undefined for state " + state + ".");

		// determining the action
		int[] actions = choices[state];
		int c = 0;
		for (int i = 0; i < actions.length; i += 2)
			if (actions[i] >= memory)
				c = actions[i + 1];
			else
				break;

		Distribution dist = new Distribution();
		dist.add(c, 1);

		return dist;
	}

	@Override
	public void reset()
	{
		memory = bound;
	}

	@Override
	public int getMemorySize()
	{
		return bound;
	}

	@Override
	public Object getCurrentMemoryElement()
	{
		return memory;
	}

	@Override
	public void setMemory(Object memory) throws InvalidStrategyStateException
	{
		if (memory instanceof Integer) {
			this.memory = (Integer) memory;
		} else {
			throw new InvalidStrategyStateException("Memory has to integer for this strategy.");
		}
	}

	@Override
	public String getDescription()
	{
		String desc = "";
		desc += "Finite memory deterministic strategy\n";
		desc += "Size of memory: " + bound + "\n";
		desc += "Size of next move function: " + chSize + " \n";
		//desc += "Memory state: " + memory;
		return desc;
	}

	@Override
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

	}

	public static void main(String[] args) throws InvalidStrategyStateException
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
	}

	/**
	 * 
	 * @return
	 */
	@Override
	public String getInfo()
	{
		return info;
	}

	/**
	 * 
	 * @return
	 */
	@Override
	public String getType()
	{
		return "Finite memory strategy";
	}

	/**
	 * 
	 * @param info
	 */
	@Override
	public void setInfo(String info)
	{
		this.info = info;
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
		// checking for supported model types
		if (model.getClass().equals(MDPSimple.class)) {
			return this.buildProductMDPSimple((MDPSimple) model);
		}
		if (model.getClass().equals(MDPSparse.class)) {
			return this.buildProductMDPSparse((MDPSparse) model);
		}
		if (model.getClass().equals(STPGExplicit.class)) {
			return this.buildProductSTPGExplicit((STPGExplicit) model);
		}
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
		for (int j = bound; j >= 1; j--) {
			// setting memory
			this.memory = j;
			for (int i = 0; i < oldStates.size(); i++) {
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
							newDistr.add(oldStates.size() * (bound - j + j == 1 ? 0 : 1) + succ, distr.get(succ));

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
							newDistr.add(oldStates.size() * (bound - j + j == 1 ? 0 : 1) + succ, distr.get(succ));

						// adding the choice
						smg.addChoice(oldStates.size() * (bound - j) + i, newDistr);
					}
				}
			}
		}

		// setting initial state for the game
		smg.addInitialState(0);

		return smg;
	}

	private Model buildProductSTPGExplicit(STPGExplicit model)
	{
		// construct a new STPG of size ModelSize * MemorySize
		STPGExplicit stpg = new STPGExplicit(model.getStatesList().size() * bound);
		int n = stpg.getNumStates();

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
				stpg.setPlayer(newStates.size() - 1, model.getPlayer(i));
			}

		// setting the states list to STPG
		stpg.setStatesList(newStates);

		// adding choices for the product STPG

		// adding transitions to the state with the next memory element
		Distribution distr, newDistr;
		for (int j = bound; j >= 1; j--) {
			// setting memory
			this.memory = j;
			for (int i = 0; i < oldStates.size(); i++) {
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
							newDistr.add(oldStates.size() * (bound - j + j == 1 ? 0 : 1) + succ, distr.get(succ));

						// adding the choice
						stpg.addChoice(oldStates.size() * (bound - j) + i, newDistr);

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
							newDistr.add(oldStates.size() * (bound - j + j == 1 ? 0 : 1) + succ, distr.get(succ));

						// adding the choice
						stpg.addChoice(oldStates.size() * (bound - j) + i, newDistr);
					}
				}
			}
		}

		// setting initial state for the game
		stpg.addInitialState(0);

		return stpg;
	}

	/**
	 * 
	 * @param model
	 * @return
	 */
	private Model buildProductMDPSparse(MDPSparse model)
	{
		return new MDPSparse(buildProductMDPSimple(new MDPSimple(model)));
	}

	/**
	 * 
	 * @param model
	 * @return
	 */
	private MDPSimple buildProductMDPSimple(MDPSimple model)
	{
		// construct a new MDP of size ModelSize * MemorySize
		MDPSimple mdp = new MDPSimple(model.getStatesList().size() * bound);
		int n = mdp.getNumStates();

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
			for (int i = 0; i < oldStates.size(); i++)
				newStates.add(new State(oldStates.get(i), mem[j]));

		// setting the states list to MDP
		mdp.setStatesList(newStates);

		// adding choices for the product MDP

		// adding transitions to the state with the next memory element
		Distribution distr, newDistr;
		for (int j = bound; j >= 1; j--) {
			// setting memory
			this.memory = j;
			for (int i = 0; i < oldStates.size(); i++) {
				// retrieving choice chosen by the optimal strategy
				try {
					distr = model.getChoice(i, this.getNextMove(i).keySet().iterator().next());

					// create a new distribution for the product
					newDistr = new Distribution();
					for (Integer succ : distr.keySet())
						// adding transition to the state with the memory
						// element one larger smaller than the current one (j)
						// except for the case where j==1, when we add
						// transition to the same
						newDistr.add(oldStates.size() * (bound - j + j == 1 ? 0 : 1) + succ, distr.get(succ));

					// adding the choice
					mdp.addChoice(oldStates.size() * (bound - j) + i, newDistr);

				} catch (InvalidStrategyStateException error) {
					// TODO Auto-generated catch block
					error.printStackTrace();
				}
			}
		}

		// setting initial state for the MDP
		mdp.addInitialState(0);

		return mdp;
	}

	@Override
	public int getInitialStateOfTheProduct(int s)
	{
		return bound;
	}

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

	public int getChoiceIndex(int state, int mode)
	{
		int[] actions = choices[state];
		int c = -1;
		for (int i = 0; i < actions.length; i += 2)
			if (actions[i] >= memory)
				c = actions[i + 1];
			else
				break;
		
		return c;
	}
	
	public Object getChoiceAction(int state, int mode)
	{
		int c = getChoiceIndex(state, mode);
		return c >= 0 ? model.getAction(state, c) : null;
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
	public void exportInducedModel(PrismLog out)
	{
		// TODO Auto-generated method stub
		
	}

	@Override
	public void exportDotFile(PrismLog out)
	{
		// TODO Auto-generated method stub
		
	};

	@Override
	public void restrictStrategyToReachableStates() throws PrismException
	{
		// TODO Auto-generated method stub
		throw new PrismException("Reach option is not supported for this strategy type");
	}

	@Override
	public void exportStratToFile(File file, StrategyExportType exportType)
	{
		// TODO Auto-generated method stub
		
	}

	@Override
	public HashMap<String, Double> getNextAction(int state) throws InvalidStrategyStateException
	{
		// TODO Auto-generated method stub
		return null;
	}
}
