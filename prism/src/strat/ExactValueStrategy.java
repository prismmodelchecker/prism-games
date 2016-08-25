package strat;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Scanner;

import parser.State;
import parser.Values;
import prism.PrismException;
import prism.PrismLog;
import prism.Prism.StrategyExportType;
import explicit.Distribution;
import explicit.Model;
import explicit.SMG;
import explicit.STPG;
import explicit.STPGExplicit;

public class ExactValueStrategy implements Strategy
{
	// model info
	protected Values lastConstants;
	
	// strategy info
	protected String info = "No information available.";
	protected double[] minValues;
	protected double[] maxValues;

	// target values
	protected double initTargetValue;

	// storing last state
	protected int lastState;
	// storing the choice indicator
	protected boolean min;
	// last randomisation probability
	protected double probMin;

	// number of states
	protected int n;

	// initial distribution function [probMin_1, probMin_0,...]
	protected double[] initialDistributionFunction;

	// strategy memory update function [state -> [action,min -> [state -> probMin], action,max -> [state -> probMin]]]
	protected Map<Integer, Double>[][] memoryUpdateFunction;

	// strategy next move function [choiceMin_1,choiceMax_1, choiceMin_2,...]
	protected Distribution[] nextMoveFunction;

	// stats
	protected int memorySize;
	protected int initSize;
	protected int updateSize;
	protected int nextSize;

	/**
	 * 
	 * Creates a ExactValueStrategy.
	 *
	 * @param minStrat minimising strategy
	 * @param minValues expected values for states for min strategy
	 * @param maxStrat maximising strategy
	 * @param maxValues expected value for states for max strategy
	 * @param targetValue value to be achieved by the strategy
	 * @param model the model to provide info about players and transitions
	 */
	public ExactValueStrategy(Strategy minStrat, double[] minValues, Strategy maxStrat, double[] maxValues,
			double targetValue, STPG model)
	{
		this.minValues = minValues;
		this.maxValues = maxValues;
		this.initTargetValue = targetValue;
		lastState = -1;
		n = model.getNumStates();
		memorySize = 2 * n;
		initSize = n;
		updateSize = 0;
		nextSize = 2 * n;

		// create strategy 
		// computing initial distribution function
		initialDistributionFunction = new double[n];
		for (int s = 0; s < n; s++)
			initialDistributionFunction[s] = (maxValues[s] - initTargetValue) / (maxValues[s] - minValues[s]);

		// computing memory update function
		memoryUpdateFunction = new Map[n][];
		for (int s = 0; s < n; s++) {
			memoryUpdateFunction[s] = new Map[2 * model.getNumChoices(s)];
			for (int c = 0; c < model.getNumChoices(s); c++) {
				memoryUpdateFunction[s][2 * c] = new HashMap<Integer, Double>();
				memoryUpdateFunction[s][2 * c + 1] = new HashMap<Integer, Double>();

				// for player 2 state adjusting the exptectation
				if (model.getPlayer(s) == 2) {
					// computing min and max expectations for the action
					Iterator<Entry<Integer, Double>> it = model.getTransitionsIterator(s, c);
					double max = 0, min = 0;
					Entry<Integer, Double> en;
					List<Integer> succs = new ArrayList<Integer>(model.getNumTransitions(s, c));
					while (it.hasNext()) {
						en = it.next();
						succs.add(en.getKey());
						min += minValues[en.getKey()] * en.getValue();
						max += maxValues[en.getKey()] * en.getValue();
					}
					// computing the randomisation coefficient
					double pMin = (max - minValues[s]) / (max - min);
					double pMax = (max - maxValues[s]) / (max - min);
					// filling in the map
					for (int succ : succs) {
						memoryUpdateFunction[s][2 * c].put(succ, pMin);
						memoryUpdateFunction[s][2 * c + 1].put(succ, pMax);
						updateSize += 2;
					}
				}
				// for player 1 just updating to max
				else {
					Iterator<Entry<Integer, Double>> it = model.getTransitionsIterator(s, c);
					Entry<Integer, Double> en;
					while (it.hasNext()) {
						en = it.next();
						memoryUpdateFunction[s][2 * c].put(en.getKey(), 1.0);
						memoryUpdateFunction[s][2 * c + 1].put(en.getKey(), 0.0);
						updateSize += 2;
					}

				}
			}
		}

		// computing next move function
		nextMoveFunction = new Distribution[2 * n];
		for (int s = 0; s < n; s++) {
			try {
				nextMoveFunction[2 * s] = minStrat.getNextMove(s);
				nextMoveFunction[2 * s + 1] = maxStrat.getNextMove(s);
			} catch (InvalidStrategyStateException error) {
				error.printStackTrace();
			}

		}
	}

	/**
	 * Creates a ExactValueStrategy.
	 *
	 * @param scan
	 */
	public ExactValueStrategy(Scanner scan)
	{
		// TODO Auto-generated constructor stub
		for (int i = 0; i < 3; i++)
			scan.nextLine();
		scan.next();
		// parsing target value
		initTargetValue = scan.nextDouble();
		while (!scan.nextLine().startsWith("Initial"))
			;
		List<Double> minValuesL = new LinkedList<Double>();
		List<Double> maxValuesL = new LinkedList<Double>();
		List<Double> initL = new LinkedList<Double>();
		n = 0;
		String line;
		Scanner local;

		// parsing initial distribution function
		while (!(line = scan.nextLine()).startsWith("Memory")) {
			n++;
			local = new Scanner(line);
			local.next();
			initL.add(local.nextDouble());
			minValuesL.add(local.nextDouble());
			local.next();
			maxValuesL.add(local.nextDouble());
		}
		initialDistributionFunction = new double[n];
		for (int i = 0; i < n; i++)
			initialDistributionFunction[i] = initL.get(i);
		initSize = n;
		initL = null;

		updateSize = 0;
		// parsing memory update function
		int s1;
		int a;
		int s2;
		double minVal, minValSucc;
		double probMinMin;
		double maxVal, maxValSucc;
		double probMaxMin;
		Map<Integer, double[]> minMaxValues = new HashMap<Integer, double[]>();
		List<List<Map<Integer, Double>>> memUp = new LinkedList<List<Map<Integer, Double>>>();
		while (!(line = scan.nextLine()).startsWith("Next")) {
			// parsing the information in located in 4 lines
			local = new Scanner(line);
			s1 = local.nextInt();
			a = local.nextInt();
			s2 = local.nextInt();
			minVal = local.nextDouble();
			probMinMin = local.nextDouble();
			minValSucc = local.nextDouble();
			local.close();
			local = new Scanner(scan.nextLine());
			local.nextInt();
			local.nextInt();
			local.nextInt();
			local.nextDouble();
			local.nextDouble();
			maxValSucc = local.nextDouble();
			local.close();
			local = new Scanner(scan.nextLine());
			local.nextInt();
			local.nextInt();
			local.nextInt();
			maxVal = local.nextDouble();
			probMaxMin = local.nextDouble();
			local.close();
			scan.nextLine();

			// ------------ storing the data
			if (!minMaxValues.containsKey(s1))
				minMaxValues.put(s1, new double[] { minVal, maxVal });
			if (!minMaxValues.containsKey(s2))
				minMaxValues.put(s2, new double[] { minValSucc, maxValSucc });

			if (memUp.size() <= s1) // if this is the first time, adding new list for choices
				memUp.add(new LinkedList<Map<Integer, Double>>());
			if (memUp.get(s1).size() <= 2 * a) // if this is the first choice adding new map for choices
			{
				memUp.get(s1).add(new HashMap<Integer, Double>());
				memUp.get(s1).add(new HashMap<Integer, Double>());
			}
			// adding memory updates
			memUp.get(s1).get(2 * a).put(s2, probMinMin);
			memUp.get(s1).get(2 * a + 1).put(s2, probMaxMin);
			updateSize += 2;
		}
		// storing min-max values
		this.minValues = new double[minMaxValues.size()];
		this.maxValues = new double[minMaxValues.size()];
		for (Integer v : minMaxValues.keySet()) {
			this.minValues[v] = minMaxValues.get(v)[0];
			this.maxValues[v] = minMaxValues.get(v)[1];
		}
		n = minValues.length;
		memorySize = 2 * n;
		// storing the function
		this.memoryUpdateFunction = new Map[memUp.size()][];
		for (int i = 0; i < memUp.size(); i++) {
			this.memoryUpdateFunction[i] = new Map[memUp.get(i).size()];
			for (int j = 0; j < memUp.get(i).size(); j++)
				this.memoryUpdateFunction[i][j] = memUp.get(i).get(j);
		}

		// parsing nextmove function
		nextSize = 0;
		nextMoveFunction = new Distribution[2 * n];
		Distribution dist;
		while (scan.hasNext()) {
			local = new Scanner(scan.nextLine());
			s1 = local.nextInt();
			local.nextDouble();
			s2 = local.nextInt();
			dist = new Distribution();
			dist.add(s2, 1.0);
			nextMoveFunction[2 * s1] = dist;
			local.close();
			local = new Scanner(scan.nextLine());
			local.nextInt();
			local.nextDouble();
			s2 = local.nextInt();
			dist = new Distribution();
			dist.add(s2, 1.0);
			nextMoveFunction[2 * s1 + 1] = dist;
			local.close();
			nextSize += 2;
		}

		lastState = -1;
		min = false;
		probMin = 0;
	}

	@Override
	public void init(int state) throws InvalidStrategyStateException
	{
		min = Math.random() < initialDistributionFunction[state];
		lastState = state;
		probMin = initialDistributionFunction[state];
	}

	@Override
	public void updateMemory(int action, int state) throws InvalidStrategyStateException
	{
		probMin = memoryUpdateFunction[lastState][2 * action + (min ? 0 : 1)].get(state);
		min = Math.random() < probMin;
		lastState = state;
	}

	@Override
	public Distribution getNextMove(int state) throws InvalidStrategyStateException
	{
		return nextMoveFunction[2 * state + (min ? 0 : 1)];
	}

	@Override
	public void reset()
	{
		lastState = -1;
	}

	@Override
	public void exportToFile(String file)
	{
		// Print adversary
		FileWriter out = null;
		try {
			out = new FileWriter(new File(file));

			out.write(Strategies.FORMAT_STRING_EXACT_VALUE_MD_STRAT + "\n");
			out.write("// Stochastic update strategy to achieve exact value in the game\n");
			out.write("// Format (memory update function): \n");
			out.write("Strategy:\n");
			out.write("targetValue " + initTargetValue + "\n");
			out.write("Initial distribution (stateId, minProb, minValue, maxProb, maxValue):\n");
			for (int s = 0; s < n; s++) {
				out.write(s + " " + initialDistributionFunction[s] + " " + minValues[s] + " "
						+ (1 - initialDistributionFunction[s]) + " " + maxValues[s] + "\n");
			}

			out
					.write("Memory update function (stateId, choiceId, successorId, memoryValue, probability, newMemoryValue):\n");
			for (int s = 0; s < n; s++) {
				for (int c = 0; c < memoryUpdateFunction[s].length / 2; c++) {
					for (int succ : memoryUpdateFunction[s][2 * c].keySet()) {
						out.write(s + " " + c + " " + succ + " " + minValues[s] + " "
								+ memoryUpdateFunction[s][2 * c].get(succ) + " " + minValues[succ] + "\n");
						out.write(s + " " + c + " " + succ + " " + minValues[s] + " "
								+ (1 - memoryUpdateFunction[s][2 * c].get(succ)) + " " + maxValues[succ] + "\n");
						out.write(s + " " + c + " " + succ + " " + maxValues[s] + " "
								+ memoryUpdateFunction[s][2 * c + 1].get(succ) + " " + minValues[succ] + "\n");
						out.write(s + " " + c + " " + succ + " " + maxValues[s] + " "
								+ (1 - memoryUpdateFunction[s][2 * c + 1].get(succ)) + " " + maxValues[succ] + "\n");
					}
				}
			}

			out.write("Next move function (stateId, memoryElement, choiceId):\n");
			for (int s = 0; s < n; s++) {
				out.write(s + " " + minValues[s] + " " + nextMoveFunction[2 * s].keySet().iterator().next() + "\n");
				out.write(s + " " + maxValues[s] + " " + nextMoveFunction[2 * s + 1].keySet().iterator().next() + "\n");
			}

			out.flush();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} finally {
			if (out != null)
				try {
					out.close();
				} catch (IOException e) {
					// do nothing
				}
		}
	}

	@Override
	public Model buildProduct(Model model) throws PrismException
	{
		// checking for supported model types
		if (model.getClass().equals(STPGExplicit.class)) {
			return this.buildProductSTPGExplicit((STPGExplicit) model);
		}
		if (model.getClass().equals(SMG.class)) {
			return this.buildProductSMG((SMG) model);
		}

		throw new UnsupportedOperationException("The product building is not supported for this class of models");

	}

	/**
	 *
	 * @param model
	 * @return
	 * @throws PrismException 
	 */
	private Model buildProductSMG(SMG model) throws PrismException
	{
		// construct a new STPG of size three times the original model
		SMG smg = new SMG(3 * model.getNumStates());
		smg.copyPlayerInfo(model);
		smg.copyCoalitionInfo(model);
		int n = smg.getNumStates();

		List<State> oldStates = model.getStatesList();

		// creating helper states
		State stateInit = new State(1), stateMin = new State(1), stateMax = new State(1);
		stateInit.setValue(0, 0); // state where memory is not yet initialised
		stateMin.setValue(0, 1); // state where target is minimum elem
		stateMax.setValue(0, 2); // state where target is maximum element

		// creating product state list
		List<State> newStates = new ArrayList<State>(n);
		int count = 0;
		for (int i = 0; i < oldStates.size(); i++) {
			newStates.add(new State(oldStates.get(i), stateInit));
			newStates.add(new State(oldStates.get(i), stateMin));
			newStates.add(new State(oldStates.get(i), stateMax));
			smg.setPlayer(count++, model.getPlayer(i));
			smg.setPlayer(count++, model.getPlayer(i));
			smg.setPlayer(count++, model.getPlayer(i));
		}

		// setting the states list to STPG
		smg.setStatesList(newStates);

		// adding choices for the product STPG
		// initial distributions
		int indx, minIndx, maxIndx;
		Distribution distr, newDistr;
		double p;
		for (int i = 0; i < oldStates.size(); i++) {
			indx = 3 * i;
			// build only for states from which the value is achievable
			if (minValues[i] <= initTargetValue && initTargetValue <= maxValues[i]) {
				p = initialDistributionFunction[i];
				distr = new Distribution();
				if (p != 0)
					distr.add(indx + 1, p);
				if (p != 1)
					distr.add(indx + 2, 1 - p);
				smg.addChoice(indx, distr);
			}
			//
			else {
				//Add self-loop only
				distr = new Distribution();
				distr.add(indx, 1.0);
				smg.addChoice(indx, distr);
			}
		}

		// all other states
		Distribution distrMin, distrMax;
		for (int i = 0; i < oldStates.size(); i++) {
			minIndx = 3 * i + 1;
			maxIndx = 3 * i + 2;

			if (model.getPlayer(i) == 1) {
				// for player 1 state just leaving max or min choice
				// constructing transition for min element
				newDistr = new Distribution();
				distr = model.getChoice(i, nextMoveFunction[2 * i].keySet().iterator().next());
				// create a new distribution for the product
				newDistr = new Distribution();
				for (Integer succ : distr.keySet())
					// adding transition to the state with the memory min memory element
					newDistr.add(succ * 3 + 1, distr.get(succ));
				// adding the choice
				smg.addChoice(minIndx, newDistr);

				// constructing transition for max element
				newDistr = new Distribution();
				distr = model.getChoice(i, nextMoveFunction[2 * i + 1].keySet().iterator().next());
				// create a new distribution for the product
				newDistr = new Distribution();
				for (Integer succ : distr.keySet())
					// adding transition to the state with the memory min memory element
					newDistr.add(succ * 3 + 2, distr.get(succ));
				// adding the choice
				smg.addChoice(maxIndx, newDistr);

			} else {
				// for player 2 state transforming every distribution#
				// computing the probability to choose min strategy
				for (int action = 0; action < model.getNumChoices(i); action++) {
					// constructing distributions
					distrMin = new Distribution();
					distrMax = new Distribution();
					distr = model.getChoice(i, action);
					for (Integer succ : distr.keySet()) {
						distrMin.add(succ * 3 + 1, distr.get(succ) * memoryUpdateFunction[i][2 * action].get(succ));
						distrMin.add(succ * 3 + 2, distr.get(succ)
								* (1 - memoryUpdateFunction[i][2 * action].get(succ)));

						distrMax.add(succ * 3 + 1, distr.get(succ) * memoryUpdateFunction[i][2 * action + 1].get(succ));
						distrMax.add(succ * 3 + 2, distr.get(succ)
								* (1 - memoryUpdateFunction[i][2 * action + 1].get(succ)));
					}

					smg.addChoice(minIndx, distrMin);
					smg.addChoice(maxIndx, distrMax);
				}
			}
		}

		// setting initial state for the game
		smg.addInitialState(0);

		return smg;
	}

	/**
	 *
	 * @param model
	 * @return
	 */
	private Model buildProductSTPGExplicit(STPGExplicit model)
	{

		// construct a new STPG of size three times the original model
		STPGExplicit stpg = new STPGExplicit(3 * model.getNumStates());
		int n = stpg.getNumStates();

		List<State> oldStates = model.getStatesList();

		// creating helper states
		State stateInit = new State(1), stateMin = new State(1), stateMax = new State(1);
		stateInit.setValue(0, 0); // state where memory is not yet initialised
		stateMin.setValue(0, 1); // state where target is minimum elem
		stateMax.setValue(0, 2); // state where target is maximum element

		// creating product state list
		List<State> newStates = new ArrayList<State>(n);
		for (int i = 0; i < oldStates.size(); i++) {
			newStates.add(new State(oldStates.get(i), stateInit));
			newStates.add(new State(oldStates.get(i), stateMin));
			newStates.add(new State(oldStates.get(i), stateMax));
		}

		// setting the states list to STPG
		stpg.setStatesList(newStates);

		// adding choices for the product STPG
		// initial distributions
		int indx, minIndx, maxIndx;
		Distribution distr, newDistr;
		double p;
		for (int i = 0; i < oldStates.size(); i++) {
			indx = 3 * i;
			// build only for states from which the value is achievable
			if (minValues[i] <= initTargetValue && initTargetValue <= maxValues[i]) {
				p = initialDistributionFunction[i];
				distr = new Distribution();
				if (p != 0)
					distr.add(indx + 1, p);
				if (p != 1)
					distr.add(indx + 2, 1 - p);
				stpg.addChoice(indx, distr);
			}
			//
			else {
				//Add self-loop only
				distr = new Distribution();
				distr.add(indx, 1.0);
				stpg.addChoice(indx, distr);
			}
		}

		// all other states
		double pmin, pmax;
		Distribution distrMin, distrMax;
		for (int i = 0; i < oldStates.size(); i++) {
			minIndx = 3 * i + 1;
			maxIndx = 3 * i + 2;

			if (model.getPlayer(i) == 1) {
				// for player 1 state just leaving max or min choice
				// constructing transition for min element
				newDistr = new Distribution();
				distr = model.getChoice(i, nextMoveFunction[2 * i].keySet().iterator().next());
				// create a new distribution for the product
				newDistr = new Distribution();
				for (Integer succ : distr.keySet())
					// adding transition to the state with the memory min memory element
					newDistr.add(succ * 3 + 1, distr.get(succ));
				// adding the choice
				stpg.addChoice(minIndx, newDistr);

				// constructing transition for max element
				newDistr = new Distribution();
				distr = model.getChoice(i, nextMoveFunction[2 * i + 1].keySet().iterator().next());
				// create a new distribution for the product
				newDistr = new Distribution();
				for (Integer succ : distr.keySet())
					// adding transition to the state with the memory min memory element
					newDistr.add(succ * 3 + 2, distr.get(succ));
				// adding the choice
				stpg.addChoice(maxIndx, newDistr);

			} else {
				// for player 2 state transforming every distribution#
				// computing the probability to choose min strategy

				for (int action = 0; action < model.getNumChoices(i); action++) {
					// constructing distributions
					distrMin = new Distribution();
					distrMax = new Distribution();
					distr = model.getChoice(i, action);
					for (Integer succ : distr.keySet()) {
						distrMin.add(succ * 3 + 1, distr.get(succ) * memoryUpdateFunction[i][2 * action].get(succ));
						distrMin.add(succ * 3 + 2, distr.get(succ)
								* (1 - memoryUpdateFunction[i][2 * action].get(succ)));

						distrMax.add(succ * 3 + 1, distr.get(succ) * memoryUpdateFunction[i][2 * action + 1].get(succ));
						distrMax.add(succ * 3 + 2, distr.get(succ)
								* (1 - memoryUpdateFunction[i][2 * action + 1].get(succ)));
					}

					stpg.addChoice(minIndx, distrMin);
					stpg.addChoice(maxIndx, distrMax);
				}

			}
		}

		// setting initial state for the game
		stpg.addInitialState(0);

		return stpg;
	}

	@Override
	public String getInfo()
	{
		return info;
	}

	@Override
	public void setInfo(String info)
	{
		this.info = info;
	}

	@Override
	public int getMemorySize()
	{
		return memorySize;
	}

	@Override
	public String getType()
	{
		return "Stochastic update strategy.";
	}

	@Override
	public Object getCurrentMemoryElement()
	{
		return new Object[] { lastState, min, probMin };
	}

	@Override
	public void setMemory(Object memory) throws InvalidStrategyStateException
	{
		if (memory instanceof Object[] && ((Object[]) memory).length == 3 && ((Object[]) memory)[0] instanceof Integer
				&& ((Object[]) memory)[1] instanceof Boolean && ((Object[]) memory)[2] instanceof Double) {
			lastState = (Integer) ((Object[]) memory)[0];
			min = (Boolean) ((Object[]) memory)[1];
			probMin = (Double) ((Object[]) memory)[2];
		} else
			throw new InvalidStrategyStateException("Memory element has to be Object array of length 3.");

	}

	private DecimalFormat df = new DecimalFormat("#.###");

	@Override
	public String getDescription()
	{
		String desc = "";
		desc += "Stochastic update strategy\n";
		desc += "Target expectation: " + initTargetValue + "\n";
		desc += "Size of memory: " + getMemorySize() + "\n";
		desc += "Size of initial dist. function: " + initSize + "\n";
		desc += "Size of memory update function: " + updateSize + "\n";
		desc += "Size of next move function: " + nextSize + "\n";
		/*desc += "Current target expectation: "
				+ (lastState < 0 ? initTargetValue : df.format(min ? minValues[lastState] : maxValues[lastState]))
				+ "\n";
		desc += "Last memory update: "
				+ (lastState < 0 ? initialDistributionFunction[0] : df.format(minValues[lastState])) + "->"
				+ df.format(probMin) + ", "
				+ (lastState < 0 ? initialDistributionFunction[0] : df.format(maxValues[lastState])) + "->"
				+ df.format((1 - probMin)) + "\n";*/

		return desc;
	}

	@Override
	public int getInitialStateOfTheProduct(int s)
	{
		return 0;
	}

	public void export(PrismLog out) {}

	@Override
	public void exportActions(PrismLog out)
	{
		// TODO Auto-generated method stub
		
	}

	@Override
	public void clear()
	{
		// TODO Auto-generated method stub
		
	}

	@Override
	public void exportIndices(PrismLog out)
	{
		// TODO Auto-generated method stub
		
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
		
	}

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
