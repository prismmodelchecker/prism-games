package strat;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.Map.Entry;

import parser.State;
import prism.PrismException;
import prism.PrismFileLog;
import prism.PrismLog;
import explicit.Distribution;
import explicit.Model;
import explicit.SMG;
import explicit.STPG;
import explicit.STPGExplicit;

public class ExactValueStrategy implements Strategy
{
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
				if (model.getPlayer(s) == STPGExplicit.PLAYER_2) {
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
		min = Math.random() < memoryUpdateFunction[lastState][2 * action + (min ? 0 : 1)].get(state);
		lastState = state;
		probMin = memoryUpdateFunction[lastState][2 * action + (min ? 0 : 1)].get(state);
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
		PrismLog out = new PrismFileLog(file);
		out.println(Strategies.FORMAT_STRING_EXACT_VALUE_MD_STRAT);
		out.print("// Stochastic update strategy to achieve exact value in the game\n");
		out.print("// Format (memory update function): \n");
		out.print("Strategy:\n");
		out.print("targetValue=" + initTargetValue + "\n");
		out.print("Initial distribution (stateId, minProb, minValue, maxProb, maxValue):\n");
		for (int s = 0; s < n; s++) {
			out.println(s + ", " + initialDistributionFunction[s] + ", " + minValues[s] + ", "
					+ (1 - initialDistributionFunction[s]) + ", " + maxValues[s]);
		}

		out
				.print("Memory update function (stateId, choiceId, successorId, memoryValue, probability, newMemoryValue):\n");
		for (int s = 0; s < n; s++) {
			for (int c = 0; c < memoryUpdateFunction[s].length / 2; c++) {
				for (int succ : memoryUpdateFunction[s][2 * c].keySet()) {
					out.println(s + ", " + c + ", " + succ + ", " + minValues[s] + ", "
							+ memoryUpdateFunction[s][2 * c].get(succ) + ", " + minValues[succ]);
					out.println(s + ", " + c + ", " + succ + ", " + minValues[s] + ", "
							+ (1 - memoryUpdateFunction[s][2 * c].get(succ)) + ", " + maxValues[succ]);
					out.println(s + ", " + c + ", " + succ + ", " + maxValues[s] + ", "
							+ memoryUpdateFunction[s][2 * c + 1].get(succ) + ", " + minValues[succ]);
					out.println(s + ", " + c + ", " + succ + ", " + maxValues[s] + ", "
							+ (1 - memoryUpdateFunction[s][2 * c + 1].get(succ)) + ", " + maxValues[succ]);
				}
			}
		}

		out.print("Next move function (stateId, memoryElement, choiceId):\n");
		for (int s = 0; s < n; s++) {
			out.println(s + ", " + minValues[s] + ", " + nextMoveFunction[2 * s].keySet().iterator().next());
			out.println(s + ", " + maxValues[s] + ", " + nextMoveFunction[2 * s + 1].keySet().iterator().next());
		}

		out.flush();
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
		smg.setPlayerMapping(new HashMap<String, Integer>(model.getPlayerMapping()));
		smg.setCoalitionInts(new HashSet<Integer>(model.getCoalition()));
		int n = smg.getNumStates();

		List<Integer> stateLabels = new ArrayList<Integer>(n);

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
			stateLabels.add(model.getPlayer(i));
			stateLabels.add(model.getPlayer(i));
			stateLabels.add(model.getPlayer(i));
		}

		// setting the states list to STPG
		smg.setStatesList(newStates);
		smg.setStateLabels(stateLabels);

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

			if (model.getPlayer(i) == STPGExplicit.PLAYER_1) {
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

			if (model.getPlayer(i) == STPGExplicit.PLAYER_1) {
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
		return new Object[] { lastState, min };
	}

	@Override
	public void setMemory(Object memory) throws InvalidStrategyStateException
	{
		if (memory instanceof Object[] && ((Object[]) memory).length == 2 && ((Object[]) memory)[0] instanceof Integer
				&& ((Object[]) memory)[1] instanceof Boolean) {
			lastState = (Integer) ((Object[]) memory)[0];
			min = (Boolean) ((Object[]) memory)[1];
		} else
			throw new InvalidStrategyStateException("Memory element has to be Object array of length 2.");

	}

	private DecimalFormat df = new DecimalFormat("#.###");

	@Override
	public String getStateDescription()
	{
		String desc = "";
		desc += "Stochastic update strategy\n";
		desc += "Target expectation: " + initTargetValue + "\n";
		desc += "Size of memory: " + getMemorySize() + "\n";
		desc += "Size of initial dist. function: " + initSize + "\n";
		desc += "Size of memory update function: " + updateSize + "\n";
		desc += "Size of next move function: " + nextSize + "\n";
		desc += "Current target expectation: "
				+ (lastState < 0 ? initTargetValue : df.format(min ? minValues[lastState] : maxValues[lastState]))
				+ "\n";
		desc += "Last memory update: "
				+ (lastState < 0 ? initialDistributionFunction[0] : df.format(minValues[lastState])) + "->"
				+ df.format(probMin) + ", "
				+ (lastState < 0 ? initialDistributionFunction[0] : df.format(maxValues[lastState])) + "->"
				+ df.format((1 - probMin)) + "\n";

		return desc;
	}

	@Override
	public int getInitialStateOfTheProduct(int s)
	{
		return 0;
	}

	//	@Override
	//	public double getExpectedValue()
	//	{
	//		// TODO Auto-generated method stub
	//		return 0;
	//	}
	//
	//	@Override
	//	public double getExpectedValue(int a, int s)
	//	{
	//		// TODO Auto-generated method stub
	//		return 0;
	//	}

}
