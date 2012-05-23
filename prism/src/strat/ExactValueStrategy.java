package strat;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
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

	// strategies to achieve optimal values
	protected Strategy minStrat, maxStrat;
	protected double[] minValues, maxValues;

	// indicator which strategy to play next
	protected boolean playMin;

	// strategy info
	protected String info = "No information available.";

	// target values
	protected double initTargetValue, currentTargetValue, probMin;

	// storing last state
	protected int lastState;

	// game model
	protected STPG game;

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
		this.minStrat = minStrat;
		this.minValues = minValues;
		this.maxStrat = maxStrat;
		this.maxValues = maxValues;
		playMin = false;
		this.initTargetValue = targetValue;
		lastState = -1;
		game = model;
	}

	@Override
	public void init(int state) throws InvalidStrategyStateException
	{

		minStrat.init(state);
		maxStrat.init(state);

		double minVal = minValues[state];
		double maxVal = maxValues[state];

		probMin = (maxVal - initTargetValue) / (maxVal - minVal);
		currentTargetValue = Math.random() < probMin ? minVal : maxVal;

		playMin = currentTargetValue == minVal;

		lastState = state;
	}

	@Override
	public void updateMemory(int action, int state) throws InvalidStrategyStateException
	{

		// computing the probability to choose min strategy
		probMin = 0;
		if (game.getPlayer(lastState) == STPGExplicit.PLAYER_1)
			probMin = playMin ? 1 : 0;
		else {
			// computing min and max expectations for the action
			Iterator<Entry<Integer, Double>> it = game.getTransitionsIterator(lastState, action);
			double max = 0, min = 0;
			Entry<Integer, Double> en;
			while (it.hasNext()) {
				en = it.next();
				min += minValues[en.getKey()] * en.getValue();
				max += maxValues[en.getKey()] * en.getValue();
			}
			// computing the randomisation coefficient
			probMin = (max - currentTargetValue) / (max - min);
		}

		minStrat.updateMemory(action, state);
		maxStrat.updateMemory(action, state);

		double minVal = minValues[state];
		double maxVal = maxValues[state];

		// determining the new current value
		currentTargetValue = Math.random() < probMin ? minVal : maxVal;
		playMin = currentTargetValue == minVal;
		lastState = state;
	}

	@Override
	public Distribution getNextMove(int state) throws InvalidStrategyStateException
	{
		return playMin ? minStrat.getNextMove(state) : maxStrat.getNextMove(state);
	}

	@Override
	public void reset()
	{
		minStrat.reset();
		maxStrat.reset();
		this.currentTargetValue = initTargetValue;
	}

	@Override
	public void exportToFile(String file)
	{
		// Print adversary
		PrismLog out = new PrismFileLog(file);
		out.print("// Stochastic update strategy to achieve exact value in the game\n");
		out.print("// Format (memory update function): \n");
		out.print("Strategy:\n");
		out.print("targetValue=" + initTargetValue + "\n");
		out.print("Initial distribution (stateId, minProb, minValue, maxProb, maxValue):\n");
		for (int s = 0; s < game.getNumStates(); s++) {
			double minVal = minValues[s];
			double maxVal = maxValues[s];
			double probMin = (maxVal - initTargetValue) / (maxVal - minVal);
			if (maxVal == minVal)
				out.println(s + ", " + 1 + ", " + minVal + ", " + 0 + ", " + maxVal);
			else
				out.println(s + ", " + probMin + ", " + minVal + ", " + (1 - probMin) + ", " + maxVal);
		}

		out
				.print("Memory update function (stateId, choiceId, successorId, memoryValue, probability, newMemoryValue):\n");
		for (int s = 0; s < game.getNumStates(); s++) {
			for (int c = 0; c < game.getNumChoices(s); c++) {
				if (game.getPlayer(s) == STPGExplicit.PLAYER_2) {
					// computing min and max expectations for the action
					Iterator<Entry<Integer, Double>> it = game.getTransitionsIterator(s, c);
					double max = 0, min = 0;
					Entry<Integer, Double> en;
					List<Integer> succs = new ArrayList<Integer>(10);
					while (it.hasNext()) {
						en = it.next();
						succs.add(en.getKey());
						min += minValues[en.getKey()] * en.getValue();
						max += maxValues[en.getKey()] * en.getValue();
					}

					// computing the randomisation coefficient
					double probMin = (max - minValues[s]) / (max - min);
					double probMax = (max - maxValues[s]) / (max - min);

					for (int succ : succs) {
						out.println(s + ", " + c + ", " + succ + ", " + minValues[s] + ", " + probMin + ", "
								+ minValues[succ]);
						out.println(s + ", " + c + ", " + succ + ", " + minValues[s] + ", " + (1 - probMin) + ", "
								+ maxValues[succ]);
						out.println(s + ", " + c + ", " + succ + ", " + maxValues[s] + ", " + probMax + ", "
								+ minValues[succ]);
						out.println(s + ", " + c + ", " + succ + ", " + maxValues[s] + ", " + (1 - probMax) + ", "
								+ maxValues[succ]);
					}
				} else {
					Iterator<Entry<Integer, Double>> it = game.getTransitionsIterator(s, c);
					Entry<Integer, Double> en;
					while (it.hasNext()) {
						en = it.next();
						out.println(s + ", " + c + ", " + en.getKey() + ", " + minValues[s] + ", " + 1 + ", "
								+ minValues[en.getKey()]);
						out.println(s + ", " + c + ", " + en.getKey() + ", " + maxValues[s] + ", " + 1 + ", "
								+ minValues[en.getKey()]);
					}

				}
			}
		}

		out.print("Next move function (stateId, memoryElement, choiceId):\n");
		for (int s = 0; s < game.getNumStates(); s++) {
			if (game.getPlayer(s) == STPGExplicit.PLAYER_1)
				try {
					out.println(s + ", " + minValues[s] + ", " + minStrat.getNextMove(s).keySet().iterator().next());
					out.println(s + ", " + maxValues[s] + ", " + maxStrat.getNextMove(s).keySet().iterator().next());
				} catch (InvalidStrategyStateException error) {
					error.printStackTrace();
				}
			else
				out.println(s + " n/a");
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
				p = (maxValues[i] - initTargetValue) / (maxValues[i] - minValues[i]);
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
		double pmin, pmax;
		Distribution distrMin, distrMax;
		for (int i = 0; i < oldStates.size(); i++) {
			minIndx = 3 * i + 1;
			maxIndx = 3 * i + 2;

			if (model.getPlayer(i) == STPGExplicit.PLAYER_1) {
				// for player 1 state just leaving max or min choice
				try {
					// constructing transition for min element
					newDistr = new Distribution();
					distr = model.getChoice(i, minStrat.getNextMove(i).keySet().iterator().next());
					// create a new distribution for the product
					newDistr = new Distribution();
					for (Integer succ : distr.keySet())
						// adding transition to the state with the memory min memory element
						newDistr.add(succ * 3 + 1, distr.get(succ));
					// adding the choice
					smg.addChoice(minIndx, newDistr);

					// constructing transition for max element
					newDistr = new Distribution();
					distr = model.getChoice(i, maxStrat.getNextMove(i).keySet().iterator().next());
					// create a new distribution for the product
					newDistr = new Distribution();
					for (Integer succ : distr.keySet())
						// adding transition to the state with the memory min memory element
						newDistr.add(succ * 3 + 2, distr.get(succ));
					// adding the choice
					smg.addChoice(maxIndx, newDistr);

				} catch (InvalidStrategyStateException error) {
					// TODO Auto-generated catch block
					error.printStackTrace();
				}

			} else {
				// for player 2 state transforming every distribution#
				// computing the probability to choose min strategy

				for (int action = 0; action < model.getNumChoices(i); action++) {
					// computing min and max expectations for the action
					Iterator<Entry<Integer, Double>> it = model.getTransitionsIterator(i, action);
					double max = 0, min = 0;
					Entry<Integer, Double> en;
					while (it.hasNext()) {
						en = it.next();
						min += minValues[en.getKey()] * en.getValue();
						max += maxValues[en.getKey()] * en.getValue();
					}

					// computing the randomisation coefficient for min and max values
					pmin = (max - minValues[i]) / (max - min);
					pmax = (max - maxValues[i]) / (max - min);

					// constructing distributions
					distrMin = new Distribution();
					distrMax = new Distribution();
					distr = model.getChoice(i, action);
					for (Integer succ : distr.keySet()) {
						distrMin.add(succ * 3 + 1, distr.get(succ) * pmin);
						distrMin.add(succ * 3 + 2, distr.get(succ) * (1 - pmin));

						distrMax.add(succ * 3 + 1, distr.get(succ) * pmax);
						distrMax.add(succ * 3 + 2, distr.get(succ) * (1 - pmax));
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
				p = (maxValues[i] - initTargetValue) / (maxValues[i] - minValues[i]);
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
				try {
					// constructing transition for min element
					newDistr = new Distribution();
					distr = model.getChoice(i, minStrat.getNextMove(i).keySet().iterator().next());
					// create a new distribution for the product
					newDistr = new Distribution();
					for (Integer succ : distr.keySet())
						// adding transition to the state with the memory min memory element
						newDistr.add(succ * 3 + 1, distr.get(succ));
					// adding the choice
					stpg.addChoice(minIndx, newDistr);

					// constructing transition for max element
					newDistr = new Distribution();
					distr = model.getChoice(i, maxStrat.getNextMove(i).keySet().iterator().next());
					// create a new distribution for the product
					newDistr = new Distribution();
					for (Integer succ : distr.keySet())
						// adding transition to the state with the memory min memory element
						newDistr.add(succ * 3 + 2, distr.get(succ));
					// adding the choice
					stpg.addChoice(maxIndx, newDistr);

				} catch (InvalidStrategyStateException error) {
					// TODO Auto-generated catch block
					error.printStackTrace();
				}

			} else {
				// for player 2 state transforming every distribution#
				// computing the probability to choose min strategy

				for (int action = 0; action < model.getNumChoices(i); action++) {
					// computing min and max expectations for the action
					Iterator<Entry<Integer, Double>> it = model.getTransitionsIterator(i, action);
					double max = 0, min = 0;
					Entry<Integer, Double> en;
					while (it.hasNext()) {
						en = it.next();
						min += minValues[en.getKey()] * en.getValue();
						max += maxValues[en.getKey()] * en.getValue();
					}

					// computing the randomisation coefficient for min and max values
					pmin = (max - minValues[i]) / (max - min);
					pmax = (max - maxValues[i]) / (max - min);

					// constructing distributions
					distrMin = new Distribution();
					distrMax = new Distribution();
					distr = model.getChoice(i, action);
					for (Integer succ : distr.keySet()) {
						distrMin.add(succ * 3 + 1, distr.get(succ) * pmin);
						distrMin.add(succ * 3 + 2, distr.get(succ) * (1 - pmin));

						distrMax.add(succ * 3 + 1, distr.get(succ) * pmax);
						distrMax.add(succ * 3 + 2, distr.get(succ) * (1 - pmax));
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
		return minValues.length + maxValues.length;
	}

	@Override
	public String getType()
	{
		return "Stochastic update strategy.";
	}

	@Override
	public Object getCurrentMemoryElement()
	{
		return new Object[] { currentTargetValue, minStrat.getCurrentMemoryElement(),
				maxStrat.getCurrentMemoryElement() };
	}

	@Override
	public void setMemory(Object memory) throws InvalidStrategyStateException
	{
		if (memory instanceof Object[] && ((Object[]) memory).length == 3 && ((Object[]) memory)[0] instanceof Double) {
			this.currentTargetValue = (Double) ((Object[]) memory)[0];
			this.minStrat.setMemory(((Object[]) memory)[1]);
			this.maxStrat.setMemory(((Object[]) memory)[2]);
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
		desc += "Size of next move function: " + getMemorySize() + "\n";
		desc += "Current target expectation: " + df.format(currentTargetValue) + "\n";
		desc += "Last memory update: " + df.format(minValues[lastState]) + "->" + df.format(probMin) + ", "
				+ df.format(maxValues[lastState]) + "->" + df.format((1 - probMin)) + "\n";

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
