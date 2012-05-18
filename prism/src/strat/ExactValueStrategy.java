package strat;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;

import prism.PrismException;
import prism.PrismFileLog;
import prism.PrismLog;
import explicit.Distribution;
import explicit.Model;
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
		// TODO Auto-generated method stub
		return null;
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
		// TODO Auto-generated method stub
		return -1;
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
