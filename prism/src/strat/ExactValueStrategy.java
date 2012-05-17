package strat;

import java.util.Iterator;
import java.util.Map.Entry;

import prism.PrismException;
import explicit.Distribution;
import explicit.Model;
import explicit.STPG;
import explicit.STPGExplicit;

public class ExactValueStrategy implements Strategy {

	// strategies to achieve optimal values
	protected Strategy minStrat, maxStrat;

	// indicator which strategy to play next
	protected boolean playMin;

	// strategy info
	protected String info = "No information available.";

	// target values
	protected double initTargetValue, currentTargetValue;

	// storing last state
	protected int lastState;

	// game model
	protected STPG game;

	/**
	 * Creates the exact value strategy
	 * 
	 * @param minStrat
	 *            strategy which minimises the value
	 * @param maxStrat
	 *            strategy which maximises the value
	 * @param targetValue
	 *            the value to be achieved by a strategy
	 */
	public ExactValueStrategy(Strategy minStrat, Strategy maxStrat,
			double targetValue, STPG model) {
		this.minStrat = minStrat;
		this.maxStrat = maxStrat;
		playMin = false;
		this.initTargetValue = targetValue;
		lastState = -1;
		game = model;
	}

	@Override
	public void init(int state) throws InvalidStrategyStateException {

		minStrat.init(state);
		maxStrat.init(state);

		double minVal = minStrat.getExpectedValue();
		double maxVal = maxStrat.getExpectedValue();

		currentTargetValue = Math.random() < (maxVal - initTargetValue)
				/ (maxVal - minVal) ? minVal : maxVal;

		playMin = currentTargetValue == minVal;

		lastState = state;
	}

	@Override
	public void updateMemory(int action, int state)
			throws InvalidStrategyStateException {

		// computing the probability to choose min strategy
		double probMin = 0;
		if (game.getPlayer(lastState) == STPGExplicit.PLAYER_1)
			probMin = playMin ? 1 : 0;
		else {
			// computing min and max expectations for the action
			Iterator<Entry<Integer, Double>> it = game.getTransitionsIterator(
					lastState, action);
			double max = 0, min = 0;
			Entry<Integer, Double> en;
			while ((en = it.next()) != null) {
				min += minStrat.getExpectedValue(action, en.getKey())
						* en.getValue();
				max += maxStrat.getExpectedValue(action, en.getKey())
						* en.getValue();
			}
			// computing the randomisation coefficient
			probMin = (max - currentTargetValue) / (max - min);
		}

		minStrat.updateMemory(action, state);
		maxStrat.updateMemory(action, state);

		double minVal = minStrat.getExpectedValue();
		double maxVal = maxStrat.getExpectedValue();

		// determining the new current value
		currentTargetValue = Math.random() < probMin ? minVal : maxVal;
		playMin = currentTargetValue == minVal;
	}

	@Override
	public Distribution getNextMove(int state)
			throws InvalidStrategyStateException {
		return playMin ? minStrat.getNextMove(state) : maxStrat
				.getNextMove(state);
	}

	@Override
	public void reset() {
		minStrat.reset();
		maxStrat.reset();
		this.currentTargetValue = initTargetValue;
	}

	@Override
	public void exportToFile(String file) {
		// TODO Auto-generated method stub

	}

	@Override
	public Model buildProduct(Model model) throws PrismException {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String getInfo() {
		return info;
	}

	@Override
	public void setInfo(String info) {
		this.info = info;
	}

	@Override
	public int getMemorySize() {
		return minStrat.getMemorySize() + maxStrat.getMemorySize();
	}

	@Override
	public String getType() {
		return "Stochastic update strategy.";
	}

	@Override
	public Object getCurrentMemoryElement() {
		return new Object[] { currentTargetValue,
				minStrat.getCurrentMemoryElement(),
				maxStrat.getCurrentMemoryElement() };
	}

	@Override
	public void setMemory(Object memory) throws InvalidStrategyStateException {
		if (memory instanceof Object[] && ((Object[]) memory).length == 3
				&& ((Object[]) memory)[0] instanceof Double) {
			this.currentTargetValue = (Double) ((Object[]) memory)[0];
			this.minStrat.setMemory(((Object[]) memory)[1]);
			this.maxStrat.setMemory(((Object[]) memory)[2]);
		} else
			throw new InvalidStrategyStateException(
					"Memory element has to be Object array of length 2.");

	}

	@Override
	public String getStateDescription() {
		return "No descrption available";
	}

	@Override
	public int getInitialStateOfTheProduct(int s) {
		// TODO Auto-generated method stub
		return -1;
	}

	@Override
	public double getExpectedValue() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double getExpectedValue(int a, int s) {
		// TODO Auto-generated method stub
		return 0;
	}

}
