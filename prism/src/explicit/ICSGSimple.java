//==============================================================================
//
//	Copyright (c) 2025-
//	Authors:
//	* Dave Parker <david.parker@cs.ox.ac.uk> (University of Oxford)
//
//------------------------------------------------------------------------------
//
//	This file is part of PRISM.
//
//	PRISM is free software; you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation; either version 2 of the License, or
//	(at your option) any later version.
//
//	PRISM is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with PRISM; if not, write to the Free Software Foundation,
//	Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
//==============================================================================

package explicit;

import common.Interval;
import parser.State;
import prism.Evaluator;
import prism.JointAction;
import prism.PlayerInfo;
import prism.PrismException;
import strat.MDStrategy;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

/**
 * Simple explicit-state representation of a (multi-player) interval concurrent stochastic game (ICSG).
 */
public class ICSGSimple<Value> extends ModelExplicitWrapper<Value> implements NondetModelSimple<Value>, IntervalModelExplicit<Value>, ICSG<Value>
{
	/**
	 * The ICSG, stored as a CSGSimple over Intervals.
	 * Also stored in {@link ModelExplicitWrapper#model} as a ModelExplicit.
	 */
	protected CSGSimple<Interval<Value>> csg;

	// Constructors

	/**
	 * Constructor: empty ICSG.
	 */
	@SuppressWarnings("unchecked")
	public ICSGSimple()
	{
		this.csg = new CSGSimple<>();
		this.model = (ModelExplicit<Value>) csg;
		createDefaultEvaluatorForCSG();
	}

	/**
	 * Constructor: new ICSG with fixed number of states.
	 */
	/*@SuppressWarnings("unchecked")
	public ICSGSimple(int numStates)
	{
		this.csg = new CSGSimple<>(numStates);
		this.model = (ModelExplicit<Value>) csg;
		createDefaultEvaluatorForCSG();
	}*/

	/**
	 * Copy constructor.
	 */
	/*@SuppressWarnings("unchecked")
	public ICSGSimple(ICSGSimple<Value> icsg)
	{
		this.csg = new CSGSimple<>(icsg.csg);
		this.model = (ModelExplicit<Value>) csg;
		createDefaultEvaluatorForCSG();
	}*/

	/**
	 * Construct an ICSG from an existing one and a state index permutation,
	 * i.e. in which state index i becomes index permut[i].
	 * Pointer to states list is NOT copied (since now wrong).
	 * Note: have to build new Distributions from scratch anyway to do this,
	 * so may as well provide this functionality as a constructor.
	 */
	@SuppressWarnings("unchecked")
	public ICSGSimple(ICSGSimple<Value> icsg, int permut[])
	{
		this.csg = new CSGSimple<>(icsg.csg, permut);
		this.model = (ModelExplicit<Value>) csg;
		createDefaultEvaluatorForCSG();
	}

	/**
	 * Add a default (double interval) evaluator to the CSG
	 */
	@SuppressWarnings("unchecked")
	private void createDefaultEvaluatorForCSG()
	{
		((ICSGSimple<Double>) this).setIntervalEvaluator(Evaluator.forDoubleInterval());
	}

	// Mutators (for ModelSimple)

	@Override
	public void clearState(int s)
	{
		csg.clearState(s);
	}

	@Override
	public int addState()
	{
		return csg.addState();
	}

	@Override
	public void addStates(int numToAdd)
	{
		csg.addStates(numToAdd);
	}

	// Mutators (for IntervalModelExplicit)

	@Override
	public void setIntervalEvaluator(Evaluator<Interval<Value>> eval)
	{
		csg.setEvaluator(eval);
	}

	// Mutators (other)

	@Override
	public void setPlayerNames(List<String> playerNames)
	{
		csg.setPlayerNames(playerNames);
	}

	/**
	 * Set the list of all action labels
	 */
	public void setActions(List<Object> actions)
	{
		csg.setActions(actions);
	}

	/**
	 * Add a choice (uncertain distribution {@code udistr}) to state {@code s} (which must exist).
	 * Returns the index of the (newly added) distribution.
	 * Returns -1 in case of error.
	 */
	public int addChoice(int s, Distribution<Interval<Value>> udistr)
	{
		return csg.addChoice(s, udistr);
	}

	/**
	 * Add a choice (uncertain distribution {@code udistr}) labelled with {@code action} to state {@code s} (which must exist).
	 * Returns the index of the (newly added) distribution.
	 * Returns -1 in case of error.
	 */
	public int addActionLabelledChoice(int s, Distribution<Interval<Value>> udistr, Object action)
	{
		return csg.addActionLabelledChoice(s, udistr, action);
	}

	/**
	 * Add a choice (distribution {@code distr}) to state {@code s} (which must exist).
	 * Behaves the same as {@link MDPSimple#addActionLabelledChoice(int, Distribution, Object)},
	 * but {@code indexes} is an array storing the (1-indexed) index for the action
	 * performed by each player in this transition, and -1 indicates that the player idles.
	 * A representation of this is stored as a {@link JointAction} (accessible via e.g.
	 * {@link #getAction(int, int)}), whereas the array of indices can be accessed via
	 * {@link #getIndexes(int, int)}.
	 */
	public int addActionLabelledChoice(int s, Distribution<Interval<Value>> distr, int[] indexes)
	{
		return csg.addActionLabelledChoice(s, distr, indexes);
	}

	/**
	 * Set the action label for choice i in some state s.
	 */
	public void setAction(int s, int i, Object action)
	{
		csg.setAction(s, i, action);
	}

	/**
	 * Delimit the intervals for probabilities for the ith choice (distribution) for state s.
	 * i.e., trim the bounds of the intervals such that at least one
	 * possible distribution takes each of the extremal values.
	 * @param s The index of the state to delimit
	 * @param i The index of the choice to delimit
	 */
	public void delimit(int s, int i)
	{
		IntervalUtils.delimit(csg.trans.get(s).get(i), getEvaluator());
	}

	// Accessors (for NondetModel)

	@Override
	public int getNumChoices(int s)
	{
		return csg.getNumChoices(s);
	}

	@Override
	public Object getAction(int s, int i)
	{
		return csg.getAction(s, i);
	}

	@Override
	public boolean allSuccessorsInSet(int s, int i, BitSet set)
	{
		return csg.allSuccessorsInSet(s, i, set);
	}

	@Override
	public boolean someSuccessorsInSet(int s, int i, BitSet set)
	{
		return csg.someSuccessorsInSet(s, i, set);
	}

	@Override
	public Iterator<Integer> getSuccessorsIterator(final int s, final int i)
	{
		return csg.getSuccessorsIterator(s, i);
	}

	@Override
	public SuccessorsIterator getSuccessors(final int s, final int i)
	{
		return csg.getSuccessors(s, i);
	}

	@Override
	public int getNumTransitions(int s, int i)
	{
		return csg.getNumTransitions(s, i);
	}

	@Override
	public Model<Value> constructInducedModel(MDStrategy<Value> strat)
	{
		throw new UnsupportedOperationException("Not yet implemented");
	}

	// Accessors (for UCSG)

	@Override
	public void checkLowerBoundsArePositive() throws PrismException
	{
		Evaluator<Interval<Value>> eval = csg.getEvaluator();
		int numStates = getNumStates();
		for (int s = 0; s < numStates; s++) {
			int numChoices = getNumChoices(s);
			for (int j = 0; j < numChoices; j++) {
				Iterator<Map.Entry<Integer, Interval<Value>>> iter = getIntervalTransitionsIterator(s, j);
				while (iter.hasNext()) {
					Map.Entry<Integer, Interval<Value>> e = iter.next();
					// NB: we phrase the check as an operation on intervals, rather than
					// accessing the lower bound directly, to make use of the evaluator
					if (!eval.gt(e.getValue(), eval.zero())) {
						List<State> sl = getStatesList();
						String state = sl == null ? "" + s : sl.get(s).toString();
						throw new PrismException("Transition probability has lower bound of 0 in state " + state);
					}
				}
			}
		}
	}
	@Override
	public double mvMultUncSingle(int s, int k, double vect[], MinMax minMax)
	{
		@SuppressWarnings("unchecked")
		DoubleIntervalDistribution did = IntervalUtils.extractDoubleIntervalDistribution(((ICSG<Double>) this).getIntervalTransitionsIterator(s, k), getNumTransitions(s, k));
		return IDTMC.mvMultUncSingle(did, vect, minMax);
	}

	// Accessors (for PlayerInfoOwner)

	@Override
	public PlayerInfo getPlayerInfo()
	{
		return csg.getPlayerInfo();
	}

	// Accessors (for IntervalModel)

	@Override
	public Evaluator<Interval<Value>> getIntervalEvaluator()
	{
		return csg.getEvaluator();
	}

	@Override
	public CSG<Interval<Value>> getIntervalModel()
	{
		return csg;
	}

	// Accessors (for ICSG)

	@Override
	public Iterator<Map.Entry<Integer, Interval<Value>>> getIntervalTransitionsIterator(int s, int i)
	{
		return csg.getTransitionsIterator(s, i);
	}

	@Override
	public UncType getUncType()
	{
		return UncType.Adv;
	}

	/**
	 * Returns arg min/max_P { sum_j P(s,j)*vect[j] }
	 */
	public Iterator<Map.Entry<Integer, Double>> getDoubleTransitionsIterator(int s, int t, double val[]) {
		{
			// Collect transitions
			MinMax minMax = this.getUncType().toMinMax();
			List<Integer> indices = new ArrayList<>();
			List<Double> lowers = new ArrayList<>();
			List<Double> uppers = new ArrayList<>();
			Iterator<Map.Entry<Integer, Interval<Double>>> iter = ((ICSG<Double>) this).getIntervalModel().getTransitionsIterator(s, t);
			while (iter.hasNext()) {
				Map.Entry<Integer, Interval<Double>> e = iter.next();
				indices.add(e.getKey());
				lowers.add(e.getValue().getLower());
				uppers.add(e.getValue().getUpper());
			}
			int size = indices.size();

			// Trivial case: singleton interval [1.0,1.0]
			if (size == 1 && lowers.get(0) == 1.0 && uppers.get(0) == 1.0) {
				Map<Integer, Double> singleton = new HashMap<>();
				singleton.put(indices.get(0), 1.0);
				return singleton.entrySet().iterator();
			}

			// Sort indices by vect values
			List<Integer> order = new ArrayList<>();
			for (int i = 0; i < size; i++) order.add(i);
			if (minMax.isMaxUnc()) {
				order.sort((o1, o2) -> -Double.compare(val[indices.get(o1)], val[indices.get(o2)]));
			} else {
				order.sort((o1, o2) -> Double.compare(val[indices.get(o1)], val[indices.get(o2)]));
			}

			// Build the extreme distribution
			Map<Integer, Double> dist = new HashMap<>();
			double totP = 1.0;
			for (int i = 0; i < size; i++) {
				dist.put(indices.get(i), lowers.get(i));
				totP -= lowers.get(i);
			}
			for (int i = 0; i < size; i++) {
				int j = order.get(i);
				double delta = uppers.get(j) - lowers.get(j);
				double add = Math.min(delta, totP);
				dist.put(indices.get(j), dist.get(indices.get(j)) + add);
				totP -= add;
				if (totP <= 0) break;
			}
			return dist.entrySet().iterator();
		}
	}
}
