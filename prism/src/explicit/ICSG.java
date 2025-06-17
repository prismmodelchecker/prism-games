package explicit;

import common.Interval;
import common.IterableStateSet;
import explicit.rewards.MDPRewards;
import parser.State;
import prism.Evaluator;
import prism.ModelType;
import prism.PrismException;

import java.util.*;

public interface ICSG<Value> extends CSG<Interval<Value>> {
    // Accessors (for Model) - default implementations

    @Override
    public default ModelType getModelType()
    {
        return ModelType.ICSG;
    }

    public enum UncType {
        Adv("Adversarial"),
        Ctrl("Controllable");

        private final String fullName;

        UncType(final String fullName) {
            this.fullName = fullName;
        }

        /**
         * Get the full name, in words, of this model type.
         */
        public String fullName()
        {
            return fullName;
        }

        /**
         * Get the PRISM keyword for this model type.
         */
        public String keyword()
        {
            return this.name().toLowerCase();
        }

        public MinMax toMinMax(boolean min) {
            return switch (this) {
                case Adv -> MinMax.minMin(min, !min).setMinUnc(!min);
                case Ctrl -> MinMax.minMin(min, !min).setMinUnc(min);
            };
        }
    }

    public UncType getUncType();
    public void setUncType(UncType type);

    /**
     * Checks that transition probability interval lower bounds are positive
     * and throws an exception if any are not.
     */
    public default void checkLowerBoundsArePositive() throws PrismException
    {
        Evaluator<Interval<Value>> eval = getEvaluator();
        int numStates = getNumStates();
        for (int s = 0; s < numStates; s++) {
            int numChoices = getNumChoices(s);
            for (int j = 0; j < numChoices; j++) {
                Iterator<Map.Entry<Integer, Interval<Value>>> iter = getTransitionsIterator(s, j);
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

    /**
     * Do a single row of matrix-vector multiplication for a specific choice k
     * i.e. return min/max_P { sum_j P(s,k,j)*vect[j] }
     * @param s State (row) index
     * @param k Choice index
     * @param vect Vector to multiply by
     * @param minMax Min/max uncertainty (via isMinUnc/isMaxUnc)
     */
    public default double mvMultUncSingle(int s, int k, double vect[], MinMax minMax)
    {
        @SuppressWarnings("unchecked")
        DoubleIntervalDistribution did = IntervalUtils.extractDoubleIntervalDistribution(((ICSG<Double>) this).getTransitionsIterator(s, k), getNumTransitions(s, k));
        return IDTMC.mvMultUncSingle(did, vect, minMax);
    }

    /**
     * Do a single row of matrix-vector multiplication for a specific choice k
     * i.e. return min/max_P { rew(s) + rew_k(s) + sum_j P(s,k,j)*vect[j] }
     * @param s State (row) index
     * @param k Choice index
     * @param vect Vector to multiply by
     * @param mdpRewards The rewards (MDP rewards)
     * @param minMax Min/max uncertainty (via isMinUnc/isMaxUnc)
     */
    public default double mvMultRewUncSingle(int s, int k, double vect[], MDPRewards<Double> mdpRewards, MinMax minMax)
    {
        double d = mdpRewards.getStateReward(s);
        d += mdpRewards.getTransitionReward(s, k);
        d += mvMultUncSingle(s, k, vect, minMax);
        return d;
    }
}
