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
import prism.ModelType;
import prism.PlayerInfoOwner;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

public interface ICSG<Value> extends IMDP<Value>, PlayerInfoOwner
{
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

        /**
         * Convert this uncertainty type to a MinMax object, assuming Player 1 is maximising, and Player 2 is minimising.
         * @return A MinMax object representing the uncertainty type.
         */
        public MinMax toMinMax() {
			switch (this) {
				case Adv:
					return MinMax.minMin(false, true).setMinUnc(true);
				case Ctrl:
					return MinMax.minMin(false, true).setMinUnc(false);
				default:
					throw new IllegalArgumentException();
			}
        }

        /**
         * Convert this uncertainty type to a MinMax object.
         * @param min If nature is minimising (true) or maximising (false) value.
         * @return A MinMax object representing the uncertainty type.
         */
        public MinMax toMinMax(boolean min) {
			switch (this) {
				case Adv:
					return MinMax.minMin(min, !min).setMinUnc(!min);
				case Ctrl:
					return MinMax.minMin(min, !min).setMinUnc(min);
				default:
					throw new IllegalArgumentException();
			}
        }
    }

    public UncType getUncType();

//    /**
//     * Do a single row of matrix-vector multiplication for a specific choice k
//     * i.e. return min/max_P { sum_j P(s,k,j)*vect[j] }
//     * @param s State (row) index
//     * @param k Choice index
//     * @param vect Vector to multiply by
//     * @param minMax Min/max uncertainty (via isMinUnc/isMaxUnc)
//     */
//    public default double mvMultUncSingle(int s, int k, double vect[], MinMax minMax)
//    {
//        @SuppressWarnings("unchecked")
//        DoubleIntervalDistribution did = IntervalUtils.extractDoubleIntervalDistribution(((ICSG<Double>) this).getTransitionsIterator(s, k), getNumTransitions(s, k));
//        return IDTMC.mvMultUncSingle(did, vect, minMax);
//    }
//
//    /**
//     * Do a single row of matrix-vector multiplication for a specific choice k
//     * i.e. return min/max_P { rew(s) + rew_k(s) + sum_j P(s,k,j)*vect[j] }
//     * @param s State (row) index
//     * @param k Choice index
//     * @param vect Vector to multiply by
//     * @param mdpRewards The rewards (MDP rewards)
//     * @param minMax Min/max uncertainty (via isMinUnc/isMaxUnc)
//     */
//    public default double mvMultRewUncSingle(int s, int k, double vect[], MDPRewards<Double> mdpRewards, MinMax minMax)
//    {
//        double d = mdpRewards.getStateReward(s);
//        d += mdpRewards.getTransitionReward(s, k);
//        d += mvMultUncSingle(s, k, vect, minMax);
//        return d;
//    }

    @Override
    CSG<Interval<Value>> getIntervalModel();

}
