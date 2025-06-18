package explicit;

import parser.ast.Coalition;
import prism.PrismComponent;
import prism.PrismException;
import strat.CSGStrategy;
import strat.ICSGStrategy;
import strat.Strategy;

import java.util.*;

public class ICSGModelChecker extends CSGModelChecker {
    protected CSGModelChecker mcCSG = null;

    /**
     * Create a new IMDPModelChecker, inherit basic state from parent (unless null).
     */
    public ICSGModelChecker(PrismComponent parent) throws PrismException
    {
        super(parent);
        mcCSG = new CSGModelChecker(this);
        mcCSG.inheritSettings(this);
    }



    @Override
    protected Strategy<?> getStrategy(CSG<?> icsg, List<List<List<Map<BitSet, Double>>>> lstrat, BitSet no, BitSet yes, BitSet inf, CSGStrategy.CSGStrategyType type) {
        return new ICSGStrategy((ICSG<Double>) icsg, lstrat, no, yes, inf, type);
    }

    public ModelCheckerResult computeReachProbs(ICSG<Double> icsg, BitSet target, MinMax minMax, int bound, Coalition coalition) throws PrismException {
        icsg.checkLowerBoundsArePositive();
        icsg.checkForDeadlocks(target);
        return super.computeReachProbs(icsg, target, minMax.isMin1(), minMax.isMin2(), bound, coalition);
    }

    public ModelCheckerResult computeUntilProbs(ICSG<Double> icsg, BitSet remain, BitSet target, int bound, MinMax minmax)
            throws PrismException {
        icsg.checkLowerBoundsArePositive();
        icsg.checkForDeadlocks(target);
        return super.computeUntilProbs(icsg, remain, target, bound, minmax.isMin1(), minmax.isMin2(), minmax.getCoalition());
    }

    public ModelCheckerResult computeUntilProbs(ICSG<Double> icsg, BitSet remain, BitSet target, MinMax minmax) throws PrismException {
        return computeUntilProbs(icsg, remain, target, maxIters, minmax);
    }

    public ModelCheckerResult computeBoundedUntilProbs(ICSG<Double> icsg, BitSet remain, BitSet target, int k, MinMax minmax)
            throws PrismException
    {
        return computeUntilProbs(icsg, remain, target, k, minmax.isMin1(), minmax.isMin2(), minmax.getCoalition());
    }

}
