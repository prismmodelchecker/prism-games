package explicit;

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


}
