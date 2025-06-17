package explicit;

import common.Interval;
import explicit.rewards.CSGRewards;
import lpsolve.LpSolve;
import prism.PrismComponent;
import prism.PrismException;
import prism.PrismUtils;
import strat.CSGStrategy;

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


    /**
     * Returns arg min/max_P { sum_j P(s,j)*vect[j] }
     */
    @Override
    public Iterator<Map.Entry<Integer, Double>> getTransitionsIterator(CSG<?> csg, int s, int t, double val[], boolean min) {
        // Collect transitions
        MinMax minMax = ((ICSG<Double>) csg).getUncType().toMinMax(min);
        List<Integer> indices = new ArrayList<>();
        List<Double> lowers = new ArrayList<>();
        List<Double> uppers = new ArrayList<>();
        Iterator<Map.Entry<Integer, Interval<Double>>> iter = ((CSG<Interval<Double>>) csg).getTransitionsIterator(s, t);
        while (iter.hasNext()) {
            Map.Entry<Integer, Interval<Double>> e = iter.next();
            indices.add(e.getKey());
            lowers.add(((Number) e.getValue().getLower()).doubleValue());
            uppers.add(((Number) e.getValue().getUpper()).doubleValue());
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
