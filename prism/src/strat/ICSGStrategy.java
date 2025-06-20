package strat;

import explicit.ICSG;
import common.Interval;

import java.util.BitSet;
import java.util.List;
import java.util.Map;

public class ICSGStrategy extends CSGStrategy<Interval<Double>> implements Strategy<Interval<Double>> {
    public ICSGStrategy(ICSG<Double> icsg, List<List<List<Map<BitSet, Double>>>> csgchoices, BitSet no, BitSet yes, BitSet inf, CSGStrategyType type) {
        super(icsg.getIntervalModel(), csgchoices, no, yes, inf, type);
    }
}