package explicit;

import common.Interval;
import prism.Evaluator;

import java.util.Iterator;
import java.util.Map;

public class ICSGSimple<Value> extends CSGSimple<Interval<Value>> implements ICSG<Value> {
    protected UncType uncType = UncType.Adv; // Default uncertainty type
    // Constructors

    /**
     * Constructor: empty ICSG.
     */
    public ICSGSimple() {
        super();
        createDefaultEvaluator();
    }

    public ICSGSimple(UncType uncType) {
        super();
        createDefaultEvaluator();
        this.uncType = uncType;
    }

    public UncType getUncType() {
        return uncType;
    }

    public void setUncType(UncType uncType) {
        this.uncType = uncType;
    }


    /**
     * Construct an ICSG from an existing one and a state index permutation,
     * i.e. in which state index i becomes index permut[i].
     * Pointer to states list is NOT copied (since now wrong).
     * Note: have to build new Distributions from scratch anyway to do this,
     * so may as well provide this functionality as a constructor.
     */
    public ICSGSimple(ICSGSimple<Value> icsg, int permut[]) {
        super(icsg, permut);
        createDefaultEvaluator();
    }




    /**
     * Delimit the intervals for probabilities for the ith choice (distribution) for state s.
     * i.e., trim the bounds of the intervals such that at least one
     * possible distribution takes each of the extremal values.
     * @param s The index of the state to delimit
     * @param i The index of the choice to delimit
     * @param evalChil An evaluator for the interval's child type (Value)
     */
    public void delimit(int s, int i, Evaluator<Value> evalChil)
    {
        IntervalUtils.delimit(trans.get(s).get(i), evalChil);
    }

    /**
     * Create a default Evaluator for double intervals (default for ModelExplicit is for doubles)
     */
    @SuppressWarnings("unchecked")
    private void createDefaultEvaluator()
    {
        ((ICSGSimple<Double>) this).setEvaluator(Evaluator.forDoubleInterval());
    }
}