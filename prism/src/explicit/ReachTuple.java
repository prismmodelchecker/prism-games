package explicit;

import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.BitSet;
import java.util.List;

/*3
 * Class represents the tuple of lower bounds of reachability probabilities
 * 
 * @author Aistis Simaitis
 * 
 */
public class ReachTuple {

	private static DecimalFormat df = new DecimalFormat("0.0000");
	private static boolean roundingOn = false;
	private static boolean roundDown = false && !roundingOn;
	private static double eps = 10e-6;

	private double[] values;

	/**
	 * Constructs a reachability tuple with 0 lower bounds
	 * 
	 * @param size
	 *            size of the tuple
	 */
	public ReachTuple(int size) {
		values = new double[size];
	}

	/**
	 * Creates a tuple with the specified lower bounds having the length of the
	 * array
	 * 
	 * @param bounds
	 */
	public ReachTuple(double[] bounds) {
		values = bounds.clone();
	}

	/**
	 * Creates a new tuple by scaling the given one by the given factor
	 * 
	 * @param tuple
	 * @param w
	 */
	public ReachTuple(ReachTuple tuple, double w) {
		values = tuple.getValues().clone();
		for (int i = 0; i < values.length; i++)
			values[i] *= w;
	}

	/**
	 * Creates the new tuple by adding two tuples together with the specified
	 * ratios
	 * 
	 * @param rt1
	 * @param rt2
	 */
	public ReachTuple(ReachTuple t1, double r1, ReachTuple t2, double r2) {
		if (t1.size() != t2.size())
			throw new IllegalArgumentException(
					"Cannot add two tuples of different sizes");
		values = new double[t1.size()];
		for (int i = 0; i < values.length; i++)
			values[i] = r1 * t1.getValue(i) + r2 * t2.getValue(i);
	}

	/**
	 * Creates a new tuple consisting of minimum/maximum elements of both tuples
	 * 
	 * @param t
	 * @param t1
	 * @param min
	 */
	public ReachTuple(ReachTuple t1, ReachTuple t2, boolean min) {
		if (t1.size() != t2.size())
			throw new IllegalArgumentException(
					"Cannot add two tuples of different sizes");
		values = new double[t1.size()];
		if (min)
			for (int i = 0; i < values.length; i++)
				values[i] = Math.min(t1.getValue(i), t2.getValue(i));
		else
			for (int i = 0; i < values.length; i++)
				values[i] = Math.max(t1.getValue(i), t2.getValue(i));
	}

	/**
	 * Addes the values of the given tuple to the current one
	 * 
	 * @param tuple
	 */
	public void add(ReachTuple tuple) {
		if (tuple.size() != size())
			throw new IllegalArgumentException(
					"Cannot add tuples of different sizes");

		for (int i = 0; i < values.length; i++)
			values[i] += tuple.getValue(i);
	}

	/**
	 * Returns values of the tuple
	 * 
	 * @return
	 */
	public double[] getValues() {
		return values;
	}

	/**
	 * Returns a bound for a specified target
	 * 
	 * @param i
	 * @return
	 */
	public double getValue(int i) {
		if (i < values.length)
			return values[i];
		else
			throw new IllegalArgumentException(
					"Index exceeds the number of targets");
	}

	/**
	 * Sets the specified value
	 */
	public void setValue(int i, double v) {
		if (i < values.length)
			values[i] = v;
		else
			throw new IllegalArgumentException(
					"Index exceeds the number of targets");
	}

	/**
	 * Checks whether the tuple is contained in the convex combinations of the
	 * two given tuples
	 * 
	 * @param t1
	 *            tuple 1
	 * @param t2
	 *            tuple 2
	 * @return
	 */
	public boolean isContained(ReachTuple t1, ReachTuple t2) {

		if (t1.size() != t2.size() || t1.size() != size())
			return false;

		if (roundDown) {
			roundDown();
			t1.roundDown();
			t2.roundDown();
		}

		double epsilon = roundingOn ? getEpsilon() : 0;

		double lower = 0;
		double upper = 1;
		double hi, lo;
		for (int i = 0; i < size(); i++) {
			if (t1.getValue(i) > t2.getValue(i)) {
				if (epsilon < getValue(i) - t1.getValue(i))
					return false;
				hi = 1;
				lo = ((getValue(i) - t2.getValue(i) - (getValue(i) == t2
						.getValue(i) ? 0 : epsilon)) / (t1.getValue(i) - t2
						.getValue(i)));
			} else if (t1.getValue(i) < t2.getValue(i)) {
				if (epsilon < getValue(i) - t2.getValue(i))
					return false;
				lo = 0;
				hi = ((getValue(i) - t2.getValue(i) - (getValue(i) == t2
						.getValue(i) ? 0 : epsilon)) / (t1.getValue(i) - t2
						.getValue(i)));
			} else {
				if (epsilon < getValue(i) - t2.getValue(i))
					return false;
				hi = 1;
				lo = 0;
			}
			if (lo > lower)
				lower = lo;
			if (hi < upper && hi >= 0)
				upper = hi;

			if (upper < lower)
				return false;
		}
		return true;
	}

	/**
	 * Checks whether the tuple contains a given tuple
	 * 
	 * @param t
	 *            the tuple
	 * @return
	 */
	public boolean contains(ReachTuple t) {
		if (t.size() != size())
			return false;

		if (roundDown) {
			roundDown();
			t.roundDown();
		}

		double epsilon = roundingOn ? getEpsilon() : 0;

		for (int i = 0; i < size(); i++)
			if (t.getValue(i) > getValue(i) + epsilon) {
				// System.out.println(this + " does not contain " + t);
				return false;
			}
		// System.out.println(this + " contains " + t);
		return true;
	}

	public void roundDown() {
		for (int i = 0; i < getValues().length; i++) {
			values[i] = Math.floor(values[i] / eps) * eps;
		}
	}

	// round probabilistically
	private double getEpsilon() {
		// if(1==1) return eps;
		if (Math.random() > -1)
			return eps;
		else
			return 0.0;
	}

	public int size() {
		return getValues().length;
	}

	public boolean isZero() {
		for (int i = 0; i < getValues().length; i++)
			if (getValue(i) != 0)
				return false;
		return true;
	}

	public void perturbateUp() {
		double sum = 0;
		for (int i = 0; i < values.length; i++)
			sum += values[i];
		if(sum == 1) return;

		for (int i = 0; i < values.length; i++)
			if(values[i] == 0)
				values[i] = 0.1;
			//			if (values[i] < 1 - 10e-1)
//				values[i] += 10e-1;
//			else
//				values[i] = 1;
	}

	public void perturbateDown() {
//		double sum = 0;
//		for (int i = 0; i < values.length; i++)
//			sum += values[i];
//		if(sum == 1) return;
		for (int i = 0; i < values.length; i++)
			if (values[i]!=1 && values[i] > 10e-2)
				values[i] -= 10e-2;
//			else
//				values[i] = 0;
	}

	@Override
	public boolean equals(Object o) {
		if (!(o instanceof ReachTuple))
			return false;
		return Arrays.equals(getValues(), ((ReachTuple) o).getValues());
	}

	@Override
	public int hashCode() {
		return Arrays.hashCode(getValues());
	}

	@Override
	public String toString() {
		String ret = "[ ";
		for (int i = 0; i < getValues().length; i++)
			ret += df.format(getValue(i)) + " ";
		ret += "]";
		return ret;
		// return Arrays.toString(getValues());
	}

	private static void testIsContained() {
		ReachTuple t1, t2, t3, t4, t5, t6;

		t1 = new ReachTuple(new double[] { 0.0, 0.8, 0.0 });
		t2 = new ReachTuple(new double[] { 0.508, 0.0, 0.48000000000000004 });
		;
		t3 = new ReachTuple(new double[] { 0.4, 0.16000000000000003, 0.4 });
		;
		t4 = new ReachTuple(new double[] { 0.49600000000000005,
				0.006400000000000002, 0.49600000000000005 });
		;
		t5 = new ReachTuple(new double[] { 0.48000000000000004,
				0.03200000000000001, 0.48000000000000004 });
		;
		t6 = new ReachTuple(new double[] { 0.54, 0.0, 0.4 });

		System.out.println(t3.isContained(t1, t4));
	}

	public static void testContains() {
		ReachTuple t1, t2;
		t1 = new ReachTuple(new double[] { 0.499968, 5.1200000000000025E-5,
				0.499968 });
		t2 = new ReachTuple(new double[] { 0.5000128, 0.0, 0.499968 });

		System.out.println(t1.contains(t2) + " " + t2.contains(t1));
	}

	public static void main(String[] args) {
		testIsContained();
		testContains();

		ReachTuple x = new ReachTuple(new double[] { 0.1, 0.7, 0.1 });
		ReachTuple y = new ReachTuple(new double[] { 0.5, 0.2, 0.9 });
		ReachTuple z = new ReachTuple(new double[] { 0.3, 0.3, 0.5 });

		z.isContained(x, y);

		double[] arr1 = { 1.0, 0.5 };
		System.out.println(Arrays.toString(arr1));

		ReachTuple rt = new ReachTuple(arr1);

		rt.setValue(0, 0.333);
		System.out.println(Arrays.toString(arr1));
		System.out.println(Arrays.toString(rt.getValues()));

	}
}
