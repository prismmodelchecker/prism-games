//==============================================================================
//	
//	Copyright (c) 2026-
//	Authors:
//	* Gabriel Santos <gabriel.santos@cs.ox.ac.uk> (University of Oxford)
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

package soplex;

/**
 * Minimal JNI wrapper for the SoPlex LP solver (floating point).
 * <p>
 * One instance holds one SoPlex object that can be reused for many LPs via {@link #clear()}.
 * Not thread-safe: use one instance per thread. Call {@link #free()} (or close) when done.
 * <p>
 * Build the model by adding columns (variables) and rows (constraints); rows refer to columns by index.
 * <p>
 * Defaults are tuned for many tiny LPs: SoPlex's presolver (simplifier) and timer are switched off
 * ({@link #setTimeLimit(double)} switches the timer back on). Scaling is SoPlex's default (bi-equilibrium);
 * see {@link #setScaler(int)}.
 */
public class SoPlex implements AutoCloseable
{
	/** SoPlex status codes (SPxSolverBase::Status) */
	public static final int ERROR = -15, ABORT_TIME = -7, ABORT_ITER = -6, SINGULAR = -4, NO_PROBLEM = -3, UNKNOWN = 0, OPTIMAL = 1, UNBOUNDED = 2,
			INFEASIBLE = 3, INF_OR_UNBD = 4, OPTIMAL_UNSCALED_VIOLATIONS = 5;

	/** Scaling methods (SoPlex::SCALER values) */
	public static final int SCALER_OFF = 0, SCALER_UNIEQUI = 1, SCALER_BIEQUI = 2, SCALER_GEO1 = 3, SCALER_GEO8 = 4, SCALER_LEASTSQ = 5, SCALER_GEOEQUI = 6;

	/** Value used for +/- infinite bounds */
	public static final double INFINITY = 1e100;

	/** Codes for {@link #setRealParam(int, double)} (SoPlex::RealParam): primal and dual feasibility tolerances (default 1e-6) */
	public static final int FEASTOL = 0, OPTTOL = 1;

	private static boolean loaded = false;

	private long ptr;

	/** Load the native library (soplexj). Called automatically by the constructor. */
	public static synchronized void loadLibrary()
	{
		if (!loaded) {
			System.loadLibrary("soplexj");
			loaded = true;
		}
	}

	public SoPlex()
	{
		loadLibrary();
		ptr = create();
		if (ptr == 0)
			throw new IllegalStateException("Could not create SoPlex object");
	}

	/** Remove all rows and columns (keeps parameters). */
	public void clear()
	{
		clear(ptr);
	}

	public void setMaximise(boolean maximise)
	{
		setObjSense(ptr, maximise);
	}

	/** Verbosity 0 (errors only, default here) to 5. */
	public void setVerbosity(int level)
	{
		setVerbosity(ptr, level);
	}

	public void setTimeLimit(double seconds)
	{
		setTimeLimit(ptr, seconds);
	}

	/** Set the scaling method (one of the SCALER_* constants; SoPlex's default is SCALER_BIEQUI). */
	public void setScaler(int scaler)
	{
		setScaler(ptr, scaler);
	}

	/** Scaler constant for a (case-insensitive) name: off, uniequi, biequi, geo1, geo8, leastsq, geoequi. */
	public static int scalerFromName(String name)
	{
		switch (name.toLowerCase()) {
		case "off": return SCALER_OFF;
		case "uniequi": return SCALER_UNIEQUI;
		case "biequi": return SCALER_BIEQUI;
		case "geo1": return SCALER_GEO1;
		case "geo8": return SCALER_GEO8;
		case "leastsq": return SCALER_LEASTSQ;
		case "geoequi": return SCALER_GEOEQUI;
		default: throw new IllegalArgumentException("Unknown SoPlex scaler \"" + name + "\"");
		}
	}

	/** Set a SoPlex integer parameter by its code (SoPlex::IntParam). */
	public void setIntParam(int code, int value)
	{
		setIntParam(ptr, code, value);
	}

	/** Set a SoPlex real parameter by its code (SoPlex::RealParam). */
	public void setRealParam(int code, double value)
	{
		setRealParam(ptr, code, value);
	}

	/** Set a SoPlex boolean parameter by its code (SoPlex::BoolParam). */
	public void setBoolParam(int code, boolean value)
	{
		setBoolParam(ptr, code, value);
	}

	/** Add a column (variable) with no nonzeros in existing rows; returns its index. */
	public int addCol(double obj, double lb, double ub)
	{
		return addCol(ptr, obj, lb, ub);
	}

	/** Add a row lhs <= sum_k vals[k] * x[idx[k]] <= rhs; returns its index. */
	public int addRow(int[] idx, double[] vals, double lhs, double rhs)
	{
		if (idx.length != vals.length)
			throw new IllegalArgumentException("index/value length mismatch");
		return addRow(ptr, idx, vals, lhs, rhs);
	}

	/**
	 * Append columns (empty column vectors; add rows afterwards) in one call.
	 * @return index of the first new column
	 */
	public int addCols(double[] obj, double[] lb, double[] ub)
	{
		if (obj.length != lb.length || obj.length != ub.length)
			throw new IllegalArgumentException("column array length mismatch");
		return addCols(ptr, obj, lb, ub);
	}

	/** Replace the objective (one coefficient per column). The basis is kept, so the next optimize() warm-starts. */
	public void changeObj(double[] obj)
	{
		if (obj.length != numCols(ptr))
			throw new IllegalArgumentException("objective length must equal the number of columns");
		changeObj(ptr, obj);
	}

	/** Change the bounds lhs <= row i <= rhs. */
	public void changeRowRange(int i, double lhs, double rhs)
	{
		changeRowRange(ptr, i, lhs, rhs);
	}

	/** Change the bounds lb <= column i <= ub. */
	public void changeColBounds(int i, double lb, double ub)
	{
		changeColBounds(ptr, i, lb, ub);
	}

	public int optimize()
	{
		return optimize(ptr);
	}

	public int getStatus()
	{
		return getStatus(ptr);
	}

	public double getObjValue()
	{
		return getObjValue(ptr);
	}

	public int numRows()
	{
		return numRows(ptr);
	}

	public int numCols()
	{
		return numCols(ptr);
	}

	public double[] getPrimal()
	{
		double[] x = new double[numCols(ptr)];
		getPrimal(ptr, x);
		return x;
	}

	public double[] getDual()
	{
		double[] y = new double[numRows(ptr)];
		getDual(ptr, y);
		return y;
	}

	public double[] getRedCost()
	{
		double[] r = new double[numCols(ptr)];
		getRedCost(ptr, r);
		return r;
	}

	public double getSolvingTime()
	{
		return getSolvingTime(ptr);
	}

	/**
	 * Solve a zero-sum matrix game A (m x n, row-major, row player maximises) in one native call,
	 * reusing this object. If {@code min} is false, the LP is for the row player
	 * (max v s.t. sum_i A[i][j] x_i >= v for all j, x a distribution) and {@code strategyOut} (length m) receives x;
	 * if true, for the column player (min w s.t. sum_j A[i][j] y_j <= w for all i) and {@code strategyOut} (length n) receives y.
	 * The value v is bounded to [vlb, vub] (use -INFINITY/INFINITY for none).
	 *
	 * @return the game value, or NaN if not solved to optimality (see {@link #getStatus()})
	 */
	public double matrixGame(double[] A, int m, int n, boolean min, double vlb, double vub, double[] strategyOut)
	{
		if (A.length != m * n)
			throw new IllegalArgumentException("matrix size mismatch");
		if (strategyOut != null && strategyOut.length != (min ? n : m))
			throw new IllegalArgumentException("strategy array has wrong length");
		return matrixGame(ptr, A, m, n, min, vlb, vub, strategyOut);
	}

	/** Release the native object. */
	public void free()
	{
		if (ptr != 0) {
			free(ptr);
			ptr = 0;
		}
	}

	@Override
	public void close()
	{
		free();
	}

	// native methods
	private static native long create();
	private static native void free(long p);
	private static native void clear(long p);
	private static native void setObjSense(long p, boolean max);
	private static native void setVerbosity(long p, int level);
	private static native void setTimeLimit(long p, double seconds);
	private static native void setScaler(long p, int scaler);
	private static native void setIntParam(long p, int code, int value);
	private static native void setRealParam(long p, int code, double value);
	private static native void setBoolParam(long p, int code, boolean value);
	private static native int addCol(long p, double obj, double lb, double ub);
	private static native int addRow(long p, int[] idx, double[] vals, double lhs, double rhs);
	private static native int addCols(long p, double[] obj, double[] lb, double[] ub);
	private static native void changeObj(long p, double[] obj);
	private static native void changeRowRange(long p, int i, double lhs, double rhs);
	private static native void changeColBounds(long p, int i, double lb, double ub);
	private static native void getRedCost(long p, double[] r);
	private static native int optimize(long p);
	private static native int getStatus(long p);
	private static native double getObjValue(long p);
	private static native int numRows(long p);
	private static native int numCols(long p);
	private static native void getPrimal(long p, double[] x);
	private static native void getDual(long p, double[] y);
	private static native double getSolvingTime(long p);
	private static native double matrixGame(long p, double[] A, int m, int n, boolean min, double vlb, double vub, double[] strategyOut);
}
