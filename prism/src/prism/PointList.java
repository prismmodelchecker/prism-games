//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Clemens Wiltsche <clemens.wiltsche@cs.ox.ac.uk> (University of Oxford)
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

package prism;

import java.util.Collections;
import java.util.Comparator;
import java.util.ArrayList;
import java.util.List;
import java.math.BigInteger;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.Arrays;

import org.apache.commons.math3.fraction.BigFraction;

import parma_polyhedra_library.Variables_Set;
import parma_polyhedra_library.Relation_Symbol;
import parma_polyhedra_library.Polyhedron;
import parma_polyhedra_library.Generator;
import parma_polyhedra_library.Variable;
import parma_polyhedra_library.Generator_Type;
import parma_polyhedra_library.Generator_System;
import parma_polyhedra_library.Coefficient;
import parma_polyhedra_library.Linear_Expression_Coefficient;
import parma_polyhedra_library.Linear_Expression;
import parma_polyhedra_library.C_Polyhedron;
import parma_polyhedra_library.Constraint;
import parma_polyhedra_library.Linear_Expression_Variable;
import parma_polyhedra_library.Linear_Expression_Times;
import parma_polyhedra_library.Partial_Function;

import parser.Values;
import parser.ast.Expression;
import parser.ast.RelOp;
import parser.ast.ExpressionVar;
import parser.type.TypeVoid;
import parser.ast.ExpressionConstant;
import parser.ast.ExpressionProb;
import parser.ast.ExpressionReward;
import explicit.MultiParameters;
import explicit.PPLSupport;
import explicit.Pareto;
import userinterface.graph.Slice;
import explicit.SMGModelChecker;

public class PointList extends TileList
{
	private static final boolean problem = false; // for debugging
	
	private static final long denominator = 10000000000L; // accuracy 10^-10 for slicing
	private static final double acc = (double) 1.0/((double) denominator);
	
	protected ArrayList<Point> boundary_points = null; // boundary points
	private ArrayList<List<Point>> poly_points = null; // union of convex polyhedra, each stored as lists of points

	// cutoff values of the polyhedra (not displayed beyond this)
	protected double[] cutoff_min;
	protected double[] cutoff_max;
	private static final double f1 = 5.0; // factor for cutting off Pareto sets
	private static final double f2 = 4.0; // factor for ignoring border points
	private static final double f3 = 3.0; // factor for display range

	// data storage
	private Pareto pareto = null; // store full-dimension Pareto set
	private List<Expression> expressions; // one expression per dimension
	private List<Double> bounds; // one bound per dimension
	protected int full_dim = 0; // dimension of full polyhedra

	// mapped dimension index of X and Y axes
	private int indexX;
	private int indexY;

	/**
	 * Accessors.
	 **/

	public double getCutoffMin(int i)
	{
		return cutoff_min[i] > 0.0 ? cutoff_min[i] / f3 : cutoff_min[i] * f3;
	}

	public double getCutoffMax(int i)
	{
		return cutoff_max[i] > 0.0 ? cutoff_max[i] * f3 : cutoff_max[i] / f3;
	}

	public double getMidpoint(int i)
	{
		if (bounds != null)
			return i < bounds.size() ? bounds.get(i) : 0.0;
		else
			return (cutoff_max[i] + cutoff_min[i]) / 2.0;
	}

	public String getXLabel()
	{
		return getDim2Var(expressions.get(indexX), false);
	}

	public String getYLabel()
	{
		return getDim2Var(expressions.get(indexY), false);
	}

	public String getLabel(int i)
	{
		return getDim2Var(expressions.get(i), false);
	}

	public int getFullDim()
	{
		return full_dim;
	}

	public List<List<Point>> getPolyPoints()
	{
		return poly_points;
	}

	/**
	 * Constructors.
	 **/

	public PointList()
	{
		super(2, null, 0);
	}
	
	/**
	 * A new list of points to display a Pareto set.
	 *
	 * @param pareto The Pareto set to be displayed: a union of convex sets
	 * @param expressions One expression per dimension of the Pareto set
	 * @param bounds The bounds, used as default slice
	 **/
	public PointList(Pareto pareto, List<Expression> expressions, List<Double> bounds) throws PrismException
	{
		// initialise as TileList to be registered correctly in the static TileList.storedTileLists for displaying
		super(2, null, 0);

		// Check that Pareto set is not empty
		if (pareto == null || pareto.size() == 0 || pareto.get(0) == null)
			throw new PrismException("Pareto set cannot be empty.");

		// here the dimensions that have identical bounds are folded. Also sets expressions and Pareto set.
		foldAndSet(pareto, expressions);

		// calculate the maximum extent of the Pareto set in each dimension, and cut Pareto set
		calculateCutoffs(this.pareto);
		cutPareto(this.pareto);

		// initialise storage for boundary points and polyhedra
		this.boundary_points = new ArrayList<Point>();
		this.poly_points = new ArrayList<List<Point>>();
		this.bounds = bounds;
	}

	/**
	 * Given a list of expressions, return a list of labels,
	 * so that all dimensions with the same label (excluding null)
	 * map to the same dimension.
	 *
	 **/
	private List<String> getDims2Vars(List<Expression> expressions)
	{
		List<String> dim2var = new ArrayList<String>();
		for (Expression expression : expressions)
			dim2var.add(getDim2Var(expression, true));
		return dim2var;
	}

	/**
	 * If padWithNull is false, then the the names of the expressions
	 * that are not bounded by a variable or constant are filled in.
	 **/
	private String getDim2Var(Expression expression, boolean padWithNull)
	{
		Expression var = null;
		if (expression instanceof ExpressionProb) {
		        ExpressionProb ep = (ExpressionProb) expression;
			if(ep.getExpression() instanceof ExpressionReward)
			        var = ((ExpressionReward) ep.getExpression()).getReward();
			else
			        var = ((ExpressionProb) expression).getProb();
		} else if (expression instanceof ExpressionReward) {
			var = ((ExpressionReward) expression).getReward();
		}

		if (var instanceof ExpressionVar) {
			String v_n = ((ExpressionVar) var).getName();
			int v_i = ((ExpressionVar) var).getIndex();
			return (v_n == null ? String.valueOf(v_i) : v_n);
		} else if (var instanceof ExpressionConstant) {
			return ((ExpressionConstant) var).getName();
		} else {
			return padWithNull ? null : expression.toString();
		}

	}

	
	/**
	 * Folds the Pareto set and expressions, and puts them into the class.
	 **/
	private void foldAndSet(Pareto pareto, List<Expression> expressions)
	{
		Pareto folded_pareto = new Pareto(pareto); // deep copy, will be modified
		List<Expression> folded_expressions = new ArrayList<Expression>(); // will be filled later

		// fold all dimensions together that have the same bound, because this means that they should be intersected
		List<String> dim2var = getDims2Vars(expressions); // associate a variable ID to expressions with variables

		// apply equality constraints to only have relevant information left
		Map<String, Integer> folded = new HashMap<String, Integer>(); // already a folding dimension assigned
		Variables_Set dimensions_to_remove = new Variables_Set(); // dimensions that will be removed
		int i = 0;
		for (String var : dim2var) {
			if (var != null) { // all dimensions that have the same var are now constrained together
				if (!folded.containsKey(var)) { // no folding dimension assigned yet
					folded.put(var, i); // assign this as the folding dimension, no folding required yet at this point
					folded_expressions.add(expressions.get(i));
				} else {
					// equality constraint
					Constraint con = new Constraint(new Linear_Expression_Variable(new Variable(i)), Relation_Symbol.EQUAL, new Linear_Expression_Variable(
							new Variable(folded.get(var))));
					// apply to every polyhedron
					for (Polyhedron p : folded_pareto.getSets())
						p.add_constraint(con);

					// this dimension will moreover be removed later
					dimensions_to_remove.add(new Variable(i));
				}
			} else {
				folded_expressions.add(expressions.get(i));
			}
			i++;
		}

		// now need to project away the irrelevant dimensions
		for (Polyhedron p : folded_pareto.getSets())
			p.remove_space_dimensions(dimensions_to_remove);

		// if dimension reduced to less than two, add space dimensions
		full_dim = folded_expressions.size();
		if (full_dim < 2) {
			for (Polyhedron p : folded_pareto.getSets())
				p.add_space_dimensions_and_embed(2 - full_dim);
			for (i = full_dim; i < 2; i++)
				folded_expressions.add(new ExpressionConstant("-", TypeVoid.getInstance()));
		}

		// now set Pareto set and expressions
		this.pareto = folded_pareto;
		this.expressions = folded_expressions;
		this.full_dim = folded_expressions.size();
	}

	/**
	 * Slice the Pareto set according to the dimension-value pairs specified in slice,
	 * and fill in the point list.
	 **/
	public void updateParetoSet(Slice slice) throws PrismException
	{
		// check if slice specifies enough dimensions to be sliced or projected
		if (slice.fullSize() != full_dim - 2)
			throw new PrismException(String.format("Need to specify at least %d dimensions for slicing or projecting. Currently only %d specified",
					(full_dim - 2), slice.fullSize()));

		// deep copy of Pareto set, since polyhedra need to be preserved for later slicing
		Pareto Ps = new Pareto(pareto);

		// SLICING
		Variables_Set sliced_vars = new Variables_Set();
		for (Integer i : slice.keySet()) {
			// dimensions to be sliced away
			sliced_vars.add(new Variable(i));
			// establish a new constraint for each dimension-value pair in the slice
			Linear_Expression lhs = new Linear_Expression_Times(new Coefficient(denominator), new Variable(i));
			BigInteger numi = BigInteger.valueOf((long) (slice.get(i) * (double) (denominator)));
			Linear_Expression rhs = new Linear_Expression_Coefficient(new Coefficient(numi));
			Constraint con = new Constraint(lhs, Relation_Symbol.EQUAL, rhs); // slice by constraining with equality
			for (Polyhedron p : Ps.getSets())
				p.add_constraint(con); // add slicing constraint to each polyhedron
		}
		// complete slicing by projecting onto the unaffected space dimensions, i.e. take away the zeroed dimensions
		for (Polyhedron p : Ps.getSets())
			p.remove_space_dimensions(sliced_vars);

		// PROJECTING - need to be careful to match indices after slicing
		Variables_Set project = new Variables_Set();
		int project_index = 0;
		for (int i = 0; i < full_dim; i++) {
			if (slice.getProject().contains(i)) // dimension projected
				project.add(new Variable(project_index));
			if (!slice.containsKey(i)) // dimension has already been sliced
				project_index++;
		}
		// complete projecting by projecting onto the unaffected dimensions
		for (Polyhedron p : Ps.getSets())
			p.remove_space_dimensions(project);

		// now fill in the points to be drawn
		fillPointList(slice, Ps);
	}

	/**
	 * Fill the list of points with the data from the Polyhedra, after they have been sliced and projected.
	 **/
	void fillPointList(Slice slice, Pareto Ps) throws PrismException
	{
		// clear the currently stored points
		boundary_points.clear();
		poly_points.clear();

		// set point list properties and add to be displayed
		boolean doneWithX = false;
		for (int j = 0; j < full_dim; j++) { // go through all dimensions
			if (!slice.removed(j)) { // dimension has not been removed
				if (!doneWithX) { // X not yet assigned and dimension not removed
					indexX = j;
					doneWithX = true;
				} else { // X already assigned and dimension not removed
					indexY = j;
					break; // done with Y
				}
			}
		}
		// if swapping dimensions, reverse order of the list
		if (slice.swapdimensions) {
			int tempX = indexX;
			indexX = indexY;
			indexY = tempX;
		}

		// now add two-dimensional points for displaying the sets
		for (int k = 0; k < Ps.size(); k++) {
			ArrayList<Point> poly_point = new ArrayList<Point>();
			for (Generator g : Ps.get(k).minimized_generators()) {
				if (g.type() == Generator_Type.RAY || g.type() == Generator_Type.LINE)
					continue; // ignore rays and lines

				// evaluate the numerator
				Map<Integer, BigInteger> num = PPLSupport.getCoefficients(g.linear_expression());
				// evaluate the denominator
				BigInteger den = g.divisor().getBigInteger();

				// get the point represented by the generator
				double[] x = new double[2]; // default initialisation to zero
				for (Entry<Integer, BigInteger> jval : num.entrySet()) { // fill x with values for dimensions specified in num
					int offset = slice.swapdimensions ? 1 - jval.getKey() : jval.getKey(); // take care of swapped dimensions
					x[offset] = (new BigFraction(jval.getValue(), den)).doubleValue();
				}
				poly_point.add(new Point(x));
			}
			poly_points.add(poly_point);
		}

		// sort points for drawing as polyhedra
		sortPoints();
	}

	public static void addStoredPointList(String name, PointList pointList)
	{
		// some edits in pointList itself
		pointList.name = name;

		// add point list statically to TileList - also fill up the stored formulae to keep indices of TileLists intact
		addStoredFormula(null);
		addStoredFormulaX(null);
		addStoredFormulaY(null);
		storedTileLists.add(pointList);
	}

	@Override
	public String toString()
	{
		// TODO
		return "Point List";
	}

	/**
	 * Calculates some reasonable maximum extend of the Pareto set in each dimension.
	 * The result is stored in the global cutoff_min and cutoff_max fields.
	 **/
	private void calculateCutoffs(Pareto pareto)
	{
		List<Polyhedron> flats = new ArrayList<Polyhedron>();
		List<Polyhedron> contributors = new ArrayList<Polyhedron>();

		go_through_polyhedra: for (Polyhedron p : pareto.getSets()) {
			for (Generator g : p.generators()) { // note, rays are ignored
				if (g.type() == Generator_Type.LINE) {
					flats.add(p); // the polyhedron is a halfspace
					continue go_through_polyhedra;
				}
			}
			// fall through if not halfspace
			contributors.add(p);
		}
		// now mutually intersect the halfspaces
		for (int i1 = 0; i1 < flats.size(); i1++) {
			Polyhedron p1 = flats.get(i1);
			for (int i2 = i1 + 1; i2 < flats.size(); i2++) {
				Polyhedron p2 = new C_Polyhedron(flats.get(i2).generators()); // deep copy
				// now intersect
				p2.intersection_assign(p1);
				// and add to contributors
				contributors.add(p2);
			}
		}
		// finally evaluate the ultimate reward in each dimension
		cutoff_min = new double[full_dim];
		cutoff_max = new double[full_dim];
		for (int i = 0; i < full_dim; i++) {
			cutoff_max[i] = Double.NEGATIVE_INFINITY;
			cutoff_min[i] = Double.POSITIVE_INFINITY;
		}

		for (Polyhedron p : contributors) {
			for (Generator g : p.generators()) {
				if (g.type() == Generator_Type.POINT) {
					Map<Variable, BigInteger> num = new HashMap<Variable, BigInteger>();
					PPLSupport.getCoefficientsFromLinearExpression(g.linear_expression(), false, BigInteger.ONE, num);
					BigInteger den = g.divisor().getBigInteger();
					for (int i = 0; i < full_dim; i++) {
						boolean foundvalue = false;
						look_for_variable: for (Variable j : num.keySet()) {
							if (j != null && i == j.id()) {
								double val_i = new BigFraction(num.get(j), den).doubleValue();
								foundvalue = true;
								if (cutoff_max[i] < val_i)
									cutoff_max[i] = val_i; // maximise
								if (cutoff_min[i] > val_i)
									cutoff_min[i] = val_i; // minimise
								break look_for_variable;
							}
						}
						if (!foundvalue) { // if no value registered, then this automatically corresponds to zero
							if (cutoff_max[i] < 0.0)
								cutoff_max[i] = 0.0; // maximise
							if (cutoff_min[i] > 0.0)
								cutoff_min[i] = 0.0; // minimise
						}
					}
				}
			}
		}

		// if some dimension is zeros, make it visible by adding some default spread
		for (int i = 0; i < full_dim; i++) {
			if (cutoff_min[i] == Double.POSITIVE_INFINITY) {
				cutoff_min[i] = -10.0;
			}
			if (cutoff_max[i] == Double.NEGATIVE_INFINITY) {
				cutoff_max[i] = 10.0;
			}
			double min_diff = 2.0;
			double diff = cutoff_max[i] - cutoff_min[i];
			if (diff < min_diff) {
				cutoff_max[i] += (min_diff - diff) / 2.0;
				cutoff_min[i] -= (min_diff - diff) / 2.0;
			}
		}
	}

	/**
	 * Cuts the Pareto set in by the values given in cutoff_min and cutoff_max for the respective dimensions.
	 * Modifies the Pareto set.
	 **/
	private void cutPareto(Pareto Vs)
	{
		// for displaying the Pareto set (only for initial state)
		for (Polyhedron q : Vs.getSets()) {
			for (int i = 0; i < full_dim; i++) {
				Linear_Expression lhs = new Linear_Expression_Variable(new Variable(i));
				Linear_Expression rhs1 = new Linear_Expression_Coefficient(new Coefficient((long) Math.ceil(cutoff_max[i] > 0.0 ? cutoff_max[i] * f1
						: cutoff_max[i] / f1)));
				Linear_Expression rhs2 = new Linear_Expression_Coefficient(new Coefficient((long) Math.floor(cutoff_min[i] > 0.0 ? cutoff_min[i] / f1
						: cutoff_min[i] * f1)));
				Constraint con1 = new Constraint(lhs, Relation_Symbol.LESS_OR_EQUAL, rhs1);
				Constraint con2 = new Constraint(lhs, Relation_Symbol.GREATER_OR_EQUAL, rhs2);
				q.add_constraint(con1);
				q.add_constraint(con2);
			}
		}
	}

	/**
	 * Sorts the points to be drawn as polyhedra. To do so, take a Jarvis march around each convex Polyhedron and put the points
	 * back in the point list in that order.
	 **/
	private void sortPoints() throws PrismException
	{
		double x_max = Double.MIN_VALUE;
		double x_min = Double.MAX_VALUE;
		double y_max = Double.MIN_VALUE;
		double y_min = Double.MAX_VALUE;

		// non-convex union of all the polyhedra of poly_points
		for (List<Point> ps : poly_points) {
			int size = ps.size();
			// do a Jarvis march around ps
			if (problem)
				System.out.printf("\nJarvis march around: [");
			for (Point p : ps) {
				// rounding points
				p.setCoord(0, ((double) Math.round(p.getCoord(0) * ((double)denominator))) * acc);
				p.setCoord(1, ((double) Math.round(p.getCoord(1) * ((double)denominator))) * acc);
				if (problem)
					System.out.printf("%s, ", p.toString());

				// determine min and max of non-border points
				if (p.getCoord(0) < (cutoff_max[indexX] > 0.0 ? cutoff_max[indexX] * f2 : cutoff_max[indexX] / f2)
						&& p.getCoord(0) > (cutoff_min[indexX] > 0.0 ? cutoff_min[indexX] / f2 : cutoff_min[indexX] * f2)
						&& p.getCoord(1) < (cutoff_max[indexY] > 0.0 ? cutoff_max[indexY] * f2 : cutoff_max[indexY] / f2)
						&& p.getCoord(1) > (cutoff_min[indexY] > 0.0 ? cutoff_min[indexY] / f2 : cutoff_min[indexY] * f2)) {
					// reset extreme values
					if (x_max < p.getCoord(0))
						x_max = p.getCoord(0);
					if (x_min > p.getCoord(0))
						x_min = p.getCoord(0);
					if (y_max < p.getCoord(1))
						y_max = p.getCoord(1);
					if (y_min > p.getCoord(1))
						y_min = p.getCoord(1);
				}
			}
			if (problem)
				System.out.printf("]\n");

			if (size <= 1)
				continue;

			// find start point in bottom left corner
			Point start = ps.get(0);
			for (Point p : ps) {
				if ((p.getCoord(1) == start.getCoord(1) && p.getCoord(0) < start.getCoord(0)) || p.getCoord(1) < start.getCoord(1)) {
					start = p;
					if (problem)
						System.out.printf("startpoint: %s\n", start);
				}
			}

			Point p = start;
			List<Point> qs = new ArrayList<Point>();
			do {
				qs.add(p);
				if (p != start)
					ps.remove(p); // no need to check again
				Point e = null;
				for (Point q : ps) {
					if (e == null) {
						e = q;
						continue; // select first point
					}

					double angle1 = Math.atan2(p.getCoord(1) - e.getCoord(1), p.getCoord(0) - e.getCoord(0)) - Math.PI;
					if (angle1 < 0.0)
						angle1 += Math.PI * 2;
					double angle2 = Math.atan2(p.getCoord(1) - q.getCoord(1), p.getCoord(0) - q.getCoord(0)) - Math.PI;
					if (angle2 < 0.0)
						angle2 += Math.PI * 2;
					if (problem)
						System.out.printf("(p: %s, e: %s, q: %s) angle1: %f, angle2: %f\n", p, e, q, angle1 == Double.MAX_VALUE ? -1.0 : angle1,
								angle2 == Double.MAX_VALUE ? -1.0 : angle2);
					if (!p.equals(q) && (p.equals(e) || angle1 >= angle2)) {
						e = q;
					}
				}
				p = e;
			} while (!p.equals(start));

			// substitute ps in poly_points by qs, which is now Jarvis-march sorted
			ps.clear();
			ps.addAll(qs);

			if (problem)
				System.out.printf("JM: %s\n", ps);
		}

		// if bounds were not initialised, do so now
		if (x_max == Double.MIN_VALUE)
			x_max = (cutoff_max[indexX] > 0.0 ? cutoff_max[indexX] * f3 : cutoff_max[indexX] / f3);
		if (x_min == Double.MAX_VALUE)
			x_min = (cutoff_min[indexX] > 0.0 ? cutoff_min[indexX] / f3 : cutoff_min[indexX] * f3);
		if (y_max == Double.MIN_VALUE)
			y_max = (cutoff_max[indexY] > 0.0 ? cutoff_max[indexY] * f3 : cutoff_max[indexY] / f3);
		if (y_min == Double.MAX_VALUE)
			y_min = (cutoff_min[indexY] > 0.0 ? cutoff_min[indexY] / f3 : cutoff_min[indexY] * f3);

		// if some dimension has zero spread, make it visible by adding some default spread
		if (Math.abs(x_max - x_min) < 0.1) {
			x_max += 1.0;
			x_min -= 1.0;
		}
		if (Math.abs(y_max - y_min) < 0.1) {
			y_max += 1.0;
			y_min -= 1.0;
		}

		// add extreme points to set boundary of the graph
		boundary_points.add(new prism.Point(new double[] { x_min, y_min }));
		boundary_points.add(new prism.Point(new double[] { x_min, y_max }));
		boundary_points.add(new prism.Point(new double[] { x_max, y_max }));
		boundary_points.add(new prism.Point(new double[] { x_max, y_min }));
	}

	@Override
	public List<Point> getPoints()
	{
		return boundary_points;
	}

	@Override
	public void addNewPoint(Point point) throws PrismException
	{
		this.boundary_points.add(point);
	}
}
