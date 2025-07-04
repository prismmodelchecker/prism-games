//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford)
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

import java.io.File;
import java.util.BitSet;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.PrimitiveIterator;
import java.util.Set;
import java.util.TreeMap;
import java.util.function.IntPredicate;

import common.IterableStateSet;
import common.IteratorTools;
import io.DotExporter;
import io.ModelExportOptions;
import io.PrismExplicitExporter;
import parser.State;
import parser.Values;
import parser.VarList;
import prism.Evaluator;
import prism.Prism;
import prism.PrismException;
import prism.PrismFileLog;
import prism.PrismLog;
import prism.PrismUtils;

import static prism.PrismSettings.DEFAULT_EXPORT_MODEL_PRECISION;

/**
 * Interface for (abstract) classes that provide (read-only) access to an explicit-state model.
 * This is generic, where probabilities/rates/etc. are of type {@code Value}.
 */
public interface Model<Value> extends prism.Model<Value>
{
	// Accessors

	/**
	 * Get iterator over initial state list.
	 */
	Iterable<Integer> getInitialStates();

	/**
	 * Get the index of the first initial state
	 * (i.e. the one with the lowest index).
	 * Returns -1 if there are no initial states.
	 */
	int getFirstInitialState();

	/**
	 * Check whether a state is an initial state.
	 */
	boolean isInitialState(int i);

	/**
	 * Get the number of states that are/were deadlocks.
	 * (Such states may have been fixed at build-time by adding self-loops)
	 */
	int getNumDeadlockStates();

	/**
	 * Get iterator over states that are/were deadlocks.
	 * (Such states may have been fixed at build-time by adding self-loops)
	 */
	Iterable<Integer> getDeadlockStates();
	
	/**
	 * Get list of states that are/were deadlocks.
	 * (Such states may have been fixed at build-time by adding self-loops)
	 */
	StateValues getDeadlockStatesList();
	
	/**
	 * Get the index of the first state that is/was a deadlock.
	 * (i.e. the one with the lowest index).
	 * Returns -1 if there are no initial states.
	 */
	int getFirstDeadlockState();

	/**
	 * Check whether a state is/was deadlock.
	 * (Such states may have been fixed at build-time by adding self-loops)
	 */
	boolean isDeadlockState(int i);
	
	/**
	 * Get access to a list of states (optionally stored).
	 */
	List<State> getStatesList();

	/** Get access to the VarList (optionally stored) */
	VarList getVarList();

	/**
	 * Get access to a list of constant values (optionally stored).
	 */
	Values getConstantValues();
	
	/**
	 * Get the states that satisfy a label in this model (optionally stored).
	 * Returns null if there is no label of this name.
	 */
	BitSet getLabelStates(String name);
	
	/**
	 * Get the labels that are (optionally) stored.
	 * Returns an empty set if there are no labels.
	 */
	Set<String> getLabels();

	/**
	 * Returns true if a label with the given name is attached to this model
	 */
	boolean hasLabel(String name);

	/**
	 * Get the mapping from labels that are (optionally) stored
	 * to the sets of states that satisfy them.
	 */
	default Map<String, BitSet> getLabelToStatesMap()
	{
		// Default implementation creates a new map on demand
		Map<String, BitSet> labels = new TreeMap<>();
		for (String name : getLabels()) {
			labels.put(name, getLabelStates(name));
		}
		return labels;
	}
	
	@Override
	default int getNumTransitions()
	{
		int numStates = getNumStates();
		int numTransitions = 0;
		for (int s = 0; s < numStates; s++) {
			numTransitions += getNumTransitions(s);
		}
		return numTransitions;
	}

	/**
	 * Get the number of transitions from state s.
	 */
	default int getNumTransitions(int s)
	{
		return Math.toIntExact(IteratorTools.count(getSuccessorsIterator(s)));
	}

	/**
	 * Get the number of transitions leaving a set of states.
	 * <br>
	 * Default implementation: Iterator over the states and sum the result of getNumTransitions(s).
	 * @param states The set of states, specified by an OfInt iterator
	 * @return the number of transitions
	 */
	default long getNumTransitions(PrimitiveIterator.OfInt states)
	{
		long count = 0;
		while (states.hasNext()) {
			int s = states.nextInt();
			count += getNumTransitions(s);
		}
		return count;
	}

	/**
	 * Get an iterator over the successors of state s.
	 * Default implementation via the SuccessorsIterator returned
	 * from {@code getSuccessors}, ensuring that there are no
	 * duplicates.
	 */
	default Iterator<Integer> getSuccessorsIterator(int s)
	{
		SuccessorsIterator successors = getSuccessors(s);
		return successors.distinct();
	}

	/**
	 * Get a SuccessorsIterator for state s.
	 */
	SuccessorsIterator getSuccessors(int s);

	/**
	 * Returns true if state s2 is a successor of state s1.
	 */
	default boolean isSuccessor(int s1, int s2)
	{
		// the code for this method is equivalent to the following stream expression,
		// but kept explicit for performance
		//
		// return getSuccessors(s1).stream().anyMatch(
		//           (t) -> {return t == s2;}
		// );

		SuccessorsIterator it = getSuccessors(s1);
		while (it.hasNext()) {
			int t = it.nextInt();
			if (t == s2)
				return true;
		}
		return false;
	}

	/**
	 * Check if all the successor states of a state are in a set.
	 * @param s The state to check
	 * @param set The set to test for inclusion
	 */
	default boolean allSuccessorsInSet(int s, BitSet set)
	{
		return allSuccessorsMatch(s, set::get);
	}

	/**
	 * Check if any successor states of a state are in a set.
	 * @param s The state to check
	 * @param set The set to test for inclusion
	 */
	default boolean someSuccessorsInSet(int s, BitSet set)
	{
		return someSuccessorsMatch(s, set::get);
	}

	/**
	 * Check if all the successor states of a state match the predicate.
	 * @param s The state to check
	 * @param p the predicate
	 */
	default boolean allSuccessorsMatch(int s, IntPredicate p)
	{
		// the code for this method is equivalent to the following stream expression,
		// but kept explicit for performance
		//
		// return getSuccessors(s).stream().allMatch(p);

		SuccessorsIterator it = getSuccessors(s);
		while (it.hasNext()) {
			int t = it.nextInt();
			if (!p.test(t))
				return false;
		}
		return true;
	}

	/**
	 * Check if any successor states of a state match the predicate.
	 * @param s The state to check
	 * @param p the predicate
	 */
	default boolean someSuccessorsMatch(int s, IntPredicate p)
	{
		// the code for this method is equivalent to the following stream expression,
		// but kept explicit for performance
		//
		// return getSuccessors(s).stream().anyMatch(p);

		SuccessorsIterator it = getSuccessors(s);
		while (it.hasNext()) {
			int t = it.nextInt();
			if (p.test(t))
				return true;
		}
		return false;
	}

	/**
	 * Perform a single step of precomputation algorithm Prob0 for a single state,
	 * i.e., for the state {@code s} returns true iff there is a transition from
	 * {@code s} to a state in {@code u}.
	 * <br>
	 * <i>Default implementation</i>: Iterates using {@code getSuccessors()} and performs the check.
	 * @param s The state in question
	 * @param u Set of states {@code u}
	 * @return true iff there is a transition from s to a state in u
	 */
	default boolean prob0step(int s, BitSet u)
	{
		for (SuccessorsIterator succ = getSuccessors(s); succ.hasNext(); ) {
			int t = succ.nextInt();
			if (u.get(t))
				return true;
		}
		return false;
	}

	/**
	 * Perform a single step of precomputation algorithm Prob0, i.e., for states i in {@code subset},
	 * set bit i of {@code result} iff there is a transition to a state in {@code u}.
	 * <br>
	 * <i>Default implementation</i>: Iterate over {@code subset} and use {@code prob0step(s,u)}
	 * to determine result for {@code s}.
	 * @param subset Only compute for these states
	 * @param u Set of states {@code u}
	 * @param result Store results here
	 */
	default void prob0step(BitSet subset, BitSet u, BitSet result)
	{
		for (PrimitiveIterator.OfInt it = new IterableStateSet(subset, getNumStates()).iterator(); it.hasNext();) {
			int s = it.nextInt();
			result.set(s, prob0step(s,u));
		}
	}

	/**
	 * Perform a single step of precomputation algorithm Prob1 for a single state,
	 * i.e., for states s return true iff there is a transition to a state in
	 * {@code v} and all transitions go to states in {@code u}.
	 * @param s The state in question
	 * @param u Set of states {@code u}
	 * @param v Set of states {@code v}
	 * @return true iff there is a transition from s to a state in v and all transitions go to u.
	 */
	default boolean prob1step(int s, BitSet u, BitSet v)
	{
		boolean allTransitionsToU = true;
		boolean hasTransitionToV = false;
		for (SuccessorsIterator succ = getSuccessors(s); succ.hasNext(); ) {
			int t = succ.nextInt();
			if (!u.get(t)) {
				allTransitionsToU = false;
				// early abort, as overall result is false
				break;
			}
			hasTransitionToV = hasTransitionToV || v.get(t);
		}
		return (allTransitionsToU && hasTransitionToV);
	}

	/**
	 * Perform a single step of precomputation algorithm Prob1, i.e., for states i in {@code subset},
	 * set bit i of {@code result} iff there is a transition to a state in {@code v} and all transitions go to states in {@code u}.
	 * @param subset Only compute for these states
	 * @param u Set of states {@code u}
	 * @param v Set of states {@code v}
	 * @param result Store results here
	 */
	default void prob1step(BitSet subset, BitSet u, BitSet v, BitSet result)
	{
		for (PrimitiveIterator.OfInt it = new IterableStateSet(subset, getNumStates()).iterator(); it.hasNext();) {
			int s = it.nextInt();
			result.set(s, prob1step(s,u,v));
		}
	}

	/**
	 * Find all deadlock states and store this information in the model.
	 * If requested (if fix=true) and if needed (i.e. for DTMCs/CTMCs),
	 * fix deadlocks by adding self-loops in these states.
	 * The set of deadlocks (before any possible fixing) can be obtained from {@link #getDeadlockStates()}.
	 * @throws PrismException if the model is unable to fix deadlocks because it is non-mutable.
	 */
	void findDeadlocks(boolean fix) throws PrismException;

	/**
	 * Checks for deadlocks and throws an exception if any exist.
	 */
	void checkForDeadlocks() throws PrismException;

	/**
	 * Checks for deadlocks and throws an exception if any exist.
	 * States in 'except' (If non-null) are excluded from the check.
	 */
	void checkForDeadlocks(BitSet except) throws PrismException;

	// Export methods (explicit files)

	/**
	 * Export to explicit format readable by PRISM (i.e. a .tra file, etc.).
	 */
	default void exportToPrismExplicit(String baseFilename) throws PrismException
	{
		exportToPrismExplicit(baseFilename, DEFAULT_EXPORT_MODEL_PRECISION);
	}

	/**
	 * Export to explicit format readable by PRISM (i.e. a .tra file, etc.).
	 * @param precision number of significant digits >= 1
	 */
	default void exportToPrismExplicit(String baseFilename, int precision) throws PrismException
	{
		// Default implementation - just output .tra file
		// (some models might override this)
		exportToPrismExplicitTra(baseFilename + ".tra", precision);
	}

	/**
	 * Export transition matrix to explicit format readable by PRISM (i.e. a .tra file).
	 */
	default void exportToPrismExplicitTra(PrismLog out, ModelExportOptions exportOptions) throws PrismException
	{
		new PrismExplicitExporter<Value>(exportOptions).exportTransitions(this, out);
	}

	/**
	 * Export transition matrix to explicit format readable by PRISM (i.e. a .tra file).
	 * @param precision number of significant digits >= 1
	 */
	default void exportToPrismExplicitTra(PrismLog out, int precision) throws PrismException
	{
		exportToPrismExplicitTra(out, new ModelExportOptions().setModelPrecision(precision));
	}

	/**
	 * Export transition matrix to explicit format readable by PRISM (i.e. a .tra file).
	 */
	default void exportToPrismExplicitTra(PrismLog out) throws PrismException
	{
		exportToPrismExplicitTra(out, new ModelExportOptions());
	}

	/**
	 * Export transition matrix to explicit format readable by PRISM (i.e. a .tra file).
	 */
	default void exportToPrismExplicitTra(String filename) throws PrismException
	{
		try (PrismFileLog out = PrismFileLog.create(filename)) {
			exportToPrismExplicitTra(out);
		}
	}

	/**
	 * Export transition matrix to explicit format readable by PRISM (i.e. a .tra file).
	 * @param precision number of significant digits >= 1
	 */
	default void exportToPrismExplicitTra(String filename, int precision) throws PrismException
	{
		try (PrismFileLog out = PrismFileLog.create(filename)) {
			exportToPrismExplicitTra(out, precision);
		}
	}

	/**
	 * Export transition matrix to explicit format readable by PRISM (i.e. a .tra file).
	 */
	default void exportToPrismExplicitTra(File file) throws PrismException
	{
		try (PrismFileLog out = PrismFileLog.create(file.getPath())) {
			exportToPrismExplicitTra(out);
		}
	}

	/**
	 * Export transition matrix to explicit format readable by PRISM (i.e. a .tra file).
	 * @param precision number of significant digits >= 1
	 */
	default void exportToPrismExplicitTra(File file, int precision) throws PrismException
	{
		try (PrismFileLog out = PrismFileLog.create(file.getPath())) {
			exportToPrismExplicitTra(out, precision);
		}
	}

	// Export methods (dot files)

	/**
	 * Export to a dot file, highlighting states in 'mark'.
	 * @param out PrismLog to export to
	 * @param exportOptions Options for export
	 * @param decorators Any Dot decorators to add (ignored if null)
	 */
	default void exportToDotFile(PrismLog out, ModelExportOptions exportOptions, Iterable<explicit.graphviz.Decorator> decorators) throws PrismException
	{
		new DotExporter<Value>(exportOptions).exportModel(this, out, decorators);
	}

	/**
	 * Export to a dot file, highlighting states in 'mark'.
	 * @param out PrismLog to export to
	 * @param exportOptions Options for export
	 */
	default void exportToDotFile(PrismLog out, ModelExportOptions exportOptions) throws PrismException
	{
		new DotExporter<Value>(exportOptions).exportModel(this, out, null);
	}

	/**
	 * Export to a dot file.
	 * @param out PrismLog to export to
	 */
	default void exportToDotFile(PrismLog out) throws PrismException
	{
		exportToDotFile(out, new ModelExportOptions());
	}

	/**
	 * Export to a dot file.
	 * @param out PrismLog to export to
	 * @param precision number of significant digits >= 1
	 */
	default void exportToDotFile(PrismLog out, int precision) throws PrismException
	{
		exportToDotFile(out, new ModelExportOptions().setModelPrecision(precision));
	}

	/**
	 * Export to a dot file, highlighting states in 'mark'.
	 * @param out PrismLog to export to
	 * @param mark States to highlight (ignored if null)
	 */
	default void exportToDotFile(PrismLog out, BitSet mark) throws PrismException
	{
		Iterable<explicit.graphviz.Decorator> decorators = (mark == null) ? null : Collections.singleton(new explicit.graphviz.MarkStateSetDecorator(mark));
		exportToDotFile(out, new ModelExportOptions(), decorators);
	}

	/**
	 * Export to a dot file, highlighting states in 'mark'.
	 * @param out PrismLog to export to
	 * @param mark States to highlight (ignored if null)
	 * @param precision number of significant digits >= 1
	 */
	default void exportToDotFile(PrismLog out, BitSet mark, int precision) throws PrismException
	{
		Iterable<explicit.graphviz.Decorator> decorators = (mark == null) ? null : Collections.singleton(new explicit.graphviz.MarkStateSetDecorator(mark));
		exportToDotFile(out, new ModelExportOptions().setModelPrecision(precision), decorators);
	}

	/**
	 * Export to a dot file, highlighting states in 'mark'.
	 * @param out PrismLog to export to
	 * @param mark States to highlight (ignored if null)
	 * @param showStates Show state info on nodes?
	 */
	default void exportToDotFile(PrismLog out, BitSet mark, boolean showStates) throws PrismException
	{
		Iterable<explicit.graphviz.Decorator> decorators = (mark == null) ? null : Collections.singleton(new explicit.graphviz.MarkStateSetDecorator(mark));
		exportToDotFile(out, new ModelExportOptions().setShowStates(showStates), decorators);
	}

	/**
	 * Export to a dot file, highlighting states in 'mark'.
	 * @param out PrismLog to export to
	 * @param mark States to highlight (ignored if null)
	 * @param showStates Show state info on nodes?
	 * @param precision number of significant digits >= 1
	 */
	default void exportToDotFile(PrismLog out, BitSet mark, boolean showStates, int precision) throws PrismException
	{
		Iterable<explicit.graphviz.Decorator> decorators = (mark == null) ? null : Collections.singleton(new explicit.graphviz.MarkStateSetDecorator(mark));
		exportToDotFile(out, new ModelExportOptions().setShowStates(showStates).setModelPrecision(precision), decorators);
	}

	/**
	 * Export to a dot file, decorating states and transitions with the provided decorators
	 * @param out PrismLog to export to
	 */
	default void exportToDotFile(PrismLog out, Iterable<explicit.graphviz.Decorator> decorators) throws PrismException
	{
		exportToDotFile(out, new ModelExportOptions(), decorators);
	}

	/**
	 * Export to a dot file, decorating states and transitions with the provided decorators
	 * @param out PrismLog to export to
	 * @param precision number of significant digits >= 1
	 */
	default void exportToDotFile(PrismLog out, Iterable<explicit.graphviz.Decorator> decorators, int precision) throws PrismException
	{
		exportToDotFile(out, new ModelExportOptions().setModelPrecision(precision), decorators);
	}

	/**
	 * Export to a dot file.
	 * @param filename Name of file to export to
	 */
	default void exportToDotFile(String filename) throws PrismException
	{
		try (PrismFileLog out = PrismFileLog.create(filename)) {
			exportToDotFile(out);
		}
	}

	/**
	 * Export to a dot file.
	 * @param filename Name of file to export to
	 * @param precision number of significant digits >= 1
	 */
	default void exportToDotFile(String filename, int precision) throws PrismException
	{
		try (PrismFileLog out = PrismFileLog.create(filename)) {
			exportToDotFile(out, precision);
		}
	}

	/**
	 * Export to a dot file, highlighting states in 'mark'.
	 * @param filename Name of file to export to
	 * @param mark States to highlight (ignored if null)
	 */
	default void exportToDotFile(String filename, BitSet mark) throws PrismException
	{
		try (PrismFileLog out = PrismFileLog.create(filename)) {
			exportToDotFile(out, mark);
		}
	}

	/**
	 * Export to a dot file, highlighting states in 'mark'.
	 * @param filename Name of file to export to
	 * @param mark States to highlight (ignored if null)
	 * @param precision number of significant digits >= 1
	 */
	default void exportToDotFile(String filename, BitSet mark, int precision) throws PrismException
	{
		try (PrismFileLog out = PrismFileLog.create(filename)) {
			exportToDotFile(out, mark, precision);
		}
	}

	/**
	 * Export to a dot file, decorating states and transitions with the provided decorators
	 * @param filename Name of the file to export to
	 */
	default void exportToDotFile(String filename, Iterable<explicit.graphviz.Decorator> decorators) throws PrismException
	{
		try (PrismFileLog out = PrismFileLog.create(filename)) {
			exportToDotFile(out, decorators);
		}
	}

	/**
	 * Export to a dot file, decorating states and transitions with the provided decorators
	 * @param filename Name of the file to export to
	 * @param precision number of significant digits >= 1
	 */
	default void exportToDotFile(String filename, Iterable<explicit.graphviz.Decorator> decorators, int precision) throws PrismException
	{
		try (PrismFileLog out = PrismFileLog.create(filename)) {
			exportToDotFile(out, decorators, precision);
		}
	}

	@Deprecated
	default void exportTransitionsToDotFile(int i, PrismLog out, Iterable<explicit.graphviz.Decorator> decorators)
	{
		// This does not need to be implemented (it will be ignored)
	}

	@Deprecated
	default void exportTransitionsToDotFile(int i, PrismLog out, Iterable<explicit.graphviz.Decorator> decorators, int precision)
	{
		// This no longer needs to be implemented (it will be ignored)
	}

	/**
	 * Export to a equivalent PRISM language model description.
	 */
	default void exportToPrismLanguage(String filename) throws PrismException
	{
		exportToPrismLanguage(filename, DEFAULT_EXPORT_MODEL_PRECISION);
	}

	/**
	 * Export to a equivalent PRISM language model description.
	 * @param precision number of significant digits >= 1
	 */
	void exportToPrismLanguage(String filename, int precision) throws PrismException;
	
	/**
	 * Export states list.
	 */
	default void exportStates(VarList varList, PrismLog out, ModelExportOptions exportOptions) throws PrismException
	{
		new PrismExplicitExporter<Value>(exportOptions).exportStates(this, varList, out);
	}

	/**
	 * @deprecated
	 * Export states list.
	 */
	@Deprecated
	default void exportStates(int exportType, VarList varList, PrismLog out) throws PrismException
	{
		exportStates(varList, out, Prism.convertExportType(exportType));
	}

	/**
	 * Report info/stats about the model as a string.
	 */
	default String infoString()
	{
		final int numStates = getNumStates();
		String s = "";
		if (getNumPlayers() > 1) {
			s += ", " + getNumPlayers() + " players";
		}
		s += numStates + " states (" + getNumInitialStates() + " initial)";
		if (this instanceof PartiallyObservableModel) {
			s += ", " + ((PartiallyObservableModel<?>) this).getNumObservations() + " observables";
			s += ", " + ((PartiallyObservableModel<?>) this).getNumUnobservations() + " unobservables";
		}
		s += ", " + getNumTransitions() + " transitions";
		if (this instanceof NondetModel) {
			s += ", " + ((NondetModel<?>) this).getNumChoices() + " choices";
			s += ", dist max/avg = " + ((NondetModel<?>) this).getMaxNumChoices() + "/" + PrismUtils.formatDouble2dp(((double) ((NondetModel<?>) this).getNumChoices()) / numStates);
		}
		return s;
	}

	/**
	 * Report info/stats about the model, tabulated, as a string.
	 */
	default String infoStringTable()
	{
		final int numStates = getNumStates();
		String s = "";
		if (getNumPlayers() > 1) {
			s += "Players:     " + getNumPlayers() + "\n";
		}
		s += "States:      " + numStates + " (" + getNumInitialStates() + " initial)\n";
		if (this instanceof PartiallyObservableModel) {
			s += "Obs/unobs:   " + ((PartiallyObservableModel<?>) this).getNumObservations() + "/" + ((PartiallyObservableModel<?>) this).getNumUnobservations() + "\n";
		}
		s += "Transitions: " + getNumTransitions() + "\n";
		if (this instanceof NondetModel) {
			s += "Choices:     " + ((NondetModel<?>) this).getNumChoices() + "\n";
			s += "Max/avg:     " + ((NondetModel<?>) this).getMaxNumChoices() + "/" + PrismUtils.formatDouble2dp(((double) ((NondetModel<?>) this).getNumChoices()) / numStates) + "\n";
		}
		return s;
	}

	/** Has this model a stored PredecessorRelation? */
	boolean hasStoredPredecessorRelation();

	/**
	 * If there is a PredecessorRelation stored for this model, return that.
	 * Otherwise, create one and return that. If {@code storeIfNew},
	 * store it for later use.
	 *
	 * @param parent a PrismComponent (for obtaining the log)
	 * @param storeIfNew if the predecessor relation is newly created, store it
	 */
	PredecessorRelation getPredecessorRelation(prism.PrismComponent parent, boolean storeIfNew);

	/** Clear any stored predecessor relation, e.g., because the model was modified */
	void clearPredecessorRelation();

	/**
	 * Get an Evaluator for the values stored in this Model for probabilities etc.
	 * This is needed, for example, to compute probability sums, check for equality to 0/1, etc.
	 * A default implementation provides an evaluator for the (usual) case when Value is Double.
	 */
	@SuppressWarnings("unchecked")
	default Evaluator<Value> getEvaluator()
	{
		return (Evaluator<Value>) Evaluator.forDouble();
	}
}
