// ==============================================================================
//	
// Copyright (c) 2002-
// Authors:
// * Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford)
// * Clemens Wiltsche <clemens.wiltsche@cs.ox.ac.uk> (University of Oxford)
//	
// ------------------------------------------------------------------------------
//	
// This file is part of PRISM.
//	
// PRISM is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//	
// PRISM is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//	
// You should have received a copy of the GNU General Public License
// along with PRISM; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//	
// ==============================================================================

package explicit;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.apache.commons.math3.fraction.BigFraction;

import parma_polyhedra_library.C_Polyhedron;
import parma_polyhedra_library.Coefficient;
import parma_polyhedra_library.Constraint;
import parma_polyhedra_library.Generator;
import parma_polyhedra_library.Generator_System;
import parma_polyhedra_library.Generator_Type;
import parma_polyhedra_library.Linear_Expression;
import parma_polyhedra_library.Linear_Expression_Coefficient;
import parma_polyhedra_library.Linear_Expression_Times;
import parma_polyhedra_library.Polyhedron;
import parma_polyhedra_library.Relation_Symbol;
import parma_polyhedra_library.Variable;
import parser.State;
import parser.ast.Coalition;
import parser.ast.ModulesFile;
import prism.ModelType;
import prism.PrismException;
import prism.PrismLog;
import prism.PrismUtils;
import explicit.rewards.SMGRewards;

/**
 * Simple explicit-state representation of a (turn-based) stochastic multi-player game (SMG).
 * States can be labelled arbitrarily with player 1..n, player 0 has a special
 * purpose of scheduling the moves of other players
 */
public class SMG extends STPGExplicit implements STPG
{
	// NB: We re-use the existing stateOwners list in the superclass to assign states to players

	// Name, used for easier association in compositional engine
	protected String name = "";

	// A definition of the players in the game, i.e., the (integer) index and name of each one.
	// It is stored as a mapping from indices to names.
	// Indices can be arbitrary and do not need to be contiguous.
	// Names are optional and can be null or "" if undefined,
	// but all normally defined names should be unique.
	// The size of this map defines the number of players in the game.
	protected Map<Integer, String> playerNames;
	
	// Optionally, a mapping from player indices to 1 or 2,
	// as induced by a coalition of players (those in the coalition
	// are mapped to 1, those who are not are mapped to 2).
	// The mapped index of the player with original index {@code p}
	// is stored as {@code coalitionPlayerMap[p]}, regardless of the
	// range of player indices actually used.
	protected int[] coalitionPlayerMap;

	// for compositional strategy synthesis need additional state information
	protected List<Integer> controlled_by; // maps (global state) -> component index
	protected List<List<Integer>> local_state; // maps (global state, component) -> local state
	protected List<List<List<Integer>>> l_action; // maps (global state, global action, component) -> local action

        public final static String TAU = "tau";

	// Constructors

	/**
	 * Constructor: empty SMG.
	 */
	public SMG()
	{
		super();
		name = null;
		playerNames = new HashMap<Integer, String>();
		coalitionPlayerMap = null;
	}

	/**
	 * Constructor: empty, named SMG.
	 */
	public SMG(String name)
	{
		this();
		this.name = name;
		playerNames = new HashMap<Integer, String>();
		coalitionPlayerMap = null;
	}

	/**
	 * Constructor: new SMG with fixed number of states.
	 */
	public SMG(int numStates)
	{
		super(numStates);
		name = null;
		playerNames = new HashMap<Integer, String>();
		coalitionPlayerMap = null;
	}

	/**
	 * Construct an SMG from an existing one and a state index permutation,
	 * i.e. in which state index i becomes index permut[i].
	 * Player and coalition info is also copied across.
	 */
	public SMG(SMG smg, int permut[])
	{
		super(smg, permut);
		name = smg.name;
		playerNames = new HashMap<Integer, String>(smg.playerNames);
		coalitionPlayerMap = smg.coalitionPlayerMap == null ? null : smg.coalitionPlayerMap.clone(); 
		// force recomputation of the strategy-related items
		l_action = null;
		local_state = null;
		controlled_by = null;
	}

	/**
	 * Copy constructor
	 */
	public SMG(SMG smg)
	{
		super(smg);
		name = smg.name;
		playerNames = new HashMap<Integer, String>(smg.playerNames);
		coalitionPlayerMap = smg.coalitionPlayerMap == null ? null : smg.coalitionPlayerMap.clone(); 
		controlled_by = new ArrayList<Integer>(smg.getControlledBy());
		local_state = new ArrayList<List<Integer>>(smg.getLocalState());
	}

	/**
	 * Remove state by disabling all accessing transitions.
	 * 
	 * @param s State to be removed
	 * @param close_P1 Close game for Player 1 (coalition players)
	 * @param close_P2 Close game for Player 2 (non-coalition players)
	 * 
	 * @return set of newly disabled states
	 */
	public BitSet disableState(int s, boolean close_P1, boolean close_P2)
	{
		throw new UnsupportedOperationException();
		
		// Disabled for now because support disabling transitions has been removed
		// during refactoring of explicit engine model classes
		
		/*
		int numStates = trans.size();
		BitSet disabled = new BitSet(numStates);
		disabled.set(s);
		
		for(int t = 0; t < numStates; t++) {
			List<Distribution> trans_t = trans.get(t);
			int num_trans_t = trans_t.size();
			int disabled_trans_t = 0;
			for(int c = 0; c < num_trans_t; c++) {
				if(trans_t.get(c).getSupport().contains(s)) { // some transition t -> s
					// disable this transition
					disableChoice(t, c);					
				}
				// count how many transitions are disabled;
				if(disabledChoices.get(t) != null && disabledChoices.get(t).get(c))
					disabled_trans_t++;
			}
			
			// test if state became disabled, depending on whether the game should be closed for a type of players
			if(getPlayer(t) == 1 && close_P1 && (num_trans_t == disabled_trans_t) |
					getPlayer(t) == 2 && close_P2 && disabled_trans_t > 0) {
				disabled.or(disableState(t, close_P1, close_P2)); // recursively disable states, and add to return set
			}
		}
		
		return disabled;
		*/
	}
	
	public String getName()
	{
		return this.name;
	}

	/** 
	 * converts this to normal form
	 * needs to correspond to StochasticUpdateStrategy.java:toNormalForm
	 * if the game is already in normal form, then nothing is done
	 * if the game is not in normal form, but consists of several components already, an exception is thrown
	 **/
	public void toNormalForm() throws PrismException
	{
	        if (isNormalForm())
			return; // already in normal form
		if (getNumComponents() > 1)
			throw new PrismException("normal form conversion only possible in atomic games");
		int gameSize = numStates;
		for (int s = 0; s < gameSize; s++) {
			if (stateOwners.get(s) == 1) {
				// P1 states are split
				// attach new P1 state at end of list
				addState(1);
				statesList.add(statesList.get(s));
				// take all outgoing moves of s and attach them to the new state
				for (Distribution d : trans.get(s)) {
					trans.get(numStates - 1).add(d);
				}
				actions.set(numStates - 1, new ArrayList<Object>());
				for (Object a : actions.get(s)) {
					actions.get(numStates - 1).add(a);
				}
				// the original state becomes a P2 state
				stateOwners.set(s, 2);
				// clear original transitions
				trans.get(s).clear();
				actions.get(s).clear();
				// add tau-Dirac transition
				Distribution d = new Distribution();
				d.add(numStates - 1, 1.0);
				addActionLabelledChoice(s, d, TAU);
			}
		}
		// force recomputation of the strategy-related items
		l_action = null;
		getLAction();
		local_state = null;
		getLocalState();
		controlled_by = null;
		getControlledBy();
	}

	// Mutators

	/**
	 * Tests whether the game is in normal form
	 */
	public boolean isNormalForm() throws PrismException
	{
		for (int s = 0; s < initialStates.size(); s++) {
		        if (stateOwners.get(initialStates.get(s)) != 2)
				return false; // initial state not P2
		}

		for (int s = 0; s < numStates; s++) { // for every state
		        if(actions.get(s) == null) { // fill in TAUs
			    actions.set(s, new ArrayList<Object>());
			    for(int c = 0; c < getNumChoices(s); c++)
				actions.get(s).add(TAU);
		        }
			for (int c = 0; c < getNumChoices(s); c++) { // for every successor
				Distribution d = trans.get(s).get(c);
				if(actions.get(s) == null) {
				    throw new PrismException("Transitions must be labelled");
				} else if (actions.get(s).get(c) == null || actions.get(s).get(c).equals(TAU)) {
				    actions.get(s).set(c, TAU);
				        if (stateOwners.get(s) != 2)
						return false; // tau-transition does not originate from P2
					if (d.size() != 1 || !PrismUtils.doublesAreEqual(d.sum(), 1.0))
						return false; // not Dirac
					for (Integer t : d.getSupport()) {
					        if (stateOwners.get(t) != 1)
							return false; // tau-transition does not lead to P1
					}
				} else {
					for (Integer t : d.getSupport()) {
					        if (stateOwners.get(t) == 1)
							return false; // non-tau-transition assigns nonzero probability to P1
					}
				}
			}
		}

		// fall through only if all conditions for normal form are met
		return true;
	}

	/**
	 * Product constructor.
	 **/
    public SMG(ModulesFile mf, List<ModulesFile> mfs, List<SMG> smgs, PrismLog mainLog, boolean[] cancel_computation) throws PrismException
	{
	    this(mf, mfs, smgs, mainLog, false, cancel_computation); // no compatibility check is default
	}

    public SMG(ModulesFile mf, List<ModulesFile> mfs, List<SMG> smgs, PrismLog mainLog, boolean checkCompatibility, boolean[] cancel_computation) throws PrismException
	{
		int n = smgs.size();
		int estSize = 1;

		// check normal form and estimate size
		for (SMG smg : smgs) {
		        if (!smg.isNormalForm())
				throw new PrismException("Game " + smgs.indexOf(smg) + " not in normal form");
			estSize *= smg.statesList.size();
		}

		// initialise name
		this.name = "";
		boolean first = true;
		for (SMG smg : smgs) {
			if (!first)
				this.name += " || ";
			first = false;
			this.name += smg.name;
		}

		// initialise with empty state space
		super.initialise(0);

		// initial distribution ... one choice
		State initial = new State(n);
		for (int i = 0; i < n; i++)
			initial.varValues[i] = smgs.get(i).getFirstInitialState();

		// player info: 2 unnamed players
		playerNames = new HashMap<Integer, String>();
		playerNames.put(1, null);
		playerNames.put(2, null);
		
		// the coalition: player 1 plays by itself in coalition
		coalitionPlayerMap = new int[3];
		coalitionPlayerMap[0] = -1;
		coalitionPlayerMap[1] = 1;
		coalitionPlayerMap[2] = 2;
		
		// state space and transition relation inductively defined
		// initial state is in the product state space
		statesList = new ArrayList<State>();
		statesList.add(initial);
		initialStates.add(statesList.indexOf(initial));
		stateOwners = new ArrayList<Integer>();
		stateOwners.add(2); // normal form initial state is always P2
		actions = new ArrayList<List<Object>>();
		actions.add(new ArrayList<Object>()); // empty actions for initial state
		trans.add(new ArrayList<Distribution>()); // empty distribution list for initial state

		// temporary transitions to be converted into actual transitions later
		List<List<List<Distribution>>> tempTrans = new ArrayList<List<List<Distribution>>>();
		tempTrans.add(new ArrayList<List<Distribution>>()); // empty list for initial state

		// inductive step
		if(estSize >= 1000)
		    mainLog.print("Progress:\n");
		int pointer = 0; // points to the state we're operating on
		while (pointer < statesList.size()) {
			// get state to explore
			State s = statesList.get(pointer);
			addTransition(smgs, s, s, new ArrayList<Distribution>(), null, 0, n, tempTrans);
			pointer++;

			// log progress
			if(pointer % 1000 == 0 && pointer > 0) {
			    mainLog.print(String.format("(%d/~%d)", pointer, estSize));
			    if (cancel_computation[0]) {
				cancel_computation[0] = false; // reset
				throw new PrismException("Computation cancelled");
			    }
			}
			if(pointer % 10000 == 0 && pointer > 0)
			    mainLog.print("\n");
		}
		if(estSize >= 1000)
		    mainLog.print("\n");

		// convert temporary transitions into proper transitions
		for (int s = 0; s < tempTrans.size(); s++) { // for each state
			for (List<Distribution> ds : tempTrans.get(s)) { // for each outgoing move
				Distribution pd = new Distribution(ds, statesList, statesList.get(s), 1.0, 0, n);
				trans.get(s).add(pd);
			}
		}

		// remove transitions to states deadlocked on TAU
		for (int s = 0; s < trans.size(); s++) { // for each state
			for (int c = 0; c < trans.get(s).size(); c++) {
			    if (TAU.equals(actions.get(s).get(c))) {
					Distribution d = trans.get(s).get(c);
					go_through_distr: for (Integer t : d.keySet()) { // go through all successors of the distribution
						if (trans.get(t).size() == 0) { // deadlocked
							// remove the transition from s to t
							trans.get(s).remove(c);
							tempTrans.get(s).remove(c);
							actions.get(s).remove(c);
							c--; // since removed current index, need to look at it again.
							break go_through_distr;
						}
					}
				}
			}
		}

		// remove unreachable states
		search_for_unreachable: while (true) {
			BitSet unreachable = new BitSet(trans.size());
			// test if reachable
			for (int s = 0; s < trans.size(); s++) { // for each state
				for (int c = 0; c < trans.get(s).size(); c++) { // for each outgoing transition
					Distribution d = trans.get(s).get(c);
					for (Integer t : d.keySet()) {
						unreachable.set(t); // set state to be reachable
					}
				}
			}
			unreachable.flip(0, trans.size());
			unreachable.clear(0); // first state is reachable
			int states_removed = 0;
			int u = unreachable.nextSetBit(0);
			if (u < 0)
				break search_for_unreachable;
			while (u > 0) {
				int v = u - states_removed; // adjust for having previously removed states
				// remove unreachable state
				for (int s = 0; s < trans.size(); s++) {
					for (int c = 0; c < trans.get(s).size(); c++) {
						Distribution d = trans.get(s).get(c);
						Distribution new_d = new Distribution();
						for (Integer t : d.keySet()) {
							if (t < v) {
								new_d.set(t, d.get(t)); // no change
							} else if (t > v) {
								new_d.set(t - 1, d.get(t)); // shift down
							} else {
								throw new PrismException("Trying to remove the reachable state " + v);
							}
						}
						trans.get(s).set(c, new_d); // update transition
					}
				}
				statesList.remove(v);
				stateOwners.remove(v);
				actions.remove(v);
				trans.remove(v);
				tempTrans.remove(v);
				states_removed++;

				u = unreachable.nextSetBit(u + 1);
			}
		}

		// set number of states
		this.numStates = statesList.size();

		// fill in: global states, local actions, controlled by, variables
		// the modules file mf contains the global order of variables
		// the modules files mfs contain the local orders of variables
		this.controlled_by = new ArrayList<Integer>(statesList.size());

		// keep track of subcomponents when composing several times
		List<Integer> cumulativeComponents = new ArrayList<Integer>(smgs.size());
		cumulativeComponents.add(0);
		for (int i = 0; i < smgs.size(); i++)
			cumulativeComponents.add(cumulativeComponents.get(cumulativeComponents.size() - 1) + smgs.get(i).getNumComponents());
		int tot_c = cumulativeComponents.get(smgs.size());

		this.local_state = new ArrayList<List<Integer>>(statesList.size());
		this.l_action = new ArrayList<List<List<Integer>>>(statesList.size());

		for (int s = 0; s < statesList.size(); s++) {
			State state = statesList.get(s);
			local_state.add(new ArrayList<Integer>(tot_c));
			l_action.add(new ArrayList<List<Integer>>(trans.get(s).size())); // one for each global action
			for (int ga = 0; ga < trans.get(s).size(); ga++) { // for each global action
				l_action.get(s).add(new ArrayList<Integer>());
			}
			Map<String, Object> newVars = new HashMap<String, Object>(mf.getNumVars());
			for (int i = 0; i < smgs.size(); i++) {
				int si = (int) (state.varValues[i]); // state of component i
				// LOCAL STATES
				for (int c = cumulativeComponents.get(i); c < cumulativeComponents.get(i + 1); c++) {
					local_state.get(s).add(smgs.get(i).getLocalState().get(si).get(c - cumulativeComponents.get(i)));
				}
				// LOCAL ACTION
				for (int c = cumulativeComponents.get(i); c < cumulativeComponents.get(i + 1); c++) {
					for (int ga = 0; ga < trans.get(s).size(); ga++) { // for each global action
						Object g_act = actions.get(s).get(ga);
						if (i < tempTrans.get(s).get(ga).size()) { // if the action is assigned
							// local distribution of component i corresponding to global action ga
							Distribution ldist = tempTrans.get(s).get(ga).get(i);
							// now need local index of action ga at component i
							// this means that both the distribution
							// and the action name need to correspond
							int la = 0;
							search_for_la: for (la = 0; la < smgs.get(i).trans.get(si).size(); la++) {
								// note: ldist may be null
								if (ldist != null && smgs.get(i).trans.get(si).get(la).equals(ldist)) { // distributions agree
								    Object m_act = smgs.get(i).actions.get(si).get(la);
								    if ((g_act == null && m_act == null) || (g_act!=null && g_act.equals(m_act))) { // action names agree
										break search_for_la;
									}
								}
							}
							if (la == smgs.get(i).trans.get(si).size()) { // fallen through loop - no local action found
								la = -1; // indicates not found
							}
							if (la >= 0) {
								// local index of action ga at component c
								int lac = smgs.get(i).getLAction().get(si).get(la).get(c - cumulativeComponents.get(i));
								// local action la corresponding to global state s, global action ga and component c
								l_action.get(s).get(ga).add(lac);
							} else {
								l_action.get(s).get(ga).add(-1);
							}
						} else {
							l_action.get(s).get(ga).add(-1);
						}
					}
				}
				// CONTROLLED BY
				if (controlled_by.size() <= s) { // not added yet
					if (smgs.get(i).stateOwners.get(si).intValue() == 1) {
						// note: there is only one component that is player 1
						controlled_by.add(cumulativeComponents.get(i) + smgs.get(i).controlledBy(si));
					} else if (stateOwners.get(s).intValue() != 1) {
						controlled_by.add(-1); // -1 default for P2
					}
				}
				// VARIABLES
				State state_i = smgs.get(i).statesList.get(si);
				for (int j = 0; j < state_i.varValues.length; j++) {
					newVars.put(mfs.get(i).getVarName(j), state_i.varValues[j]);
				}
			}
			List<Object> newVarValues = new ArrayList<Object>();
			for (int j = 0; j < mf.getNumVars(); j++) {
				Object value_j = newVars.get(mf.getVarName(j));
				if (value_j != null)
					newVarValues.add(value_j);
			}
			state.varValues = newVarValues.toArray();
		}

		// set number of transitions (count like PRISM does)
		numTransitions = 0;
		for (int s = 0; s < trans.size(); s++) // for each state
		    for(int j = 0; j < trans.get(s).size(); j++) // for each outgoing move
			numTransitions += trans.get(s).get(j).size();

		// check compatibility
		if (checkCompatibility) {
			for (int s = 0; s < trans.size(); s++) { // for each state
			    numTransitions += trans.get(s).size();
				if (stateOwners.get(s) == 1) { // if a player 1 state
					int cb = controlled_by.get(s); // controlling component		    
					// records for which local actions at the controlling component a transition has been added already
					Set<Integer> local_trans_present = new HashSet<Integer>();
					for (int j = 0; j < trans.get(s).size(); j++) { // for each outgoing move
						int la = l_action.get(s).get(j).get(cb); // local action at controlling component
						if (local_trans_present.contains(la))
							throw new PrismException("Games not compatible: more than one transition enabled in state " + statesList.get(s).toString());
						else
							local_trans_present.add(la);
					}
				}
			}
		}
	}

	public boolean movePresent(Object action, List<Object> as, List<Distribution> distrs, List<List<Distribution>> dss)
	{
		outer: for (int j = 0; j < dss.size(); j++) {
			// now compare pairwise
			if (dss.get(j).size() != distrs.size()) { // distributions disagree on size
				continue outer;
			}
			Object asj = as.get(j);
			if ((action==null && asj != null) || (action!=null && !action.equals(asj))) { // actions don't agree
				continue outer;
			}
			// same size and same action
			inner: for (int i = 0; i < dss.get(j).size(); i++) {
				if (dss.get(j).get(i) == null) {
					if (distrs.get(i) != null) {
						continue outer; // not equal
					} else { // equal on this element
						continue inner;
					}
				}
				if (!dss.get(j).get(i).equals(distrs.get(i)))
					continue outer;
			}
			// fall through only if distributions are equal
			return true;
		}
		// fall through only if no distributions are equal
		return false;
	}

	// returns true if a transition has been added,
	// i.e. more states have become reachable
	private boolean addTransition(List<SMG> smgs, State s, State t, List<Distribution> distrs, Object action, int i, int n,
			List<List<List<Distribution>>> tempTrans) throws PrismException
	{
		boolean result = false;
		if (i == n) { // no more games to be looked at
			if (action != null) { // check (C1)
				// in state t is now the state to be added
				// sanity check: evaluate player and check (C4)
				boolean p1 = false;
				for (int ii = 0; ii < n; ii++) {
					if (smgs.get(ii).stateOwners.get((Integer) t.varValues[ii]) == 1) {
						// check (C4)
						if (distrs.get(ii) == null) {
							return false;
						}
						if (p1) {
							throw new PrismException(String.format("Condition (C5) violated: trying to compose two PLAYER 1 states: %s", t)); // already a P1 state present
						}
						p1 = true;
					}
				}

				// TRANSITIONS - NEED TO BE ADDED EVEN IF STATE ISN'T
				// note that t may not be present in the list, so it will be added
				int s_index = statesList.indexOf(s);
				if (!movePresent(action, actions.get(s_index), distrs, tempTrans.get(s_index))) { // not already present at s
					tempTrans.get(s_index).add(distrs);
					actions.get(s_index).add(action);
				}

				// ADD STATE
				// check whether t is duplicate (note that index n holds extra info, which is ignored)
				boolean duplicate = false;
				int t_index = -1; // index of duplicate state if found
				search_for_duplicate: for (State r : statesList) {
					if (r.compareTo(t) == 0) {
						duplicate = true;
						t_index = statesList.indexOf(r);
						break search_for_duplicate;
					}
				}
				if (!duplicate) { // state is unique and hence a new state needs to be added
					statesList.add(t);
					t_index = statesList.indexOf(t);

					stateOwners.add(p1 ? 1 : 2); // which player
					tempTrans.add(new ArrayList<List<Distribution>>()); // add new empty temporary transition list for state t
					trans.add(new ArrayList<Distribution>()); // add new empty transition list for state t
					actions.add(new ArrayList<Object>()); // add new empty action list for state t
				}
				return true;
			}
		} else { // assert: i < n ... more games to be looked at
			List<Distribution> ti = smgs.get(i).trans.get((Integer) s.varValues[i]);
			List<Object> ci = smgs.get(i).actions.get((Integer) s.varValues[i]);
			for (int j = 0; j < ti.size(); j++) { // go through all action-distribution pairs of game i
				// get action distribution pair j for game i
				Distribution di = ti.get(j);
				Object ai = ci.get(j);
				// tau-transitions are added directly
				if (TAU.equals(ai)) { // tau-transition
					// pass on if not null, and if this game doesn't block it
					if (action != null) {
						boolean in_alphabet = false;
						evaluate_action_alphabet: for (List<Object> actions_i : smgs.get(i).actions) {
							if (actions_i != null && actions_i.contains(action)) { // action in action alphabet of game i
								in_alphabet = true;
								break evaluate_action_alphabet;
							}
						}
						if (in_alphabet) { // blocked but in alphabet ... no transition, (C2) violated
							// just continue
						} else { // blocked and not in alphabet ... (C3) comes in
							List<Distribution> new_distrs = new ArrayList<Distribution>(distrs);
							new_distrs.add(null);
							result |= addTransition(smgs, s, t, new_distrs, action, i + 1, n, tempTrans);
						}
					}
					// or deal with tau
					if (di.size() > 1) {
						throw new PrismException("Game not in normal form, tau-transition not Dirac");
					}
					go_through_successors: for (Integer ri : di.getSupport()) { // go through successors of distribution
						State t_copy = new State(s);
						t_copy.varValues[i] = ri; // only this component moves

						// evaluate player
						boolean p1 = false;
						for (int ii = 0; ii < n; ii++) {
							if (smgs.get(ii).stateOwners.get((Integer) t_copy.varValues[ii]) == 1) {
								if (p1) {
									continue go_through_successors; // don't need to add states with several components controlled by P1 entered on a tau transition
								}
								p1 = true;
							}
						}
						List<Distribution> new_distrs = new ArrayList<Distribution>();
						for (int nulls = 0; nulls < distrs.size(); nulls++) {
							new_distrs.add(null);
						}
						new_distrs.add(di);
						result |= addTransition(smgs, s, t_copy, new_distrs, ai, n, n, tempTrans); // add t_copy
					}
				}
				// deal with other cases
				if (action == null) { // action still not specified
					// either pass
					if (smgs.get(i).stateOwners.get(((Integer) s.varValues[i])).intValue() != 1) { // P1 must be involved, (C4)
						List<Distribution> new_distrs = new ArrayList<Distribution>(distrs);
						new_distrs.add(null);
						result |= addTransition(smgs, s, s, new_distrs, action, i + 1, n, tempTrans);
					}
					// or fix action to ai and restart i
					if (!TAU.equals(ai)) {
						result |= addTransition(smgs, s, s, new ArrayList<Distribution>(), ai, 0, n, tempTrans);
					}
				} else if (action.equals(ai)) { // synchronised
					for (Integer ri : di.getSupport()) { // go through successors of distribution
						State t_copy = new State(t);
						t_copy.varValues[i] = ri; // involve this component (C2), (C1)
						List<Distribution> new_distrs = new ArrayList<Distribution>(distrs);
						new_distrs.add(di);
						result |= addTransition(smgs, s, t_copy, new_distrs, action, i + 1, n, tempTrans);
					}
				} else { // (ai != action) blocked
					boolean in_alphabet = false;
					evaluate_action_alphabet: for (List<Object> actions_i : smgs.get(i).actions) {
						if (actions_i != null && actions_i.contains(action)) { // action in action alphabet of game i
							in_alphabet = true;
							break evaluate_action_alphabet;
						}
					}
					if (in_alphabet) { // blocked but in alphabet ... no transition, (C2) violated
						// just continue
					} else { // blocked and not in alphabet ... (C3) comes in
						// stay put in component i
						List<Distribution> new_distrs = new ArrayList<Distribution>(distrs);
						new_distrs.add(null);
						result |= addTransition(smgs, s, t, new_distrs, action, i + 1, n, tempTrans);
					}
				}
			}
		}
		return result;
	}

	/**
	 * Add a new (player 1) state and return its index.
	 */
	@Override
	public int addState()
	{
		return addState(0);
	}

	/**
	 * Add multiple new (player 1) states.
	 */
	@Override
	public void addStates(int numToAdd)
	{
		for (int i = 0; i < numToAdd; i++) {
			addState();
		}
	}

	/**
	 * Set the info about players, i.e., the (integer) index and name of each one.
	 * This is given as a mapping from indices to names.
	 * Indices can be arbitrary and do not need to be contiguous.
	 * Names are optional and can be null or "" if undefined,
	 * but all normally defined names should be unique.
	 */
	public void setPlayerInfo(Map<Integer, String> playerNames)
	{
		this.playerNames = new HashMap<Integer, String>(playerNames);
	}

	/**
	 * Copy player info from another SMG.
	 */
	public void copyPlayerInfo(SMG smg)
	{
		setPlayerInfo(smg.playerNames);
	}

	/**
	 * Set a coalition of players for this SMG
	 * (which effectively makes it an STPG with player 1 representing the coalition and 2 the rest).
	 * Pass null to remove any coalition info from this SMG.
	 * 
	 * @param coalition Coalition info object 
	 */
	public void setCoalition(Coalition coalition) throws PrismException
	{
		// Clear info if coalition is null
		if (coalition == null) {
			coalitionPlayerMap = null;
			return;
		}
		
		// This is first time we really need the {@code playerNames} info,
		// so if it has not been set, we tried to create it based on {@code stateOwners} 
		if (playerNames.isEmpty()) {
			for (int i = 0; i < numStates; i++) {
				int p = stateOwners.get(i);
				if (!playerNames.containsKey(p)) {
					playerNames.put(p, null);
				}
			}
		}
		
		// Find max player index
		int maxIndex = 0;
		for (int index : playerNames.keySet()) {
			maxIndex = Math.max(maxIndex, index);
		}
		
		// Construct mapping
		coalitionPlayerMap = new int[maxIndex + 1];
		for (int i = 0; i < maxIndex + 1; i++) {
			coalitionPlayerMap[i] = -1;
		}
		for (Entry<Integer, String> entry : playerNames.entrySet()) {
			int playerIndex = entry.getKey();
			boolean inCoalition = coalition.isPlayerIndexInCoalition(playerIndex, playerNames);
			// In coalition => player 1; not in coalition (or undefined) => player 2
			coalitionPlayerMap[playerIndex] = inCoalition ? 1 : 2;
		}
	}

	/**
	 * Copy coalition info from another SMG.
	 */
	public void copyCoalitionInfo(SMG smg)
	{
		coalitionPlayerMap = smg.coalitionPlayerMap == null ? null : smg.coalitionPlayerMap.clone();
	}

	/**
	 * Makes a half-deep (up to one reference level) copy of itself
	 */
	public SMG clone()
	{
		SMG smg = new SMG();
		smg.copyFrom(this);
		smg.actions = new ArrayList<List<Object>>(this.actions);
		smg.allowDupes = this.allowDupes;
		smg.maxNumDistrs = this.maxNumDistrs;
		smg.maxNumDistrsOk = this.maxNumDistrsOk;
		smg.numDistrs = this.numDistrs;
		smg.numTransitions = this.numTransitions;
		smg.stateOwners = new ArrayList<Integer>(this.stateOwners);
		smg.trans = new ArrayList<List<Distribution>>(this.trans);
		smg.controlled_by = new ArrayList<Integer>(this.getControlledBy());
		return smg;
	}

	// Accessors (for Model)

	@Override
	public ModelType getModelType()
	{
		return ModelType.SMG;
	}

	// Accessors (for STPG/SMG)

	//@Override
	public int getNumPlayers()
	{
		return playerNames.size();
	}
	
	@Override
	public int getPlayer(int s)
	{
		int playerIndex = stateOwners.get(s); 
		if (coalitionPlayerMap == null) {
			// No coalition: just return index
			return playerIndex;
		} else {
			// Coalition defined: look up if player 1 or 2
			// (note: undefined players are mapped to player 2)
			return playerIndex == -1 ? 2 : coalitionPlayerMap[playerIndex];
		}
	}

	@Override
	public void exportToDotFile(PrismLog out, BitSet mark, boolean states)
	{
		BitSet players = new BitSet(numStates);
		for (int i = 0; i < numStates; i++) {
			if (stateOwners.get(i) == 1) {
				players.set(i);
			}
		}
		super.exportToDotFile(out, mark, players, states);
	}

	/**
	 * take X^k and apply F(X^k)(s) for each state (cf. MFCS'13 and QEST'13)
	 *
	 * arguments:
	 * @param gaussSeidel Gauss Seidel update allowed
	 * @param rounding rounding enabled
	 * @param union_with_previous take union with previous Pareto set
	 * @param cut cut off everything that is strictly above the negative orthant (used for energy objectives)
	 * @param M maximum bound on Pareto sets (quantity is positive)
	 */
	public Pareto[] pMultiObjective(Pareto[] Xk, List<SMGRewards> rewards, boolean gaussSeidel,
					long baseline_accuracy, double[] biggest_reward,
					boolean[] cancel_computation,
					List<Pareto>[] stochasticStates, boolean rounding,
					boolean union_with_previous, boolean cut, long M)
	    throws PrismException
	{
		Pareto[] result = new Pareto[Xk.length];
		Pareto[] Yk = gaussSeidel ? null : new Pareto[Xk.length]; // if Gauss-Seidel, no memory allocation required
		System.arraycopy(Xk, 0, gaussSeidel ? result : Yk, 0, Xk.length); // if Gauss-Seidel, update result in-place
		// iterate for each state separately
		for (int s = 0; s < numStates; s++) {
			// first, check if cancelled
			if (cancel_computation[0]) {
				cancel_computation[0] = false; // reset
				throw new PrismException("Computation cancelled");
			}
			// initialize the polyhedra for the stochastic states of s
			List<Pareto> distPolys = new ArrayList<Pareto>(trans.get(s).size());
			// apply F to (X^k)(s)
			//double t0 = (double)System.nanoTime();
			result[s] = pMultiObjectiveSingle(s, gaussSeidel ? result : Yk, rewards, baseline_accuracy, biggest_reward, distPolys, rounding,
							  union_with_previous, cut, M);
			//System.out.printf("total: %f s\n", ((double) (System.nanoTime() - t0)) / 1e9);
			// store stochastic states if requested (by the reference being non-null)
			if (stochasticStates != null)
				stochasticStates[s] = distPolys;
		}

		// return X^{k+1}
		return result;
	}

    private Polyhedron round(Generator_System ngs, long baseline_accuracy, double[] biggest_reward, boolean energy_objective) throws PrismException
	{
		int n = biggest_reward.length;
		// accuracy
		long[] accuracy = new long[n];
		for (int i = 0; i < n; i++) {
		        long tmp_a = energy_objective ? baseline_accuracy : ((long) (((double) baseline_accuracy) / biggest_reward[i]));
			// sanity check to prevent overflow
			accuracy[i] = tmp_a < Long.MAX_VALUE && tmp_a > 0 ? tmp_a : Long.MAX_VALUE;
		}

		Generator_System new_ngs = new Generator_System();
		for (Generator ng : ngs) {
			if (ng.type() == Generator_Type.POINT) {
				Linear_Expression le = ng.linear_expression();
				Coefficient c = ng.divisor();
				Map<Variable, BigInteger> map = new HashMap<Variable, BigInteger>();
				PPLSupport.getCoefficientsFromLinearExpression(le, false, BigInteger.ONE, map);

				// new denominator at baseline accuracy
				Coefficient new_c = new Coefficient(BigInteger.valueOf(baseline_accuracy));

				// new linear expression
				Linear_Expression new_le;
				if (map.containsKey(null)) { // there is a coefficient without a variable
					if (map.get(null).compareTo(BigInteger.ZERO) != 0)
						throw new PrismException("Exception in Polyhedron presentation.");
					new_le = new Linear_Expression_Coefficient(new Coefficient(map.get(null)));
				} else {
					new_le = new Linear_Expression_Coefficient(new Coefficient(BigInteger.ZERO));
				}
				for (Variable k : map.keySet()) {
					if (k != null) {
						BigFraction round_test = new BigFraction(map.get(k), c.getBigInteger());
						long rounded = ((long) (Math.floor(round_test.doubleValue() * accuracy[(int)k.id()]) * baseline_accuracy / ((double) accuracy[(int)k.id()])));
						new_le = new_le.sum(new Linear_Expression_Times(new Coefficient(rounded), k));
					}
				}
				new_ngs.add(Generator.point(new_le, new_c));
			} else if (ng.type() == Generator_Type.RAY) {
				new_ngs.add(ng);
			}
		}

		Polyhedron result = new C_Polyhedron(new_ngs);
		// add zero dimensions if rounding deleted them
		if (result.space_dimension() != n) {
			result.add_space_dimensions_and_project(n - result.space_dimension());
		}

		return result;
	}

    protected Pareto stochasticState(int s, Distribution distr, int d, Pareto[] Xk, List<SMGRewards> rewards, double[] extra_rewards, boolean cut, long M)
			throws PrismException
	{
		int n = rewards.size();

		// the successors of the distribution d
		ArrayList<Integer> states = new ArrayList<Integer>(distr.keySet());
		int b = states.size();

		Pareto cp = null;
		if (b == 0) {
			throw new PrismException("Distribution " + s + ", " + d + " has no successors.");
		} else if (b == 1) {
			// distribution assigns 1 to first successor
			cp = Xk[states.get(0)];
		} else { // need to compute Minkowski sum
			Linear_Expression lhs, rhs;

			// first need to make sure probabilities add to one
			BigFraction[] probs = new BigFraction[b];
			BigFraction residual = BigFraction.ONE;
			int supdim = 0;
			for (Integer t : states) {
				BigFraction prob = new BigFraction(distr.get(t));
				probs[supdim] = prob;
				residual = residual.subtract(prob);
				supdim++;
			}
			probs[0] = probs[0].add(residual); // just add residual to first probability

			supdim = 0;
			for (Integer t : states) {
				C_Polyhedron p = new C_Polyhedron((C_Polyhedron) Xk[t].get()); // deep copy!
				p.add_space_dimensions_and_embed(b);
				BigFraction prob = probs[supdim];
				for (int i = 0; i < b; i++) {
					if (i == supdim) {
						lhs = new Linear_Expression_Times(new Coefficient(prob.getNumerator()), new Variable(n + i));
						rhs = new Linear_Expression_Coefficient(new Coefficient(prob.getDenominator()));
					} else {
						lhs = new Linear_Expression_Times(new Coefficient(BigInteger.ONE), new Variable(n + i));
						rhs = new Linear_Expression_Coefficient(new Coefficient(BigInteger.ZERO));
					}
					Constraint c = new Constraint(lhs, Relation_Symbol.EQUAL, rhs);
					p.add_constraint(c);
				}
				if (cp == null)
					cp = new Pareto(p);
				else
					cp.get().upper_bound_assign(p);

				supdim++;
			}

			supdim = 0;
			for (Integer t : states) {
				lhs = new Linear_Expression_Times(new Coefficient(BigInteger.ONE), new Variable(n + supdim));
				rhs = new Linear_Expression_Coefficient(new Coefficient(BigInteger.ONE));
				Constraint c = new Constraint(lhs, Relation_Symbol.EQUAL, rhs);
				cp.get().add_constraint(c);
				supdim++;
			}

			// project away the unneccessary dimensions
			//double t0 = (double)System.nanoTime();
			cp.get().remove_higher_space_dimensions(n);
			//System.out.printf("rhsd: %f s\n", ((double) (System.nanoTime() - t0)) / 1e9);
			// now in cp have the Minkowski sum for that particular distribution d
		}

		// add rewards
		Polyhedron Yk1 = PPLSupport.add_rewards(cp.get(), s, d, rewards, extra_rewards);

		// cut everything but the negative orthant bounded by -M
		if (cut) PPLSupport.cutBox(Yk1, M);

		// add transition rewards and return polyhedron
		return new Pareto(Yk1);
	}

	// distPolys will hold the polyhedra of the stochastic states
	private Pareto pMultiObjectiveSingle(int s, Pareto[] Xk, List<SMGRewards> rewards, long baseline_accuracy, double[] biggest_reward, List<Pareto> distPolys,
					     boolean rounding, boolean union_with_previous, boolean cut, long M) throws PrismException
	{
		int n = rewards.size();

		// distributions of successor states
		List<Distribution> dists = trans.get(s);

		// ------------------------------------------------------------------------------
		// STOCHASTIC STATE OPERATIONS

		// step through all stochastic successors of s
		int d = 0;
		for (Distribution distr : dists) {
			// add polyhedron to the list of polyhedra in the successors of s
		        distPolys.add(stochasticState(s, distr, d, Xk, rewards, null, cut, M));
			d++;
		}

		// ------------------------------------------------------------------------------
		// PLAYER ONE AND PLAYER TWO OPERATIONS

		// Xk1s holds the polyhedron of state s
		// need deep copy here because want to retain Minkowski sums
		Polyhedron Xk1s;
		if (distPolys.size() > 0) {
			if (getPlayer(s) == 1) {
			        // Player 1
			        Xk1s = new C_Polyhedron(distPolys.get(0).get().generators());
				int cp_start = 0;
			        get_first_cp: for (cp_start = 0; cp_start < distPolys.size(); cp_start++) {
				    if(!distPolys.get(cp_start).get().is_empty()) {
					Xk1s = new C_Polyhedron(distPolys.get(cp_start).get().generators());
					break get_first_cp;
				    }
				}
				if(distPolys.get(0).get().space_dimension() > Xk1s.space_dimension())
				        Xk1s.add_space_dimensions_and_project(distPolys.get(0).get().space_dimension() - Xk1s.space_dimension());

				for (int cp_i = cp_start+1; cp_i < distPolys.size(); cp_i++) {
				    if(!distPolys.get(cp_i).get().is_empty())
					Xk1s.upper_bound_assign(distPolys.get(cp_i).get());
				}
			} else {
				// Player 2
 			        Xk1s = new C_Polyhedron(distPolys.get(0).get().generators());
				Xk1s.add_space_dimensions_and_project(distPolys.get(0).get().space_dimension() - Xk1s.space_dimension());
				for (int cp_i = 1; cp_i < distPolys.size(); cp_i++) {
				    if(!Xk1s.is_empty() && !distPolys.get(cp_i).get().is_empty())
					Xk1s.intersection_assign(distPolys.get(cp_i).get());
				    else if(!Xk1s.is_empty()) // now the other polyhedron must be empty
					Xk1s = new C_Polyhedron(distPolys.get(cp_i).get().generators());
				}
			}
		} else { // deadlock
		        Xk1s = Xk[s].get(); //new C_Polyhedron(new Generator_System()); // empty
		}

		// ------------------------------------------------------------------------------
		// ADD STATE REWARDS
		Xk1s = PPLSupport.add_rewards(Xk1s, s, Integer.MIN_VALUE, rewards, null);

		// ------------------------------------------------------------------------------
		// ROUNDING (if required)
		if (rounding) Xk1s = round(Xk1s.generators(), baseline_accuracy, biggest_reward, cut);

		// ------------------------------------------------------------------------------
		// CLEAN UP: UNION WITH PREVIOUS RESULT OR CUT, MINIMIZE REPRESENTATION, DIMENSIONALITY

		// union with previous result (after rounding)
		if (rounding && union_with_previous) Xk1s.upper_bound_assign(Xk[s].get());
		// cut everything but the negative orthant bounded by -M
		if (cut) PPLSupport.cutBox(Xk1s, M);

		// minimize representation
		Xk1s = new C_Polyhedron(Xk1s.minimized_generators());

		// add zero dimensions if minimization deleted them
		if (Xk1s.space_dimension() != n)
			Xk1s.add_space_dimensions_and_project(n - Xk1s.space_dimension());

		return new Pareto(Xk1s);
	}

	// Standard methods

	@Override
	public String toString()
	{
		int i, j, n;
		Object o;
		String s = "";
		s = "[ ";
		for (i = 0; i < numStates; i++) {
			if (i > 0)
				s += ", ";
			if (statesList.size() > i)
				s += i + "(P-" + stateOwners.get(i) + " " + statesList.get(i) + "): ";
			else
				s += i + "(P-" + stateOwners.get(i) + "): ";
			s += "[";
			n = getNumChoices(i);
			for (j = 0; j < n; j++) {
				if (j > 0)
					s += ",";
				o = getAction(i, j);
				if (o != null)
					s += o + ":";
				s += trans.get(i).get(j);
			}
			s += "]";
		}
		s += " ]\n";
		return s;
	}
	// which component controls the current state
	public int controlledBy(int i)
	{
		return getControlledBy().get(i);
	}

	public List<Integer> getControlledBy()
	{
		if (controlled_by == null) { // fill with default
			controlled_by = new ArrayList<Integer>(numStates);
			for (int i = 0; i < numStates; i++) {
				if (stateOwners.get(i).intValue() == 1) {
					controlled_by.add(0); // P1 is component zero by default
				} else {
					controlled_by.add(-1); // P2 is irrelevant
				}
			}
		}
		return controlled_by;
	}

	// returns the number of components
	public int getNumComponents()
	{
		return getLocalState().get(0).size();
	}

	// mapping from global to local action
	public List<List<List<Integer>>> getLAction()
	{
		if (l_action == null) { // fill with default
			l_action = new ArrayList<List<List<Integer>>>(numStates);
			for (int i = 0; i < numStates; i++) {
				l_action.add(new ArrayList<List<Integer>>(trans.get(i).size()));
				for (int a = 0; a < trans.get(i).size(); a++) {
					l_action.get(i).add(new ArrayList<Integer>(1));
					l_action.get(i).get(a).add(a); // identity mapping is default
				}
			}
		}
		return l_action;
	}

	// mapping from global states to local state
	public List<List<Integer>> getLocalState()
	{
		if (local_state == null) { // fill with default
			local_state = new ArrayList<List<Integer>>(numStates);
			for (int i = 0; i < numStates; i++) {
				List<Integer> gs = new ArrayList<Integer>(1); // one component
				gs.add(i); // identity mapping is default
				local_state.add(gs);
			}
		}
		return local_state;
	}

	public boolean deadlocksAllowed()
	{
		return false;
	}

}
