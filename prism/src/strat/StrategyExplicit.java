//==============================================================================
//	
//	Copyright (c) 2022-
//	Authors:
//	* Dave Parker <d.a.parker@cs.bham.ac.uk> (University of Birmingham)
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

package strat;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.List;
import java.util.function.Consumer;

import explicit.InducedStateMap;
import explicit.Model;
import explicit.NondetModel;
import explicit.TurnBasedGame;
import io.ModelExportFormat;
import prism.Evaluator;
import prism.PrismException;
import prism.PrismLog;

/**
 * Base class for implementations of Strategy associated with an explicit engine model.
 */
public abstract class StrategyExplicit<Value> extends StrategyWithStates<Value>
{
	// Model associated with the strategy
	protected NondetModel<Value> model;

	// For a game, the players (0-indexed) that the strategy decides for, the coalitions
	// they are partitioned into, and the states in which they make the choice
	// (all null if the strategy decides for all players of the model, as one coalition;
	// the states are also null if the choices in all states are made by these players)
	protected BitSet players = null;
	protected List<BitSet> coalitions = null;
	protected BitSet playerStates = null;

	public StrategyExplicit(NondetModel<Value> model)
	{
		this.model = model;
		if (model.getStatesList() != null) {
			setStateLookUp(state -> model.getStatesList().indexOf(state));
		} else {
			setStateLookUp(state -> -1);
		}
	}

	/**
	 * For a game, specify the players that this strategy decides for, as a list of coalitions
	 * (see {@link #coalitions()}); the players are their union. In the states of other players,
	 * its decision is {@link #OTHER_PLAYER}: the strategy has no decision there (and any choice
	 * information stored for them is not used). For a turn-based game, these are the
	 * states owned by other players; states owned by no player keep their decisions.
	 * Otherwise (e.g., for a model with one player), the strategy decides in all states.
	 * @param coalitions The coalitions of players (0-indexed) that this strategy decides for
	 */
	public void setPlayers(List<BitSet> coalitions)
	{
		this.coalitions = new ArrayList<>();
		this.players = new BitSet();
		for (BitSet coalition : coalitions) {
			this.coalitions.add((BitSet) coalition.clone());
			this.players.or(coalition);
		}
		// For a turn-based game, find the states owned by these players (or by no player)
		playerStates = null;
		if (model instanceof TurnBasedGame) {
			TurnBasedGame game = (TurnBasedGame) model;
			playerStates = new BitSet();
			int n = model.getNumStates();
			for (int s = 0; s < n; s++) {
				int p = game.getPlayerIgnoringCoalition(s);
				if (p == -1 || players.get(p)) {
					playerStates.set(s);
				}
			}
		}
	}

	/**
	 * Is the choice in state s made by one of the players that this strategy decides for?
	 * If not, the strategy's decision is {@link #OTHER_PLAYER}.
	 */
	protected boolean isPlayerState(int s)
	{
		return playerStates == null || playerStates.get(s);
	}

	@Override
	public BitSet players()
	{
		if (players == null) {
			// By default, all players of the model (or just player 1, index 0, if it has none)
			BitSet all = new BitSet();
			all.set(0, Math.max(1, getNumModelPlayers()));
			return all;
		}
		// Return a copy: BitSets are mutable, and a caller modifying the stored set would change
		// the strategy, and make it inconsistent with its coalitions and player states
		return (BitSet) players.clone();
	}

	@Override
	public List<BitSet> coalitions()
	{
		if (coalitions == null) {
			return new ArrayList<>(Collections.singletonList(players()));
		}
		// Return copies, as for players(), since BitSets are mutable
		List<BitSet> copies = new ArrayList<>(coalitions.size());
		for (BitSet coalition : coalitions) {
			copies.add((BitSet) coalition.clone());
		}
		return copies;
	}

	@Override
	public NondetModel<Value> getModel()
	{
		return model;
	}

	@Override
	public Evaluator<Value> getEvaluator()
	{
		return model.getEvaluator();
	}

	@Override
	public Model<Value> constructInducedModel(StrategyExportOptions options) throws PrismException
	{
		return constructInducedModel(options, null);
	}

	/**
	 * Construct the model induced by this strategy.
	 * If {@code stateMapConsumer} is non-null, the map from states of the induced model
	 * to states of the original model (and strategy memory) is also stored
	 * and then passed to it.
	 */
	public abstract Model<Value> constructInducedModel(StrategyExportOptions options, Consumer<InducedStateMap> stateMapConsumer) throws PrismException;

	@Override
	public void exportInducedModel(PrismLog out, StrategyExportOptions options) throws PrismException
	{
		constructInducedModel(options).export(out, options.getInducedModelExportOptions());
	}

	@Override
	public void exportDotFile(PrismLog out, StrategyExportOptions options) throws PrismException
	{
		// Usually, this has been set up so we can treat it an induced model export
		// in Dot format; if so, just delegate to the induced model export method
		if (options.getInducedModelExportOptions() != null && options.getInducedModelExportOptions().getFormat() == ModelExportFormat.DOT) {
			exportInducedModel(out, options);
		}
		// Otherwise call Dot export directly, passing a few options
		else {
			constructInducedModel(options).exportToDotFile(out, null, options.getShowStates(), options.getModelPrecision());
		}
	}
}
