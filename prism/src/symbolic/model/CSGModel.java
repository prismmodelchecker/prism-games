//==============================================================================
//
//	Copyright (c) 2002-
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

package symbolic.model;

import java.io.File;
import java.io.FileNotFoundException;

import io.ModelExportOptions;
import jdd.JDDVars;
import jdd.JDDNode;
import parser.VarList;
import prism.ModelType;
import prism.PlayerInfo;
import prism.PlayerInfoOwner;
import prism.PrismException;
import prism.PrismNotSupportedException;

/**
 * Class for symbolic (BDD/MTBDD-based) representation of a CSG (concurrent stochastic game),
 * as built by {@code symbolic.build.CSG2MTBDD}.
 *
 * <p>Unlike {@link GamesModel} (SMG's one-hot "whose turn is it" turn-based encoding), a CSG's
 * nondeterminism dimension is a per-player block of action variables ({@code Act_p}): every
 * player chooses simultaneously, every step, so there is no single one-hot "current player"
 * cube to expose the way {@link GamesModel} does. {@link #actDDVars} stores these per-player
 * blocks instead; {@link NondetModel#getAllDDNondetVars()} (inherited) is their union, exactly
 * as {@code CSG2MTBDD.allActDDVars()} computes it during the build.
 *
 * <p>The transition matrix this class is constructed with is expected to already be gated by
 * {@code Enabled} -- i.e. zero outside the joint action tuples actually legally available at a
 * state -- rather than the raw, everywhere-defined {@code T} that the module-identity
 * fallback construction produces internally (it makes every module's own row sum to 1
 * unconditionally, purely as a build-time convenience for the per-module product; {@code T}
 * alone is not the model's real transition function). Folding {@code Enabled} into {@code trans}
 * up front, rather than keeping it as a separate stored BDD, is what lets this class inherit
 * {@link NondetModel}'s standard reachability/deadlock/model-checking machinery (which knows
 * nothing about {@code Enabled}) unmodified -- using the raw, unrestricted {@code T} directly
 * here would let the module-identity fallback's rows leak into reachability as if they were real,
 * legally-available transitions. See {@code CSG2MTBDD.translate()}, the only intended way to
 * construct one of these.
 *
 * <p>No reward-structure support yet -- {@link #getModelType()} models are always constructed
 * with zero reward structs. No qualitative or game-theoretic model-checking algorithms
 * (sure/almost/limit, rPATL) exist yet either -- see {@code symbolic.comp.CSGModelChecker}'s
 * own doc.
 */
public class CSGModel extends NondetModel implements PlayerInfoOwner
{
	/** Per-player Act_p DD variable blocks. Their union is allDDNondetVars. */
	private JDDVars[] actDDVars;

	/** Player names, for {@link PlayerInfoOwner} / rPATL-style {@code <<...>>} property syntax. */
	protected PlayerInfo playerInfo;

	public CSGModel(JDDNode trans, JDDNode start, JDDVars allDDRowVars, JDDVars allDDColVars, JDDVars allDDNondetVars, ModelVariablesDD modelVariables,
					VarList varList, JDDVars[] varDDRowVars, JDDVars[] varDDColVars, JDDVars[] actDDVars, PlayerInfo playerInfo)
	{
		super(trans, start, allDDRowVars, allDDColVars, allDDNondetVars, modelVariables, varList, varDDRowVars, varDDColVars);
		this.actDDVars = actDDVars;
		this.playerInfo = playerInfo;
	}

	// Accessors (for Model)

	@Override
	public ModelType getModelType()
	{
		return ModelType.CSG;
	}

	@Override
	public void clear()
	{
		super.clear();
		JDDVars.derefAllArray(actDDVars);
	}

	@Override
	public void exportToFile(int exportType, boolean explicit, File file, int precision) throws FileNotFoundException, PrismException
	{
		throw new PrismNotSupportedException("Symbolic engine does not support export of " + getModelType() + "s");
	}

	@Override
	public void exportToFile(File file, ModelExportOptions exportOptions) throws FileNotFoundException, PrismException
	{
		throw new PrismNotSupportedException("Symbolic engine does not support export of " + getModelType() + "s");
	}

	// Accessors (for PlayerInfoOwner)

	@Override
	public PlayerInfo getPlayerInfo()
	{
		return playerInfo;
	}

	// Accessors (for CSG)

	/** DD variables encoding player {@code p}'s own action, {@code Act_p}. */
	public JDDVars getActDDVars(int p)
	{
		return actDDVars[p];
	}
}
