package symbolic.build;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import jdd.JDD;
import jdd.JDDNode;
import jdd.JDDVars;
import mtbdd.PrismMTBDD;
import parser.VarList;
import parser.Values;
import parser.ast.CSGPlayerActions;
import parser.ast.Command;
import parser.ast.Expression;
import parser.ast.Module;
import parser.ast.ModulesFile;
import parser.ast.Update;
import parser.ast.Updates;
import prism.ModelType;
import prism.PlayerInfo;
import prism.Prism;
import prism.PrismException;
import prism.PrismLangException;
import prism.PrismNativeComponent;
import prism.PrismNotSupportedException;
import prism.PrismUtils;
import symbolic.comp.CSGStateModelChecker;
import symbolic.model.CSGModel;
import symbolic.model.Model;
import symbolic.model.ModelVariablesDD;

/**
 * Symbolic (ADD/BDD) builder for CSG models.
 *
 * <p><b>Not a subclass of {@link Modules2MTBDD}.</b> The existing {@code *2MTBDD} builders in
 * this package ({@link Modules2MTBDD}, {@link ExplicitFiles2MTBDD}, {@link
 * ModelGenerator2MTBDD}) are siblings, each extending {@link PrismNativeComponent} directly,
 * not each other, and {@code Modules2MTBDD} has no extension points built for subclassing. It
 * also has no useful CSG behaviour to inherit: its nondeterminism combination
 * ({@code combineCommandsNondet}) and SMG's one-hot turn encoding ({@code buildDDGame}) are
 * both incompatible with CSG's simultaneous per-player action structure.
 * What genuinely is reused is called out inline below: {@link ModelVariablesDD} for variable
 * allocation, and {@link CSGStateModelChecker} (a small subclass of {@code StateModelChecker},
 * itself reused as-is) for guard/expression translation.
 *
 * <h2>Scope of this first slice</h2>
 * Builds the state and {@code Act_p} ADD variables, the per-command action-list BDD
 * (see {@link #buildActionListBDD}), the transition ADD {@code T} for player-owned and
 * independent modules, the {@code Enabled} BDD (idle-action support) that
 * gates which joint action tuples are legal, and reachability restricted to it.
 * <b>Deliberately not yet implemented</b>, and not silently assumed correct -- flagged at the
 * relevant point below, not just here:
 * <ul>
 * <li>reward structures.</li>
 * </ul>
 * Idle-action handling and the module-identity fallback ARE
 * implemented -- see {@link #translateModule} and {@link #buildEnabled}. So is local
 * rectangularity (including the precedence/"shadowing" resolution it turned out
 * {@code translateModule} itself needed -- see {@link #buildEffectiveMatchSets}), global
 * completeness (the global, multi-command companion to local rectangularity, once the
 * row-sum formulation stopped being a useful proxy for it -- see
 * {@link #checkGlobalCompleteness}) and primed-variable dependency acyclicity
 * ({@link #checkPrimedVarAcyclicity}). None of these are restricted to models like
 * {@code unilateral.prism}/{@code idle_forced.prism} that happen not to need them.
 *
 * <h2>CUDD reference-counting lifecycle</h2>
 * CUDD nodes are refcounted, not garbage-collected, so this class follows the same
 * two-tier discipline {@link Modules2MTBDD#translate} and {@code symbolic.model.ModelSymbolic}
 * already use elsewhere in this codebase (worth restating explicitly, since getting this
 * wrong is exactly what silently exhausts memory across a batch of builds rather than
 * failing loudly on the first one): DDs that exist only as build-time scaffolding (the
 * per-variable identity/range ADDs, {@code range} itself, the per-action enabled-action BDDs)
 * are derefed automatically and unconditionally at the end of {@link #build}, success or
 * failure, exactly like {@code Modules2MTBDD.cleanup()} and its
 * {@code expr2mtbdd.clearDummyModel()} call (this class makes the identical call, for the
 * identical reason -- {@link CSGStateModelChecker}'s abbreviated constructor, like
 * {@code Modules2MTBDD}'s own use of it, creates a throwaway dummy {@code ProbModel} purely
 * to have somewhere to hang expression translation, and that dummy model's DDs are never
 * otherwise released). The final artifacts -- {@link #getTrans}, {@link #getEnabled},
 * {@link #getStart}, {@link #getReach}, and the DD variables needed to interpret them -- are
 * NOT touched by that automatic cleanup, since build() can't know when the caller is done
 * with them; call {@link #clear} explicitly once finished with this object, the same
 * "build, use, then clear()" lifecycle {@code ModelSymbolic.clear()} follows. This matters
 * most exactly where it's easy to overlook: a loop that builds many models in one
 * JVM/CUDD session (e.g. batch experiments) will leak every one of them for the life of
 * that session if {@link #clear} is never called.
 */
public class CSG2MTBDD extends PrismNativeComponent
{
	private final ModulesFile modulesFile;

	// Basic model info
	private VarList varList;
	private int numVars;
	private int numModules;
	private int numPlayers;
	private Values constantValues;

	// State variable allocation
	private ModelVariablesDD modelVariables;
	private JDDVars allDDRowVars;
	private JDDVars allDDColVars;
	private JDDVars[] varDDRowVars;
	private JDDVars[] varDDColVars;
	private JDDNode[] varIdentities;
	private JDDNode[] varColRangeDDs;
	private JDDNode range;

	// Action ownership and coding
	private CSGPlayerActions playerActions;
	/** Idle (⊥_p) is always code 0 within its player's Act_p block, regardless of use. */
	private static final int IDLE_CODE = 0;
	/** Per player: map from (globally unique) action index to its 0-based code within that player's Act_p block. */
	private Map<Integer, Integer>[] playerActionCode;
	/** Per player: how many named (non-idle) actions that player owns. */
	private int[] playerNumNamedActions;

	// Act_p variable allocation
	private JDDVars[] actDDVars;

	// Expression translation
	private CSGStateModelChecker expr2mtbdd;

	// Scratch state for translateUpdate/translateUpdateElement (used by the identity-padding step)
	private boolean[] varsUsed;

	// Enabled actions: per (globally unique, *primary-labelled*) action index,
	// the state-only BDD of "some command owned by this action's player, primarily labelled
	// with it, has a satisfied guard" -- i.e. the symbolic analogue of Updater's moves[p]
	// construction (built purely from each active command's *primary* synch index, per
	// module -- see CSGPlayerActions/Updater.calculateEnabledCommands).
	private Map<Integer, JDDNode> enabledByActionIndex;

	// Global completeness: per (globally unique, *primary-labelled*) action
	// index, the state x joint-action BDD "PrimaryCovered_a" -- OR, across EVERY module in the
	// model (unlike enabledByActionIndex, not restricted to player-owned ones -- an independent
	// module may also carry a command primarily labelled with a player's own action, and its
	// coverage counts too), of that command's *effective* match set (shadowing-resolved
	// Guard_κ ∧ A_κ). See checkGlobalCompleteness.
	private Map<Integer, JDDNode> primaryCoverageByActionIndex;

	// Result
	private JDDNode trans;
	private JDDNode enabled;
	private JDDNode start;
	private JDDNode reach;

	public CSG2MTBDD(Prism prism, ModulesFile mf) throws PrismException
	{
		super(prism);
		if (mf.getModelType() != ModelType.CSG) {
			throw new PrismException("CSG2MTBDD only supports CSG models (got " + mf.getModelType() + ")");
		}
		this.modulesFile = mf;
	}

	/**
	 * Build the model. After this returns, {@link #getTrans()}, {@link #getReach()} and
	 * friends are available. See the class doc for what is and isn't implemented yet, and for
	 * this class's CUDD reference-counting lifecycle in particular -- in short, call
	 * {@link #clear()} once done with the built model.
	 */
	public void build() throws PrismException
	{
		try {
			varList = modulesFile.createVarList();
			if (modulesFile.containsUnboundedVariables()) {
				throw new PrismNotSupportedException("Cannot build a model that contains a variable with unbounded range");
			}
			numVars = varList.getNumVars();
			constantValues = modulesFile.getConstantValues();
			numModules = modulesFile.getNumModules();
			numPlayers = modulesFile.getNumPlayers();

			// Ownership, shared with the explicit engine via CSGPlayerActions -- see that
			// class's own doc for why this must not be a second, independently-maintained copy.
			playerActions = new CSGPlayerActions(modulesFile);
			assignActionCodes();

			allocateDDVars();
			buildIdentitiesAndRanges();

			expr2mtbdd = new CSGStateModelChecker(prism, varList, allDDRowVars, varDDRowVars, varDDColVars, constantValues);

			// Sequencing: local rectangularity (inline in translateModules, below) and
			// primed-variable acyclicity need only per-command guards/action-lists, not the
			// assembled T -- so they run before anything else, ahead of reachability pruning.
			checkPrimedVarAcyclicity();

			buildEnabledActions();
			enabled = buildEnabled();
			translateModules();
			checkGlobalCompleteness();
			buildInitialStates();
			buildReachability();
		} finally {
			// Unconditional (success or failure), mirroring Modules2MTBDD.cleanup() -- see the
			// class doc's "CUDD reference-counting lifecycle" section for why this is split
			// from clear().
			cleanupBuildScaffolding();
		}
	}

	/**
	 * Deref DDs that were only ever build-time scaffolding: {@link #varIdentities},
	 * {@link #varColRangeDDs}, {@link #range} and {@link #enabledByActionIndex}'s values are not
	 * part of what any getter exposes, and {@link CSGStateModelChecker}'s dummy model
	 * ({@code expr2mtbdd.clearDummyModel()}, exactly {@code Modules2MTBDD}'s own call) is pure
	 * translation scaffolding with no life beyond the build. Called from {@link #build}'s own
	 * {@code finally} block, not by {@link #clear}, since these need releasing whether or not
	 * the build actually succeeded, unlike the final-artifact DDs {@link #clear} handles.
	 */
	private void cleanupBuildScaffolding()
	{
		if (varIdentities != null) {
			JDD.DerefArrayNonNull(varIdentities);
		}
		if (varColRangeDDs != null) {
			JDD.DerefArrayNonNull(varColRangeDDs);
		}
		JDD.DerefNonNull(range);
		if (enabledByActionIndex != null) {
			for (JDDNode dd : enabledByActionIndex.values()) {
				JDD.DerefNonNull(dd);
			}
			enabledByActionIndex.clear();
		}
		if (primaryCoverageByActionIndex != null) {
			for (JDDNode dd : primaryCoverageByActionIndex.values()) {
				JDD.DerefNonNull(dd);
			}
			primaryCoverageByActionIndex.clear();
		}
		if (expr2mtbdd != null) {
			expr2mtbdd.clearDummyModel();
		}
	}

	/**
	 * Release every CUDD reference the built model itself holds: {@link #getTrans()},
	 * {@link #getEnabled()}, {@link #getStart()}, {@link #getReach()}, and every DD variable
	 * this builder allocated ({@code Act_p} included). Call this once, when completely done
	 * with this object -- e.g. after reading off {@link #getNumReachableStates()}/
	 * {@link #getNumTransitions()} -- mirroring {@code symbolic.model.ModelSymbolic}'s own
	 * "build, use, then clear()" lifecycle. Not called automatically by {@link #build()} (see
	 * the class doc): only the caller knows when it is actually finished with the result.
	 * Essential for anything that builds many {@code CSG2MTBDD} instances within one
	 * JVM/CUDD session (e.g. batch experiments) -- CUDD nodes are refcounted, not
	 * garbage-collected, so an instance that is simply dropped without calling this leaks
	 * every DD it built for the life of that CUDD manager.
	 */
	public void clear()
	{
		JDD.DerefNonNull(trans);
		JDD.DerefNonNull(enabled);
		JDD.DerefNonNull(start);
		JDD.DerefNonNull(reach);
		if (varDDRowVars != null) {
			JDDVars.derefAllArray(varDDRowVars);
		}
		if (varDDColVars != null) {
			JDDVars.derefAllArray(varDDColVars);
		}
		if (actDDVars != null) {
			JDDVars.derefAllArray(actDDVars);
		}
		if (allDDRowVars != null) {
			allDDRowVars.derefAll();
		}
		if (allDDColVars != null) {
			allDDColVars.derefAll();
		}
		if (modelVariables != null) {
			modelVariables.clear();
		}
	}

	/**
	 * Higher-level entry point, mirroring {@code Modules2MTBDD.translate()}'s contract: call
	 * {@link #build()} and package the result as a proper {@link symbolic.model.Model}
	 * ({@link CSGModel}) usable by the rest of PRISM's symbolic-engine infrastructure (the
	 * standard {@code -mtbdd} build dispatch in {@code Prism.doBuildModel}, and, via
	 * {@code StateModelChecker.createModelChecker}, {@code CSGModelChecker}), rather than this
	 * class's own raw {@link #getTrans()}/{@link #getEnabled()}/{@link #getReach()}/
	 * {@link #getStart()} accessors.
	 *
	 * <p>{@code Enabled} is folded directly into the model's transition matrix
	 * ({@code T x Enabled}, zero outside legally-available joint action tuples) rather than kept
	 * as a separate stored BDD -- see {@link CSGModel}'s own doc for why: it is what lets the
	 * returned model reuse {@code NondetModel}'s stock reachability/deadlock machinery
	 * unmodified, with no CSG-specific override needed. This build's own already-Enabled-restricted
	 * {@link #reach} (from {@link #buildReachability}) is reused directly via
	 * {@code CSGModel.setReach} rather than recomputed -- it is already exactly what the
	 * model's generic {@code doReachability()} would (redundantly) derive from the gated trans.
	 *
	 * <p>No reward structures yet -- the returned model always has zero
	 * reward structs.
	 *
	 * <p><b>Ownership note:</b> after calling this, do not also call {@link #build()},
	 * {@link #clear()}, or the raw {@code get*} accessors on this builder -- every DD this
	 * builder held is transferred to the returned {@link CSGModel}, which owns their lifecycle
	 * from here (release via the model's own {@code clear()}, not this class's). On failure
	 * (either {@link #build()} itself throwing, or one of the wrapping steps below), this
	 * releases whatever was already built instead of leaking it -- mirroring
	 * {@code Modules2MTBDD.translate()}'s own {@code catch (Exception e) { if (model != null)
	 * model.clear(); ...}} pattern: {@link #clear()} (this class's own) if the failure happened
	 * before a {@link CSGModel} existed to own the DDs, the model's own {@code clear()} if it
	 * happened after.
	 */
	public Model translate() throws PrismException
	{
		CSGModel model = null;
		try {
			build();

			JDDNode gatedTrans = JDD.Apply(JDD.TIMES, trans, enabled);
			trans = null;
			enabled = null;

			PlayerInfo playerInfo = new PlayerInfo();
			for (int p = 0; p < numPlayers; p++) {
				playerInfo.addPlayer(modulesFile.getPlayer(p).getName());
			}

			model = new CSGModel(gatedTrans, start, allDDRowVars, allDDColVars, allActDDVars(), modelVariables,
					varList, varDDRowVars, varDDColVars, actDDVars, playerInfo);
			start = null;
			allDDRowVars = null;
			allDDColVars = null;
			varDDRowVars = null;
			varDDColVars = null;
			modelVariables = null;
			actDDVars = null;

			model.setSynchs(new ArrayList<String>(modulesFile.getSynchs()));
			model.setRewards(new JDDNode[0], new JDDNode[0], new String[0]);
			model.setConstantValues(constantValues);

			model.setReach(reach);
			reach = null;

			model.filterReachableStates();
			model.findDeadlocks(prism.getFixDeadlocks());

			return model;
		} catch (PrismException e) {
			if (model != null) {
				model.clear();
			} else {
				clear();
			}
			throw e;
		}
	}

	// -------------------------------------------------------------------------------
	// Action ownership and canonical codes
	// -------------------------------------------------------------------------------

	@SuppressWarnings("unchecked")
	private void assignActionCodes()
	{
		playerActionCode = new Map[numPlayers];
		playerNumNamedActions = new int[numPlayers];
		for (int p = 0; p < numPlayers; p++) {
			BitSet owned = playerActions.getActionsForPlayer(p);
			playerActionCode[p] = new HashMap<>();
			// Codes are assigned in increasing action-index order, which (action indices
			// being allocated in file-declaration order by the parser) coincides with
			// first-declared order -- deterministic and reproducible either way.
			int code = IDLE_CODE + 1;
			for (int idx = owned.nextSetBit(0); idx >= 0; idx = owned.nextSetBit(idx + 1)) {
				playerActionCode[p].put(idx, code);
				code++;
			}
			playerNumNamedActions[p] = owned.cardinality();
		}
	}

	/** n_p: player p's action count, including the reserved idle slot. */
	private int numActionsForPlayer(int p)
	{
		return playerNumNamedActions[p] + 1;
	}

	// -------------------------------------------------------------------------------
	// Variable allocation (state vars + Act_p)
	// -------------------------------------------------------------------------------

	private void allocateDDVars()
	{
		modelVariables = new ModelVariablesDD();

		varDDRowVars = new JDDVars[numVars];
		varDDColVars = new JDDVars[numVars];
		for (int i = 0; i < numVars; i++) {
			varDDRowVars[i] = new JDDVars();
			varDDColVars[i] = new JDDVars();
		}
		// State variables: interleaved row/col, in declaration order. Variable-ordering
		// EXPERIMENTATION -- including whether Act_p belongs before, after or interleaved with
		// these -- is explicitly deferred; this is a reasonable starting point, not a
		// claimed-optimal default. Act_p is NOT required to sit after every state variable
		// here, since nothing downstream of this slice needs that (no hybrid solver).
		for (int i = 0; i < numVars; i++) {
			int n = varList.getRangeLogTwo(i);
			for (int j = 0; j < n; j++) {
				varDDRowVars[i].addVar(modelVariables.allocateVariable(varList.getName(i) + "." + j));
				varDDColVars[i].addVar(modelVariables.allocateVariable(varList.getName(i) + "'." + j));
			}
		}
		allDDRowVars = new JDDVars();
		allDDColVars = new JDDVars();
		for (int i = 0; i < numVars; i++) {
			allDDRowVars.mergeVarsFrom(varDDRowVars[i]);
			allDDColVars.mergeVarsFrom(varDDColVars[i]);
		}

		// Act_p: one block per player, sized to n_p = (named actions) + 1 for ⊥_p.
		// No row/col pairing -- an action is chosen once per transition step, not carried
		// across it.
		actDDVars = new JDDVars[numPlayers];
		for (int p = 0; p < numPlayers; p++) {
			actDDVars[p] = new JDDVars();
			int nBits = (int) Math.ceil(PrismUtils.log2(numActionsForPlayer(p)));
			String playerName = modulesFile.getPlayer(p).getName();
			for (int j = 0; j < nBits; j++) {
				actDDVars[p].addVar(modelVariables.allocateVariable(playerName + ".a" + j));
			}
		}
	}

	private JDDVars allActDDVars()
	{
		JDDVars all = new JDDVars();
		for (int p = 0; p < numPlayers; p++) {
			all.mergeVarsFrom(actDDVars[p]);
		}
		return all;
	}

	// -------------------------------------------------------------------------------
	// Per-action-slot equality BDD
	// -------------------------------------------------------------------------------

	/** BDD testing whether player p's Act_p block encodes exactly {@code code}. */
	private JDDNode isAction(int p, int code)
	{
		int nBits = actDDVars[p].n();
		JDDNode decoder = JDD.Constant(0);
		for (int val = 0; val < (1 << nBits); val++) {
			decoder = JDD.SetVectorElement(decoder, actDDVars[p], val, val);
		}
		JDDNode result = JDD.Apply(JDD.EQUALS, decoder, JDD.Constant(code));
		return result;
	}

	// -------------------------------------------------------------------------------
	// Assembling the action-list BDD for a command
	// -------------------------------------------------------------------------------

	/**
	 * The action-list BDD: conjunction of per-action equality tests for a command's action list. Players not
	 * mentioned in the action list are don't-cares, and an unlabelled command ([]) has action
	 * list [0], which is excluded here, giving the constant-true BDD -- it synchronises with
	 * every tuple.
	 */
	private JDDNode buildActionListBDD(Command command) throws PrismLangException
	{
		JDDNode actionListDD = JDD.Constant(1);
		for (int idx : command.getSynchIndices()) {
			if (idx == 0) {
				continue;
			}
			Integer p = playerActions.getPlayerForActionIndex(idx);
			if (p == null) {
				// CSGPlayerActions.validateIndependentModules should already have caught
				// this at construction time; defensive only.
				throw new PrismLangException("Action " + modulesFile.getSynch(idx - 1) + " is not associated to any player", command);
			}
			int code = playerActionCode[p].get(idx);
			actionListDD = JDD.Apply(JDD.TIMES, actionListDD, isAction(p, code));
		}
		return actionListDD;
	}

	// -------------------------------------------------------------------------------
	// Per-command transition, per-module moduleTrans, and their product into T
	// -------------------------------------------------------------------------------

	private JDDNode translateExpression(Expression e) throws PrismException
	{
		return expr2mtbdd.checkExpressionDD(e, JDD.ONE.copy());
	}

	/**
	 * The state x action-tuple set where {@code command} fires: its guard AND action-list BDD,
	 * before folding in its update. Reused both to build the per-command transition (below) and, summed
	 * across a module's commands, to detect when the module doesn't react at all (see
	 * {@code moduleIdentity} in {@link #translateModule}).
	 */
	private JDDNode translateCommandGuard(Command command) throws PrismException
	{
		JDDNode guardDD = translateExpression(command.getGuard());
		guardDD = JDD.Times(guardDD, range.copy());
		JDDNode actionListDD = buildActionListBDD(command);
		guardDD = JDD.Times(guardDD, actionListDD);
		return guardDD;
	}

	/** The per-command transition: guard AND action-list BDD, x sum of weighted updates. */
	private JDDNode translateCommand(int m, Command command, JDDNode guardDD) throws PrismException
	{
		if (guardDD.equals(JDD.ZERO)) {
			return guardDD;
		}
		JDDNode upDD = translateUpdates(m, command.getUpdates(), guardDD.copy());
		return JDD.Times(upDD, guardDD);
	}

	private JDDNode translateUpdates(int m, Updates u, JDDNode guard) throws PrismException
	{
		JDDNode dd = JDD.Constant(0);
		int n = u.getNumUpdates();
		for (int i = 0; i < n; i++) {
			JDDNode udd = translateUpdate(m, u.getUpdate(i), guard.copy());
			Expression p = u.getProbability(i);
			JDDNode pdd = (p == null) ? JDD.Constant(1.0) : translateExpression(p);
			udd = JDD.Times(udd, pdd);
			dd = JDD.Plus(dd, udd);
		}
		JDD.Deref(guard);
		return dd;
	}

	private JDDNode translateUpdate(int m, Update c, JDDNode guard) throws PrismException
	{
		for (int i = 0; i < numVars; i++) {
			varsUsed[i] = false;
		}
		JDDNode dd = JDD.Constant(1);
		int n = c.getNumElements();
		for (int i = 0; i < n; i++) {
			JDDNode udd = translateUpdateElement(c, i, guard.copy());
			dd = JDD.Times(dd, udd);
		}
		// Identity-padding: any variable of this module (or global) not mentioned in
		// this update is assumed unchanged.
		for (int i = 0; i < numVars; i++) {
			if ((varList.getModule(i) == m || varList.getModule(i) == -1) && !varsUsed[i]) {
				dd = JDD.Times(dd, varIdentities[i].copy());
			}
		}
		JDD.Deref(guard);
		return dd;
	}

	/**
	 * Translate a single update element x'=e. e may reference primed variables -- this falls
	 * straight out of using {@link CSGStateModelChecker} for translateExpression,
	 * no special-casing needed here.
	 */
	private JDDNode translateUpdateElement(Update c, int i, JDDNode guard) throws PrismException
	{
		String s = c.getVar(i);
		int v = varList.getIndex(s);
		if (v == -1) {
			throw new PrismLangException("Unknown variable \"" + s + "\" in update", c.getVarIdent(i));
		}
		varsUsed[v] = true;
		int l = varList.getLow(v);
		int h = varList.getHigh(v);

		JDDNode tmp1 = JDD.Constant(0);
		for (int j = l; j <= h; j++) {
			tmp1 = JDD.SetVectorElement(tmp1, varDDColVars[v], j - l, j);
		}
		JDDNode tmp2 = translateExpression(c.getExpression(i));
		tmp2 = JDD.Times(tmp2, guard.copy());
		JDDNode cl = JDD.Apply(JDD.EQUALS, tmp1, tmp2);
		cl = JDD.Times(cl, guard);
		cl = JDD.Times(cl, varColRangeDDs[v].copy());
		cl = JDD.Times(cl, range.copy());
		return cl;
	}

	/**
	 * {@code moduleTrans}: sum of a module's active commands' per-command transitions (using each command's
	 * precedence-resolved <i>effective</i> match set, not its raw one -- see
	 * {@link #buildEffectiveMatchSets}), plus the module-identity fallback -- if NO command of this
	 * module reacts to a given {@code (s,ā)} (whichever player owns it, if any, chose an
	 * action -- or was forced into idle -- that this module has no matching command for), the
	 * module contributes the identity transition on its own variables instead of a zero row.
	 * This is what keeps every module's own row summing to 1 regardless of {@code Enabled} --
	 * unconditional on Enabled is deliberate: it is only *within* Enabled that this needs to hold,
	 * and Enabled is what excludes the combinations where "no command reacts" would otherwise be
	 * wrong (e.g. a player choosing an action that isn't actually available). This is not the
	 * same thing as the global completeness *check* ({@link #checkGlobalCompleteness}) --
	 * this is unconditional construction, not validation; it cannot mask the specific bug that
	 * check is for (two commands sharing a primary label failing to jointly cover the
	 * co-action space), only the (fine, expected) case of a module simply not being party to
	 * the tuple at all.
	 */
	private JDDNode translateModule(int m, Module module) throws PrismException
	{
		int n = module.getNumCommands();
		Command[] commands = new Command[n];
		JDDNode[] matchSets = new JDDNode[n];
		for (int i = 0; i < n; i++) {
			commands[i] = module.getCommand(i);
			matchSets[i] = translateCommandGuard(commands[i]);
		}
		// buildEffectiveMatchSets consumes matchSets.
		JDDNode[] effective = buildEffectiveMatchSets(commands, matchSets);
		try {
			checkLocalRectangularity(module, commands, effective);
		} catch (PrismLangException e) {
			JDD.DerefArrayNonNull(effective);
			throw e;
		}

		JDDNode moduleTrans = JDD.Constant(0);
		JDDNode moduleActive = JDD.Constant(0);
		for (int i = 0; i < n; i++) {
			JDDNode eff = effective[i];
			if (eff.equals(JDD.ZERO)) {
				JDD.Deref(eff);
				continue;
			}
			moduleActive = JDD.Or(moduleActive, eff.copy());
			recordPrimaryCoverage(commands[i], eff.copy());
			JDDNode commandTrans = translateCommand(m, commands[i], eff);
			moduleTrans = JDD.Plus(moduleTrans, commandTrans);
		}
		JDDNode moduleIdentity = buildModuleIdentity(m);
		JDDNode moduleInactive = JDD.Not(moduleActive);
		moduleTrans = JDD.Plus(moduleTrans, JDD.Apply(JDD.TIMES, moduleInactive, moduleIdentity));
		return moduleTrans;
	}

	/**
	 * Does command {@code higher} take precedence over (shadow) command {@code lower} within
	 * the same module, wherever both would otherwise match the same {@code (s,ā)}?
	 * Mirrors {@code Updater.calculateTransitionsCSG}'s conflict-pruning pass exactly: a
	 * labelled command always shadows an unlabelled one (an unlabelled command has no owning
	 * player to resolve the tie, so it yields), and among commands sharing the same primary
	 * label, the one with more synch indices (the more specific co-action requirement)
	 * shadows one with fewer. Two labelled commands with <i>different</i> primary labels never
	 * shadow each other -- Updater's own pruning doesn't resolve that case either, so an
	 * overlap there is a genuine ambiguity for {@link #checkLocalRectangularity} to catch, not
	 * something to silently resolve here.
	 */
	private boolean shadows(Command higher, Command lower)
	{
		int higherPrimary = higher.getSynchIndices().get(0);
		int lowerPrimary = lower.getSynchIndices().get(0);
		if (lowerPrimary == 0) {
			return higherPrimary != 0;
		}
		if (higherPrimary == 0) {
			return false;
		}
		return higherPrimary == lowerPrimary && higher.getSynchIndices().size() > lower.getSynchIndices().size();
	}

	/**
	 * Resolve within-module overlap by precedence. Each command's raw match set
	 * ({@code Guard_κ ∧ A_κ}, {@code matchSets[i]}) has every *shadowing* command's own raw
	 * match set subtracted out -- mirroring {@code Updater}'s
	 * {@code expansions.get(m).get(j).remove(prod)} pass exactly, which prunes against each
	 * shadowing command's ORIGINAL expansion set, not a recursively-pruned one, so this uses
	 * {@code matchSets} (not the effective sets being built) as the subtrahend throughout.
	 * Consumes {@code matchSets} entirely -- callers must not use it afterwards.
	 */
	private JDDNode[] buildEffectiveMatchSets(Command[] commands, JDDNode[] matchSets)
	{
		int n = commands.length;
		JDDNode[] effective = new JDDNode[n];
		for (int j = 0; j < n; j++) {
			JDDNode eff = matchSets[j].copy();
			for (int i = 0; i < n; i++) {
				if (i != j && shadows(commands[i], commands[j])) {
					eff = JDD.Apply(JDD.TIMES, eff, JDD.Not(matchSets[i].copy()));
				}
			}
			effective[j] = eff;
		}
		JDD.DerefArrayNonNull(matchSets);
		return effective;
	}

	/**
	 * After precedence resolution ({@link #buildEffectiveMatchSets}), verify no
	 * two distinct commands of this module can still both match the same {@code (s,ā)} -- the
	 * symbolic analogue of {@code Updater}'s "Module ... has multiple active commands for
	 * action ..." hard error, checked over the whole state x action space in one pass instead
	 * of state-by-state. Also covers, as a special case falling out of the same mechanism
	 * (neither ever shadows the other), two unlabelled commands in the same independent
	 * module with overlapping guards.
	 */
	private void checkLocalRectangularity(Module module, Command[] commands, JDDNode[] effective) throws PrismLangException
	{
		int n = commands.length;
		for (int i = 0; i < n; i++) {
			for (int j = i + 1; j < n; j++) {
				JDDNode overlap = JDD.And(effective[i].copy(), effective[j].copy());
				boolean bad = !overlap.equals(JDD.ZERO);
				JDD.Deref(overlap);
				if (bad) {
					throw new PrismLangException("Module " + module.getName() + " has commands " + describeCommand(commands[i])
							+ " and " + describeCommand(commands[j]) + " that can both be active for the same state and joint action"
							+ " -- local rectangularity violated."
							+ " (No witness state extracted yet.)");
				}
			}
		}
	}

	private String describeCommand(Command c)
	{
		return c.isUnlabelled() ? "<unlabelled>" : c.getSynchs().toString();
	}

	/**
	 * Fold {@code command}'s own (shadowing-resolved) effective match set into
	 * its *primary* label's running coverage, across every module in the model -- see
	 * {@link #primaryCoverageByActionIndex}'s field doc and {@link #checkGlobalCompleteness}.
	 * Unlabelled commands (primary index 0) contribute nothing here: they aren't any player's
	 * action and so have no primary-label coverage obligation to discharge (same reasoning as
	 * local rectangularity's own treatment of this command category). Consumes {@code effCopy}.
	 */
	private void recordPrimaryCoverage(Command command, JDDNode effCopy)
	{
		int label = command.getSynchIndices().get(0);
		if (label == 0) {
			JDD.Deref(effCopy);
			return;
		}
		JDDNode existing = primaryCoverageByActionIndex.get(label);
		primaryCoverageByActionIndex.put(label, existing == null ? effCopy : JDD.Or(existing, effCopy));
	}

	/**
	 * Global completeness: the *global*, multi-command companion to local rectangularity, since
	 * the module-identity fallback makes the row-sum formulation alone insufficient (every module's own row sums to 1
	 * unconditionally, whether or not it actually reacts). What can still go wrong, and what
	 * this actually catches: player p chooses (owns) some named action a as part of a
	 * legally-available joint tuple ā ({@link #enabled} already confirms *some* primary-a
	 * command's guard holds somewhere), but across every module in the model -- not just p's
	 * own, since an independent module may also carry a command primarily labelled a -- no
	 * command whose primary label is a actually has both its guard AND its full (possibly
	 * multi-label) action list {@code A_κ} satisfied by this specific ā. That is exactly
	 * {@code Updater.calculateTransitionsCSG}'s "Missing specification for action product [...]"
	 * pass, generalised from one concrete state to the whole state x action space: {@code tmp}
	 * there starts as the joint tuple's own action-index bits and gets a bit cleared whenever
	 * some active command's full expansion (guard + action list, not guard alone) covers the
	 * tuple for that bit's action; anything left uncleared is missing specification. Here,
	 * {@link #primaryCoverageByActionIndex} (already accumulated from the same
	 * shadowing-resolved effective match sets {@link #translateModule} builds, via
	 * {@link #recordPrimaryCoverage}) plays the role of "some command's expansion covers ā",
	 * and {@code IsAction_p(code(a)) ∧ Enabled} plays the role of "ā chooses a and is otherwise
	 * legal".
	 */
	private void checkGlobalCompleteness() throws PrismLangException
	{
		List<Integer> uncovered = new ArrayList<Integer>();
		for (Map.Entry<Integer, JDDNode> ce : primaryCoverageByActionIndex.entrySet()) {
			int label = ce.getKey();
			JDDNode covered = ce.getValue();
			Integer p = playerActions.getPlayerForActionIndex(label);
			if (p == null) {
				// CSGPlayerActions.validateIndependentModules should already have caught this at
				// construction time (a primary label used anywhere must already be owned); defensive only.
				continue;
			}
			int code = playerActionCode[p].get(label);
			JDDNode chosen = isAction(p, code);
			JDDNode badLabel = JDD.Apply(JDD.TIMES, chosen, enabled.copy());
			badLabel = JDD.Apply(JDD.TIMES, badLabel, JDD.Not(covered.copy()));
			boolean isZero = badLabel.equals(JDD.ZERO);
			JDD.Deref(badLabel);
			if (!isZero) {
				uncovered.add(label);
			}
		}
		if (!uncovered.isEmpty()) {
			StringBuilder names = new StringBuilder();
			for (int i = 0; i < uncovered.size(); i++) {
				if (i > 0) {
					names.append(", ");
				}
				names.append(modulesFile.getSynch(uncovered.get(i) - 1));
			}
			throw new PrismLangException("Missing specification for action(s) [" + names + "]"
					+ " -- some legally-available joint action tuple choosing this action has no"
					+ " command, in any module, whose primary label is this action and whose guard and full action"
					+ " list actually cover that tuple."
					+ " (No witness state extracted yet.)");
		}
	}

	/**
	 * The module-identity fallback: the identity transition on every variable module {@code m}
	 * could conceivably touch -- its own variables plus globals, exactly the set
	 * {@link #translateUpdate}'s own identity-padding already treats as "this module's to
	 * account for". Used only when the module has no active command at all for a given
	 * {@code (s,ā)}.
	 */
	private JDDNode buildModuleIdentity(int m)
	{
		JDDNode id = JDD.Constant(1);
		for (int i = 0; i < numVars; i++) {
			if (varList.getModule(i) == m || varList.getModule(i) == -1) {
				id = JDD.Times(id, varIdentities[i].copy());
			}
		}
		return id;
	}

	/** T = product of all modules' moduleTrans. */
	private void translateModules() throws PrismException
	{
		varsUsed = new boolean[numVars];
		primaryCoverageByActionIndex = new HashMap<Integer, JDDNode>();
		trans = JDD.Constant(1);
		for (int m = 0; m < numModules; m++) {
			JDDNode moduleTrans = translateModule(m, modulesFile.getModule(m));
			trans = JDD.Apply(JDD.TIMES, trans, moduleTrans);
		}
	}

	// -------------------------------------------------------------------------------
	// Per-player enabled actions and the Enabled BDD gating legal joint action tuples
	// -------------------------------------------------------------------------------

	/**
	 * Populates {@link #enabledByActionIndex}: for every synchronising-action index that is
	 * some player's own *primary* command label (i.e. every index {@code CSGPlayerActions}
	 * assigns an owner to), the state-only BDD of "some command across that player's own
	 * modules, primarily labelled with this action, has a satisfied guard" -- exactly
	 * Updater's {@code moves[p]} construction, generalised to the whole state space at once.
	 */
	private void buildEnabledActions() throws PrismException
	{
		enabledByActionIndex = new HashMap<Integer, JDDNode>();
		for (int m = 0; m < numModules; m++) {
			if (playerActions.getPlayerForModule(m) == -1) {
				// Independent modules aren't a player's strategic choice, so they don't
				// participate in anyone's enabled-action set -- see buildActionListBDD's own
				// comment on unlabelled/[] commands for the companion fact.
				continue;
			}
			Module module = modulesFile.getModule(m);
			int n = module.getNumCommands();
			for (int i = 0; i < n; i++) {
				Command command = module.getCommand(i);
				// Player-owned commands are always labelled (CSGPlayerActions.computeOwnership
				// already enforces this), so get(0) is always a real action index here.
				int primary = command.getSynchIndices().get(0);
				JDDNode guard = translateExpression(command.getGuard());
				guard = JDD.Times(guard, range.copy());
				JDDNode existing = enabledByActionIndex.get(primary);
				enabledByActionIndex.put(primary, existing == null ? guard : JDD.Or(existing, guard));
			}
		}
	}

	/** Does player p have some enabled (guard-true) named action at s? */
	private JDDNode enabledForPlayer(int p)
	{
		JDDNode enabled = JDD.Constant(0);
		BitSet owned = playerActions.getActionsForPlayer(p);
		for (int idx = owned.nextSetBit(0); idx >= 0; idx = owned.nextSetBit(idx + 1)) {
			JDDNode a = enabledByActionIndex.get(idx);
			enabled = JDD.Or(enabled, a == null ? JDD.Constant(0) : a.copy());
		}
		return enabled;
	}

	/**
	 * LegalAct_p(s, Act_p): the code player p's Act_p block holds is one they could actually
	 * have chosen at s -- idle iff no action is enabled for p, otherwise that specific named action's own
	 * guard. Codes outside {@code 0..n_p-1} (allocation padding) are legal nowhere,
	 * with no special-casing needed: they simply never match any disjunct below.
	 */
	private JDDNode legalActForPlayer(int p)
	{
		JDDNode noEnabledAction = JDD.Not(enabledForPlayer(p));
		JDDNode legal = JDD.Apply(JDD.TIMES, isAction(p, IDLE_CODE), noEnabledAction);
		BitSet owned = playerActions.getActionsForPlayer(p);
		for (int idx = owned.nextSetBit(0); idx >= 0; idx = owned.nextSetBit(idx + 1)) {
			int code = playerActionCode[p].get(idx);
			JDDNode a = enabledByActionIndex.get(idx);
			JDDNode term = JDD.Apply(JDD.TIMES, isAction(p, code), a == null ? JDD.Constant(0) : a.copy());
			legal = JDD.Or(legal, term);
		}
		return legal;
	}

	/** Enabled(s, ā) = AND_p LegalAct_p(s, ā_p): ā is a joint tuple every player could actually choose. */
	private JDDNode buildEnabled()
	{
		JDDNode result = JDD.Constant(1);
		for (int p = 0; p < numPlayers; p++) {
			result = JDD.Apply(JDD.TIMES, result, legalActForPlayer(p));
		}
		return result;
	}

	// -------------------------------------------------------------------------------
	// Primed-variable dependency acyclicity check
	// -------------------------------------------------------------------------------

	/**
	 * Reject a model where some reachable {@code (s,ā)} gives the primed-variable-in-RHS
	 * updates (e.g. {@code s1'=c'?1:0}) no well-defined resolution -- a cycle among the
	 * *currently active* update definitions. Summary of what this builds:
	 * <ol>
	 * <li>{@code Vars} = every variable ever targeted by a primed-referencing update anywhere
	 *     in the model. Empty is the common case -- mirrors {@code ChoiceListFlexi}'s own
	 *     {@code getContainsPrimes} fast path, just moved to build time and made whole-model,
	 *     so this returns immediately without allocating anything.</li>
	 * <li>A small auxiliary {@code Var}/{@code Var'} ADD-variable block indexes {@code Vars}
	 *     (unrelated to the model's own {@code S}/{@code S'}/{@code Act_p} blocks) -- plus a
	 *     transient {@code Var''} block, needed only for the relational-composition step
	 *     below and derefed before this method returns.</li>
	 * <li>{@code Edge(s,ā,v_i,v_j)}: disjunction, over every (command κ, target variable x)
	 *     pair with static dependency set {@code D = Deps(κ,x)} (the primed-variable names
	 *     {@code κ}'s own update for {@code x} syntactically references, via the existing
	 *     engine-agnostic {@code Expression.getPrimedVars()}), of
	 *     {@code (Γ_κ ∧ A_κ) x [Var=code(x)] x (OR_{y in D} [Var'=code(y)])}.</li>
	 *     <li>Transitive closure of {@code Edge} over the {@code Var}/{@code Var'} dimension,
	 *     via the same {@code ThereExists}/{@code PermuteVariables} relational-composition
	 *     idiom the reachability fixpoint already uses -- just over this tiny dimension
	 *     instead of the state space, bounded by {@code |Vars|} iterations.</li>
	 * <li>{@code Bad_v = R(s,ā,v,v)} for each {@code v} -- a genuine cycle through {@code v}.
	 *     Checked per-variable (not just {@code ∃v. R(s,ā,v,v)}) so the diagnostic can name
	 *     the actual cyclic variables, matching {@code ChoiceListFlexi}'s own "Cyclic updates
	 *     with variables [...]" runtime error.</li>
	 * </ol>
	 * Not subsumed by global completeness: a tautological cycle ({@code x'=y' ∧ y'=x'}) can
	 * leave the model underdetermined without reliably violating row-sum-to-1, so this stays
	 * a dedicated check.
	 */
	private void checkPrimedVarAcyclicity() throws PrismException
	{
		List<Integer> varsList = new ArrayList<Integer>();
		Map<Command, Map<Integer, Set<Integer>>> depsByCommand = new HashMap<Command, Map<Integer, Set<Integer>>>();

		for (int m = 0; m < numModules; m++) {
			Module module = modulesFile.getModule(m);
			int nc = module.getNumCommands();
			for (int c = 0; c < nc; c++) {
				Command command = module.getCommand(c);
				Updates ups = command.getUpdates();
				int nu = ups.getNumUpdates();
				Map<Integer, Set<Integer>> depsForCommand = null;
				for (int u = 0; u < nu; u++) {
					Update up = ups.getUpdate(u);
					int ne = up.getNumElements();
					for (int e = 0; e < ne; e++) {
						int v = varList.getIndex(up.getVar(e));
						if (v == -1) {
							continue; // an unknown update target is translateUpdateElement's error to raise, not this check's
						}
						List<String> primedNames = up.getExpression(e).getPrimedVars();
						if (primedNames.isEmpty()) {
							continue;
						}
						if (!varsList.contains(v)) {
							varsList.add(v);
						}
						if (depsForCommand == null) {
							depsForCommand = new HashMap<Integer, Set<Integer>>();
							depsByCommand.put(command, depsForCommand);
						}
						Set<Integer> deps = depsForCommand.get(v);
						if (deps == null) {
							deps = new HashSet<Integer>();
							depsForCommand.put(v, deps);
						}
						for (String depName : primedNames) {
							int dv = varList.getIndex(depName);
							if (dv != -1) {
								deps.add(dv);
								if (!varsList.contains(dv)) {
									varsList.add(dv);
								}
							}
						}
					}
				}
			}
		}

		if (varsList.isEmpty()) {
			return;
		}

		Map<Integer, Integer> varCode = new HashMap<Integer, Integer>();
		for (int i = 0; i < varsList.size(); i++) {
			varCode.put(varsList.get(i), i);
		}
		int nBits = Math.max((int) Math.ceil(PrismUtils.log2(varsList.size())), 1);
		JDDVars rowVars = new JDDVars();
		JDDVars colVars = new JDDVars();
		JDDVars tmpVars = new JDDVars();
		for (int b = 0; b < nBits; b++) {
			rowVars.addVar(modelVariables.allocateVariable("_primeDep.v" + b));
			colVars.addVar(modelVariables.allocateVariable("_primeDep.v'" + b));
			tmpVars.addVar(modelVariables.allocateVariable("_primeDep.v''" + b));
		}

		try {
			JDDNode edge = JDD.Constant(0);
			for (Map.Entry<Command, Map<Integer, Set<Integer>>> ce : depsByCommand.entrySet()) {
				JDDNode matchSet = translateCommandGuard(ce.getKey());
				if (matchSet.equals(JDD.ZERO)) {
					JDD.Deref(matchSet);
					continue;
				}
				for (Map.Entry<Integer, Set<Integer>> ve : ce.getValue().entrySet()) {
					JDDNode depsDD = JDD.Constant(0);
					for (int y : ve.getValue()) {
						depsDD = JDD.Or(depsDD, indexEquals(colVars, varCode.get(y)));
					}
					JDDNode term = JDD.Apply(JDD.TIMES, matchSet.copy(), indexEquals(rowVars, varCode.get(ve.getKey())));
					term = JDD.Apply(JDD.TIMES, term, depsDD);
					edge = JDD.Or(edge, term);
				}
				JDD.Deref(matchSet);
			}

			JDDNode r = edge.copy();
			for (int iter = 0; iter < varsList.size(); iter++) {
				JDDNode rRenamed = JDD.PermuteVariables(r.copy(), colVars, tmpVars);
				JDDNode edgeRenamed = JDD.PermuteVariables(edge.copy(), rowVars, tmpVars);
				JDDNode composed = JDD.ThereExists(JDD.And(rRenamed, edgeRenamed), tmpVars);
				JDDNode next = JDD.Or(r.copy(), composed);
				if (next.equals(r)) {
					JDD.Deref(next);
					break;
				}
				JDD.Deref(r);
				r = next;
			}

			List<Integer> cyclic = new ArrayList<Integer>();
			for (int v : varsList) {
				int code = varCode.get(v);
				JDDNode pinned = JDD.Apply(JDD.TIMES, indexEquals(rowVars, code), indexEquals(colVars, code));
				JDDNode badV = JDD.And(r.copy(), pinned);
				boolean isZero = badV.equals(JDD.ZERO);
				JDD.Deref(badV);
				if (!isZero) {
					cyclic.add(v);
				}
			}
			JDD.Deref(edge);
			JDD.Deref(r);

			if (!cyclic.isEmpty()) {
				StringBuilder names = new StringBuilder();
				for (int i = 0; i < cyclic.size(); i++) {
					if (i > 0) {
						names.append(", ");
					}
					names.append(varList.getName(cyclic.get(i)));
				}
				throw new PrismLangException("Cyclic primed-variable update dependency among variables [" + names
						+ "] -- some reachable (state, joint action) combination has no well-defined"
						+ " resolution for these updates. (No witness state extracted yet.)");
			}
		} finally {
			rowVars.derefAll();
			colVars.derefAll();
			tmpVars.derefAll();
		}
	}

	/** Small decoder-equality BDD over an auxiliary index block -- same pattern as {@link #isAction}. */
	private JDDNode indexEquals(JDDVars vars, int code)
	{
		int nBits = vars.n();
		JDDNode decoder = JDD.Constant(0);
		for (int val = 0; val < (1 << nBits); val++) {
			decoder = JDD.SetVectorElement(decoder, vars, val, val);
		}
		return JDD.Apply(JDD.EQUALS, decoder, JDD.Constant(code));
	}

	// -------------------------------------------------------------------------------
	// Identities/ranges (same construction as Modules2MTBDD.sortIdentities/sortRanges,
	// at variable granularity only -- module-level identities are not needed by this slice)
	// -------------------------------------------------------------------------------

	private void buildIdentitiesAndRanges()
	{
		varIdentities = new JDDNode[numVars];
		for (int i = 0; i < numVars; i++) {
			JDDNode id = JDD.Constant(0);
			for (int j = 0; j < varList.getRange(i); j++) {
				id = JDD.SetMatrixElement(id, varDDRowVars[i], varDDColVars[i], j, j, 1);
			}
			varIdentities[i] = id;
		}
		range = JDD.Constant(1);
		varColRangeDDs = new JDDNode[numVars];
		for (int i = 0; i < numVars; i++) {
			JDDNode varRangeDD = JDD.SumAbstract(varIdentities[i].copy(), varDDColVars[i]);
			varColRangeDDs[i] = JDD.SumAbstract(varIdentities[i].copy(), varDDRowVars[i]);
			range = JDD.Apply(JDD.TIMES, range, varRangeDD);
		}
	}

	// -------------------------------------------------------------------------------
	// Initial states (same construction as Modules2MTBDD.buildInitialStates)
	// -------------------------------------------------------------------------------

	private void buildInitialStates() throws PrismException
	{
		if (modulesFile.getInitialStates() != null) {
			start = translateExpression(modulesFile.getInitialStates());
			start = JDD.And(start, range.copy());
			if (start.equals(JDD.ZERO)) {
				throw new PrismLangException("No initial states: \"init\" construct evaluates to false", modulesFile.getInitialStates());
			}
		} else {
			start = JDD.Constant(1);
			for (int i = 0; i < numVars; i++) {
				Object startObj = modulesFile.getVarDeclaration(i).getStartOrDefault().evaluate(constantValues);
				int startInt = varList.encodeToInt(i, startObj);
				JDDNode tmp = JDD.SetVectorElement(JDD.Constant(0), varDDRowVars[i], startInt, 1);
				start = JDD.And(start, tmp);
			}
		}
	}

	// -------------------------------------------------------------------------------
	// Reachability, via the existing PrismMTBDD.Reachability primitive (the same one
	// ModelSymbolic.doReachability() uses) rather than a hand-rolled fixpoint.
	// -------------------------------------------------------------------------------

	private void buildReachability() throws PrismException
	{
		JDDVars actVars = allActDDVars();
		// T_det: collapse the action dimension to "is it possible to transition",
		// restricted to Enabled -- otherwise a joint tuple no player could actually have
		// chosen (e.g. a non-idle action at a state where it isn't guard-enabled) would still
		// explore successors via the module-identity fallback, which is harmless for that
		// tuple in isolation but would wrongly validate the tuple as "possible" here.
		JDDNode transRestricted = JDD.Apply(JDD.TIMES, trans.copy(), enabled.copy());
		JDDNode transDet = JDD.ThereExists(JDD.GreaterThan(transRestricted, 0), actVars);
		// actVars is allActDDVars()'s own fresh copies (mergeVarsFrom/copyVarsFrom), read-only
		// to ThereExists (which -- like every other JDDVars consumer in this class -- does not
		// itself deref the variable list it's given) -- ours to release once done with it.
		actVars.derefAll();
		// PrismMTBDD.Reachability is documented "[ REFS: result, DEREFS: none ]" -- it derefs
		// neither argument (see NondetModel.doReachability()'s own use of it for the established
		// pattern: an explicit Deref of a throwaway trans-like argument right after the call,
		// and no copy of `start`, since -- like here -- it's a field that must stay valid
		// afterward, not a throwaway). transDet is throwaway (needed only to feed this call);
		// start is this class's own field, still needed by translate() afterward, so passed
		// directly rather than copied -- Reachability leaves it untouched either way.
		reach = PrismMTBDD.Reachability(transDet, allDDRowVars, allDDColVars, start);
		JDD.Deref(transDet);
	}

	// -------------------------------------------------------------------------------
	// Accessors, including the state/transition counts used for cross-checking against the
	// explicit engine (an oracle-style comparison).
	// -------------------------------------------------------------------------------

	public JDDNode getTrans()
	{
		return trans;
	}

	/** Enabled(s, ā): the BDD of joint action tuples legally available at each state. */
	public JDDNode getEnabled()
	{
		return enabled;
	}

	public JDDNode getReach()
	{
		return reach;
	}

	public JDDNode getStart()
	{
		return start;
	}

	public double getNumReachableStates()
	{
		return JDD.GetNumMinterms(reach, allDDRowVars.n());
	}

	/**
	 * Number of (state, action-tuple) pairs with a nonzero-probability transition, restricted
	 * to reachable states and to {@link #getEnabled() legally available} tuples -- matching what
	 * the explicit engine ever builds a choice for.
	 */
	public double getNumTransitions()
	{
		JDDNode transReachable = JDD.Times(trans.copy(), reach.copy());
		transReachable = JDD.Times(transReachable, enabled.copy());
		// Bit count only, no DD node needed -- avoids allocating (and having to deref) a fresh
		// allActDDVars() copy just to read its size.
		int numActDDVars = 0;
		for (int p = 0; p < numPlayers; p++) {
			numActDDVars += actDDVars[p].n();
		}
		int numDDVars = allDDRowVars.n() + allDDColVars.n() + numActDDVars;
		double count = JDD.GetNumMinterms(JDD.GreaterThan(transReachable, 0), numDDVars);
		return count;
	}
}
