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

package symbolic.comp;

import java.util.BitSet;
import java.util.List;

import jdd.JDD;
import jdd.JDDNode;
import jdd.JDDVars;
import parser.ast.Coalition;
import parser.ast.Expression;
import parser.ast.ExpressionStrategy;
import parser.ast.ExpressionStrategyQual;
import parser.ast.ExpressionTemporal;
import parser.ast.ExpressionUnaryOp;
import parser.ast.PropertiesFile;
import parser.type.TypeBool;
import prism.Prism;
import prism.PrismException;
import prism.PrismNotSupportedException;
import symbolic.model.CSGModel;
import symbolic.model.Model;
import symbolic.states.StateValues;
import symbolic.states.StateValuesMTBDD;

/**
 * Model checker for CSGs (concurrent stochastic games), symbolic engine.
 *
 * <p><b>Still not a full CSG model checker.</b> rPATL-style {@code <<coalition>> P/R} property
 * syntax remains separate, not-yet-done work.
 * What this class now provides, as direct Java entry points over a coalition fixed via
 * {@link #setCoalition}, mirroring {@code explicit.CSGModelChecker}'s own methods of the same
 * name (that class is the ground truth these were derived from, including the "legality
 * guard convention" that makes them
 * correct, re-derived per primitive rather than assumed from a sibling each time):
 * <ul>
 * <li><b>sure mode</b>: {@link #G}/{@link #SF}/{@link #SFG}/{@link #SGF}, built from the shared
 *     {@link #pre1};</li>
 * <li><b>almost-sure mode</b>: {@link #AF}/{@link #AFG}/{@link #AGF}, built from the shared
 *     {@link #apreXY}/{@link #apreXYZ};</li>
 * <li><b>limit-sure mode</b>: {@link #LF}/{@link #LFG}/{@link #LGF}, built from the shared
 *     {@link #lpreXY}/{@link #lpreXYZ} -- {@code G} is reused as-is for both
 *     almost and limit mode (safety provably coincides across all three qualitative modes).</li>
 * </ul>
 * All ten qualitative operators are now implemented, and reachable from real property syntax:
 * {@code <<coalition>>}/{@code [[coalition]]
 * sure/almost/limit [ G psi | F psi | F G psi | G F psi ]} for a Boolean state formula
 * {@code psi}, via {@link #checkExpression}'s override recognising {@link ExpressionStrategy}
 * -- mirroring {@code symbolic.comp.GamesModelChecker}'s own precedent for the *shape* of that
 * override (this is the pattern actually used elsewhere in this package for recognising
 * {@code <<>>}/{@code [[]]}, not {@code ProbModelChecker}, which has no such logic on the
 * symbolic side at all), and {@code explicit.CSGModelChecker.checkExpressionStrategyQual} for
 * the mode/shape dispatch logic itself. Not supported: rPATL {@code P}/{@code R} strategy
 * queries, and the "combo" (G F & F G) non-atomic shape the explicit engine handles via
 * Rabin-chain machinery this class doesn't have -- both rejected with a clear
 * {@link PrismNotSupportedException} rather than silently mishandled.
 *
 * <p>Boolean state formulas/labels/filters still work via {@link StateModelChecker}'s own
 * generic {@code checkExpression} dispatch, which this class's own override falls back to for
 * anything that isn't an {@link ExpressionStrategy}.
 *
 * <p>Deliberately does NOT extend {@link NonProbModelChecker}: that class's CTL {@code E}/
 * {@code A} operators quantify over the <i>whole</i> joint transition relation with no
 * per-player distinction ("does some/every combined choice of every player lead to..."),
 * which is not a meaningful semantics for a genuinely multi-player, simultaneous-move game --
 * it would silently answer a different question than the one a CSG property is actually
 * asking, rather than admitting it can't yet. Anything beyond the generic base dispatch and
 * this class's own methods (temporal operators, {@code ExpressionProb}, rPATL-style
 * {@code ExpressionStrategy} operands) instead falls through to
 * {@link StateModelChecker#checkExpression}'s own honest {@code "Couldn't check <class>"}
 * fallback.
 */
public class CSGModelChecker extends StateModelChecker
{
	private CSGModel csgModel;

	// Coalition context: mutable, set once per top-level query via
	// setCoalition, consulted implicitly by pre1/G/SF/SFG/SGF -- mirrors
	// explicit.CSGModelChecker's own buildCoalitions()-then-implicit-instance-field design
	// (that class's coalitionIndexes/actionIndexes, consulted deep inside its own fixpoint
	// loops, rather than threaded as an explicit parameter to every operator call).
	private JDDVars actC;
	private JDDVars actNC;
	private JDDNode enabled;
	private JDDNode legalC;
	private JDDNode legalNC;

	public CSGModelChecker(Prism prism, Model m, PropertiesFile pf) throws PrismException
	{
		super(prism, m, pf);
		this.csgModel = (CSGModel) m;
	}

	// -------------------------------------------------------------------------------
	// Property-syntax wiring: <<coalition>>/[[coalition]] sure/almost/limit
	// [ G/F/F G/G F psi ], dispatching to the ten methods below.
	// -------------------------------------------------------------------------------

	@Override
	public StateValues checkExpression(Expression expr, JDDNode statesOfInterest) throws PrismException
	{
		StateValues res;
		if (expr instanceof ExpressionStrategy) {
			res = checkExpressionStrategy((ExpressionStrategy) expr, statesOfInterest);
		} else {
			res = super.checkExpression(expr, statesOfInterest);
		}
		if (res instanceof StateValuesMTBDD && reach != null) {
			res.filter(reach);
		}
		return res;
	}

	/**
	 * Recognise {@code <<coalition>>}/{@code [[coalition]]} and delegate to
	 * {@link #checkExpressionStrategyQual} -- mirrors {@code symbolic.comp.GamesModelChecker
	 * .checkExpressionStrategy}'s own shape (that is the existing precedent for this kind of
	 * override on the symbolic side, not {@code ProbModelChecker}, which has none). Only the
	 * qualitative operator ({@code sure}/{@code almost}/{@code limit}) is supported here --
	 * rPATL-style {@code P}/{@code R} strategy queries are not (still separate,
	 * not-yet-done work). Consumes {@code statesOfInterest} (via {@link #checkExpressionDD} or
	 * the eventual {@link #checkExpressionStrategyQual} call, on every path).
	 */
	protected StateValues checkExpressionStrategy(ExpressionStrategy expr, JDDNode statesOfInterest) throws PrismException
	{
		boolean forAll = !expr.isThereExists();
		if (expr.getNumCoalitions() > 1) {
			JDD.Deref(statesOfInterest);
			throw new PrismNotSupportedException("The " + expr.getOperatorString() + " operator can only contain one coalition");
		}
		Coalition coalition = expr.getCoalition();
		if (coalition == null) {
			JDD.Deref(statesOfInterest);
			throw new PrismNotSupportedException("The " + expr.getOperatorString() + " operator requires an explicit coalition for a CSG");
		}
		List<Expression> exprs = expr.getOperands();
		if (exprs.size() != 1) {
			JDD.Deref(statesOfInterest);
			throw new PrismNotSupportedException("Cannot currently check strategy operators with lists of expressions");
		}
		Expression exprSub = exprs.get(0);
		if (!(exprSub instanceof ExpressionStrategyQual)) {
			JDD.Deref(statesOfInterest);
			throw new PrismNotSupportedException("Only the qualitative (sure/almost/limit) strategy operator is currently "
					+ "supported for CSGs by the symbolic engine");
		}
		return checkExpressionStrategyQual((ExpressionStrategyQual) exprSub, forAll, coalition, statesOfInterest);
	}

	/**
	 * Compute the (Boolean) result of a qualitative strategy operator, e.g.
	 * {@code sure/almost/limit [ G psi | F psi | F G psi | G F psi ]}, for a single coalition
	 * against its complement. Mirrors {@code explicit.CSGModelChecker.checkExpressionStrategyQual}
	 * almost exactly: same shape-classification via
	 * {@code ExpressionTemporal.isGlobally}/{@code isFinally}/{@code isGloballyFinally}/
	 * {@code isFinallyGlobally}, same {@code (mode, shape)} dispatch table -- restricted here to
	 * the four atomic shapes; the "combo" (G F & F G) shape the explicit engine also supports
	 * via Rabin-chain machinery (not implemented symbolically) is rejected explicitly rather
	 * than silently mishandled.
	 *
	 * <p>{@code forAll} ({@code [[C]]}) plays the adversarial role: per
	 * {@code explicit.CSGModelChecker}'s own convention ("[[C]] Op[phi] is definitionally 'the
	 * complement of C plays the achieving role for phi instead'"), this flips which players are
	 * passed to {@link #setCoalition} -- {@code setCoalition}'s argument always means "the
	 * achieving/maximising side," so {@code [[C]]} sets the coalition to {@code C}'s complement
	 * rather than {@code C} itself.
	 */
	protected StateValues checkExpressionStrategyQual(ExpressionStrategyQual expr, boolean forAll, Coalition coalition, JDDNode statesOfInterest)
			throws PrismException
	{
		Expression phi = expr.getExpression();
		while (Expression.isParenth(phi)) {
			phi = ((ExpressionUnaryOp) phi).getOperand();
		}

		ExpressionStrategyQual.QualMode mode = expr.getMode();

		String shape;
		Expression psi;
		if (ExpressionTemporal.isGloballyFinally(phi)) {
			shape = "G F";
			psi = ((ExpressionTemporal) ((ExpressionTemporal) phi).getOperand2()).getOperand2();
		} else if (ExpressionTemporal.isFinallyGlobally(phi)) {
			shape = "F G";
			psi = ((ExpressionTemporal) ((ExpressionTemporal) phi).getOperand2()).getOperand2();
		} else if (ExpressionTemporal.isGlobally(phi)) {
			shape = "G";
			psi = ((ExpressionTemporal) phi).getOperand2();
		} else if (ExpressionTemporal.isFinally(phi)) {
			shape = "F";
			psi = ((ExpressionTemporal) phi).getOperand2();
		} else {
			JDD.Deref(statesOfInterest);
			throw new PrismNotSupportedException("The " + mode + " operator currently supports only G, F, G F and F G "
					+ "path formulas for CSGs -- got " + phi);
		}

		if (!(psi.getType() instanceof TypeBool)) {
			JDD.Deref(statesOfInterest);
			throw new PrismException("The " + mode + " operator's path formula must wrap a Boolean state formula, not " + psi);
		}

		List<String> playerNames = csgModel.getPlayerNames();
		BitSet players = new BitSet();
		for (int p = 0; p < playerNames.size(); p++) {
			if (coalition.isPlayerIndexInCoalition(p, playerNames)) {
				players.set(p);
			}
		}
		if (forAll) {
			players.flip(0, playerNames.size());
		}
		setCoalition(players);

		JDDNode target = checkExpressionDD(psi, csgModel.getReach().copy());

		JDDNode win;
		if (shape.equals("G")) {
			win = G(target);
		} else if (shape.equals("F") && mode == ExpressionStrategyQual.QualMode.SURE) {
			win = SF(target);
		} else if (shape.equals("F") && mode == ExpressionStrategyQual.QualMode.ALMOST) {
			win = AF(target);
		} else if (shape.equals("F") && mode == ExpressionStrategyQual.QualMode.LIMIT) {
			win = LF(target);
		} else if (shape.equals("F G") && mode == ExpressionStrategyQual.QualMode.SURE) {
			win = SFG(target);
		} else if (shape.equals("F G") && mode == ExpressionStrategyQual.QualMode.ALMOST) {
			win = AFG(target);
		} else if (shape.equals("F G") && mode == ExpressionStrategyQual.QualMode.LIMIT) {
			win = LFG(target);
		} else if (shape.equals("G F") && mode == ExpressionStrategyQual.QualMode.SURE) {
			win = SGF(target);
		} else if (shape.equals("G F") && mode == ExpressionStrategyQual.QualMode.ALMOST) {
			win = AGF(target);
		} else {
			win = LGF(target);
		}

		clearCoalition();
		JDD.Deref(statesOfInterest);
		return new StateValuesMTBDD(win, csgModel);
	}

	/**
	 * Fix the coalition C (0-indexed player indices set in {@code players}) that every
	 * subsequent {@link #pre1}/{@link #G}/{@link #SF}/{@link #SFG}/{@link #SGF} call on this
	 * checker instance uses, until the next call to this method -- see the class doc for why
	 * this is a mutable, set-once-per-formula field rather than a parameter
	 * threaded through every operator call.
	 *
	 * <p>Derives {@link #actC}/{@link #actNC} (the coalition's own / the complement's
	 * {@code Act_p} variable blocks, merged via the same copy-safe
	 * {@code JDDVars.mergeVarsFrom} pattern {@code CSG2MTBDD.allActDDVars()} already uses) and
	 * {@link #enabled}/{@link #legalC}/{@link #legalNC} ({@code Enabled}, and its projection
	 * onto each side) from the model's own {@code trans}. Degenerate coalitions (all-players or
	 * empty) need no special-casing: merging zero players' blocks into {@code actC} or
	 * {@code actNC} just leaves that side an empty {@code JDDVars}, and quantifying over an
	 * empty variable block is a no-op.
	 */
	public void setCoalition(BitSet players)
	{
		clearCoalition();
		actC = new JDDVars();
		actNC = new JDDVars();
		int numPlayers = csgModel.getPlayerNames().size();
		for (int p = 0; p < numPlayers; p++) {
			if (players.get(p)) {
				actC.mergeVarsFrom(csgModel.getActDDVars(p));
			} else {
				actNC.mergeVarsFrom(csgModel.getActDDVars(p));
			}
		}
		JDDNode trans01 = JDD.GreaterThan(csgModel.getTrans().copy(), 0);
		enabled = JDD.ThereExists(trans01, csgModel.getAllDDColVars());
		legalC = JDD.ThereExists(enabled.copy(), actNC);
		legalNC = JDD.ThereExists(enabled.copy(), actC);
	}

	/**
	 * Release the current coalition context ({@link #actC}/{@link #actNC}/{@link #enabled}/
	 * {@link #legalC}/{@link #legalNC}), if any -- a no-op if {@link #setCoalition} was never
	 * called. Called both from {@link #setCoalition} itself (releasing the *previous* context
	 * before rebuilding, for a checker instance that gets reused across coalitions -- e.g.
	 * {@code CSGQualitativeValidationTest}) and from {@link #checkExpressionStrategyQual} once
	 * the operator's result has been computed and the context is no longer needed.
	 *
	 * <p>The latter call site is not optional cleanliness: {@code Prism.modelCheck} creates a
	 * <i>fresh</i> {@code CSGModelChecker} instance for every single property checked (confirmed
	 * against {@code Prism.java}'s own model-checking driver), so a checker used via property
	 * syntax only ever calls {@link #setCoalition} once before being discarded -- the "release
	 * the previous context" branch here never fires for that path, and without this explicit
	 * end-of-query release, every property check leaks its own coalition context permanently.
	 * Caught empirically as a real, reproducible leak under {@code -dddebug}:
	 * checking N properties against the same model left exactly N sets of coalition-context
	 * nodes with nonzero references at shutdown.
	 */
	private void clearCoalition()
	{
		if (actC != null) {
			actC.derefAll();
			actNC.derefAll();
			JDD.Deref(enabled);
			JDD.Deref(legalC);
			JDD.Deref(legalNC);
			actC = null;
			actNC = null;
			enabled = null;
			legalC = null;
			legalNC = null;
		}
	}

	private void requireCoalition()
	{
		if (actC == null) {
			throw new IllegalStateException("setCoalition(...) must be called before pre1/G/SF/SFG/SGF");
		}
	}

	/**
	 * Sure-mode one-step controllable-predecessor operator (the {@code Pre1} operator): the
	 * states from which the coalition fixed via {@link #setCoalition} has a legal joint move
	 * such that, against every legal response from the complement, every possible successor is
	 * in {@code x}. Mirrors {@code explicit.CSGModelChecker.pre1} exactly, generalised to the
	 * whole state space at once via {@code ForAll}/{@code ThereExists}/{@code Implies} instead
	 * of per-state matrix games.
	 *
	 * <p><b>The {@link #legalC} guard on the outer existential is mandatory, not an
	 * optimisation</b> (the "legality guard convention"): existential/universal
	 * quantification here ranges over the <i>entire</i> {@code Act_p} bit-vector space,
	 * including codes that aren't actually available at a given state, unlike the explicit
	 * engine's rows/columns, which only ever range over states' actually-enumerated joint
	 * choices. Without the {@code legalC} guard, an illegal coalition move {@code āC} would
	 * make the inner "for all legal complement responses" universally quantified formula
	 * vacuously true (no legal response exists to falsify it against), wrongly treating an
	 * unplayable move as a winning one. Consumes {@code x}.
	 */
	public JDDNode pre1(JDDNode x)
	{
		requireCoalition();
		// "for every legal complement response aNC, every successor of (s,aC,aNC) is in x" -- per (s,aC)
		JDDNode guarded = JDD.ForAll(JDD.Implies(enabled.copy(), subsetOf(x)), actNC);
		// "some legal coalition move aC achieves the above" -- per s
		return JDD.ThereExists(JDD.And(legalC.copy(), guarded), actC);
	}

	/**
	 * "Every positive-probability successor of {@code (s,āC,āNC)} is in {@code y}" -- the
	 * {@code Distribution.isSubsetOf} analogue, shared by {@link #pre1} and {@link #A}. A
	 * {@code (S,ActC,ActNC)}-shaped BDD. Consumes {@code y}.
	 */
	private JDDNode subsetOf(JDDNode y)
	{
		JDDNode yCol = JDD.PermuteVariables(y, csgModel.getAllDDRowVars(), csgModel.getAllDDColVars());
		JDDNode trans01 = JDD.GreaterThan(csgModel.getTrans().copy(), 0);
		return JDD.ForAll(JDD.Implies(trans01, yCol), csgModel.getAllDDColVars());
	}

	/**
	 * "Some positive-probability successor of {@code (s,āC,āNC)} is in {@code x}" -- the
	 * {@code Distribution.containsOneOf} analogue, used by {@link #B}. A
	 * {@code (S,ActC,ActNC)}-shaped BDD. Consumes {@code x}.
	 */
	private JDDNode hitsSome(JDDNode x)
	{
		JDDNode xCol = JDD.PermuteVariables(x, csgModel.getAllDDRowVars(), csgModel.getAllDDColVars());
		JDDNode trans01 = JDD.GreaterThan(csgModel.getTrans().copy(), 0);
		return JDD.ThereExists(JDD.And(trans01, xCol), csgModel.getAllDDColVars());
	}

	/**
	 * {@code A(V2,Y)(s,āC)}: coalition move {@code āC} is "good" w.r.t.
	 * {@code (V2,Y)} -- every legal complement response either keeps the resulting
	 * distribution's support entirely within {@code Y}, or is already excluded via
	 * {@code V2}. Ground truth: {@code explicit.CSGModelChecker.A}.
	 *
	 * <p>Unlike {@link #pre1}'s outer existential, this needs no separate {@link #legalC}
	 * guard: {@code A}'s own quantifier ranges over {@link #actNC} alone, guarded by the
	 * <i>joint</i> {@link #enabled} -- for an illegal {@code āC}, {@code enabled} is false for
	 * every {@code āNC}, so this formula reads (spuriously) true there. That spurious bit is
	 * harmless: every caller of {@code A} (via {@link #B}) always re-conjoins with
	 * {@code enabled} at its own point of use, which washes it out before it can matter -- see
	 * {@link #B}. Consumes {@code v2} and {@code y}.
	 */
	private JDDNode A(JDDNode v2, JDDNode y)
	{
		return JDD.ForAll(JDD.Implies(enabled.copy(), JDD.Or(subsetOf(y), v2)), actNC);
	}

	/**
	 * {@code B(V1,X)(s,āNC)}: complement move {@code āNC} is "threatened" --
	 * some coalition move {@code āC} in {@code V1} has a positive-probability chance of
	 * reaching {@code X} against it. Ground truth: {@code explicit.CSGModelChecker.B}. The
	 * existential over {@link #actC} is guarded directly by the joint {@link #enabled} as a
	 * plain conjunct -- sufficient on its own (no separate {@link #legalC} needed) since
	 * {@code enabled} already excludes illegal {@code āC} regardless of what {@code v1} claims
	 * for it (see {@link #A}'s doc). Consumes {@code v1} and {@code x}.
	 */
	private JDDNode B(JDDNode v1, JDDNode x)
	{
		return JDD.ThereExists(JDD.And(v1, JDD.And(enabled.copy(), hitsSome(x))), actC);
	}

	/**
	 * {@code apreXY(X,Y)(s)}: {@code Apre1(Y,X)}, the almost-sure one-step operator: the
	 * coalition can, against every legal complement response, guarantee
	 * either staying safely within {@code Y} while also threatening to reach {@code X} --
	 * using, in general, a different {@code Y}-safe row per complement response (this is what
	 * makes the guarantee <i>almost</i>-sure rather than sure: repeating the game lets the
	 * coalition eventually pick whichever safe row threatens whatever the complement actually
	 * played). Ground truth: {@code explicit.CSGModelChecker.apreXY}.
	 *
	 * <p>The outer {@code ForAll} needs its own {@link #legalNC} guard: by this point
	 * {@code enabled} has already been quantified away inside {@link #B}, so the vacuous-truth
	 * trap resurfaces here exactly as it does for {@link #pre1}'s outer {@link #legalC} guard
	 * (the "legality guard convention"). Consumes {@code x} and {@code y}.
	 */
	public JDDNode apreXY(JDDNode x, JDDNode y)
	{
		requireCoalition();
		JDDNode goodRows = A(JDD.Constant(0), y);
		JDDNode threatenedCols = B(goodRows, x);
		return JDD.ForAll(JDD.Implies(legalNC.copy(), threatenedCols), actNC);
	}

	/**
	 * {@code apreXYZ(X,Y,Z)(s)}: the almost-sure co-Buchi one-step operator (AFpre1).
	 * Unlike {@link #apreXY}, this needs its own inner fixpoint over the
	 * coalition's own action dimension, parameterised by state -- ground truth:
	 * {@code explicit.CSGModelChecker.apreXYZ}'s own per-state {@code v} while-loop,
	 * generalised here to an {@code (S,ActC)}-shaped BDD fixpoint covering every state at
	 * once. Still just the same "iterate a BDD to a fixed point" idiom used throughout this
	 * class (e.g. {@link #G}/{@link #SF}), not a fundamentally new technique -- no extra
	 * scratch variable copies are needed the way a genuine relational-composition fixpoint
	 * (like reachability's row/col alternation) would need, since each round recomputes
	 * directly from the current {@code v} via {@link #A}/{@link #B}.
	 *
	 * <p>{@code v_0 = legalC} (all legal coalition moves, matching the explicit engine's
	 * {@code v.set(0, mdist.size())} -- "all enumerated rows"); {@code v_(k+1) = A(0,z) AND
	 * A(B(v_k,x),y) AND legalC}, monotonically shrinking (a nu-fixpoint over the row
	 * dimension) until it stabilises. The result, {@code ThereExists(v_final, actC)}, is
	 * exactly "{@code !v.isEmpty()}" per state. Consumes {@code x}, {@code y} and {@code z}.
	 *
	 * <p><b>The extra {@code legalC} conjunct on every round is mandatory, not defensive
	 * belt-and-braces</b> -- a second instance of the "legality guard
	 * convention", distinct from the ones already accounted for inside {@link #A}/{@link #B}
	 * themselves. {@link #A}'s raw output reads (spuriously) true at an illegal {@code āC}
	 * (see {@link #A}'s own doc) and normally gets washed out by {@link #B}'s own {@code enabled}
	 * conjunct wherever {@code A}'s result feeds into {@code B} -- which is every use in
	 * {@link #apreXY}, but NOT here: {@code a1 = A(0,z)} and {@code a2 = A(B(v,x),y)} are both
	 * consumed <i>directly</i> into {@code next}, with no intervening {@code B} call to wash
	 * either one. Without the extra conjunct, those spurious illegal-{@code āC} bits creep
	 * into {@code v} on the very first round (overwriting the clean {@code v_0 = legalC}) and
	 * corrupt the final {@code ThereExists(v_final, actC)} into reading true at states that
	 * have no genuine legal witness at all -- caught empirically as a real bug during
	 * validation via a direct {@code A}/{@code B}/{@code apreXYZ} minterm
	 * trace on {@code skirmish.prism}, not spotted by inspection alone.
	 */
	public JDDNode apreXYZ(JDDNode x, JDDNode y, JDDNode z)
	{
		requireCoalition();
		JDDNode v = legalC.copy();
		while (true) {
			JDDNode a1 = A(JDD.Constant(0), z.copy());
			JDDNode a2 = A(B(v.copy(), x.copy()), y.copy());
			JDDNode next = JDD.And(a1, JDD.And(a2, legalC.copy()));
			if (next.equals(v)) {
				JDD.Deref(next);
				break;
			}
			JDD.Deref(v);
			v = next;
		}
		JDD.Deref(x);
		JDD.Deref(y);
		JDD.Deref(z);
		return JDD.ThereExists(v, actC);
	}

	/**
	 * {@code lpreXY(X,Y)(s)}: {@code Lpre1(Y,X)}, the limit-sure one-step operator. Unlike
	 * {@link #apreXY}, this needs its own inner fixpoint over the coalition's own
	 * action dimension -- ground truth: {@code explicit.CSGModelChecker.lpreXY}'s per-state
	 * {@code w} while-loop, generalised to an {@code (S,ActC)}-shaped BDD fixpoint as usual.
	 * {@code w_0 = 0}; {@code w_(k+1) = A(B(w_k,X),Y)}, until it stabilises; result =
	 * {@code ForAll(actNC, Implies(legalNC, B(w_final,X)))}.
	 *
	 * <p><b>No extra {@code legalC} conjunct needed here, unlike {@link #apreXYZ}</b> -- the
	 * legality-guard placement has to be re-derived per primitive, not
	 * assumed from a sibling. Every use of {@code A}'s raw output in this method (both
	 * {@code w} feeding the next round's {@code B(w,X)}, and the final {@code B(w_final,X)})
	 * is immediately re-consumed by another {@link #B} call, whose own {@code enabled} conjunct
	 * washes out {@code A}'s illegal-{@code āC} vacuous-truth artifact regardless of what
	 * {@code w} itself contains there (see {@link #B}'s doc) -- unlike {@link #apreXYZ}, where
	 * {@code A}'s raw output was combined directly with no intervening {@code B}. Consumes
	 * {@code x} and {@code y}.
	 */
	public JDDNode lpreXY(JDDNode x, JDDNode y)
	{
		requireCoalition();
		JDDNode w = JDD.Constant(0);
		while (true) {
			JDDNode next = A(B(w.copy(), x.copy()), y.copy());
			if (next.equals(w)) {
				JDD.Deref(next);
				break;
			}
			JDD.Deref(w);
			w = next;
		}
		JDDNode result = JDD.ForAll(JDD.Implies(legalNC.copy(), B(w, x)), actNC);
		JDD.Deref(y);
		return result;
	}

	/**
	 * {@code lpreXYZ(X,Y,Z)(s)}: the limit-sure co-Buchi one-step operator.
	 * Ground truth: {@code explicit.CSGModelChecker.lpreXYZ}'s doubly-nested {@code v}/{@code w}
	 * while-loop -- an outer {@code v} fixpoint (starting at {@link #legalC}, matching
	 * {@link #apreXYZ}'s own {@code v_0}), each round of which runs a full inner {@code w}
	 * fixpoint (starting at {@code 0} every outer round) using an {@code ay = A(B(v,X),Y)} term
	 * held fixed for that whole inner fixpoint: {@code w_(k+1) = A(B(w_k,X),Z) AND ay}, then
	 * {@code v <- w_final} once the inner fixpoint stabilises, until <i>that</i> stabilises too.
	 * Result: {@code ThereExists(v_final, actC)}.
	 *
	 * <p><b>The extra {@code legalC} conjunct on the inner {@code sol} is mandatory</b>, exactly
	 * {@link #apreXYZ}'s fix, re-derived fresh here rather than assumed: {@code sol = az AND ay}
	 * (both raw {@link #A} outputs) becomes {@code w} directly, with no intervening {@link #B}
	 * call, and {@code w} eventually becomes {@code v} -- whose <i>final</i> value feeds the
	 * same raw, unwashed {@code ThereExists(v, actC)} that {@link #apreXYZ} does. Without the
	 * guard, the same illegal-{@code āC} contamination would corrupt the result the same way.
	 * Consumes {@code x}, {@code y} and {@code z}.
	 */
	public JDDNode lpreXYZ(JDDNode x, JDDNode y, JDDNode z)
	{
		requireCoalition();
		JDDNode v = legalC.copy();
		while (true) {
			JDDNode ay = A(B(v.copy(), x.copy()), y.copy());
			JDDNode w = JDD.Constant(0);
			while (true) {
				JDDNode az = A(B(w.copy(), x.copy()), z.copy());
				JDDNode sol = JDD.And(az, JDD.And(ay.copy(), legalC.copy()));
				if (sol.equals(w)) {
					JDD.Deref(sol);
					break;
				}
				JDD.Deref(w);
				w = sol;
			}
			JDD.Deref(ay);
			if (w.equals(v)) {
				JDD.Deref(w);
				break;
			}
			JDD.Deref(v);
			v = w;
		}
		JDD.Deref(x);
		JDD.Deref(y);
		JDD.Deref(z);
		return JDD.ThereExists(v, actC);
	}

	/**
	 * Sure-mode safety: {@code nu X. (b AND pre1(X))}. Mirrors
	 * {@code explicit.CSGModelChecker.G} exactly. Consumes {@code b}.
	 */
	public JDDNode G(JDDNode b)
	{
		JDDNode x = JDD.Constant(1);
		while (true) {
			JDDNode next = JDD.And(b.copy(), pre1(x.copy()));
			if (next.equals(x)) {
				JDD.Deref(next);
				break;
			}
			JDD.Deref(x);
			x = next;
		}
		JDD.Deref(b);
		return x;
	}

	/**
	 * Sure-mode reachability: {@code mu X. (b OR pre1(X))}. Mirrors
	 * {@code explicit.CSGModelChecker.SF} exactly. Consumes {@code b}.
	 */
	public JDDNode SF(JDDNode b)
	{
		JDDNode x = JDD.Constant(0);
		while (true) {
			JDDNode next = JDD.Or(b.copy(), pre1(x.copy()));
			if (next.equals(x)) {
				JDD.Deref(next);
				break;
			}
			JDD.Deref(x);
			x = next;
		}
		JDD.Deref(b);
		return x;
	}

	/**
	 * Sure-mode co-Buchi (F G b): reach the maximal b-invariant, then stay. {@code SF(G(b))}
	 * exactly, mirroring {@code explicit.CSGModelChecker.SFG}'s own decomposition -- see that
	 * method's doc comment for why this decomposition (rather than a coupled fixpoint) is sound
	 * specifically in sure mode. Consumes {@code b}.
	 */
	public JDDNode SFG(JDDNode b)
	{
		return SF(G(b));
	}

	/**
	 * Sure-mode Buchi (G F b): visit b infinitely often, forced with certainty. Mirrors
	 * {@code explicit.CSGModelChecker.SGF}'s coupled {@code nu Y. mu X. [(b AND pre1(Y)) OR
	 * (not b AND pre1(X))]} exactly. Consumes {@code b}.
	 */
	public JDDNode SGF(JDDNode b)
	{
		JDDNode notB = JDD.Not(b.copy());
		JDDNode y = JDD.Constant(1);
		while (true) {
			JDDNode x = JDD.Constant(0);
			while (true) {
				JDDNode fromB = JDD.And(b.copy(), pre1(y.copy()));
				JDDNode fromNotB = JDD.And(notB.copy(), pre1(x.copy()));
				JDDNode next = JDD.Or(fromB, fromNotB);
				if (next.equals(x)) {
					JDD.Deref(next);
					break;
				}
				JDD.Deref(x);
				x = next;
			}
			if (x.equals(y)) {
				JDD.Deref(x);
				break;
			}
			JDD.Deref(y);
			y = x;
		}
		JDD.Deref(b);
		JDD.Deref(notB);
		return y;
	}

	/**
	 * Almost-sure reachability: {@code nu Y. mu X. (b OR apreXY(X,Y))}. Mirrors
	 * {@code explicit.CSGModelChecker.AF} exactly. Consumes {@code b}.
	 */
	public JDDNode AF(JDDNode b)
	{
		JDDNode y = JDD.Constant(1);
		while (true) {
			JDDNode x = JDD.Constant(0);
			while (true) {
				JDDNode next = JDD.Or(b.copy(), apreXY(x.copy(), y.copy()));
				if (next.equals(x)) {
					JDD.Deref(next);
					break;
				}
				JDD.Deref(x);
				x = next;
			}
			if (x.equals(y)) {
				JDD.Deref(x);
				break;
			}
			JDD.Deref(y);
			y = x;
		}
		JDD.Deref(b);
		return y;
	}

	/**
	 * Almost-sure co-Buchi (F G b): eventually stay in b forever, forced almost-surely.
	 * Mirrors {@code explicit.CSGModelChecker.AFG}'s triple-nested fixpoint exactly:
	 * {@code Z = nu Z. mu X. nu Y. [(b AND apreXYZ(X,Y,Z)) OR (not b AND apreXY(X,Z))]}.
	 * Consumes {@code b}.
	 */
	public JDDNode AFG(JDDNode b)
	{
		JDDNode notB = JDD.Not(b.copy());
		JDDNode z = JDD.Constant(1);
		while (true) {
			JDDNode x = JDD.Constant(0);
			while (true) {
				JDDNode y = JDD.Constant(1);
				while (true) {
					JDDNode fromB = JDD.And(b.copy(), apreXYZ(x.copy(), y.copy(), z.copy()));
					JDDNode fromNotB = JDD.And(notB.copy(), apreXY(x.copy(), z.copy()));
					JDDNode next = JDD.Or(fromB, fromNotB);
					if (next.equals(y)) {
						JDD.Deref(next);
						break;
					}
					JDD.Deref(y);
					y = next;
				}
				if (y.equals(x)) {
					JDD.Deref(y);
					break;
				}
				JDD.Deref(x);
				x = y;
			}
			if (x.equals(z)) {
				JDD.Deref(x);
				break;
			}
			JDD.Deref(z);
			z = x;
		}
		JDD.Deref(b);
		JDD.Deref(notB);
		return z;
	}

	/**
	 * Almost-sure Buchi (G F b): visit b infinitely often, forced almost-surely. Mirrors
	 * {@code explicit.CSGModelChecker.AGF}'s coupled {@code nu Y. mu X. [(not b AND
	 * apreXY(X,Y)) OR (b AND pre1(Y))]} exactly -- the same {@link #pre1}/{@link #apreXY}
	 * substitution relative to {@link #SGF} that the explicit engine's own AGF makes relative
	 * to SGF. Consumes {@code b}.
	 */
	public JDDNode AGF(JDDNode b)
	{
		JDDNode notB = JDD.Not(b.copy());
		JDDNode y = JDD.Constant(1);
		while (true) {
			JDDNode x = JDD.Constant(0);
			while (true) {
				JDDNode fromNotB = JDD.And(notB.copy(), apreXY(x.copy(), y.copy()));
				JDDNode fromB = JDD.And(b.copy(), pre1(y.copy()));
				JDDNode next = JDD.Or(fromNotB, fromB);
				if (next.equals(x)) {
					JDD.Deref(next);
					break;
				}
				JDD.Deref(x);
				x = next;
			}
			if (x.equals(y)) {
				JDD.Deref(x);
				break;
			}
			JDD.Deref(y);
			y = x;
		}
		JDD.Deref(b);
		JDD.Deref(notB);
		return y;
	}

	/**
	 * Limit-sure reachability: {@code nu Y. mu X. (b OR lpreXY(X,Y))}. Mirrors
	 * {@code explicit.CSGModelChecker.LF} exactly. Consumes {@code b}.
	 */
	public JDDNode LF(JDDNode b)
	{
		JDDNode y = JDD.Constant(1);
		while (true) {
			JDDNode x = JDD.Constant(0);
			while (true) {
				JDDNode next = JDD.Or(b.copy(), lpreXY(x.copy(), y.copy()));
				if (next.equals(x)) {
					JDD.Deref(next);
					break;
				}
				JDD.Deref(x);
				x = next;
			}
			if (x.equals(y)) {
				JDD.Deref(x);
				break;
			}
			JDD.Deref(y);
			y = x;
		}
		JDD.Deref(b);
		return y;
	}

	/**
	 * Limit-sure co-Buchi (F G b): eventually stay in b forever, forced in the limit. Mirrors
	 * {@code explicit.CSGModelChecker.LFG}'s triple-nested fixpoint exactly:
	 * {@code Z = nu Z. mu X. nu Y. [(b AND lpreXYZ(X,Y,Z)) OR (not b AND lpreXY(X,Z))]}.
	 * Consumes {@code b}.
	 */
	public JDDNode LFG(JDDNode b)
	{
		JDDNode notB = JDD.Not(b.copy());
		JDDNode z = JDD.Constant(1);
		while (true) {
			JDDNode x = JDD.Constant(0);
			while (true) {
				JDDNode y = JDD.Constant(1);
				while (true) {
					JDDNode fromB = JDD.And(b.copy(), lpreXYZ(x.copy(), y.copy(), z.copy()));
					JDDNode fromNotB = JDD.And(notB.copy(), lpreXY(x.copy(), z.copy()));
					JDDNode next = JDD.Or(fromB, fromNotB);
					if (next.equals(y)) {
						JDD.Deref(next);
						break;
					}
					JDD.Deref(y);
					y = next;
				}
				if (y.equals(x)) {
					JDD.Deref(y);
					break;
				}
				JDD.Deref(x);
				x = y;
			}
			if (x.equals(z)) {
				JDD.Deref(x);
				break;
			}
			JDD.Deref(z);
			z = x;
		}
		JDD.Deref(b);
		JDD.Deref(notB);
		return z;
	}

	/**
	 * Limit-sure Buchi (G F b): visit b infinitely often, forced in the limit. Mirrors
	 * {@code explicit.CSGModelChecker.LGF}'s coupled {@code nu Y. mu X. [(not b AND
	 * lpreXY(X,Y)) OR (b AND pre1(Y))]} exactly. Consumes {@code b}.
	 */
	public JDDNode LGF(JDDNode b)
	{
		JDDNode notB = JDD.Not(b.copy());
		JDDNode y = JDD.Constant(1);
		while (true) {
			JDDNode x = JDD.Constant(0);
			while (true) {
				JDDNode fromNotB = JDD.And(notB.copy(), lpreXY(x.copy(), y.copy()));
				JDDNode fromB = JDD.And(b.copy(), pre1(y.copy()));
				JDDNode next = JDD.Or(fromNotB, fromB);
				if (next.equals(x)) {
					JDD.Deref(next);
					break;
				}
				JDD.Deref(x);
				x = next;
			}
			if (x.equals(y)) {
				JDD.Deref(x);
				break;
			}
			JDD.Deref(y);
			y = x;
		}
		JDD.Deref(b);
		JDD.Deref(notB);
		return y;
	}
}
