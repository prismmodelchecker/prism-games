package parser.ast;

import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

import prism.PrismLangException;

/**
 * Determines, for a CSG (concurrent stochastic game) {@link ModulesFile}, which player owns
 * each module and each synchronising action, and validates the CSG-specific well-formedness
 * rules that this ownership assignment depends on:
 * <ul>
 * <li>every command in a player-owned module must be labelled (carry an action);</li>
 * <li>action ownership is disjoint across players (an action index claimed by one player's
 *     module cannot also be claimed by another player's module) -- though a single player's
 *     own action may legitimately be claimed by more than one of that player's own modules,
 *     since a player's modules must react jointly to that player's actions;</li>
 * <li>independent (non-player) modules may only reference actions already owned by some
 *     player.</li>
 * </ul>
 *
 * This is pure AST-level analysis with no dependency on any particular model-checking
 * engine (explicit-state or symbolic), so it is meant to be shared rather than
 * re-implemented per engine. It does not consult the explicit {@code [action]} item some
 * player declarations may carry (see {@link ModulesFile#getPlayerForAction}) -- for CSGs,
 * ownership is always inferred from which player-owned module's commands use a given action
 * label, never from an explicit declaration; that mechanism exists for SMGs, which can leave
 * ownership of a synchronising action genuinely ambiguous from module structure alone, a
 * situation CSG's disjointness requirement rules out by construction.
 *
 * <br>Ported out of {@code simulator.Updater#initialiseCSG()}, which previously computed
 * this independently; both now delegate to this single implementation so that the
 * explicit-state engine and the (future) symbolic CSG builder cannot silently disagree about
 * action ownership or well-formedness.
 */
public class CSGPlayerActions
{
	private final ModulesFile modulesFile;
	private final int numPlayers;
	private final int numModules;

	/** For each module index, the (0-indexed) owning player, or -1 if the module is independent. */
	private final int[] playersIndexes;
	/** For each player index, the {@link BitSet} of synchronising-action indices that player owns. */
	private final BitSet[] playersActionsIndexes;
	/** Map from synchronising-action index to its (0-indexed) owning player index. */
	private final Map<Integer, Integer> actionIndexPlayerMap;

	/**
	 * Compute and validate action ownership for {@code modulesFile}.
	 * @throws PrismLangException if the model violates a CSG action-ownership well-formedness rule
	 */
	public CSGPlayerActions(ModulesFile modulesFile) throws PrismLangException
	{
		this.modulesFile = modulesFile;
		this.numPlayers = modulesFile.getNumPlayers();
		this.numModules = modulesFile.getNumModules();
		this.playersIndexes = new int[numModules];
		this.playersActionsIndexes = new BitSet[numPlayers];
		this.actionIndexPlayerMap = new HashMap<Integer, Integer>();
		Arrays.fill(playersIndexes, -1);
		for (int p = 0; p < numPlayers; p++) {
			playersActionsIndexes[p] = new BitSet();
		}
		if (numPlayers > 0) {
			computeOwnership();
			validateIndependentModules();
		}
	}

	/**
	 * For each player-owned module, assign every command's (primary) action index to that
	 * player, checking that player-owned commands are labelled and that action ownership is
	 * disjoint across (but not necessarily within) players.
	 */
	private void computeOwnership() throws PrismLangException
	{
		BitSet seen = new BitSet();
		for (int m = 0; m < numModules; m++) {
			playersIndexes[m] = modulesFile.getPlayerForModule(modulesFile.getModuleName(m));
			if (playersIndexes[m] == -1) {
				continue;
			}
			Module module = modulesFile.getModule(m);
			int numCommands = module.getNumCommands();
			for (int c = 0; c < numCommands; c++) {
				Command command = module.getCommand(c);
				if (command.isUnlabelled()) {
					throw new PrismLangException("Commands in a player-owned module cannot be unlabelled", command);
				}
				int index = command.getSynchIndices().get(0);
				if (!seen.get(index) || playersActionsIndexes[playersIndexes[m]].get(index)) {
					seen.set(index);
					playersActionsIndexes[playersIndexes[m]].set(index);
					actionIndexPlayerMap.put(index, playersIndexes[m]);
				} else {
					throw new PrismLangException("Action " + command.getSynchs().get(0) + " of module "
							+ module.getName()
							+ " had already been associated to a different module. Action sets must be disjoint");
				}
			}
		}
	}

	/**
	 * Independent (non-player) modules may only synchronise on actions some player already owns.
	 */
	private void validateIndependentModules() throws PrismLangException
	{
		for (int m = 0; m < numModules; m++) {
			if (playersIndexes[m] != -1) {
				continue;
			}
			Module module = modulesFile.getModule(m);
			int numCommands = module.getNumCommands();
			for (int c = 0; c < numCommands; c++) {
				for (int i : module.getCommand(c).getSynchIndices()) {
					if (!actionIndexPlayerMap.containsKey(i) && i != 0) {
						throw new PrismLangException("Label \"" + modulesFile.getSynch(i - 1) + "\" of command " + c
								+ " of module " + modulesFile.getModuleName(m) + " is not associated to any player."
								+ " Independent modules can only synchronize on players' actions.");
					}
				}
			}
		}
	}

	/**
	 * The (0-indexed) player owning module {@code m}, or -1 if {@code m} is independent.
	 */
	public int getPlayerForModule(int m)
	{
		return playersIndexes[m];
	}

	/**
	 * The {@link BitSet} of synchronising-action indices owned by player {@code p}.
	 */
	public BitSet getActionsForPlayer(int p)
	{
		return playersActionsIndexes[p];
	}

	/**
	 * The (0-indexed) player owning synchronising-action index {@code actionIndex}, or
	 * {@code null} if unowned (which can only happen for index 0, the unlabelled/idle case).
	 */
	public Integer getPlayerForActionIndex(int actionIndex)
	{
		return actionIndexPlayerMap.get(actionIndex);
	}

	/**
	 * For each module index (0 to {@code modulesFile.getNumModules() - 1}), the owning
	 * player, or -1 if independent. Returned array is a live reference, not a copy, to match
	 * the field it replaces in {@code simulator.Updater}; callers should not mutate it.
	 */
	public int[] getPlayersIndexes()
	{
		return playersIndexes;
	}

	/**
	 * For each player index, the {@link BitSet} of synchronising-action indices that player
	 * owns. Returned array is a live reference, not a copy, to match the field it replaces in
	 * {@code simulator.Updater}; callers should not mutate it.
	 */
	public BitSet[] getPlayersActionsIndexes()
	{
		return playersActionsIndexes;
	}

	/**
	 * Map from synchronising-action index to its owning player index, unmodifiable.
	 */
	public Map<Integer, Integer> getActionIndexPlayerMap()
	{
		return Collections.unmodifiableMap(actionIndexPlayerMap);
	}
}
