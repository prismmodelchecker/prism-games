package simulator;

import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

import parser.State;
import parser.Values;
import parser.VarList;
import parser.ast.Expression;
import parser.ast.LabelList;
import parser.ast.ModulesFile;
import parser.ast.Player;
import parser.ast.RewardStruct;
import parser.type.Type;
import prism.DefaultModelGenerator;
import prism.ModelType;
import prism.PrismComponent;
import prism.PrismException;
import prism.PrismLangException;

public class ModulesFileModelGenerator extends DefaultModelGenerator
{
	// Parent PrismComponent (logs, settings etc.)
	protected PrismComponent parent;
	
	// PRISM model info
	/** The original modules file (might have unresolved constants) */
	private ModulesFile originalModulesFile;
	/** The modules file used for generating (has no unresolved constants after {@code initialise}) */
	private ModulesFile modulesFile;
	private ModelType modelType;
	private Values mfConstants;
	private VarList varList;
	private LabelList labelList;
	
	// Wwhether to use the "compositional" system construction (for SMGSs)
	private boolean compositional = false;

	// Model exploration info
	
	// State currently being explored
	private State exploreState;
	// Updater object for model
	protected Updater updater;
	// List of currently available transitions
	protected TransitionList transitionList;
	// Has the transition list been built? 
	protected boolean transitionListBuilt;
	
	/**
	 * Build a ModulesFileModelGenerator for a particular PRISM model, represented by a ModuleFile instance.
	 * @param modulesFile The PRISM model
	 */
	public ModulesFileModelGenerator(ModulesFile modulesFile) throws PrismException
	{
		this(modulesFile, null);
	}
	
	/**
	 * Build a ModulesFileModelGenerator for a particular PRISM model, represented by a ModuleFile instance.
	 * @param modulesFile The PRISM model
	 */
	public ModulesFileModelGenerator(ModulesFile modulesFile, PrismComponent parent) throws PrismException
	{
		this.parent = parent;
		
		// No support for PTAs yet
		if (modulesFile.getModelType() == ModelType.PTA) {
			throw new PrismException("Sorry - the simulator does not currently support PTAs");
		}
		// No support for system...endsystem yet
		if (modulesFile.getSystemDefn() != null) {
			throw new PrismException("Sorry - the simulator does not currently handle the system...endsystem construct");
		}
		
		// Store basic model info
		this.modulesFile = modulesFile;
		this.originalModulesFile = modulesFile;
		modelType = modulesFile.getModelType();
		
		// If there are no constants to define, go ahead and initialise;
		// Otherwise, setSomeUndefinedConstants needs to be called when the values are available  
		mfConstants = modulesFile.getConstantValues();
		if (mfConstants != null) {
			initialise();
		}
	}
	
	/**
	 * (Re-)Initialise the class ready for model exploration
	 * (can only be done once any constants needed have been provided)
	 */
	private void initialise() throws PrismLangException
	{
		// Evaluate constants on (a copy) of the modules file, insert constant values and optimize arithmetic expressions
		modulesFile = (ModulesFile) modulesFile.deepCopy().replaceConstants(mfConstants).simplify();

		// Get info
		varList = modulesFile.createVarList();
		labelList = modulesFile.getLabelList();
		
		// Create data structures for exploring model
		updater = new Updater(modulesFile, varList, parent);
		transitionList = new TransitionList();
		transitionListBuilt = false;
	}
	
	/**
	 * Get access to the ModulesFile being used to generate the model.
	 * @return
	 */
	public ModulesFile getModulesFile()
	{
		return modulesFile;
	}
	
	/**
	 * Set whether to use the "compositional" system construction (for SMGSs) 
	 */
	public void setCompositional(boolean compositional)
	{
		this.compositional = compositional;
	}
	
	// Model info
	
	@Override
	public ModelType getModelType()
	{
		return modelType;
	}
	
	@Override
	public void setSomeUndefinedConstants(Values someValues) throws PrismException
	{
		// We start again with a copy of the original modules file
		// and set the constants in the copy.
		// As {@code initialise()} can replace references to constants
		// with the concrete values in modulesFile, this ensures that we
		// start again at a place where references to constants have not
		// yet been replaced.
		modulesFile = (ModulesFile) originalModulesFile.deepCopy();
		modulesFile.setSomeUndefinedConstants(someValues);
		mfConstants = modulesFile.getConstantValues();
		initialise();
	}
	
	@Override
	public Values getConstantValues()
	{
		return mfConstants;
	}
	
	@Override
	public boolean containsUnboundedVariables()
	{
		return modulesFile.containsUnboundedVariables();
	}
	
	@Override
	public int getNumVars()
	{
		return modulesFile.getNumVars();
	}
	
	@Override
	public List<String> getVarNames()
	{
		return modulesFile.getVarNames();
	}

	@Override
	public List<Type> getVarTypes()
	{
		return modulesFile.getVarTypes();
	}

	@Override
	public int getNumLabels()
	{
		return labelList.size();	
	}
	
	@Override
	public String getLabelName(int i) throws PrismException
	{
		return labelList.getLabelName(i);
	}
	
	@Override
	public int getLabelIndex(String label)
	{
		return labelList.getLabelIndex(label);
	}
	
	@Override
	public int getNumPlayers()
	{
		return modulesFile.getNumPlayers();
	}

	@Override
	public Player getPlayer(int i)
	{
		return modulesFile.getPlayer(i);
	}
	
	//@Override
	public VarList getVarList()
	{
		return varList;
	}
	
	@Override
	public int getNumRewardStructs()
	{
		return modulesFile.getNumRewardStructs();
	}
	
	@Override
	public int getRewardStructIndex(String name)
	{
		return modulesFile.getRewardStructIndex(name);
	}
	
	@Override
	public RewardStruct getRewardStruct(int i)
	{
		return modulesFile.getRewardStruct(i);
	}

	@Override
	public boolean hasSingleInitialState() throws PrismException
	{
		return modulesFile.getInitialStates() == null;
	}
	
	@Override
	public State getInitialState() throws PrismException
	{
		if (modulesFile.getInitialStates() == null) {
			return modulesFile.getDefaultInitialState();
		} else {
			// Inefficient but probably won't be called
			return getInitialStates().get(0);
		}
	}
	
	@Override
	public List<State> getInitialStates() throws PrismException
	{
		List<State> initStates = new ArrayList<State>();
		// Easy (normal) case: just one initial state
		if (modulesFile.getInitialStates() == null) {
			State state = modulesFile.getDefaultInitialState();
			initStates.add(state);
		}
		// Otherwise, there may be multiple initial states
		// For now, we handle this is in a very inefficient way
		else {
			Expression init = modulesFile.getInitialStates();
			List<State> allPossStates = varList.getAllStates();
			for (State possState : allPossStates) {
				if (init.evaluateBoolean(modulesFile.getConstantValues(), possState)) {
					initStates.add(possState);
				}
			}
		}
		return initStates;
	}

	@Override
	public void exploreState(State exploreState) throws PrismException
	{
		this.exploreState = exploreState;
		transitionListBuilt = false;
	}
	
	@Override
	public State getExploreState()
	{
		return exploreState;
	}
	
	@Override
	public int getNumChoices() throws PrismException
	{
		return getTransitionList().getNumChoices();
	}

	@Override
	public int getNumTransitions() throws PrismException
	{
		return getTransitionList().getNumTransitions();
	}

	@Override
	public int getNumTransitions(int index) throws PrismException
	{
		return getTransitionList().getChoice(index).size();
	}

	@Override
	public String getTransitionAction(int index, int offset) throws PrismException
	{
		TransitionList transitions = getTransitionList();
		int a = transitions.getTransitionModuleOrActionIndex(transitions.getTotalIndexOfTransition(index, offset));
		return a < 0 ? null : modulesFile.getSynch(a - 1);
	}

	//@Override
	public String getTransitionAction(int index) throws PrismException
	{
		int a = getTransitionList().getTransitionModuleOrActionIndex(index);
		return a < 0 ? null : modulesFile.getSynch(a - 1);
	}

	public int getTransitionModuleOrActionIndex(int i, int offset) throws PrismException
	{
		TransitionList transitions = getTransitionList();
		return transitions.getTransitionModuleOrActionIndex(transitions.getTotalIndexOfTransition(i, offset));
	}
	
	public String getTransitionModuleOrAction(int i, int offset) throws PrismException
	{
		TransitionList transitions = getTransitionList();
		return transitions.getTransitionModuleOrAction(transitions.getTotalIndexOfTransition(i, offset));
	}
	
	@Override
	public double getTransitionProbability(int index, int offset) throws PrismException
	{
		TransitionList transitions = getTransitionList();
		return transitions.getChoice(index).getProbability(offset);
	}

	//@Override
	public double getTransitionProbability(int index) throws PrismException
	{
		TransitionList transitions = getTransitionList();
		return transitions.getTransitionProbability(index);
	}

	@Override
	public State computeTransitionTarget(int index, int offset) throws PrismException
	{
		return getTransitionList().getChoice(index).computeTarget(offset, exploreState);
	}

	//@Override
	public State computeTransitionTarget(int index) throws PrismException
	{
		return getTransitionList().computeTransitionTarget(index, exploreState);
	}
	
	@Override
    public int getPlayerNumberForChoice(int i) throws PrismException
	{
    		// Normal game
		if (!compositional) {
			String modAct = getTransitionModuleOrAction(i, 0);
			// NB: 0-indexed to 1-indexed
			return modulesFile.getPlayerForModuleOrAction(modAct) + 1;
		}
		// Compositional game
		else {
			// Get index of module/action (all transitions within choice have same action, so use offset 0)
			int modAct = getTransitionModuleOrActionIndex(i, 0);

			// Synchronous action
			int player = -1;
			if (modAct > 0) {
				String actionName = modulesFile.getSynch(modAct - 1);
				// get inputs and outputs of subsystem
				Vector<String> inputs = new Vector<String>();
				Vector<String> outputs = new Vector<String>();
				for (int n = 0; n < modulesFile.getNumModules(); n++) {
					inputs.addAll(modulesFile.getModule(n).getAllInputActions());
					outputs.addAll(modulesFile.getModule(n).getAllOutputActions());
				}
				if (actionName == null || "".equals(actionName)) {
					player = 2; // tau controlled by player 2
				} else if (inputs.contains(actionName)) {
					player = 2; // inputs controlled by player 2
				} else if (outputs.contains(actionName)) {
					player = 1; // outputs controlled by player 1
				} else {
					throw new PrismException("Action \"" + actionName + "\" is not assigned to any player");
				}
			}
			// Asynchronous action
			else {
				player = 2; // tau controlled by player 2
			}
			return player;
		}
	}

	//@Override
	public void calculateStateRewards(State state, double[] store) throws PrismLangException
	{
		updater.calculateStateRewards(state, store);
	}
	
	//@Override
	public void getRandomInitialState(RandomNumberGenerator rng, State initialState) throws PrismException
	{
		if (modulesFile.getInitialStates() == null) {
			initialState.copy(modulesFile.getDefaultInitialState());
		} else {
			throw new PrismException("Random choice of multiple initial states not yet supported");
		}
	}

	@Override
	public boolean isLabelTrue(int i) throws PrismException
	{
		Expression expr = labelList.getLabel(i);
		return expr.evaluateBoolean(exploreState);
	}
	
	@Override
	public double getStateReward(int index, State state) throws PrismException
	{
		RewardStruct rewStr = modulesFile.getRewardStruct(index);
		int n = rewStr.getNumItems();
		double d = 0;
		for (int i = 0; i < n; i++) {
			Expression guard = rewStr.getStates(i);
			if (guard.evaluateBoolean(modulesFile.getConstantValues(), state)) {
				double rew = rewStr.getReward(i).evaluateDouble(modulesFile.getConstantValues(), state);
				if (Double.isNaN(rew))
					throw new PrismLangException("Reward structure evaluates to NaN at state " + state, rewStr.getReward(i));
				d += rew;
			}
		}
		return d;
	}
	
	// Local utility methods
	
	/**
	 * Returns the current list of available transitions, generating it first if this has not yet been done.
	 */
	private TransitionList getTransitionList() throws PrismException
	{
		// Compute the current transition list, if required
		if (!transitionListBuilt) {
			updater.calculateTransitions(exploreState, transitionList);
			transitionListBuilt = true;
		}
		return transitionList;
	}

	@Override
	public VarList createVarList()
	{
		return varList;
	}
}
