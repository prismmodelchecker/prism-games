//==============================================================================
//	
//	Copyright (c) 2014-
//	Authors:
//	* Clemens Wiltsche <clemens.wiltsche@cs.ox.ac.uk> (University of Oxford)
//	* Dave Parker <d.a.parker@cs.bham.ac.uk> (University of Birmingham/Oxford)
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

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import parma_polyhedra_library.C_Polyhedron;
import parma_polyhedra_library.Partial_Function;
import parma_polyhedra_library.Polyhedron;
import parser.Values;
import parser.ast.Expression;
import parser.ast.ExpressionFunc;
import parser.ast.ExpressionLiteral;
import parser.ast.LabelList;
import parser.ast.ModulesFile;
import parser.ast.PropertiesFile;
import parser.ast.SystemBrackets;
import parser.ast.SystemDefn;
import parser.ast.SystemFullParallel;
import parser.ast.SystemModule;
import parser.ast.SystemParallel;
import parser.ast.SystemReference;
import prism.PointList;
import prism.PrismComponent;
import prism.PrismException;
import prism.PrismSettings;
import prism.Result;
import simulator.SimulatorEngine;
import strat.StochasticUpdateStrategyProduct;
import strat.Strategy;

/**
 * Model checker for compositional methods on SMGs
 */
public class CompositionalSMGModelChecker extends PrismComponent
{
	// Model file
	private ModulesFile modulesFile;
	// Properties file
	private PropertiesFile propertiesFile;
	// Simulator engine
	private SimulatorEngine engine;
	// Constants from model
	private Values constantValues;
	// Labels from the model
	private LabelList labelListModel;
	// Labels from the property file
	private LabelList labelListProp;

	// the SMG model checker for the subtasks
	private SMGModelChecker mc;

	// options
	private boolean computePareto = true; // computing Pareto set, or doing actual model checking?
        private boolean checkCompatibility = false; // check compatibility of subsystems
	private boolean[] cancel_computation = new boolean[1]; // false by default

	// the subsystems, their properties, and the return values of the separate subtasks
	private List<ModulesFile> subsystemModulesFiles; // one modules File per subsystem
	private List<SMG> subsystems; // the explicit models, one per subsystem
	private List<Pareto> subsystemParetos; // one Pareto set for the first initial state of each subsystem
        private List<MultiParameters> subsystemParams; // one parameter object per query for displaying the Pareto set
	private List<Strategy> subsystemStrategies; // on computed strategy per subsystem

        // only simulate strategies
        private boolean simulateStrategiesOnly = false;

	// the subsytem results composed together
	protected Strategy strategy = null; // the product strategy
	protected SMG currentModelExpl; // the composed model

	/**
	 * Constructor.
	 */
	public CompositionalSMGModelChecker(PrismComponent parent, ModulesFile modulesFile, PropertiesFile propertiesFile, SimulatorEngine engine)
			throws PrismException
	{
		super(parent);
		this.modulesFile = modulesFile;
		this.propertiesFile = propertiesFile;
		this.engine = engine;

		// Get combined constant values from model/properties
		constantValues = new Values();
		constantValues.addValues(modulesFile.getConstantValues());
		if (propertiesFile != null)
			constantValues.addValues(propertiesFile.getConstantValues());
		this.labelListModel = modulesFile.getLabelList();
		this.labelListProp = propertiesFile.getLabelList();

		// initialise the model checker
		mc = new SMGModelChecker(parent);
	}

	public void setCheckCompatibility(boolean b)
	{
	        checkCompatibility = b;
	}

	public void setComputeParetoSet(boolean computePareto)
	{
		this.computePareto = computePareto;
	}

	public void setCancel(boolean[] cancel_computation)
	{
		this.cancel_computation = cancel_computation;
	}


        public Strategy getStrategy()
        {
	        return strategy;
        }

	/**
	 * Model check a property.
	 */
	public Result check(Expression expr) throws PrismException
	{
		Result res;
		String resultString;
		long timer;

		// Starting model checking
		timer = System.currentTimeMillis();

		// Do model checking
		res = checkExpression(expr);

		// Model checking complete
		timer = System.currentTimeMillis() - timer;
		mainLog.println("\nModel checking completed in " + (timer / 1000.0) + " secs.");

		// Print result to log
		resultString = "Result";
		if (!("Result".equals(expr.getResultName())))
			resultString += " (" + expr.getResultName().toLowerCase() + ")";
		resultString += ": " + res;
		mainLog.print("\n" + resultString + "\n");

		// Return result
		return res;
	}

	private Result checkExpression(Expression expr) throws PrismException
	{
		Result res;

		// Just support the comp() operator
		if (expr instanceof ExpressionFunc && ((ExpressionFunc) expr).getNameCode() == ExpressionFunc.COMP) {
		        res = checkCompositionalExpression((ExpressionFunc) expr);
		} else
			throw new PrismException("Compositional model checking just works with the comp() function");

		return res;
	}

	private Result checkCompositionalExpression(ExpressionFunc expr) throws PrismException
	{
		mainLog.print("//////////////////////////////////////////////////////////////////////////////\n");
	        mainLog.print("//                   STARTING COMPOSITIONAL MODEL CHECKING                  //\n");
	        mainLog.print("//////////////////////////////////////////////////////////////////////////////\n\n");

		Result result = new Result(true);

		// set subsystem model checker options
		mc.setComputeParetoSet(computePareto);
		mc.setCancel(cancel_computation);

		// build the subsystems
		mainLog.print("Building Model ... \n");
		buildSubsystems(false);
		int numSubsystems = subsystems.size();

		if(expr.getNumOperands() < numSubsystems)
		        throw new PrismException("Compositional property does not have a specification for each subsystem.");
		else if(expr.getNumOperands() == numSubsystems + 1)
		        // rebuild full model to be used for the last specification
		        buildSubsystems(true);
		else if (expr.getNumOperands() > numSubsystems + 1)
		        mainLog.printWarning("Compositional property has too many arguments, the superfluous ones are ignored.");
		

		// check each model separately	    
		subsystemStrategies = new ArrayList<Strategy>(numSubsystems); // a strategy for each subsystem
		subsystemParetos = new ArrayList<Pareto>(numSubsystems); // a list of polyhedra for each subsystem
		subsystemParams = new ArrayList<MultiParameters>(numSubsystems); // a parameter for each subsystem
		List<Result> subsystemResults = new ArrayList<Result>(numSubsystems); // one result per subsystem

		for (int i = 0; i < numSubsystems; i++) {
			mc.setModulesFileAndPropertiesFile(subsystemModulesFiles.get(i), propertiesFile);
			Result result_i = mc.check(subsystems.get(i), expr.getOperand(i));
			subsystemResults.add(result_i);
			subsystemStrategies.add(result_i.getStrategy());
			subsystemParetos.add(mc.pareto_set);
			subsystemParams.add(mc.parsed_params);
		}

		mainLog.print("//////////////////////////////////////////////////////////////////////////////\n");
		mainLog.print("//                   COMPLETED COMPOSITIONAL MODEL CHECKING                 //\n");
		mainLog.print("//////////////////////////////////////////////////////////////////////////////\n\n");

		if (computePareto) { // Pareto set computation
		        showCompositionalParetoSet(expr, currentModelExpl.getName());
		        return result;
		} else if (settings.getBoolean(PrismSettings.PRISM_GENERATE_STRATEGY)) { // construct product strategy
			strategy = new StochasticUpdateStrategyProduct(subsystemStrategies);
			//strategy.setComposition(currentModelExpl);
			strategy.setInfo("Property: " + expr.toString() + "\n" + "Type: " + strategy.getType() + "\nMemory size: " + strategy.getMemorySize());
			// put stratetegy in result
			result.setStrategy(strategy);
		}

		for(Result r : subsystemResults) {
		    if(r.getResult() instanceof Boolean) {
			result.setResult(((Boolean) result.getResult()) && ((Boolean)r.getResult()));
		    } else
			throw new PrismException("Unexpected result type");
		}

		return result;		
	}

	/**
	 * Builds the subsystems one by one.
	 **/
	private void buildSubsystems(boolean buildFullModel) throws PrismException
	{
	    ConstructModel constructModel = new ConstructModel(this, engine);
	    constructModel.setCheckCompatibility(checkCompatibility);
	    subsystems = new ArrayList<SMG>();
	    subsystemModulesFiles = new ArrayList<ModulesFile>();
   	    currentModelExpl = (SMG) constructModel.constructSMGModelCompositionally(modulesFile, false, true, true, subsystems, subsystemModulesFiles, buildFullModel, cancel_computation);

	    currentModelExpl.findDeadlocks(false); // do not fix deadlocks in composition
	    checkForDeadlocksExpl(currentModelExpl); // check for deadlocks - if found, abort
	}


 	/**
	* check for deadlocks in an explicit model - same as in prism.Prism without fixing deadlocks
	**/
	private void checkForDeadlocksExpl(Model modelExpl) throws PrismException
	{
		StateValues deadlocks = modelExpl.getDeadlockStatesList();
		int numDeadlocks = modelExpl.getNumDeadlockStates();
		if (numDeadlocks > 0) {
			if (!(modelExpl == null))
				mainLog.print(modelExpl.infoStringTable());
			mainLog.print("\n" + numDeadlocks + " deadlock states found");
			if (!getVerbose() && numDeadlocks > 10) {
				mainLog.print(". The first 10 are below. Use verbose mode to view them all.\n");
				deadlocks.print(mainLog, 10);
			} else {
				mainLog.print(":\n");
				deadlocks.print(mainLog);
				mainLog.print(String.format("deadlocked model: %s\n", modelExpl));
			}
			mainLog.print("\nTip: Use the \"fix deadlocks\" option to automatically add self-loops in deadlock states.\n");

			throw new PrismException("Model contains " + numDeadlocks + " deadlock state" + (numDeadlocks > 1 ? "s" : ""));
		}
	}


	

	/**
	* Prepares the Pareto sets for displaying by mapping them to the appropriate dimensions
	* and intersecting them with one another. The folding is then done within the PointList class.
	**/
	private void showCompositionalParetoSet(ExpressionFunc expr, String name) throws PrismException
	{
		int numSubsystems = subsystems.size();

		// all expressions and bounds in the compositional expression
		List<Expression> expressionList = new ArrayList<Expression>();
		List<Double> boundsList = new ArrayList<Double>();
		for(MultiParameters params : subsystemParams) {
		    expressionList.addAll(params.expressions);
		    boundsList.addAll(params.bounds);
		}
		int dim = expressionList.size(); // full dimension over all subsystems

		// build list of lifted Pareto sets
		List<Pareto> paretos = new ArrayList<Pareto>();

		// for each subsystem, lift the Pareto set into its own space, will fold later
		int base_dim = 0; // dimension up to last subsystem
		for (int i = 0; i < numSubsystems; i++) {
			Pareto sPs = new Pareto(subsystemParetos.get(i)); // deep copy
			paretos.add(sPs);
			int sdim = (int) sPs.getDimension();

			// build partial function to map dimensions
			Partial_Function pfunc = new Partial_Function();
			for (int j = 0; j < sdim; j++)
				pfunc.insert(j, base_dim + j);
			for (int j = 0; j < base_dim; j++)
				// fill up to get consistent mapping
				pfunc.insert(sdim + j, j);

			for (Polyhedron p : sPs.getSets()) {
				// fill in full dimensions
				p.add_space_dimensions_and_embed(dim - p.space_dimension());
				// map polyhedron to the right dimensions
				p.map_space_dimensions(pfunc);
				// fill in full dimensions if removed
				p.add_space_dimensions_and_embed(dim - p.space_dimension());
			}
			base_dim += sdim; // dimension up to this subsystem
		}

		// CROSS-INTERSECT
		Pareto new_pareto = null;
		int i = 0;
		for (Pareto p : paretos) {
			if (p.size() == 0)
				continue;
			if (new_pareto == null) {
				new_pareto = p;
				continue; // first non-empty Pareto set found
			} else {
				Pareto temp_pareto = new Pareto();
				// intersect new_pareto and p
				middle: for (Polyhedron p1 : new_pareto.getSets()) {
					if (p1.is_empty())
						continue middle;
					for (Polyhedron p2 : p.getSets()) {
						Polyhedron q = new C_Polyhedron((C_Polyhedron) p1); // deep copy
						i++;
						q.intersection_assign(p2);
						temp_pareto.add(q); // add p1 intersect p2 to temp_pareto
					}
				}
				new_pareto = temp_pareto;
			}
		}

		// now build a point list
		PointList pointList = new PointList(new_pareto, expressionList, boundsList);
		// register point list to be displayed
		pointList.addStoredPointList(name, pointList);

	}

	public boolean getVerbose()
	{
		return settings.getBoolean(PrismSettings.PRISM_VERBOSE);
	}
}
