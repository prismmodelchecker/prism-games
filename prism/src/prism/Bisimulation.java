package prism;

import java.io.*;
import java.util.Arrays;
import java.util.HashMap;

import parser.ast.*;

public class Bisimulation {

	public static void main(String[] args) {
		new Bisimulation().go(args);
	}
	
	public void go(String args[]) {
		try {
			//Initialize
			PrismLog mainLog = new PrismFileLog("stdout");
			Prism prism = new Prism(mainLog);
			prism.initialise();
			prism.setEngine(Prism.EXPLICIT);
			
			//Parse model
			ModulesFile modulesFile = prism.parseModelFile(new File(args[0]));
			prism.loadPRISMModel(modulesFile);
			
			//Get model
			prism.buildModel();
			explicit.Model model = prism.getBuiltModelExplicit();
			explicit.MDPSparse mdp = (explicit.MDPSparse) model;
			
			//Compute action matrices
			int numStates = mdp.getNumStates();
			HashMap<String, double[][]> actionMatrices = new HashMap();
			for(int i = 0; i < numStates; i++) {
				
			}
	
			
			System.out.println(mdp.getAction(0, 6));
		}
		catch(FileNotFoundException e) {
			System.out.println("Error: " + e);
			System.exit(1);
		}
		catch(PrismException e) {
			System.out.println("Error: " + e);
			System.exit(1);
		}
	}
	
	
}
