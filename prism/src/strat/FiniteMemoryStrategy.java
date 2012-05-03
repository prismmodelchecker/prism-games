package strat;

import explicit.Distribution;
import explicit.Model;

public class FiniteMemoryStrategy implements Strategy {

	@Override
	public void init(int state) throws InvalidStrategyStateException {
		// TODO Auto-generated method stub

	}

	@Override
	public void updateMemory(int action, int state)
			throws InvalidStrategyStateException {
		// TODO Auto-generated method stub

	}

	@Override
	public Distribution getNextMove(int state)
			throws InvalidStrategyStateException {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void reset() {
		// TODO Auto-generated method stub

	}

	@Override
	public void exportToFile(String file) {
		// TODO Auto-generated method stub

	}

	@Override
	public Model buildProduct(Model model) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String getInfo() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void setInfo(String info) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public int getMemorySize() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public String getType() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Object getCurrentMemoryElement() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void setMemory(Object memory) throws InvalidStrategyStateException {
		// TODO Auto-generated method stub
		
	}

}
