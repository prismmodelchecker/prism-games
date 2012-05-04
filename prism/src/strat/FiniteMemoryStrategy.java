package strat;

import explicit.Distribution;
import explicit.Model;

public abstract class FiniteMemoryStrategy implements Strategy
{
	private String info = null;

	@Override
	abstract public void init(int state) throws InvalidStrategyStateException;

	@Override
	abstract public void updateMemory(int action, int state) throws InvalidStrategyStateException;

	@Override
	abstract public Distribution getNextMove(int state) throws InvalidStrategyStateException;

	@Override
	abstract public void reset();

	@Override
	public void exportToFile(String file)
	{
		// TODO Auto-generated method stub
	}

	@Override
	public Model buildProduct(Model model)
	{
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String getInfo()
	{
		return info;
		}

	@Override
	public void setInfo(String info)
	{
		this.info=info;
	}

	@Override
	abstract public int getMemorySize();

	@Override
	public String getType()
	{
		return "Finite memory strategy";
	}

	@Override
	abstract public Object getCurrentMemoryElement();
	
	@Override
	abstract public void setMemory(Object memory) throws InvalidStrategyStateException;

	@Override
	abstract public String getStateDescription();
}
