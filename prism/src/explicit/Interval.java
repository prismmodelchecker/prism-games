package explicit;

import java.text.DecimalFormat;

public class Interval {
	 private static DecimalFormat df = new DecimalFormat("0.000");
		
	
	public double lhs=0;
	public double rhs=1;
	
	public Interval(double l, double r)
	{
		lhs=l;rhs=r;
	}
	
	/**
	 * Checks whether the interval int1 is contained in this interval
	 * @param int1
	 * @return
	 */
	public boolean achievableFrom(Interval int1)
	{
		return lhs <= int1.lhs && rhs >= int1.rhs;
	}
	
	/**
	 * Checks whether the interval contained in this interval
	 * can be achieved by a convex combination of int1 and int2
	 * @param int1
	 * @param int2
	 * @return
	 */
	public boolean achievableFrom(Interval int1, Interval int2)
	{
		double a = (lhs - int2.lhs) / (int1.lhs-int2.lhs);
		if(a>=0 && a<=1)
		{
			if(rhs >= a*int1.rhs + (1-a)*int2.rhs)
				return true;
		}
		
		
		double b = (rhs - int2.rhs) / (int1.rhs-int2.rhs);
		if(b>=0 && b<=1)
		{
			if(lhs <= b*int1.lhs + (1-b)*int2.lhs)
				return true;
		}
		
		return false;
	}
	
	/**
	 * HACK
	 */
	public boolean convexAchievableFrom(Interval int1, Interval int2)
	{
		
		
		double a = (lhs - int2.lhs) / (int1.lhs-int2.lhs);
		if(a>=0 && a<=1)
		{
			if(rhs >= a*int1.rhs + (1-a)*int2.rhs)
				return true;
		}
		
		
		double b = (rhs - int2.rhs) / (int1.rhs-int2.rhs);
		if(b>=0 && b<=1)
		{
			if(lhs <= b*int1.lhs + (1-b)*int2.lhs)
				return true;
		}
		
		return false;
	}
	
	public static Interval getUnion(Interval int1, Interval int2)
	{
		return new Interval(Math.min(int1.lhs, int2.lhs), Math.max(int1.rhs, int2.rhs));
	}
	
	@Override
	public String toString()
	{
		return "(" + df.format(lhs) + ","+df.format(rhs) + ")";
	}
	
	@Override
	public boolean equals(Object o)
	{
		if(o instanceof Interval && ((Interval)o).lhs==lhs && ((Interval)o).rhs==rhs)
			return true;
		else
			return false;
	}
	
	@Override
	public int hashCode()
	{
		return (int) (1000000*((lhs+rhs)+lhs*rhs));
	}
	
	public static void main(String[] args)
	{
		Interval int1, int2, int3, i;
		
		i = new Interval(0.25,0.5);
		int1 = new Interval(0,0);
		int2 = new Interval(0.4,1);
		int3 = new Interval(0.3,0.4);
		
		System.out.println(i.achievableFrom(int3));
		System.out.println(i.achievableFrom(int3, int1));
	}
	
}
