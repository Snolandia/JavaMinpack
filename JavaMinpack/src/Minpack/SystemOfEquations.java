package Minpack;

import java.util.ArrayList;
import java.util.function.Function;

public class SystemOfEquations{
	
	//Base format for the System of Equations for Minpack to use.
	
	private ArrayList<Function<double[],Double>> functions = new ArrayList<Function<double[],Double>>();
	
	public SystemOfEquations() {
		
	}
	
	public double[] evaluate(double[] x0) {
		
		int fSize = functions.size();
		double x[] = new double[fSize];
		for(int i = 0;i<fSize;i++) {
			x[i] = functions.get(i).apply(x0);
		}
		
		return x;
	}
	
	public void addFunction(Function<double[],Double> func) {
		functions.add(func);
	}
	
	public void removeFunction(int r) {
		functions.remove(r);
	}
	
	public void clearFunctions() {
		functions.clear();
	}
	
	public Function<double[],Double> getFunction(int i){
		return functions.get(i);
	}
	public int size() {
		return functions.size();
	}
	
}
