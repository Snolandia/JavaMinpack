package Minpack;

import java.util.ArrayList;
import java.util.function.Function;

public class SystemOfEquations{
	
	//Base format for the System of Equations for Minpack to use.
	
	private ArrayList<Function<double[],Double>> functions = new ArrayList<Function<double[],Double>>();
	private ArrayList<ArrayList<Function<double[],Double>>> Jacobian = new ArrayList<ArrayList<Function<double[],Double>>>(); //Just added, not fully tested
	
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
	public double[][] evaluateJacobian(double[] x0) { //Just added, not fully tested
		
		int jLength = Jacobian.size();
		int jWidth =  Jacobian.get(0).size();
		double x[][] = new double[jWidth][jLength];
		for(int i = 0;i<jLength;i++) {
			for(int j = 0;i<jWidth;i++) {
				x[i][j] = Jacobian.get(i).get(j).apply(x0);
			}
		}
		
		return x;
	}

	public void setJacobian(ArrayList<ArrayList<Function<double[],Double>>> jac){
		Jacobian = jac;

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
