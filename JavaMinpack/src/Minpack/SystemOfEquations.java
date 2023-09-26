package Minpack;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.function.Function;

import Minpack.systemInterface.func;
import Minpack.systemInterface.jacob;


public class SystemOfEquations{
	
	//Base format for the System of Equations for Minpack to use.
	
	private func functionMethod;
	private jacob jacobianMethod;
	int jWidth;
	int jHeight;
	int fSize;
	private boolean useFuncMethod = false;
	private boolean useJacMethod = false;
	
	private ArrayList<Function<double[],Double>> functions = new ArrayList<Function<double[],Double>>();
	private ArrayList<ArrayList<Function<double[],Double>>> Jacobian = new ArrayList<ArrayList<Function<double[],Double>>>(); //Just added, not fully tested
	
	public SystemOfEquations() {
		
	}
		
	public double[] evaluate(double[] x0) {
		double[] x;
		if(useFuncMethod) {
			x = functionMethod.operation(x0);
		}else {
			int fSize = functions.size();
			x = new double[fSize];
			for(int i = 0;i<fSize;i++) {
				x[i] = functions.get(i).apply(x0);
			}
		}
		return x;
	}
	public double[][] evaluateJacobian(double[] x0) { //Just added, not fully tested
		double[][] x;
		if(useJacMethod) {
			x = jacobianMethod.operation(x0);
		}else {
			int jLength = Jacobian.size();
			int jWidth =  Jacobian.get(0).size();
			x = new double[jWidth][jLength];
			for(int i = 0;i<jLength;i++) {
				for(int j = 0;i<jWidth;i++) {
					x[i][j] = Jacobian.get(i).get(j).apply(x0);
				}
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
	
	public void useFunctionMethod(func funct, int n) {
		useFuncMethod = true;
		fSize = n;
		functionMethod = funct;
	}
	public void useJacobianMethod(jacob funct, int n, int m) {
		jWidth = m;
		jHeight = n;
		useJacMethod = true;
		jacobianMethod = funct;
	}
	public void removeFunctionMethod(func funct) {
		useFuncMethod = false;
	}
	public void removeJacobianMethod(func funct) {
		useJacMethod = false;
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
		if(useFuncMethod) {
			return fSize;
		}else {
			return functions.size();
		}
	}
	
}
