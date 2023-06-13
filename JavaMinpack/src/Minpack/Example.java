package Minpack;

import java.util.function.Function;

public class Example {
	
	public static void main(String args[]) {
		example();
	}
	
	private static void example() {
		// Create a SystemOfEquations to represent the equations needing to be solved
		SystemOfEquations fcn = new SystemOfEquations();
		// Load the system up with Function<double[],Double> functions
		// Use a double[] array to represent the variables in the equation, such that it
		// starts at 0
		// 2 different methods, you can add the function directly in like this
		fcn.addFunction(x -> (Math.pow((x[0]), 2)) + (Math.pow((x[1]), 2)) + (Math.pow((x[2]), 2)) - 100); // A sphere
//											 with radius of 10, at point 0,0,0 in 3d space
		fcn.addFunction(x -> (Math.pow((x[0] - 10), 2)) + (Math.pow((x[1]), 2)) + (Math.pow((x[2]), 2)) - 100); // A sphere
//											with radius of 10, at point 10,0,0 in 3d space
		// or you can declare a function then add it, like this
		Function<double[],Double> fx = x -> (Math.pow((x[0]-10), 2)) + (Math.pow((x[1]-10), 2)) + (Math.pow((x[2]-8), 2)) - 100; // A sphere 
//											with radius of 10, at point 10,10,8 in 3d space
		fcn.addFunction(fx);
		//Now that the SystemOfEquations has been created, we need to declare an initial guess of the solution. Attempt to guess as close as possible, 
//		as a far away guess may lead to an unexpected or wanted solution. 
		//The initial guess will need to be in a double[] array of length n, n = the number of variables.
		double[] x0 = {25,25,25};
		//Next we need to set xTol, the desired tolerance of the solution, as a double.
		double xTol = 1.0E-8;
		//Now, from Minpack we will call hybrd1 to solve this for us. Note, hybrd1 required that n, the number 
			//of variables, be equal to m, the number of equations. Hybrd will return the solution as an double[] array of length n.
		double[] solution;
		solution = Minpack.hybrd1(fcn, x0, xTol);
		//Now printing the solution results in : 5.0, 1.357433228762054, 8.553208464047431
		for(int i = 0;i<solution.length;i++) {
			System.out.println(solution[i]);
		}
		//If we where to change our initial guess;
		x0 =new double[] {-25,0,-8};
		solution = Minpack.hybrd1(fcn, x0, xTol);
		//Now printing the solution results in : 5.0, 8.642566771237943, -0.5532084640474304
		for(int i = 0;i<solution.length;i++) {
			System.out.println(solution[i]);
		}
		//fin
		//More features to come in future builds
		
	}
	
	
	
}
