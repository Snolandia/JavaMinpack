package Minpack;

public class Minpack {

	//Converted from Fortran, then modified slightly. Work in progress to finish the rest 
//	of the functions as currently, only hybrd1/hybrd is fully functional
	
	//Disclaimer from original Minpack:
	
//	Minpack Copyright Notice (1999) University of Chicago.  All rights reserved
//
//	Redistribution and use in source and binary forms, with or
//	without modification, are permitted provided that the
//	following conditions are met:
//
//	1. Redistributions of source code must retain the above
//	copyright notice, this list of conditions and the following
//	disclaimer.
//
//	2. Redistributions in binary form must reproduce the above
//	copyright notice, this list of conditions and the following
//	disclaimer in the documentation and/or other materials
//	provided with the distribution.
//
//	3. The end-user documentation included with the
//	redistribution, if any, must include the following
//	acknowledgment:
//
//	   "This product includes software developed by the
//	   University of Chicago, as Operator of Argonne National
//	   Laboratory.
//
//	Alternately, this acknowledgment may appear in the software
//	itself, if and wherever such third-party acknowledgments
//	normally appear.
//
//	4. WARRANTY DISCLAIMER. THE SOFTWARE IS SUPPLIED "AS IS"
//	WITHOUT WARRANTY OF ANY KIND. THE COPYRIGHT HOLDER, THE
//	UNITED STATES, THE UNITED STATES DEPARTMENT OF ENERGY, AND
//	THEIR EMPLOYEES: (1) DISCLAIM ANY WARRANTIES, EXPRESS OR
//	IMPLIED, INCLUDING BUT NOT LIMITED TO ANY IMPLIED WARRANTIES
//	OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE
//	OR NON-INFRINGEMENT, (2) DO NOT ASSUME ANY LEGAL LIABILITY
//	OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR
//	USEFULNESS OF THE SOFTWARE, (3) DO NOT REPRESENT THAT USE OF
//	THE SOFTWARE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS, (4)
//	DO NOT WARRANT THAT THE SOFTWARE WILL FUNCTION
//	UNINTERRUPTED, THAT IT IS ERROR-FREE OR THAT ANY ERRORS WILL
//	BE CORRECTED.
//
//	5. LIMITATION OF LIABILITY. IN NO EVENT WILL THE COPYRIGHT
//	HOLDER, THE UNITED STATES, THE UNITED STATES DEPARTMENT OF
//	ENERGY, OR THEIR EMPLOYEES: BE LIABLE FOR ANY INDIRECT,
//	INCIDENTAL, CONSEQUENTIAL, SPECIAL OR PUNITIVE DAMAGES OF
//	ANY KIND OR NATURE, INCLUDING BUT NOT LIMITED TO LOSS OF
//	PROFITS OR LOSS OF DATA, FOR ANY REASON WHATSOEVER, WHETHER
//	SUCH LIABILITY IS ASSERTED ON THE BASIS OF CONTRACT, TORT
//	(INCLUDING NEGLIGENCE OR STRICT LIABILITY), OR OTHERWISE,
//	EVEN IF ANY OF SAID PARTIES HAS BEEN WARNED OF THE
//	POSSIBILITY OF SUCH LOSS OR DAMAGES.


//module minpack_module

	static double n;
	//static double epsmch = Math.ulp(n); //removing as too much precision causes issues
	static double epsmch = 0.0000000000000002220446049250313; //Using this instead due to above statement
	static String[] infoDict = {
			"Info 0 : Input parameters are incorrect.",
			"Info 1 : Relative error between two consecutive iterates is at most `xtol`.",
			"Info 2 : Iterations has exceded maxfev, the max allowed iterations",
			"Info 3 : xTol is too small. No further improvement in the approximate solution `x` is possible.",
			"Info 4 : Iteration is not making good progress, measured by the improvement from the last five jacobian evaluations.",
			"Info 5 : Tteration is not making good progress, as measured by the improvement from the last ten iterations."};

	/**
	 * <p>
	 * NOT TESTED YET
	 * This subroutine checks the gradients of m nonlinear functions in n variables,
	 * evaluated at a point x, for consistency with the functions themselves.
	 * <p>
	 * The subroutine does not perform reliably if cancellation or rounding errors
	 * cause a severe loss of significance in the evaluation of a function.
	 * Therefore, none of the components of x should be unusually small (in
	 * particular, zero) or any other value which may cause loss of significance.
	 * <p>
	 * 
	 * @param m      int -> A positive integer variable set to the number of
	 *               functions.(Length of functions array)
	 *               <p>
	 * @param n      int -> A positive integer variable set to the number of
	 *               variables.(Length of array x)
	 *               <p>
	 * @param ldfjac int -> A positive integer, not less than m, which specifies the
	 *               leading dimension of the array fjac. (i.e., The number of
	 *               functions/equations.ldfjac can just be set to m)
	 *               <p>
	 * @param mode   int -> An integer input set to 1 on the first call and 2 on the
	 *               second. The user must call chkder twice,first with mode = 1 and
	 *               then with mode = 2. mode = 1 -> x must contain the point of
	 *               evaluation. xp is set to a neighboring point. mode = 2 -> fvec
	 *               must contain the functions and the rows of fjac must contain
	 *               the gradients of the respective functions each evaluated at x,
	 *               and fvecp must contain the functions evaluated at xp. Err
	 *               contains measures of correctness of the respective gradients.
	 *               <p>
	 * @param x      double[n] -> Input array (i.e. The variables)
	 *               <p>
	 * @param fvec   double[m] -> An array of length m. When mode = 2, fvec must
	 *               contain the functions evaluated at x.
	 *               <p>
	 * @param fjac   double[ldfjac][n] -> an m by n array. When mode = 2, the rows
	 *               of fjac must contain the gradients of the respective functions
	 *               evaluated at x.
	 *               <p>
	 * @param xp     double[n] -> An array of length n. When mode = 1, xp is set to
	 *               a neighboring point of x.
	 *               <p>
	 * @param fvecp  double[m] -> An array of length m. When mode = 2, fvecp must
	 *               contain the functions evaluated at xp.
	 *               <p>
	 * @param err    double[m] -> An array of length m. When mode = 2, err contains
	 *               measures of correctness of the respective gradients. If there
	 *               is no severe loss of significance, then if err(i) is 1.0 the
	 *               i-th gradient is correct, while if err(i) is 0.0 the i-th
	 *               gradient is incorrect. For values of err between 0.0 and 1.0,
	 *               the categorization is less certain. In general, a value of
	 *               err(i) greater than 0.5 indicates that the i-th gradient is
	 *               probably correct, while a value of err(i) less than 0.5
	 *               indicates that the i-th gradient is probably incorrect.
	 * 
	 * @return Nothing is returned, rather the objects themselves are modified.
	 */
	public static void chkder(int m, int n, double[] x, double[] fvec, double[][] fjac, int ldfjac, double[] xp,
			double[] fvecp, int mode, double[] err) {
		
		//NOT TESTED YET
		
		int i, j; // Only used for loops, not important.
		double temp; // Just used as a temp variable
		double eps = Math.sqrt(epsmch); // sqrt of machine epsilon
		double factor = 100.0; // Factor, or percent of change that happens.
		double epsf = factor * epsmch; // Factor of epsmch, or percent of change that will happen.
		double epslog = Math.log10(epsf); // log 10 of epsf

		if (mode == 2) {
			err = new double[m];
			for (j = 0; j <= n; j++) {
				temp = Math.abs(x[j]);
				if (temp == 0) {
					temp = 1.0;
				}
				for (i = 0; i <= m; i++) {
					err[i] = err[i] + temp * fjac[i][j];
				}
			}
			for (i = 0; i <= m; i++) {
				temp = 1.0;
				if (fvec[i] != 0 && fvecp[i] != 0 && Math.abs(fvecp[i] - fvec[i]) >= epsf * Math.abs(fvec[i])) {
					temp = eps * Math.abs((fvecp[i] - fvec[i]) / eps - err[i])
							/ (Math.abs(fvec[i]) + Math.abs(fvecp[i]));
				}
				err[i] = 1.0;
				if (temp > epsmch && temp < eps) {
					err[i] = (Math.log10(temp) - epslog) / epslog;
				}
				if (temp >= eps) {
					err[i] = 0.0;
				}
			}
		} else if (mode == 1) {
			for (j = 0; j <= n; j++) {
				temp = eps * Math.abs(x[j]);
				if (temp == 0) {
					temp = eps;
				}
				xp[j] = x[j] + temp;
			}
		} else {
			throw new IllegalArgumentException("Invalid mode inside the Chkder method");
		}
	}

	/**
	 * <p>
	 * Given an m by n matrix (a), an n by n non-singular diagonal matrix (d), an
	 * m-vector (b), and a positive number (delta), the problem is to determine the
	 * convex combination x of the gauss-newton and scaled gradient directions that
	 * minimizes (a*x - b) in the least squares sense, subject to the restriction
	 * that the euclidean norm of d*x be at most delta.
	 * <p>
	 * This subroutine completes the solution of the problem if it is provided with
	 * the necessary information from the qr factorization of (a). That is, if a =
	 * q*r, where q has orthogonal columns and r is an upper triangular matrix, then
	 * dogleg expects the full upper triangle of r and the first n components of (q
	 * transpose)*b.
	 * 
	 * @param n     int -> A positive integer variable set to the order of r. This
	 *              should be equivalent to the length of array x.
	 *              <p>
	 * @param lr    int -> A positive integer variable not less than (n*(n+1))/2.
	 *              Can just be set as lr=(n*(n+1))/2
	 *              <p>
	 * @param delta int -> A positive integer variable which specifies an upper
	 *              bound on the euclidean norm of d*x.
	 *              <p>
	 * @param r     double[lr] -> An array of length lr which must contain the upper
	 *              triangular matrix r stored by rows.
	 *              <p>
	 * @param diag  double[n] -> An array of length n which must contain the
	 *              diagonal elements of the matrix d.
	 *              <p>
	 * @param qtb   double[n] -> An array of length n which must contain the first n
	 *              elements of the vector (q transpose)*b.
	 *              <p>
	 * @param x     double[n] -> An array of length n which contains the desired
	 *              convex combination of the gauss-newton direction and the scaled
	 *              gradient direction.
	 *              <p>
	 * @param wa1   double[n] -> A working array of length n.
	 *              <p>
	 * @param wa2   double[n] -> A working array of length n.
	 *              <p>
	 * @return Nothing is returned, rather the objects themselves are modified.
	 * 
	 */
	public static void dogleg(int n, double[] r, int lr, double[] diag, double[] qtb, double delta, double[] x,
			double[] wa1, double[] wa2) {

		int i, j, jj, jp1, k, l; // temporary ints
		double alpha, bnorm, gnorm, qnorm, sgnorm, sum, temp; // temporary doubles
//	    First, calculate the gauss-newton direction.
		jj = (n * (n + 1)) / 2;
		for (k = 0; k < n; k++) {
			j = n - (k+1);
			jp1 = j + 1;
			jj = jj - (k+1);
			l = jj + 1;
			sum = 0;
			if (n >= jp1) {
				for (i = jp1; i < n; i++) {
					sum = sum + r[l] * x[i];
					l = l + 1;
				}
			}
			temp = r[jj];
			if (temp == 0) {
				l = j;
				for (i = 0; i < j; i++) {
					if (temp < Math.abs(r[l])) {
						temp = Math.abs(r[l]);
					}
					l = l + n - 1 - i;
				}
				temp = epsmch * temp;
				if (temp == 0) {
					temp = epsmch;
				}
			}
			x[j] = (qtb[j] - sum) / temp;
		}
//		Test whether the gauss-newton direction is acceptable.
		for (j = 0; j < n; j++) {
			wa1[j] = 0;
			wa2[j] = diag[j] * x[j];
		}
		qnorm = enorm(n, wa2);
		if (qnorm > delta) {
//		The gauss-newton direction is not acceptable.
//		Next, calculate the scaled gradient direction.
			l = 0;
			for (j = 0; j < n; j++) {
				temp = qtb[j];
				for (i = j; i < n; i++) {
					wa1[i] = wa1[i] + r[l] * temp;
					l = l + 1;
				}
				wa1[j] = wa1[j] / diag[j];
			}
//			Calculate the norm of the scaled gradient and test for
//			the special case in which the scaled gradient is zero.
			gnorm = enorm(n, wa1);
			sgnorm = 0;
			alpha = delta / qnorm;
			if (gnorm != 0) {
//			Calculate the point along the scaled gradient
//			at which the quadratic is minimized.
				for (j = 0; j < n; j++) {
					wa1[j] = (wa1[j] / gnorm) / diag[j];
				}
				l = 0;
				for (j = 0; j < n; j++) {
					sum = 0;
					for (i = j; i < n; i++) {
						sum = sum + r[l] * wa1[i];
						l = l + 1;
					}
					wa2[j] = sum;
				}
				temp = enorm(n, wa2);
				sgnorm = (gnorm / temp) / temp;
//				Test whether the scaled gradient direction is acceptable.
				alpha = 0;
				if (sgnorm < delta) {
//				The scaled gradient direction is not acceptable.
//				Finally, calculate the point along the dogleg
//				at which the quadratic is minimized.
					bnorm = enorm(n, qtb);
					temp = (bnorm / gnorm) * (bnorm / qnorm) * (sgnorm / delta);
					temp = temp - (delta / qnorm) * ((sgnorm / delta) * (sgnorm / delta))
							+ Math.sqrt((temp - (delta / qnorm)) * (temp - (delta / qnorm)))
							+ (1 - ((delta / qnorm) * (delta / qnorm))) * (1 - ((sgnorm / delta) * (sgnorm / delta)));
					alpha = ((delta / qnorm) * (1 - ((sgnorm / delta) * (sgnorm / delta)))) / temp;
				}
			}
//			Form appropriate convex combination of the gauss-newton
//			direction and the scaled gradient direction.
			if (sgnorm < delta) {
				temp = (1 - alpha) * sgnorm;
			} else {
				temp = (1 - alpha) * delta;
			}
			for (j = 0; j < n; j++) {
				x[j] = temp * wa1[j] + alpha * x[j];
			}
		}
	}

	/**
	 * <p>
	 * Given an n-vector x, this function calculates the euclidean norm of x.
	 * 
	 * <p>
	 * The euclidean norm is computed by accumulating the sum of ! squares in three
	 * different sums. The sums of squares for the small and large components are
	 * scaled so that no overflows occur. Non-destructive underflows are permitted.
	 * Underflows and overflows do not occur in the computation of the un-scaled sum
	 * of squares for the intermediate components. The definitions of small,
	 * intermediate and large components depend on two constants, rdwarf and rgiant.
	 * The main restrictions on these constants are that rdwarf**2 prevents
	 * underflow and rgiant**2 prevents overflow.
	 * 
	 * @param n int -> A positive integer variable.
	 *          <p>
	 * @param x double[n] -> An array of length n.
	 *          <p>
	 * @return enorm double -> Containing the euclidean norm of x.
	 */

	public static double enorm(int n, double[] x) {

		int i; // Used for the loops
		double agiant, s1, s2, s3, xabs, x1max, x3max; // Some temporary holders
		double enorm = 0; // The returned variable.
		double rdwarf = 0.00000000000000000003834; // 3.834e-20 Unsure if this is still needed, or if this is
//													the best number for this process. Leaving as is.
		double rgiant = 13040000000000000000.0; // 1.304e19 Unsure if this is still needed, or if this is
//													the best number for this process. Leaving as is.
		s1 = 0;
		s2 = 0;
		s3 = 0;
		x1max = 0;
		x3max = 0;
		agiant = rgiant / Double.valueOf(n);
		for (i = 0; i < n; i++) {
			xabs = Math.abs(x[i]);
			if (xabs > rdwarf && xabs < agiant) {
//				Sum for intermediate components.
				s2 = s2 + (xabs * xabs);
			} else if (xabs <= rdwarf) {
//				Sum for small components.
				if (xabs <= x3max) {
					if (xabs != 0) {
						s3 = s3 + ((xabs / x3max) * (xabs / x3max));
					}
				} else {
					s3 = 1 + s3 * ((x3max / xabs) * (x3max / xabs));
					x3max = xabs;
				}
//				Sum for large components.
			} else if (xabs <= x1max) {
				s1 = s1 + ((xabs / x1max) * (xabs / x1max));
			} else {
				s1 = 1 + s1 * ((x1max / xabs) * (x1max / xabs));
				x1max = xabs;
			}
		}
//				Calculation of norm.
		if (s1 != 0) {
			enorm = x1max * Math.sqrt(s1 + (s2 / x1max) / x1max);
		} else if (s2 == 0) {
			enorm = x3max * Math.sqrt(s3);
		} else {
			if (s2 >= x3max) {
				enorm = Math.sqrt(s2 * (1 + (x3max / s2) * (x3max * s3)));
			}
			if (s2 < x3max) {
				enorm = Math.sqrt(x3max * ((s2 / x3max) + (x3max * s3)));
			}
		}
//			
		return enorm;
	}
	
	/**
	 * 
	 * <p>
	 * This subroutine computes a forward-difference approximation to the n by n
	 * jacobian matrix associated with a specified problem of n functions in n
	 * variables. If the jacobian has a banded form, then function evaluations are
	 * saved by only approximating the nonzero terms.
	 * 
	 * 
	 * procedure(func) :: fcn -> the user-supplied subroutine which calculates the
	 * functions.
	 * 
	 * @param n      int -> A positive integer variable set to the number of
	 *               functions and variables.
	 *               <p>
	 * @param ldfjac int -> A positive integer variable not less than n which
	 *               specifies the leading dimension of the array fjac. I.e. Number
	 *               of functions.
	 *               <p>
	 * @param iFlag  int -> An integer variable which can be used to terminate the
	 *               execution of fdjac1.
	 *               <p>
	 * @param ml     int -> A nonnegative integer variable which specifies the
	 *               number of sub-diagonals within the band of the jacobian matrix.
	 *               If the jacobian is not banded, set ml to at least n - 1.
	 *               <p>
	 * @param mu     int -> A nonnegative integer variable which specifies the
	 *               number of super-diagonals within the band of the jacobian
	 *               matrix. If the jacobian is not banded, set mu to at least n -
	 *               1.
	 *               <p>
	 * @param epsfcn double -> A variable used in determining a suitable step length
	 *               for the forward-difference approximation. This approximation
	 *               assumes that the relative errors in the functions are of the
	 *               order of epsfcn. If epsfcn is less than the machine precision,
	 *               it is assumed that the relative errors in the functions are of
	 *               the order of the machine precision.
	 *               <p>
	 * @param x      double[n] -> An array of length n. The variables.
	 *               <p>
	 * @param fvec   double[n] -> An array of length n which must contain the
	 *               functions evaluated at x.
	 *               <p>
	 * @param fjac   double[ldfjac][n] -> A n by n array which contains the
	 *               approximation to the jacobian matrix evaluated at x.
	 *               <p>
	 * @param wa1    double[n] -> A working array of length n.
	 *               <p>
	 * @param wa2    double[n] -> A working array of length n. If ml + mu + 1 is at
	 *               least n, then the jacobian is considered dense, and wa2 is not
	 *               referenced.
	 * @return Nothing is returned, rather the objects themselves are modified.
	 */       

	public static void fdjac1(SystemOfEquations fcn, int n, double[] x, double[] fvec, double[][] fjac, int ldfjac,
			int iflag, int ml, int mu, double epsfcn, double[] wa1, double[] wa2) {
		
		int i, j, k, msum; //Temporary variables
		double eps, h, temp; //Temporary variables

		//double epsmch = 0.0000000000000002220446049250313; // Removing this and using the global causes fdjac1 to return too small of a gradient change
		
		if (epsfcn > epsmch) {
			eps = Math.sqrt(epsfcn);
		} else {
			eps = Math.sqrt(epsmch);
		}
		msum = ml + mu + 1;
		if (msum < n) {
			// computation of banded approximate jacobian.
			for (k = 0; k < msum+1; k++) {
				for (j = k; j < n; j = +msum) {
					wa2[j] = x[j];
					h = eps * Math.abs(wa2[j]);
					if (h == 0) {
						h = eps;
					}
					x[j] = wa2[j ] + h;
				}
				wa1 = fcn.evaluate(x); // This part needs to be changed later
				for (j = k; k < n; k = +msum) {
					x[j] = wa2[j];
					h = eps * Math.abs(wa2[j]);
					if (h == 0) {
						h = eps;
					}
					for (i = 0; i < n; i++) {
						fjac[i][j] = 0;
						if (i >= j - mu && i <= j + ml) {
							fjac[i][j] = (wa1[i] - fvec[i]) / h;
						}
					}
				}
			}
		} else {
//	        Computation of dense approximate jacobian.
			for (j = 0; j < n; j++) {
				temp = x[j];
				h = eps * Math.abs(temp);
				if (h == 0) {
					h = eps;
				}
				x[j] = temp + h;
				wa1 = fcn.evaluate(x); // This part needs to be changed later
				x[j] = temp;
				for (i = 0; i < n; i++) {
					fjac[i][j] = (wa1[i] - fvec[i]) / h;
				}
			}
		}
	}

/**
 * 
 * <p>
 * NOT TESTED YET
 * This subroutine computes a forward-difference approximation to the m by n
 * jacobian matrix associated with a specified problem of m functions in n
 * variables.
 * 
 * procedure(func2) :: fcn !! the user-supplied subroutine which calculates the
 * functions.
 * <p>
 * 
 * @param m      int -> A positive integer variable set to the number of
 *               functions.
 *               <p>
 * @param n      int -> A positive integer variable set to the number of
 *               variables. n must not exceed m.
 *               <p>
 * @param ldfjac int -> A positive integer variable not less than m which
 *               specifies the leading dimension of the array fjac.
 *               <p>
 * @param iFlag  int -> An integer variable which can be used to terminate the
 *               execution of fdjac2.
 *               <p>
 * @param epsfcn double -> An variable used in determining a suitable step
 *               length for the forward-difference approximation. This
 *               approximation assumes that the relative errors in the functions
 *               are of the order of epsfcn. If epsfcn is less than the machine
 *               precision, it is assumed that the relative errors in the
 *               functions are of the order of the machine precision.
 *               <p>
 * @param x      double[n] -> An array of length n.
 *               <p>
 * @param fvec   double[m] -> An array of length m which must contain the
 *               functions evaluated at x.
 *               <p>
 * @param fjac   double[ldfjac, n] -> A m by n array which contains the
 *               approximation to the jacobian matrix evaluated at x.
 *               <p>
 * @param wa     double[m] -> A working array of length m.
 *               <p>
 * @return Nothing is returned, rather the objects themselves are modified.
 */
	public static void fdjac2(SystemOfEquations fcn,int m,int n,double[] x,double[] fvec,double[][] fjac,int Ldfjac,int iFlag,double epsfcn,double[] wa) {
		
		//NOT TESTED YET
		
		int i,j; //used for loops
		double eps, h, temp; //temporary variables
		
		if(epsfcn>epsmch) {
			eps = Math.sqrt(epsfcn);
		}else {
			eps = Math.sqrt(epsmch);
		}
		for(j = 0;j<n;j++) {
	            temp = x[j];
	            h = eps*Math.abs(temp);
	            if (h == 0) {
	            	h = eps;
	            }
	            x[j] = temp + h;
	            fcn.evaluate(x);
	            if (iFlag < 0) {
	            	throw new IllegalArgumentException("Class method fdjac2, at call to fcn, has returned an iFlag less than 0");
	            }
	            x[j] = temp;
	            for(i = 0;i<m;i++) {
	                fjac[i][j] = (wa[i] - fvec[i])/h;
	            }
		}
	}
		

/**
 * <p>The purpose of hybrd is to find a zero of a system of n nonlinear functions
 * in n variables by a modification of the powell hybrid method. the user must
 * provide a subroutine which calculates the functions. The jacobian is then
 * calculated by a forward-difference approximation.
 * 
 * procedure(func) :: fcn !! user-supplied subroutine which calculates the
 * functions
 * 
 * @param n      int -> A positive integer variable set to the number of
 *               functions and variables.
 *               <p>
 * @param maxFev int -> A positive integer variable. Termination occurs
 *               when the number of calls to `fcn` is at least `maxFev` by the
 *               end of an iteration.
 *               <p>
 * @param ml     int -> A nonnegative integer variable which specifies the
 *               number of sub-diagonals within the band of the jacobian matrix.
 *               If the jacobian is not banded, set `ml` to at least `n - 1`.
 *               <p>
 * @param mu     int -> A nonnegative integer variable which specifies the
 *               number of super-diagonals within the band of the jacobian
 *               matrix. If the jacobian is not banded, set `mu` to at least` n
 *               - 1`.
 *               <p>
 * @param mode   int -> If `mode = 1`, the variables will be scaled internally.
 *               If `mode = 2`, the scaling is specified by the input `diag`.
 *               Other values of `mode` are equivalent to `mode = 1`.
 *               <p>
 * @param nPrint int -> An integer variable that enables controlled
 *               printing of iterates if it is positive. In this case, `fcn` is
 *               called with `iFlag = 0` at the beginning of the first iteration
 *               and every `nPrint` iterations thereafter and immediately prior
 *               to return, with `x` and `fvec` available for printing. If
 *               `nPrint` is not positive, no special calls of `fcn` with `iFlag
 *               = 0` are made.
 *               <p>
 * @param info   int -> An integer variable. If the user has terminated
 *               execution, `info` is set to the (negative) value of `iflag`.
 *               see description of `fcn`. Otherwise, `info` is set as follows:
 *               <p>
 *               * ***info = 0*** improper input parameters.
 *               <p>
 *               * ***info = 1*** relative error between two consecutive
 *               iterates is at most `xTol`.
 *               <p>
 *               * ***info = 2*** number of calls to `fcn` has reached or
 *               exceeded `maxFev`.
 *               <p>
 *               * ***info = 3*** `xTol` is too small. no further improvement in
 *               the approximate solution `x` is possible.
 *               <p>
 *               * ***info = 4*** iteration is not making good progress, as
 *               measured by the improvement from the last five jacobian
 *               evaluations.
 *               <p>
 *               * ***info = 5*** iteration is not making good progress, as
 *               measured by the improvement from the last ten iterations.
 *               <p>
 * @param nFev   int -> Variable set to the number of calls to `fcn`.
 *               <p>
 * @param ldfjac int -> A positive integer variable not less than `n`
 *               which specifies the leading dimension of the array `fjac`.
 *               <p>
 * @param lr     int -> A positive integer input variable not less than
 *               `(n*(n+1))/2`.
 *               <p>
 * @param xTol   double -> A nonnegative variable. Termination occurs when
 *               the relative error between two consecutive iterates is at most
 *               `xTol`.
 *               <p>
 * @param epsfcn double -> A variable used in determining a suitable step
 *               length for the forward-difference approximation. This
 *               approximation assumes that the relative errors in the functions
 *               are of the order of `epsfcn`. If `epsfcn` is less than the
 *               machine precision, it is assumed that the relative errors in
 *               the functions are of the order of the machine precision.
 *               <p>
 * @param factor double -> A positive variable used in determining the
 *               initial step bound. This bound is set to the product of
 *               `factor` and the euclidean norm of `diag*x` if nonzero, or else
 *               to `factor` itself. In most cases factor should lie in the
 *               interval (.1,100.). -> (100.0) is a generally recommended value.
 *               <p>
 * @param x      double[n] -> Array of length n. On input `x` must contain an
 *               initial estimate of the solution vector. On output `x` contains
 *               the final estimate of the solution vector.
 *               <p>
 * @param fvec   double[n] -> An output array of length `n` which contains the
 *               functions evaluated at the output `x`.
 *               <p>
 * @param diag   double[n] -> An array of length `n`. If `mode = 1` (see below),
 *               `diag` is internally set. If `mode = 2`, `diag` must contain
 *               positive entries that serve as multiplicative scale factors for
 *               the variables.
 *               <p>
 * @param fjac   double[ldfjac][n] -> Array which contains the orthogonal matrix
 *               `q` produced by the QR factorization of the final approximate
 *               jacobian.
 *               <p>
 * @param r      double[lr] -> An array which contains the upper
 *               triangular matrix produced by the QR factorization of the final
 *               approximate jacobian, stored row-wise.
 *               <p>
 * @param qtf    double[n] -> An output array of length `n` which contains the
 *               vector `(q transpose)*fvec`.
 *               <p>
 * @param wa1    double[n] -> A working array.
 *               <p>
 * @param wa2    double[n] -> A working array.
 *               <p>
 * @param wa3    double[n] -> A working array.
 *               <p>
 * @param wa4    double[n] -> A working array.
 *               <p>
 * 
 * @return Returns a double[] that contains the solved values, x
 * 
 */
	
public static double[] hybrd(SystemOfEquations fcn, int n, double[] x, double[] fvec, double xTol, int maxFev, int ml,
			int mu, double epsfcn, double[] diag, int mode, double factor, int nPrint, int info, int nfev,
			double[][] fjac, int ldfjac, double[] r, int lr, double[] qtf, double[] wa1, double[] wa2, double[] wa3,
			double[] wa4) {

		int i, iFlag, iter, j, jm1, l, msum, ncFail, ncSuc, nSlow1, nSlow2; //Temporary Variables
		int[] iwa = new int[0]; //Temporary Variables
		boolean jeval, sing; //Temporary Variables
		double actred, delta = 0, fNorm, fNorm1, pNorm, prered, ratio, sum, temp, xNorm = 0; //Temporary Variables
		double[][] qtfHold = new double[1][n]; //Used to create a sub-array to pass

		double p1 = 0.1;
		double p5 = 0.5;
		double p001 = 0.001;
		double p0001 = 0.0001;

		info = 0;
		iFlag = 0;
		nfev = 0;
		boolean main = true;
		boolean inner = true;
		boolean outer = true;

//			
//	        main : block
		do { // main loop
//	        check the input parameters for errors.
			if (n <= 0 || xTol < 0 || maxFev <= 0 || ml < 0 || mu < 0 || factor <= 0 || ldfjac < n
					|| lr < (n * (n + 1)) / 2) {
				throw new IllegalArgumentException("Input parameters do not match as expected");
			}
			if (mode == 2) {
				for (int a = 0; a < n; a++) {
					if (diag[a] <= 0) {
						throw new IllegalArgumentException("Diagonal Error at Hybrd");
					}
				}
			}
//	        Evaluate the function at the starting point
//	        and calculate its norm.
			iFlag = 1;
			fvec = fcn.evaluate(x);
			nfev = 1;
			if (iFlag < 0) {
				throw new IllegalArgumentException("iFlag was returned as negative");
			}
			fNorm = enorm(n, fvec);
//	        Determine the number of calls to fcn needed to compute
//	        the jacobian matrix.
			if (ml + mu + 1 < n) {
				msum = ml + mu + 1;
			} else {
				msum = n + 1;
			}
//	        Initialize iteration counter and monitors.
			iter = 1;
			ncSuc = 0;
			ncFail = 0;
			nSlow1 = 0;
			nSlow2 = 0;

//	        Beginning of the outer loop.
			do { // outer loop
				jeval = true;
//	            Calculate the jacobian matrix.
				iFlag = 2;
				fdjac1(fcn, n, x, fvec, fjac, ldfjac, iFlag, ml, mu, epsfcn, wa1, wa2);
				nfev = nfev + msum;
//	            Compute the qr factorization of the jacobian.
				qrfac(n, n, fjac, ldfjac, false, iwa, 1, wa1, wa2, wa3);
//	            On the first iteration and if mode is 1, scale according
//	            to the norms of the columns of the initial jacobian.
				if (iter == 1) {
					if (mode != 2) {
						for (j = 0; j < n; j++) {
							diag[j] = wa2[j];
							if (wa2[j] == 0) {
								diag[j] = 1;
							}
						}
					}
//					On the first iteration, calculate the norm of the scaled x
//					and initialize the step bound delta.
					for (j = 0; j < n; j++) {
						wa3[j] = diag[j] * x[j];
						
					}
					xNorm = enorm(n, wa3);
					delta = factor * xNorm;
					if (delta == 0) {
						delta = factor;
					}
				}
//	            Form (q transpose)*fvec and store in qtf.
				for (i = 0; i < n; i++) {
					qtf[i] = fvec[i];
				}
				for (j = 0; j < n; j++) {
					if (fjac[j][j] != 0) {
						sum = 0;
						for (i = j; i < n; i++) {
							sum = sum + fjac[i][j] * qtf[i];
						}
						temp = -sum / fjac[j][j];
						for (i = j; i < n; i++) {
							qtf[i] = qtf[i] + fjac[i][j] * temp;
						}
					}
				}
//				Copy the triangular factor of the qr factorization into r.
				sing = false;
				for (j = 0; j < n; j++) {
					l = j;
					jm1 = j;
					if (jm1 >= 0) {
						for (i = 0; i < jm1; i++) {
							r[l] = fjac[i][j];
							l = l + n - 1 - i;
						}
					}
					r[l] = wa1[j];
					if (wa1[j] == 0) {
						sing = true;
					}
				}
//				Accumulate the orthogonal factor in fjac.
				qform(n, n, fjac, ldfjac, wa1);
//				Re-scale if necessary.
				if (mode != 2) {
					for (j = 0; j < n; j++) {
						if (diag[j] > wa2[j]) {
							diag[j] = diag[j];
						} else {
							diag[j] = wa2[j];
						}
					}
				}
				do { // Beginning of the inner loop.
//					If requested, call fcn to enable printing of iterates.
					if (nPrint > 0) {
						iFlag = 0;
						if ((iter - 1 % nPrint) == 0) {
							// add a method to do the printing
						}
					}
//					Determine the direction p.
					dogleg(n, r, lr, diag, qtf, delta, wa1, wa2, wa3);
//					Store the direction p and x + p. calculate the norm of p.
					for (j = 0; j < n; j++) {
						wa1[j] = -wa1[j];
						wa2[j] = x[j] + wa1[j];
						wa3[j] = diag[j] * wa1[j];
					}
					pNorm = enorm(n, wa3);
//					On the first iteration, adjust the initial step bound.
					if (iter == 1) {
						if (delta > pNorm) {
							delta = pNorm;
						}
					}
//					Evaluate the function at x + p and calculate its norm.
					iFlag = 1;
					wa4 = fcn.evaluate(wa2);

//					System.out.println("test " + wa2[j]);
//					for (i = 0; i < n; i++) {
//						for (j = 0; j < n; j++) {
//							System.out.println("test " + wa2[j]);
//							System.out.println("test " + wa4[j]);
//						}
//					}
//					System.exit(0);
					nfev = nfev + 1;
					if (iFlag < 0) {
						throw new IllegalArgumentException("iFlag returned from fcn as a negative");
					}
					fNorm1 = enorm(n, wa4);
//					Compute the scaled actual reduction.
					actred = -1;
					if (fNorm1 < fNorm) {
						actred = 1 - ((fNorm1 / fNorm) * (fNorm1 / fNorm));
					}
//					Compute the scaled predicted reduction.
					l = 0;
					for (i = 0; i < n; i++) {
						sum = 0;
						for (j = i; j < n; j++) {
							sum = sum + r[l] * wa1[j];
							l = l + 1;
						}
						wa3[i] = qtf[i] + sum;
					}
					temp = enorm(n, wa3);

					prered = 0;
					if (temp < fNorm) {
						prered = 1 - ((temp / fNorm) * (temp / fNorm));
					}
//					Compute the ratio of the actual to the predicted
//					reduction.
					ratio = 0;
					if (prered > 0) {
						ratio = actred / prered;
					}
//					Update the step bound.
					if (ratio >= p1) {
						ncFail = 0;
						ncSuc = ncSuc + 1;
						if (ratio >= p5 || ncSuc > 1) {
							if (delta < pNorm / p5) {
								delta = pNorm / p5;
							}
						}
						if (Math.abs(ratio - 1) <= p1) {
							delta = pNorm / p5;
						}
					} else {
						ncSuc = 0;
						ncFail = ncFail + 1;
						delta = p5 * delta;
					}
//					Test for successful iteration.
					if (ratio >= p0001) {
//						Successful iteration. update x, fvec, and their norms.
						for (j = 0; j < n; j++) {
							x[j] = wa2[j];
							wa2[j] = diag[j] * x[j];
							fvec[j] = wa4[j];
						}
						xNorm = enorm(n, wa2);
						fNorm = fNorm1;
						iter = iter + 1;
					}
//					Determine the progress of the iteration.
					nSlow1 = nSlow1 + 1;
					if (actred >= p001) {
						nSlow1 = 0;
					}
					if (jeval) {
						nSlow2 = nSlow2 + 1;
					}
					if (actred >= p1) {
						nSlow2 = 0;
					}
//					Test for convergence.
					if (delta <= xTol * xNorm || fNorm == 0) {
						info = 1;
						inner = false;
						outer = false;
						main = false;
						break;
					}
					if (info != 0) {
						throw new IllegalArgumentException(infoDict[info]);
					}
//					Tests for termination and stringent tolerances.
					if (nfev >= maxFev) {
						info = 2;
					}
					if (p1 * pNorm <= epsmch * xNorm || p1 * (p1 * delta) <= epsmch * xNorm) {
						info = 3;
					}
					if (nSlow2 == 5) {
						info = 4;
					}
					if (nSlow1 == 10) {
						info = 5;
					}
					if (info != 0) {
						for(int zx = 0;zx<n;zx++) {
							//System.out.println(x[zx]);
						}
						return x;
						//throw new IllegalArgumentException(infoDict[info]);
					}
//					Criterion for recalculating jacobian approximation
//					by forward differences.
					if (ncFail == 2) {
						break; // Breaks to outer cycle to recalculate jacobian
					}
//					Calculate the rank one modification to the jacobian
//					and update qtf if necessary.
					for (j = 0; j < n; j++) {
						sum = 0;
						for (i = 0; i < n; i++) {
							sum = sum + fjac[i][j] * wa4[i];
						}
						wa2[j] = (sum - wa3[j]) / pNorm;
						wa1[j] = diag[j] * ((diag[j] * wa1[j]) / pNorm);
						if (ratio >= p0001) {
							qtf[j] = sum;
						}
					}
//					Compute the qr factorization of the updated jacobian.
					r1updt(n, n, r, lr, wa1, wa2, wa3, sing);
					r1mpyq(n, n, fjac, ldfjac, wa2, wa3); 
					for (int z = 0; z < n; z++) {
						qtfHold[0][z] = qtf[z]; //Used to pass a smaller array to r1mpyq
					}
					qtfHold = r1mpyq(1, n, qtfHold, 1, wa2, wa3);
					for (int z = 0; z < n; z++) {
						qtf[z] = qtfHold[0][z]; //Used to pass r1mpyq results back to qtf
					}
					jeval = false; 
				} while (inner); //end of the inner loop.
			} while (outer); //end of the outer loop.
		} while (main); //end of main loop		
		return x;
	}

	public static double[] hybrd1(SystemOfEquations fcn, double[] x, double tol) {
		// !*****************************************************************************************
		// !>

		// ! the purpose of hybrd1 is to find a zero of a system of
		// ! n nonlinear functions in n variables by a modification
		// ! of the powell hybrid method. this is done by using the
		// ! more general nonlinear equation solver hybrd. the user
		// ! must provide a subroutine which calculates the functions.
		// ! the jacobian is then calculated by a forward-difference
		// ! approximation.
		//
//	    subroutine hybrd1(fcn, n, x, Fvec, Tol, Info, Wa, Lwa)
		//
//	        implicit none
		//
//	        procedure(func)                     :: fcn      !! user-supplied subroutine which calculates the functions
//	        integer, intent(in)                  :: n        !! a positive integer input variable set to the number
//	                                                    !! of functions and variables.
//	        integer, intent(out)                 :: info     !! an integer output variable. if the user has
//	                                                    !! terminated execution, info is set to the (negative)
//	                                                    !! value of `iflag`. see description of `fcn`. otherwise,
//	                                                    !! `info` is set as follows:
//	                                                    !!
//	                                                    !!  * ***info = 0*** improper input parameters.
//	                                                    !!  * ***info = 1*** algorithm estimates that the relative error
//	                                                    !!  between `x` and the solution is at most `tol`.
//	                                                    !!  * ***info = 2*** number of calls to `fcn` has reached or exceeded
//	                                                    !!  `200*(n+1)`.
//	                                                    !!  * ***info = 3*** `tol` is too small. no further improvement in
//	                                                    !!  the approximate solution `x` is possible.
//	                                                    !!  * ***info = 4*** iteration is not making good progress.
//	        real(wp), intent(in)                 :: tol      !! a nonnegative input variable. termination occurs
//	                                                    !! when the algorithm estimates that the relative error
//	                                                    !! between `x` and the solution is at most `tol`.
//	        real(wp), dimension(n), intent(inout) :: x        !! an array of length `n`. on input `x` must contain
//	                                                    !! an initial estimate of the solution vector. on output `x`
//	                                                    !! contains the final estimate of the solution vector.
//	        real(wp), dimension(n), intent(out)   :: fvec     !! an output array of length `n` which contains
//	                                                    !! the functions evaluated at the output `x`.
//	        integer, intent(in) :: Lwa !! a positive integer input variable not less than
//	                              !! (n*(3*n+13))/2.
//	        real(wp), intent(inout) :: Wa(Lwa) !! a work array of length lwa.
		//
//	        integer :: index, j, lr, maxfev, ml, mode, mu, nfev, nprint
//	        real(wp) :: epsfcn, xtol

		//
//	        reaL(wp), parameter :: factor = 100.0_wp

		//
//	        Info = 0
		//
//	        ! check the input parameters for errors.
//				
		if (fcn.size() != x.length) {
			System.out.println("not square, exiting now");
			System.exit(0);
		}
		if (tol <= 0) {
			tol = 1.490116119384766E-008;// real(wp), intent(in):: tol !! a nonnegative input variable. termination
									// occurs
//	          !! when the algorithm estimates that the relative error
//	          !! between `x` and the solution is at most `tol`.
		}
		int info = 0;
//		        integer, intent(out) :: info : an integer output variable. if the user has
//	            !! terminated execution, `info` is set to the (negative)
//	            !! value of `iflag`. see description of `fcn`. otherwise,
//	            !! `info` is set as follows:
//	            !!
//	            !!  * ***info = 0*** improper input parameters.
//	            !!  * ***info = 1*** relative error between two consecutive iterates
//	            !!    is at most `xtol`.
//	            !!  * ***info = 2*** number of calls to `fcn` has reached or exceeded
//	            !!    `maxfev`.
//	            !!  * ***info = 3*** `xtol` is too small. no further improvement in
//	            !!    the approximate solution `x` is possible.
//	            !!  * ***info = 4*** iteration is not making good progress, as
//	            !!    measured by the improvement from the last
//	            !!    five jacobian evaluations.
//	            !!  * ***info = 5*** iteration is not making good progress, as
//	            !!    measured by the improvement from the last
//	            !!    ten iterations.
		int n = x.length;
//		        integer, intent(in) :: n : a positive integer input variable set to the number
//	            !! of functions and variables.
		int lwa = (n * (3 * n + 13)) / 2;// integer, intent(in) :: Lwa !! a positive integer input variable not less
											// than (n*(3*n+13))/2.
		if (n > 0 & tol >= 0 & lwa >= (n * (3 * n + 13)) / 2) {
			// ! call hybrd.
			int maxfev = 200 * (n + 1);
//		        integer, intent(in) :: maxfev : a positive integer input variable. termination
//	          !! occurs when the number of calls to `fcn` is at least `maxfev`
//	          !! by the end of an iteration.
			double xtol = tol;
//		        real(wp), intent(in) :: xtol : a nonnegative input variable. termination
//	          !! occurs when the relative error between two consecutive
//	          !! iterates is at most `xtol`.
			int ml = n;
//		        integer, intent(in) :: ml : a nonnegative integer input variable which specifies
//	          !! the number of subdiagonals within the band of the
//	          !! jacobian matrix. if the jacobian is not banded, set
//	          !! `ml` to at least `n - 1`.
			int mu = n;
//		        integer, intent(in) :: mu: a nonnegative integer input variable which specifies
//	          !! the number of superdiagonals within the band of the
//	          !! jacobian matrix. if the jacobian is not banded, set
//	          !! `mu` to at least` n - 1`.
			double epsfcn = 0;
//	          real(wp), intent(in) :: epsfcn           !! an input variable used in determining a suitable
//	          !! step length for the forward-difference approximation. this
//	          !! approximation assumes that the relative errors in the
//	          !! functions are of the order of `epsfcn`. if `epsfcn` is less
//	          !! than the machine precision, it is assumed that the relative
//	          !! errors in the functions are of the order of the machine
//	          !! precision.
			int mode = 2;
			// mode = 1;
//		        integer, intent(in) :: mode : if `mode = 1`, the
//	          !! variables will be scaled internally. if `mode = 2`,
//	          !! the scaling is specified by the input `diag`. other
//	          !! values of `mode` are equivalent to `mode = 1`.

			int nprint = 0;
//		        integer, intent(in)  :: nprint : an integer input variable that enables controlled
//	          !! printing of iterates if it is positive. in this case,
//	          !! `fcn` is called with `iflag = 0` at the beginning of the first
//	          !! iteration and every `nprint` iterations thereafter and
//	          !! immediately prior to return, with `x` and `fvec` available
//	          !! for printing. if `nprint` is not positive, no special calls
//	          !! of `fcn` with `iflag = 0` are made.
			int lr = (n * (n + 1)) / 2;
//		        integer, intent(in) :: lr : a positive integer input variable not less than `(n*(n+1))/2`.
			int nfev = 0;
//		        integer, intent(out) :: nfev :output variable set to the number of calls to `fcn`.
			double[] fvec = new double[n];
//		        real(wp), intent(out) :: fvec(n) : an output array of length `n` which contains
//		                                            !! the functions evaluated at the output `x`.
			double[] diag = new double[n];
//		        real(wp), intent(inout) :: diag(n)       !! an array of length `n`. if `mode = 1` (see
//	          !! below), `diag` is internally set. if `mode = 2`, `diag`
//	          !! must contain positive entries that serve as
//	          !! multiplicative scale factors for the variables.
			for (int j = 0; j < n; j++) { // this is for if mode is set to 2
				diag[j] = 1;
			}

			double factor = 100.0;
//	          	real(wp), intent(in) :: factor : a positive input variable used in determining the
//	          !! initial step bound. this bound is set to the product of
//	          !! `factor` and the euclidean norm of `diag*x` if nonzero, or else
//	          !! to `factor` itself. in most cases factor should lie in the
//	          !! interval (.1,100.). 100. is a generally recommended value.
			double[][] fjac = new double[n][n];// array which contains the
//	          !! orthogonal matrix `q` produced by the QR factorization
//	          !! of the final approximate jacobian.
			int ldfjac = n;
//		        integer, intent(in):: ldfjac : a positive integer input variable not less than `n`
//	          !! which specifies the leading dimension of the array `fjac`.
			double[] r = new double[lr]; // an output array which contains the
//	          !! upper triangular matrix produced by the QR factorization
//	          !! of the final approximate jacobian, stored rowwise.
			double[] qtf = new double[n]; // an output array of length `n` which contains the vector `(q
											// transpose)*fvec`.

			double[] wa1 = new double[n];// real(wp), intent(inout) :: wa1(n) !! work array
			double[] wa2 = new double[n];// real(wp), intent(inout) :: wa2(n) !! work array
			double[] wa3 = new double[n];// real(wp), intent(inout) :: wa3(n) !! work array
			double[] wa4 = new double[n];// real(wp), intent(inout) :: wa4(n) !! work array

			x = hybrd(fcn, n, x, fvec, xtol, maxfev, ml, mu, epsfcn, diag, mode, factor, nprint, info, nfev, fjac,
					ldfjac, r, lr, qtf, wa1, wa2, wa3, wa4);

//	          hybrd(double[]fcn,int n,double[] x, double[] fvec, double Xtol,int maxfev,int ml,int mu,double epsfcn,double[] diag,int mode,
//	               double factor,int nprint,int info,int nfev,double[][] fjac,int ldfjac,double[] r,int lr
//	          	,double[] qtf,double[] wa1, double[]wa2,double[] wa3,double[] wa4)

			if (info == 5) {
				info = 4;
			}
		}

//	        if (n > 0 .and. Tol >= zero .and. Lwa >= (n*(3*n + 13))/2) then
//	            ! call hybrd.
//	            maxfev = 200*(n + 1)
//	            xtol = Tol
//	            ml = n - 1
//	            mu = n - 1
//	            epsfcn = zero
//	            mode = 2
//	            do j = 1, n
//	                Wa(j) = one
//	            end do
//	            nprint = 0
//	            lr = (n*(n + 1))/2
//	            index = 6*n + lr
//	            call hybrd(fcn, n, x, Fvec, xtol, maxfev, ml, mu, epsfcn, Wa(1), mode, &
//	                       factor, nprint, Info, nfev, Wa(index + 1), n, Wa(6*n + 1), lr, &
//	                       Wa(n + 1), Wa(2*n + 1), Wa(3*n + 1), Wa(4*n + 1), Wa(5*n + 1))
//	            if (Info == 5) Info = 4
//	        end if
		//
//	    end subroutine hybrd1
		// !*****************************************************************************************
		//
		return x;
	}

	// !*****************************************************************************************
	// !>
	// ! the purpose of hybrj is to find a zero of a system of
	// ! n nonlinear functions in n variables by a modification
	// ! of the powell hybrid method. the user must provide a
	// ! subroutine which calculates the functions and the jacobian.
	//
//	    subroutine hybrj(fcn, n, x, Fvec, Fjac, Ldfjac, Xtol, Maxfev, Diag, Mode, &
//	                     Factor, Nprint, Info, Nfev, Njev, r, Lr, Qtf, Wa1, Wa2, &
//	                     Wa3, Wa4)
	//
//	        implicit none
	//
//	        procedure(fcn_hybrj) :: fcn !! the user-supplied subroutine which
//	                                !! calculates the functions and the jacobian
//	        integer, intent(in) :: n !! a positive integer input variable set to the number
//	                            !! of functions and variables.
//	        integer, intent(in) :: Ldfjac !! a positive integer input variable not less than n
//	                                    !! which specifies the leading dimension of the array fjac.
//	        integer, intent(in) :: Maxfev !! a positive integer input variable. termination
//	                                    !! occurs when the number of calls to fcn with iflag = 1
//	                                    !! has reached maxfev.
//	        integer, intent(in) :: Mode !! an integer input variable. if mode = 1, the
//	                                !! variables will be scaled internally. if mode = 2,
//	                                !! the scaling is specified by the input diag. other
//	                                !! values of mode are equivalent to mode = 1.
//	        integer, intent(in) :: Nprint !! an integer input variable that enables controlled
//	                                    !! printing of iterates if it is positive. in this case,
//	                                    !! fcn is called with iflag = 0 at the beginning of the first
//	                                    !! iteration and every nprint iterations thereafter and
//	                                    !! immediately prior to return, with x and fvec available
//	                                    !! for printing. fvec and fjac should not be altered.
//	                                    !! if nprint is not positive, no special calls of fcn
//	                                    !! with iflag = 0 are made.
//	        integer, intent(out) :: Info !! an integer output variable. if the user has
//	                                !! terminated execution, info is set to the (negative)
//	                                !! value of iflag. see description of fcn. otherwise,
//	                                !! info is set as follows:
//	                                !!
//	                                !!  * ***info = 0***   improper input parameters.
//	                                !!  * ***info = 1***   relative error between two consecutive iterates
//	                                !!    is at most xtol.
//	                                !!  * ***info = 2***   number of calls to fcn with iflag = 1 has
//	                                !!    reached maxfev.
//	                                !!  * ***info = 3***   xtol is too small. no further improvement in
//	                                !!    the approximate solution x is possible.
//	                                !!  * ***info = 4***   iteration is not making good progress, as
//	                                !!    measured by the improvement from the last
//	                                !!    five jacobian evaluations.
//	                                !!  * ***info = 5***   iteration is not making good progress, as
//	                                !!    measured by the improvement from the last
//	                                !!    ten iterations.
//	        integer, intent(out) :: Nfev !! an integer output variable set to the number of
//	                                !! calls to fcn with iflag = 1.
//	        integer, intent(out) :: Njev !! an integer output variable set to the number of
//	                                !! calls to fcn with iflag = 2.
//	        integer, intent(in) :: Lr !! a positive integer input variable not less than
//	                                !! (n*(n+1))/2.
//	        real(wp), intent(in) :: Xtol !! a nonnegative input variable. termination
//	                                !! occurs when the relative error between two consecutive
//	                                !! iterates is at most xtol.
//	        real(wp), intent(in) :: Factor !! a positive input variable used in determining the
//	                                    !! initial step bound. this bound is set to the product of
//	                                    !! factor and the euclidean norm of diag*x if nonzero, or else
//	                                    !! to factor itself. in most cases factor should lie in the
//	                                    !! interval (.1,100.). 100. is a generally recommended value.
//	        real(wp), intent(inout) :: x(n) !! an array of length n. on input x must contain
//	                                    !! an initial estimate of the solution vector. on output x
//	                                    !! contains the final estimate of the solution vector.
//	        real(wp), intent(out) :: Fvec(n) !! an output array of length n which contains
//	                                    !! the functions evaluated at the output x.
//	        real(wp), intent(out) :: Fjac(Ldfjac, n) !! an output n by n array which contains the
//	                                            !! orthogonal matrix q produced by the qr factorization
//	                                            !! of the final approximate jacobian.
//	        real(wp), intent(inout) :: Diag(n) !! an array of length n. if mode = 1 (see
//	                                        !! below), diag is internally set. if mode = 2, diag
//	                                        !! must contain positive entries that serve as
//	                                        !! multiplicative scale factors for the variables.
//	        real(wp), intent(out) :: r(Lr) !! an output array of length lr which contains the
//	                                    !! upper triangular matrix produced by the qr factorization
//	                                    !! of the final approximate jacobian, stored rowwise.
//	        real(wp), intent(out) :: Qtf(n) !! an output array of length n which contains
//	                                    !! the vector (q transpose)*fvec.
//	        real(wp), intent(inout) :: Wa1(n) !! work array of length n.
//	        real(wp), intent(inout) :: Wa2(n) !! work array of length n.
//	        real(wp), intent(inout) :: Wa3(n) !! work array of length n.
//	        real(wp), intent(inout) :: Wa4(n) !! work array of length n.
	//
//	        integer :: i, iflag, iter, j, jm1, l, ncfail, ncsuc, nslow1, nslow2
//	        integer :: iwa(1)
//	        logical :: jeval, sing
//	        real(wp) :: actred, delta, fnorm, fnorm1, pnorm, prered, ratio, sum, temp, xnorm
	//
//	        real(wp), parameter :: p1 = 1.0e-1_wp
//	        real(wp), parameter :: p5 = 5.0e-1_wp
//	        real(wp), parameter :: p001 = 1.0e-3_wp
//	        real(wp), parameter :: p0001 = 1.0e-4_wp
		public static void hybrj() {
		}
//	        Info = 0
//	        iflag = 0
//	        Nfev = 0
//	        Njev = 0
	//
//	        main : block
	//
//	            ! check the input parameters for errors.
	//
//	            if (n <= 0 .or. Ldfjac < n .or. Xtol < zero .or. Maxfev <= 0 .or. &
//	                Factor <= zero .or. Lr < (n*(n + 1))/2) exit main
//	            if (Mode == 2) then
//	                do j = 1, n
//	                    if (Diag(j) <= zero) exit main
//	                end do
//	            end if
	//
//	            ! evaluate the function at the starting point
//	            ! and calculate its norm.
	//
//	            iflag = 1
//	            call fcn(n, x, Fvec, Fjac, Ldfjac, iflag)
//	            Nfev = 1
//	            if (iflag < 0) exit main
//	            fnorm = enorm(n, Fvec)
	//
//	            ! initialize iteration counter and monitors.
	//
//	            iter = 1
//	            ncsuc = 0
//	            ncfail = 0
//	            nslow1 = 0
//	            nslow2 = 0
	//
//	            ! beginning of the outer loop.
//	            outer : do
	//
//	                jeval = .true.
	//
//	                ! calculate the jacobian matrix.
	//
//	                iflag = 2
//	                call fcn(n, x, Fvec, Fjac, Ldfjac, iflag)
//	                Njev = Njev + 1
//	                if (iflag < 0) exit main
	//
//	                ! compute the qr factorization of the jacobian.
	//
//	                call qrfac(n, n, Fjac, Ldfjac, .false., iwa, 1, Wa1, Wa2, Wa3)
	//
//	                ! on the first iteration and if mode is 1, scale according
//	                ! to the norms of the columns of the initial jacobian.
	//
//	                if (iter == 1) then
//	                    if (Mode /= 2) then
//	                        do j = 1, n
//	                            Diag(j) = Wa2(j)
//	                            if (Wa2(j) == zero) Diag(j) = one
//	                        end do
//	                    end if
	//
//	                    ! on the first iteration, calculate the norm of the scaled x
//	                    ! and initialize the step bound delta.
	//
//	                    do j = 1, n
//	                        Wa3(j) = Diag(j)*x(j)
//	                    end do
//	                    xnorm = enorm(n, Wa3)
//	                    delta = Factor*xnorm
//	                    if (delta == zero) delta = Factor
//	                end if
	//
//	                ! form (q transpose)*fvec and store in qtf.
	//
//	                do i = 1, n
//	                    Qtf(i) = Fvec(i)
//	                end do
//	                do j = 1, n
//	                    if (Fjac(j, j) /= zero) then
//	                        sum = zero
//	                        do i = j, n
//	                            sum = sum + Fjac(i, j)*Qtf(i)
//	                        end do
//	                        temp = -sum/Fjac(j, j)
//	                        do i = j, n
//	                            Qtf(i) = Qtf(i) + Fjac(i, j)*temp
//	                        end do
//	                    end if
//	                end do
	//
//	                ! copy the triangular factor of the qr factorization into r.
	//
//	                sing = .false.
//	                do j = 1, n
//	                    l = j
//	                    jm1 = j - 1
//	                    if (jm1 >= 1) then
//	                        do i = 1, jm1
//	                            r(l) = Fjac(i, j)
//	                            l = l + n - i
//	                        end do
//	                    end if
//	                    r(l) = Wa1(j)
//	                    if (Wa1(j) == zero) sing = .true.
//	                end do
	//
//	                ! accumulate the orthogonal factor in fjac.
	//
//	                call qform(n, n, Fjac, Ldfjac, Wa1)
	//
//	                ! rescale if necessary.
	//
//	                if (Mode /= 2) then
//	                    do j = 1, n
//	                        Diag(j) = max(Diag(j), Wa2(j))
//	                    end do
//	                end if
	//
//	                ! beginning of the inner loop.
//	                inner : do
	//
//	                    ! if requested, call fcn to enable printing of iterates.
	//
//	                    if (Nprint > 0) then
//	                        iflag = 0
//	                        if (mod(iter - 1, Nprint) == 0) &
//	                            call fcn(n, x, Fvec, Fjac, Ldfjac, iflag)
//	                        if (iflag < 0) exit main
//	                    end if
	//
//	                    ! determine the direction p.
	//
//	                    call dogleg(n, r, Lr, Diag, Qtf, delta, Wa1, Wa2, Wa3)
	//
//	                    ! store the direction p and x + p. calculate the norm of p.
	//
//	                    do j = 1, n
//	                        Wa1(j) = -Wa1(j)
//	                        Wa2(j) = x(j) + Wa1(j)
//	                        Wa3(j) = Diag(j)*Wa1(j)
//	                    end do
//	                    pnorm = enorm(n, Wa3)
	//
//	                    ! on the first iteration, adjust the initial step bound.
	//
//	                    if (iter == 1) delta = min(delta, pnorm)
	//
//	                    ! evaluate the function at x + p and calculate its norm.
	//
//	                    iflag = 1
//	                    call fcn(n, Wa2, Wa4, Fjac, Ldfjac, iflag)
//	                    Nfev = Nfev + 1
//	                    if (iflag < 0) exit main
	//
//	                    fnorm1 = enorm(n, Wa4)
	//
//	                    ! compute the scaled actual reduction.
	//
//	                    actred = -one
//	                    if (fnorm1 < fnorm) actred = one - (fnorm1/fnorm)**2
	//
//	                    ! compute the scaled predicted reduction.
	//
//	                    l = 1
//	                    do i = 1, n
//	                        sum = zero
//	                        do j = i, n
//	                            sum = sum + r(l)*Wa1(j)
//	                            l = l + 1
//	                        end do
//	                        Wa3(i) = Qtf(i) + sum
//	                    end do
//	                    temp = enorm(n, Wa3)
//	                    prered = zero
//	                    if (temp < fnorm) prered = one - (temp/fnorm)**2
	//
//	                    ! compute the ratio of the actual to the predicted
//	                    ! reduction.
	//
//	                    ratio = zero
//	                    if (prered > zero) ratio = actred/prered
	//
//	                    ! update the step bound.
	//
//	                    if (ratio >= p1) then
//	                        ncfail = 0
//	                        ncsuc = ncsuc + 1
//	                        if (ratio >= p5 .or. ncsuc > 1) delta = max(delta, pnorm/p5)
//	                        if (abs(ratio - one) <= p1) delta = pnorm/p5
//	                    else
//	                        ncsuc = 0
//	                        ncfail = ncfail + 1
//	                        delta = p5*delta
//	                    end if
	//
//	                    ! test for successful iteration.
	//
//	                    if (ratio >= p0001) then
	//
//	                        ! successful iteration. update x, fvec, and their norms.
	//
//	                        do j = 1, n
//	                            x(j) = Wa2(j)
//	                            Wa2(j) = Diag(j)*x(j)
//	                            Fvec(j) = Wa4(j)
//	                        end do
//	                        xnorm = enorm(n, Wa2)
//	                        fnorm = fnorm1
//	                        iter = iter + 1
//	                    end if
	//
//	                    ! determine the progress of the iteration.
	//
//	                    nslow1 = nslow1 + 1
//	                    if (actred >= p001) nslow1 = 0
//	                    if (jeval) nslow2 = nslow2 + 1
//	                    if (actred >= p1) nslow2 = 0
	//
//	                    ! test for convergence.
	//
//	                    if (delta <= Xtol*xnorm .or. fnorm == zero) Info = 1
//	                    if (Info /= 0) exit main
	//
//	                    ! tests for termination and stringent tolerances.
	//
//	                    if (Nfev >= Maxfev) Info = 2
//	                    if (p1*max(p1*delta, pnorm) <= epsmch*xnorm) Info = 3
//	                    if (nslow2 == 5) Info = 4
//	                    if (nslow1 == 10) Info = 5
//	                    if (Info /= 0) exit main
	//
//	                    ! criterion for recalculating jacobian.
	//
//	                    if (ncfail == 2) cycle outer
	//
//	                    ! calculate the rank one modification to the jacobian
//	                    ! and update qtf if necessary.
	//
//	                    do j = 1, n
//	                        sum = zero
//	                        do i = 1, n
//	                            sum = sum + Fjac(i, j)*Wa4(i)
//	                        end do
//	                        Wa2(j) = (sum - Wa3(j))/pnorm
//	                        Wa1(j) = Diag(j)*((Diag(j)*Wa1(j))/pnorm)
//	                        if (ratio >= p0001) Qtf(j) = sum
//	                    end do
	//
//	                    ! compute the qr factorization of the updated jacobian.
	//
//	                    call r1updt(n, n, r, Lr, Wa1, Wa2, Wa3, sing)
//	                    call r1mpyq(n, n, Fjac, Ldfjac, Wa2, Wa3)
//	                    call r1mpyq(1, n, Qtf, 1, Wa2, Wa3)
	//
//	                    jeval = .false.
	//
//	                end do inner  ! end of the inner loop.
	//
//	            end do outer  ! end of the outer loop.
	//
//	        end block main
	//
//	        ! termination, either normal or user imposed.
	//
//	        if (iflag < 0) Info = iflag
//	        iflag = 0
//	        if (Nprint > 0) call fcn(n, x, Fvec, Fjac, Ldfjac, iflag)
	//
//	    end subroutine hybrj
	// !*****************************************************************************************
	//
	// !*****************************************************************************************
	// !>
	// ! the purpose of hybrj1 is to find a zero of a system of
	// ! n nonlinear functions in n variables by a modification
	// ! of the powell hybrid method. this is done by using the
	// ! more general nonlinear equation solver hybrj. the user
	// ! must provide a subroutine which calculates the functions
	// ! and the jacobian.
	//
//	    subroutine hybrj1(fcn, n, x, Fvec, Fjac, Ldfjac, Tol, Info, Wa, Lwa)
	//
//	        implicit none
	//
//	        procedure(fcn_hybrj) :: fcn !! the user-supplied subroutine which
//	                                !! calculates the functions and the jacobian
//	        integer, intent(in) :: n !! a positive integer input variable set to the number
//	                            !! of functions and variables.
//	        integer, intent(in) :: Ldfjac !! a positive integer input variable not less than n
//	                                 !! which specifies the leading dimension of the array fjac.
//	        integer, intent(out) :: Info !! an integer output variable. if the user has
//	                                !! terminated execution, info is set to the (negative)
//	                                !! value of iflag. see description of fcn. otherwise,
//	                                !! info is set as follows:
//	                                !!
//	                                !!  * ***info = 0***   improper input parameters.
//	                                !!  * ***info = 1***   algorithm estimates that the relative error
//	                                !!    between x and the solution is at most tol.
//	                                !!  * ***info = 2***   number of calls to fcn with iflag = 1 has
//	                                !!    reached 100*(n+1).
//	                                !!  * ***info = 3***   tol is too small. no further improvement in
//	                                !!    the approximate solution x is possible.
//	                                !!  * ***info = 4***   iteration is not making good progress.
//	        integer, intent(in) :: Lwa !! a positive integer input variable not less than
//	                              !! (n*(n+13))/2.
//	        real(wp), intent(in) :: Tol !! a nonnegative input variable. termination occurs
//	                                !! when the algorithm estimates that the relative error
//	                                !! between x and the solution is at most tol.
//	        real(wp), intent(inout) :: x(n) !! an array of length n. on input x must contain
//	                                    !! an initial estimate of the solution vector. on output x
//	                                    !! contains the final estimate of the solution vector.
//	        real(wp), intent(out) :: Fvec(n) !! an output array of length n which contains
//	                                    !! the functions evaluated at the output x.
//	        real(wp), intent(out) :: Fjac(Ldfjac, n) !! an output n by n array which contains the
//	                                            !! orthogonal matrix q produced by the qr factorization
//	                                            !! of the final approximate jacobian.
//	        real(wp), intent(inout) :: Wa(Lwa) !! a work array of length lwa.
	//
		public static void hyberj1() {
			
		}
//	        integer :: j, lr, maxfev, mode, nfev, njev, nprint
//	        real(wp) :: xtol
	//
//	        real(wp), parameter :: factor = 100.0_wp
	//
//	        Info = 0
	//
//	        ! check the input parameters for errors.
	//
//	        if (n > 0 .and. Ldfjac >= n .and. Tol >= zero .and. Lwa >= (n*(n + 13))/2) then
//	            ! call hybrj.
//	            maxfev = 100*(n + 1)
//	            xtol = Tol
//	            mode = 2
//	            do j = 1, n
//	                Wa(j) = one
//	            end do
//	            nprint = 0
//	            lr = (n*(n + 1))/2
//	            call hybrj(fcn, n, x, Fvec, Fjac, Ldfjac, xtol, maxfev, Wa(1), mode, &
//	                       factor, nprint, Info, nfev, njev, Wa(6*n + 1), lr, Wa(n + 1), &
//	                       Wa(2*n + 1), Wa(3*n + 1), Wa(4*n + 1), Wa(5*n + 1))
//	            if (Info == 5) Info = 4
//	        end if
	//
//	    end subroutine hybrj1
	// !*****************************************************************************************
	//
	// !*****************************************************************************************
	// !>
	// ! the purpose of lmder is to minimize the sum of the squares of
	// ! m nonlinear functions in n variables by a modification of
	// ! the levenberg-marquardt algorithm. the user must provide a
	// ! subroutine which calculates the functions and the jacobian.
	//
//	    subroutine lmder(fcn, m, n, x, Fvec, Fjac, Ldfjac, Ftol, Xtol, Gtol, Maxfev, &
//	                     Diag, Mode, Factor, Nprint, Info, Nfev, Njev, Ipvt, Qtf, &
//	                     Wa1, Wa2, Wa3, Wa4)
	//
//	        implicit none
	//
//	        procedure(fcn_lmder) :: fcn !! the user-supplied subroutine which
//	                                !! calculates the functions and the jacobian
//	        integer, intent(in) :: m !! a positive integer input variable set to the number
//	                            !! of functions.
//	        integer, intent(in) :: n !! a positive integer input variable set to the number
//	                            !! of variables. n must not exceed m.
//	        integer, intent(in) :: Ldfjac !! a positive integer input variable not less than m
//	                                 !! which specifies the leading dimension of the array fjac.
//	        integer, intent(in) :: Maxfev !! a positive integer input variable. termination
//	                                 !! occurs when the number of calls to fcn with iflag = 1
//	                                 !! has reached maxfev.
//	        integer, intent(in) :: Mode !! an integer input variable. if mode = 1, the
//	                                !! variables will be scaled internally. if mode = 2,
//	                                !! the scaling is specified by the input diag. other
//	                                !! values of mode are equivalent to mode = 1.
//	        integer, intent(in) :: Nprint !! an integer input variable that enables controlled
//	                                 !! printing of iterates if it is positive. in this case,
//	                                 !! fcn is called with iflag = 0 at the beginning of the first
//	                                 !! iteration and every nprint iterations thereafter and
//	                                 !! immediately prior to return, with x, fvec, and fjac
//	                                 !! available for printing. fvec and fjac should not be
//	                                 !! altered. if nprint is not positive, no special calls
//	                                 !! of fcn with iflag = 0 are made.
//	        integer, intent(out) :: Info !! an integer output variable. if the user has
//	                                !! terminated execution, info is set to the (negative)
//	                                !! value of iflag. see description of fcn. otherwise,
//	                                !! info is set as follows:
//	                                !!
//	                                !!  * ***info = 0***  improper input parameters.
//	                                !!  * ***info = 1***  both actual and predicted relative reductions
//	                                !!    in the sum of squares are at most ftol.
//	                                !!  * ***info = 2***  relative error between two consecutive iterates
//	                                !!    is at most xtol.
//	                                !!  * ***info = 3***  conditions for info = 1 and info = 2 both hold.
//	                                !!  * ***info = 4***  the cosine of the angle between fvec and any
//	                                !!    column of the jacobian is at most gtol in
//	                                !!    absolute value.
//	                                !!  * ***info = 5***  number of calls to fcn with iflag = 1 has
//	                                !!    reached maxfev.
//	                                !!  * ***info = 6***  ftol is too small. no further reduction in
//	                                !!    the sum of squares is possible.
//	                                !!  * ***info = 7***  xtol is too small. no further improvement in
//	                                !!    the approximate solution x is possible.
//	                                !!  * ***info = 8***  gtol is too small. fvec is orthogonal to the
//	                                !!    columns of the jacobian to machine precision.
//	        integer, intent(out) :: Nfev !! an integer output variable set to the number of
//	                                !! calls to fcn with iflag = 1.
//	        integer, intent(out) :: Njev !! an integer output variable set to the number of
//	                                !! calls to fcn with iflag = 2.
//	        integer, intent(out) :: Ipvt(n) !! an integer output array of length n. ipvt
//	                                   !! defines a permutation matrix p such that jac*p = q*r,
//	                                   !! where jac is the final calculated jacobian, q is
//	                                   !! orthogonal (not stored), and r is upper triangular
//	                                   !! with diagonal elements of nonincreasing magnitude.
//	                                   !! column j of p is column ipvt(j) of the identity matrix.
//	        real(wp), intent(in) :: Ftol !! a nonnegative input variable. termination
//	                                !! occurs when both the actual and predicted relative
//	                                !! reductions in the sum of squares are at most ftol.
//	                                !! therefore, ftol measures the relative error desired
//	                                !! in the sum of squares.
//	        real(wp), intent(in) :: Xtol !! a nonnegative input variable. termination
//	                                !! occurs when the relative error between two consecutive
//	                                !! iterates is at most xtol. therefore, xtol measures the
//	                                !! relative error desired in the approximate solution.
//	        real(wp), intent(in) :: Gtol !! a nonnegative input variable. termination
//	                                !! occurs when the cosine of the angle between fvec and
//	                                !! any column of the jacobian is at most gtol in absolute
//	                                !! value. therefore, gtol measures the orthogonality
//	                                !! desired between the function vector and the columns
//	                                !! of the jacobian.
//	        real(wp), intent(in) :: Factor !! a positive input variable used in determining the
//	                                  !! initial step bound. this bound is set to the product of
//	                                  !! factor and the euclidean norm of diag*x if nonzero, or else
//	                                  !! to factor itself. in most cases factor should lie in the
//	                                  !! interval (.1,100.).100. is a generally recommended value.
//	        real(wp), intent(inout) :: x(n) !! an array of length n. on input x must contain
//	                                   !! an initial estimate of the solution vector. on output x
//	                                   !! contains the final estimate of the solution vector.
//	        real(wp), intent(out) :: Fvec(m) !! an output array of length m which contains
//	                                    !! the functions evaluated at the output x.
//	        real(wp), intent(out) :: Fjac(Ldfjac, n) !! an output m by n array. the upper n by n submatrix
//	                                            !! of fjac contains an upper triangular matrix r with
//	                                            !! diagonal elements of nonincreasing magnitude such that
//	                                            !!```
//	                                            !!        t     t           t
//	                                            !!       p *(jac *jac)*p = r *r,
//	                                            !!```
//	                                            !! where p is a permutation matrix and jac is the final
//	                                            !! calculated jacobian. column j of p is column ipvt(j)
//	                                            !! (see below) of the identity matrix. the lower trapezoidal
//	                                            !! part of fjac contains information generated during
//	                                            !! the computation of r.
//	        real(wp), intent(inout) :: Diag(n) !! an array of length n. if mode = 1 (see
//	                                      !! below), diag is internally set. if mode = 2, diag
//	                                      !! must contain positive entries that serve as
//	                                      !! multiplicative scale factors for the variables.
//	        real(wp), intent(out) :: Qtf(n) !! an output array of length n which contains
//	                                   !! the first n elements of the vector (q transpose)*fvec.
//	        real(wp), intent(inout) :: Wa1(n) !! work array of length n.
//	        real(wp), intent(inout) :: Wa2(n) !! work array of length n.
//	        real(wp), intent(inout) :: Wa3(n) !! work array of length n.
//	        real(wp), intent(inout) :: Wa4(m) !! work array of length m.
	//
		public static void lmder() {
			
		}
//	        integer :: i, iflag, iter, j, l
//	        real(wp) :: actred, delta, dirder, fnorm, fnorm1, gnorm, par, &
//	                    pnorm, prered, ratio, sum, temp, temp1, temp2, xnorm
	//
//	        real(wp), parameter :: p1 = 1.0e-1_wp
//	        real(wp), parameter :: p5 = 5.0e-1_wp
//	        real(wp), parameter :: p25 = 2.5e-1_wp
//	        real(wp), parameter :: p75 = 7.5e-1_wp
//	        real(wp), parameter :: p0001 = 1.0e-4_wp
	//
//	        Info = 0
//	        iflag = 0
//	        Nfev = 0
//	        Njev = 0
	//
//	        main : block
	//
//	            ! check the input parameters for errors.
	//
//	            if (n > 0 .and. m >= n .and. Ldfjac >= m .and. Ftol >= zero .and. &
//	                Xtol >= zero .and. Gtol >= zero .and. Maxfev > 0 .and. &
//	                Factor > zero) then
//	                if (Mode == 2) then
//	                    do j = 1, n
//	                        if (Diag(j) <= zero) exit main
//	                    end do
//	                end if
//	            else
//	                exit main
//	            end if
	//
//	            ! evaluate the function at the starting point
//	            ! and calculate its norm.
	//
//	            iflag = 1
//	            call fcn(m, n, x, Fvec, Fjac, Ldfjac, iflag)
//	            Nfev = 1
//	            if (iflag < 0) exit main
//	            fnorm = enorm(m, Fvec)
	//
//	            ! initialize levenberg-marquardt parameter and iteration counter.
	//
//	            par = zero
//	            iter = 1
	//
//	            ! beginning of the outer loop.
	//
//	            outer : do
	//
//	                ! calculate the jacobian matrix.
	//
//	                iflag = 2
//	                call fcn(m, n, x, Fvec, Fjac, Ldfjac, iflag)
//	                Njev = Njev + 1
//	                if (iflag < 0) exit main
	//
//	                ! if requested, call fcn to enable printing of iterates.
	//
//	                if (Nprint > 0) then
//	                    iflag = 0
//	                    if (mod(iter - 1, Nprint) == 0) &
//	                        call fcn(m, n, x, Fvec, Fjac, Ldfjac, iflag)
//	                    if (iflag < 0) exit main
//	                end if
	//
//	                ! compute the qr factorization of the jacobian.
	//
//	                call qrfac(m, n, Fjac, Ldfjac, .true., Ipvt, n, Wa1, Wa2, Wa3)
	//
//	                ! on the first iteration and if mode is 1, scale according
//	                ! to the norms of the columns of the initial jacobian.
	//
//	                if (iter == 1) then
//	                    if (Mode /= 2) then
//	                        do j = 1, n
//	                            Diag(j) = Wa2(j)
//	                            if (Wa2(j) == zero) Diag(j) = one
//	                        end do
//	                    end if
	//
//	                    ! on the first iteration, calculate the norm of the scaled x
//	                    ! and initialize the step bound delta.
	//
//	                    do j = 1, n
//	                        Wa3(j) = Diag(j)*x(j)
//	                    end do
//	                    xnorm = enorm(n, Wa3)
//	                    delta = Factor*xnorm
//	                    if (delta == zero) delta = Factor
//	                end if
	//
//	                ! form (q transpose)*fvec and store the first n components in
//	                ! qtf.
	//
//	                do i = 1, m
//	                    Wa4(i) = Fvec(i)
//	                end do
//	                do j = 1, n
//	                    if (Fjac(j, j) /= zero) then
//	                        sum = zero
//	                        do i = j, m
//	                            sum = sum + Fjac(i, j)*Wa4(i)
//	                        end do
//	                        temp = -sum/Fjac(j, j)
//	                        do i = j, m
//	                            Wa4(i) = Wa4(i) + Fjac(i, j)*temp
//	                        end do
//	                    end if
//	                    Fjac(j, j) = Wa1(j)
//	                    Qtf(j) = Wa4(j)
//	                end do
	//
//	                ! compute the norm of the scaled gradient.
	//
//	                gnorm = zero
//	                if (fnorm /= zero) then
//	                    do j = 1, n
//	                        l = Ipvt(j)
//	                        if (Wa2(l) /= zero) then
//	                            sum = zero
//	                            do i = 1, j
//	                                sum = sum + Fjac(i, j)*(Qtf(i)/fnorm)
//	                            end do
//	                            gnorm = max(gnorm, abs(sum/Wa2(l)))
//	                        end if
//	                    end do
//	                end if
	//
//	                ! test for convergence of the gradient norm.
	//
//	                if (gnorm <= Gtol) Info = 4
//	                if (Info /= 0) exit main
	//
//	                ! rescale if necessary.
	//
//	                if (Mode /= 2) then
//	                    do j = 1, n
//	                        Diag(j) = max(Diag(j), Wa2(j))
//	                    end do
//	                end if
	//
//	                ! beginning of the inner loop.
//	                inner : do
	//
//	                    ! determine the levenberg-marquardt parameter.
	//
//	                    call lmpar(n, Fjac, Ldfjac, Ipvt, Diag, Qtf, delta, par, Wa1, Wa2, Wa3, Wa4)
	//
//	                    ! store the direction p and x + p. calculate the norm of p.
	//
//	                    do j = 1, n
//	                        Wa1(j) = -Wa1(j)
//	                        Wa2(j) = x(j) + Wa1(j)
//	                        Wa3(j) = Diag(j)*Wa1(j)
//	                    end do
//	                    pnorm = enorm(n, Wa3)
	//
//	                    ! on the first iteration, adjust the initial step bound.
	//
//	                    if (iter == 1) delta = min(delta, pnorm)
	//
//	                    ! evaluate the function at x + p and calculate its norm.
	//
//	                    iflag = 1
//	                    call fcn(m, n, Wa2, Wa4, Fjac, Ldfjac, iflag)
//	                    Nfev = Nfev + 1
//	                    if (iflag < 0) exit main
//	                    fnorm1 = enorm(m, Wa4)
	//
//	                    ! compute the scaled actual reduction.
	//
//	                    actred = -one
//	                    if (p1*fnorm1 < fnorm) actred = one - (fnorm1/fnorm)**2
	//
//	                    ! compute the scaled predicted reduction and
//	                    ! the scaled directional derivative.
	//
//	                    do j = 1, n
//	                        Wa3(j) = zero
//	                        l = Ipvt(j)
//	                        temp = Wa1(l)
//	                        do i = 1, j
//	                            Wa3(i) = Wa3(i) + Fjac(i, j)*temp
//	                        end do
//	                    end do
//	                    temp1 = enorm(n, Wa3)/fnorm
//	                    temp2 = (sqrt(par)*pnorm)/fnorm
//	                    prered = temp1**2 + temp2**2/p5
//	                    dirder = -(temp1**2 + temp2**2)
	//
//	                    ! compute the ratio of the actual to the predicted
//	                    ! reduction.
	//
//	                    ratio = zero
//	                    if (prered /= zero) ratio = actred/prered
	//
//	                    ! update the step bound.
	//
//	                    if (ratio <= p25) then
//	                        if (actred >= zero) temp = p5
//	                        if (actred < zero) temp = p5*dirder/(dirder + p5*actred)
//	                        if (p1*fnorm1 >= fnorm .or. temp < p1) temp = p1
//	                        delta = temp*min(delta, pnorm/p1)
//	                        par = par/temp
//	                    elseif (par == zero .or. ratio >= p75) then
//	                        delta = pnorm/p5
//	                        par = p5*par
//	                    end if
	//
//	                    ! test for successful iteration.
	//
//	                    if (ratio >= p0001) then
//	                        ! successful iteration. update x, fvec, and their norms.
//	                        do j = 1, n
//	                            x(j) = Wa2(j)
//	                            Wa2(j) = Diag(j)*x(j)
//	                        end do
//	                        do i = 1, m
//	                            Fvec(i) = Wa4(i)
//	                        end do
//	                        xnorm = enorm(n, Wa2)
//	                        fnorm = fnorm1
//	                        iter = iter + 1
//	                    end if
	//
//	                    ! tests for convergence.
//	                    if (abs(actred) <= Ftol .and. prered <= Ftol .and. p5*ratio <= one) Info = 1
//	                    if (delta <= Xtol*xnorm) Info = 2
//	                    if (abs(actred) <= Ftol .and. prered <= Ftol .and. p5*ratio <= one .and. Info == 2) Info = 3
//	                    if (Info /= 0) exit main
	//
//	                    ! tests for termination and stringent tolerances.
//	                    if (Nfev >= Maxfev) Info = 5
//	                    if (abs(actred) <= epsmch .and. prered <= epsmch .and. p5*ratio <= one) Info = 6
//	                    if (delta <= epsmch*xnorm) Info = 7
//	                    if (gnorm <= epsmch) Info = 8
//	                    if (Info /= 0) exit main
	//
//	                    if (ratio >= p0001) exit inner
	//
//	                end do inner ! end of the inner loop. repeat if iteration unsuccessful.
	//
//	            end do outer ! end of the outer loop
	//
//	        end block main
	//
//	        ! termination, either normal or user imposed.
	//
//	        if (iflag < 0) Info = iflag
//	        iflag = 0
//	        if (Nprint > 0) call fcn(m, n, x, Fvec, Fjac, Ldfjac, iflag)
	//
//	    end subroutine lmder
	// !*****************************************************************************************
	//
	// !*****************************************************************************************
	// !>
	// ! the purpose of lmder1 is to minimize the sum of the squares of
	// ! m nonlinear functions in n variables by a modification of the
	// ! levenberg-marquardt algorithm. this is done by using the more
	// ! general least-squares solver lmder. the user must provide a
	// ! subroutine which calculates the functions and the jacobian.
	//
//	    subroutine lmder1(fcn, m, n, x, Fvec, Fjac, Ldfjac, Tol, Info, Ipvt, Wa, Lwa)
//	        implicit none
	//
//	        procedure(fcn_lmder) :: fcn !! user-supplied subroutine which
//	                                    !! calculates the functions and the jacobian.
//	        integer, intent(in) :: m !! a positive integer input variable set to the number
//	                                !! of functions.
//	        integer, intent(in) :: n !! a positive integer input variable set to the number
//	                                !! of variables. n must not exceed m.
//	        integer, intent(in) :: Ldfjac !! a positive integer input variable not less than m
//	                                     !! which specifies the leading dimension of the array fjac.
//	        integer, intent(out) :: Info !! an integer output variable. if the user has
//	                                    !! terminated execution, info is set to the (negative)
//	                                    !! value of iflag. see description of fcn. otherwise,
//	                                    !! info is set as follows.
//	                                    !!
//	                                    !!  * ***info = 0***  improper input parameters.
//	                                    !!  * ***info = 1***  algorithm estimates that the relative error
//	                                    !!    in the sum of squares is at most tol.
//	                                    !!  * ***info = 2***  algorithm estimates that the relative error
//	                                    !!    between x and the solution is at most tol.
//	                                    !!  * ***info = 3***  conditions for info = 1 and info = 2 both hold.
//	                                    !!  * ***info = 4***  fvec is orthogonal to the columns of the
//	                                    !!    jacobian to machine precision.
//	                                    !!  * ***info = 5***  number of calls to fcn with iflag = 1 has
//	                                    !!    reached 100*(n+1).
//	                                    !!  * ***info = 6***  tol is too small. no further reduction in
//	                                    !!    the sum of squares is possible.
//	                                    !!  * ***info = 7***  tol is too small. no further improvement in
//	                                    !!    the approximate solution x is possible.
//	        integer, intent(in) :: Lwa !! a positive integer input variable not less than 5*n+m.
//	        integer, intent(out) :: Ipvt(n) !! an integer output array of length n. ipvt
//	                                       !! defines a permutation matrix p such that jac*p = q*r,
//	                                       !! where jac is the final calculated jacobian, q is
//	                                       !! orthogonal (not stored), and r is upper triangular
//	                                       !! with diagonal elements of nonincreasing magnitude.
//	                                       !! column j of p is column ipvt(j) of the identity matrix.
//	        real(wp), intent(in) :: Tol !! a nonnegative input variable. termination occurs
//	                                   !! when the algorithm estimates either that the relative
//	                                   !! error in the sum of squares is at most tol or that
//	                                   !! the relative error between x and the solution is at
//	                                   !! most tol.
//	        real(wp), intent(inout) :: x(n) !! an array of length n. on input x must contain
//	                                       !! an initial estimate of the solution vector. on output x
//	                                       !! contains the final estimate of the solution vector.
//	        real(wp), intent(out) :: Fvec(m) !! an output array of length m which contains
//	                                        !! the functions evaluated at the output x.
//	        real(wp), intent(out) :: Fjac(Ldfjac, n) !! an output m by n array. the upper n by n submatrix
//	                                                !! of fjac contains an upper triangular matrix r with
//	                                                !! diagonal elements of nonincreasing magnitude such that
//	                                                !!```
//	                                                !!        t     t           t
//	                                                !!       p *(jac *jac)*p = r *r,
//	                                                !!```
//	                                                !! where p is a permutation matrix and jac is the final
//	                                                !! calculated jacobian. column j of p is column ipvt(j)
//	                                                !! (see below) of the identity matrix. the lower trapezoidal
//	                                                !! part of fjac contains information generated during
//	                                                !! the computation of r.
//	        real(wp), intent(inout) :: Wa(Lwa) !! a work array of length lwa.
	//
//	        integer :: maxfev, mode, nfev, njev, nprint
//	        real(wp) :: ftol, gtol, xtol
	//
//	        real(wp), parameter :: factor = 100.0_wp
	//		
		public static void lmder1() {
			
		}
//	        Info = 0
	//
			
//	        ! check the input parameters for errors.
	//
//	        if (n > 0 .and. m >= n .and. Ldfjac >= m .and. Tol >= zero .and. &
//	            Lwa >= 5*n + m) then
//	            ! call lmder.
//	            maxfev = 100*(n + 1)
//	            ftol = Tol
//	            xtol = Tol
//	            gtol = zero
//	            mode = 1
//	            nprint = 0
//	            call lmder(fcn, m, n, x, Fvec, Fjac, Ldfjac, ftol, xtol, gtol, maxfev,   &
//	                     & Wa(1), mode, factor, nprint, Info, nfev, njev, Ipvt, Wa(n + 1)&
//	                     & , Wa(2*n + 1), Wa(3*n + 1), Wa(4*n + 1), Wa(5*n + 1))
//	            if (Info == 8) Info = 4
//	        end if
	//
//	    end subroutine lmder1
	// !*****************************************************************************************
	//
	// !*****************************************************************************************
	// !>
	// ! the purpose of lmdif is to minimize the sum of the squares of
	// ! m nonlinear functions in n variables by a modification of
	// ! the levenberg-marquardt algorithm. the user must provide a
	// ! subroutine which calculates the functions. the jacobian is
	// ! then calculated by a forward-difference approximation.
	//
//	    subroutine lmdif(fcn, m, n, x, Fvec, Ftol, Xtol, Gtol, Maxfev, Epsfcn, Diag, &
//	                     Mode, Factor, Nprint, Info, Nfev, Fjac, Ldfjac, Ipvt, &
//	                     Qtf, Wa1, Wa2, Wa3, Wa4)
//	        implicit none
	//
//	        procedure(func2) :: fcn !! the user-supplied subroutine which
//	                                !! calculates the functions.
//	        integer, intent(in) :: m !! a positive integer input variable set to the number
//	                                !! of functions.
//	        integer, intent(in) :: n !! a positive integer input variable set to the number
//	                                !! of variables. n must not exceed m.
//	        integer, intent(in) :: Maxfev !! a positive integer input variable. termination
//	                                     !! occurs when the number of calls to fcn is at least
//	                                     !! maxfev by the end of an iteration.
//	        integer, intent(in) :: Mode !! an integer input variable. if mode = 1, the
//	                                   !! variables will be scaled internally. if mode = 2,
//	                                   !! the scaling is specified by the input diag. other
//	                                   !! values of mode are equivalent to mode = 1.
//	        integer, intent(in) :: Nprint !! an integer input variable that enables controlled
//	                                     !! printing of iterates if it is positive. in this case,
//	                                     !! fcn is called with iflag = 0 at the beginning of the first
//	                                     !! iteration and every nprint iterations thereafter and
//	                                     !! immediately prior to return, with x and fvec available
//	                                     !! for printing. if nprint is not positive, no special calls
//	                                     !! of fcn with iflag = 0 are made.
//	        integer, intent(out) :: Info !! an integer output variable. if the user has
//	                                    !! terminated execution, info is set to the (negative)
//	                                    !! value of iflag. see description of fcn. otherwise,
//	                                    !! info is set as follows:
//	                                    !!
//	                                    !!  * ***info = 0***  improper input parameters.
//	                                    !!  * ***info = 1***  both actual and predicted relative reductions
//	                                    !!    in the sum of squares are at most ftol.
//	                                    !!  * ***info = 2***  relative error between two consecutive iterates
//	                                    !!    is at most xtol.
//	                                    !!  * ***info = 3***  conditions for info = 1 and info = 2 both hold.
//	                                    !!  * ***info = 4***  the cosine of the angle between fvec and any
//	                                    !!    column of the jacobian is at most gtol in
//	                                    !!    absolute value.
//	                                    !!  * ***info = 5***  number of calls to fcn has reached or
//	                                    !!    exceeded maxfev.
//	                                    !!  * ***info = 6***  ftol is too small. no further reduction in
//	                                    !!    the sum of squares is possible.
//	                                    !!  * ***info = 7***  xtol is too small. no further improvement in
//	                                    !!    the approximate solution x is possible.
//	                                    !!  * ***info = 8***  gtol is too small. fvec is orthogonal to the
//	                                    !!    columns of the jacobian to machine precision.
//	        integer, intent(out) :: Nfev !! an integer output variable set to the number of
//	                                    !! calls to fcn.
//	        integer, intent(in) :: Ldfjac !! a positive integer input variable not less than m
//	                                     !! which specifies the leading dimension of the array fjac.
//	        integer, intent(out) :: Ipvt(n) !! an integer output array of length n. ipvt
//	                                       !! defines a permutation matrix p such that jac*p = q*r,
//	                                       !! where jac is the final calculated jacobian, q is
//	                                       !! orthogonal (not stored), and r is upper triangular
//	                                       !! with diagonal elements of nonincreasing magnitude.
//	                                       !! column j of p is column ipvt(j) of the identity matrix.
//	        real(wp), intent(in) :: Ftol !! a nonnegative input variable. termination
//	                                    !! occurs when both the actual and predicted relative
//	                                    !! reductions in the sum of squares are at most ftol.
//	                                    !! therefore, ftol measures the relative error desired
//	                                    !! in the sum of squares.
//	        real(wp), intent(in) :: Xtol !! a nonnegative input variable. termination
//	                                    !! occurs when the relative error between two consecutive
//	                                    !! iterates is at most xtol. therefore, xtol measures the
//	                                    !! relative error desired in the approximate solution.
//	        real(wp), intent(in) :: Gtol !! a nonnegative input variable. termination
//	                                    !! occurs when the cosine of the angle between fvec and
//	                                    !! any column of the jacobian is at most gtol in absolute
//	                                    !! value. therefore, gtol measures the orthogonality
//	                                    !! desired between the function vector and the columns
//	                                    !! of the jacobian.
//	        real(wp), intent(in) :: Epsfcn !! an input variable used in determining a suitable
//	                                      !! step length for the forward-difference approximation. this
//	                                      !! approximation assumes that the relative errors in the
//	                                      !! functions are of the order of epsfcn. if epsfcn is less
//	                                      !! than the machine precision, it is assumed that the relative
//	                                      !! errors in the functions are of the order of the machine
//	                                      !! precision.
//	        real(wp), intent(in) :: Factor !! a positive input variable used in determining the
//	                                      !! initial step bound. this bound is set to the product of
//	                                      !! factor and the euclidean norm of diag*x if nonzero, or else
//	                                      !! to factor itself. in most cases factor should lie in the
//	                                      !! interval (.1,100.). 100. is a generally recommended value.
//	        real(wp), intent(inout) :: x(n) !!  an array of length n. on input x must contain
//	                                       !! an initial estimate of the solution vector. on output x
//	                                       !! contains the final estimate of the solution vector.
//	        real(wp), intent(out) :: Fvec(m) !! an output array of length m which contains
//	                                        !! the functions evaluated at the output x.
//	        real(wp), intent(inout) :: Diag(n) !! an array of length n. if mode = 1 (see
//	                                          !! below), diag is internally set. if mode = 2, diag
//	                                          !! must contain positive entries that serve as
//	                                          !! multiplicative scale factors for the variables.
//	        real(wp), intent(out) :: Fjac(Ldfjac, n) !! an output m by n array. the upper n by n submatrix
//	                                                !! of fjac contains an upper triangular matrix r with
//	                                                !! diagonal elements of nonincreasing magnitude such that
//	                                                !!```
//	                                                !!        t     t           t
//	                                                !!       p *(jac *jac)*p = r *r,
//	                                                !!```
//	                                                !! where p is a permutation matrix and jac is the final
//	                                                !! calculated jacobian. column j of p is column ipvt(j)
//	                                                !! (see below) of the identity matrix. the lower trapezoidal
//	                                                !! part of fjac contains information generated during
//	                                                !! the computation of r.
//	        real(wp), intent(out) :: Qtf(n) !! an output array of length n which contains
//	                                       !! the first n elements of the vector (q transpose)*fvec.
//	        real(wp), intent(inout) :: Wa1(n) !! work array of length n.
//	        real(wp), intent(inout) :: Wa2(n) !! work array of length n.
//	        real(wp), intent(inout) :: Wa3(n) !! work array of length n.
//	        real(wp), intent(inout) :: Wa4(m) !! work array of length n.
	//
//	        integer :: i, iflag, iter, j, l
//	        real(wp) :: actred, delta, dirder, fnorm, &
//	                    fnorm1, gnorm, par, pnorm, prered, &
//	                    ratio, sum, temp, temp1, temp2, xnorm
	//
//	        real(wp), parameter :: p1 = 1.0e-1_wp
//	        real(wp), parameter :: p5 = 5.0e-1_wp
//	        real(wp), parameter :: p25 = 2.5e-1_wp
//	        real(wp), parameter :: p75 = 7.5e-1_wp
//	        real(wp), parameter :: p0001 = 1.0e-4_wp
		public static void lmdif() {
			
		}
//	        Info = 0
//	        iflag = 0
//	        Nfev = 0
	//
//	        main : block
	//
//	            ! check the input parameters for errors.
	//
//	            if (n > 0 .and. m >= n .and. Ldfjac >= m .and. Ftol >= zero .and. &
//	                Xtol >= zero .and. Gtol >= zero .and. Maxfev > 0 .and. &
//	                Factor > zero) then
//	                if (Mode == 2) then
//	                    do j = 1, n
//	                        if (Diag(j) <= zero) exit main
//	                    end do
//	                end if
//	            else
//	                exit main
//	            end if
	//
//	            ! evaluate the function at the starting point
//	            ! and calculate its norm.
	//
//	            iflag = 1
//	            call fcn(m, n, x, Fvec, iflag)
//	            Nfev = 1
//	            if (iflag < 0) exit main
	//
//	            fnorm = enorm(m, Fvec)
	//
//	            ! initialize levenberg-marquardt parameter and iteration counter.
	//
//	            par = zero
//	            iter = 1
	//
//	            ! beginning of the outer loop.
	//
//	            outer : do
	//
//	                ! calculate the jacobian matrix.
	//
//	                iflag = 2
//	                call fdjac2(fcn, m, n, x, Fvec, Fjac, Ldfjac, iflag, Epsfcn, Wa4)
//	                Nfev = Nfev + n
//	                if (iflag < 0) exit main
	//
//	                ! if requested, call fcn to enable printing of iterates.
	//
//	                if (Nprint > 0) then
//	                    iflag = 0
//	                    if (mod(iter - 1, Nprint) == 0) &
//	                        call fcn(m, n, x, Fvec, iflag)
//	                    if (iflag < 0) exit main
//	                end if
	//
//	                ! compute the qr factorization of the jacobian.
	//
//	                call qrfac(m, n, Fjac, Ldfjac, .true., Ipvt, n, Wa1, Wa2, Wa3)
	//
//	                ! on the first iteration and if mode is 1, scale according
//	                ! to the norms of the columns of the initial jacobian.
	//
//	                if (iter == 1) then
//	                    if (Mode /= 2) then
//	                        do j = 1, n
//	                            Diag(j) = Wa2(j)
//	                            if (Wa2(j) == zero) Diag(j) = one
//	                        end do
//	                    end if
	//
//	                    ! on the first iteration, calculate the norm of the scaled x
//	                    ! and initialize the step bound delta.
	//
//	                    do j = 1, n
//	                        Wa3(j) = Diag(j)*x(j)
//	                    end do
//	                    xnorm = enorm(n, Wa3)
//	                    delta = Factor*xnorm
//	                    if (delta == zero) delta = Factor
//	                end if
	//
//	                ! form (q transpose)*fvec and store the first n components in
//	                ! qtf.
	//
//	                do i = 1, m
//	                    Wa4(i) = Fvec(i)
//	                end do
//	                do j = 1, n
//	                    if (Fjac(j, j) /= zero) then
//	                        sum = zero
//	                        do i = j, m
//	                            sum = sum + Fjac(i, j)*Wa4(i)
//	                        end do
//	                        temp = -sum/Fjac(j, j)
//	                        do i = j, m
//	                            Wa4(i) = Wa4(i) + Fjac(i, j)*temp
//	                        end do
//	                    end if
//	                    Fjac(j, j) = Wa1(j)
//	                    Qtf(j) = Wa4(j)
//	                end do
	//
//	                ! compute the norm of the scaled gradient.
	//
//	                gnorm = zero
//	                if (fnorm /= zero) then
//	                    do j = 1, n
//	                        l = Ipvt(j)
//	                        if (Wa2(l) /= zero) then
//	                            sum = zero
//	                            do i = 1, j
//	                                sum = sum + Fjac(i, j)*(Qtf(i)/fnorm)
//	                            end do
//	                            gnorm = max(gnorm, abs(sum/Wa2(l)))
//	                        end if
//	                    end do
//	                end if
	//
//	                ! test for convergence of the gradient norm.
	//
//	                if (gnorm <= Gtol) Info = 4
//	                if (Info /= 0) exit main
	//
//	                ! rescale if necessary.
	//
//	                if (Mode /= 2) then
//	                    do j = 1, n
//	                        Diag(j) = max(Diag(j), Wa2(j))
//	                    end do
//	                end if
	//
//	                ! beginning of the inner loop.
	//
//	                inner : do
	//
//	                    ! determine the levenberg-marquardt parameter.
	//
//	                    call lmpar(n, Fjac, Ldfjac, Ipvt, Diag, Qtf, delta, par, Wa1, &
//	                               Wa2, Wa3, Wa4)
	//
//	                    ! store the direction p and x + p. calculate the norm of p.
	//
//	                    do j = 1, n
//	                        Wa1(j) = -Wa1(j)
//	                        Wa2(j) = x(j) + Wa1(j)
//	                        Wa3(j) = Diag(j)*Wa1(j)
//	                    end do
//	                    pnorm = enorm(n, Wa3)
	//
//	                    ! on the first iteration, adjust the initial step bound.
	//
//	                    if (iter == 1) delta = min(delta, pnorm)
	//
//	                    ! evaluate the function at x + p and calculate its norm.
	//
//	                    iflag = 1
//	                    call fcn(m, n, Wa2, Wa4, iflag)
//	                    Nfev = Nfev + 1
//	                    if (iflag < 0) exit main
	//
//	                    fnorm1 = enorm(m, Wa4)
	//
//	                    ! compute the scaled actual reduction.
	//
//	                    actred = -one
//	                    if (p1*fnorm1 < fnorm) actred = one - (fnorm1/fnorm)**2
	//
//	                    ! compute the scaled predicted reduction and
//	                    ! the scaled directional derivative.
	//
//	                    do j = 1, n
//	                        Wa3(j) = zero
//	                        l = Ipvt(j)
//	                        temp = Wa1(l)
//	                        do i = 1, j
//	                            Wa3(i) = Wa3(i) + Fjac(i, j)*temp
//	                        end do
//	                    end do
//	                    temp1 = enorm(n, Wa3)/fnorm
//	                    temp2 = (sqrt(par)*pnorm)/fnorm
//	                    prered = temp1**2 + temp2**2/p5
//	                    dirder = -(temp1**2 + temp2**2)
	//
//	                    ! compute the ratio of the actual to the predicted
//	                    ! reduction.
	//
//	                    ratio = zero
//	                    if (prered /= zero) ratio = actred/prered
	//
//	                    ! update the step bound.
	//
//	                    if (ratio <= p25) then
//	                        if (actred >= zero) temp = p5
//	                        if (actred < zero) &
//	                            temp = p5*dirder/(dirder + p5*actred)
//	                        if (p1*fnorm1 >= fnorm .or. temp < p1) temp = p1
//	                        delta = temp*min(delta, pnorm/p1)
//	                        par = par/temp
//	                    elseif (par == zero .or. ratio >= p75) then
//	                        delta = pnorm/p5
//	                        par = p5*par
//	                    end if
	//
//	                    ! test for successful iteration.
	//
//	                    if (ratio >= p0001) then
	//
//	                        ! successful iteration. update x, fvec, and their norms.
	//
//	                        do j = 1, n
//	                            x(j) = Wa2(j)
//	                            Wa2(j) = Diag(j)*x(j)
//	                        end do
//	                        do i = 1, m
//	                            Fvec(i) = Wa4(i)
//	                        end do
//	                        xnorm = enorm(n, Wa2)
//	                        fnorm = fnorm1
//	                        iter = iter + 1
//	                    end if
	//
//	                    ! tests for convergence.
	//
//	                    if (abs(actred) <= Ftol .and. prered <= Ftol .and. &
//	                        p5*ratio <= one) Info = 1
//	                    if (delta <= Xtol*xnorm) Info = 2
//	                    if (abs(actred) <= Ftol .and. prered <= Ftol .and. &
//	                        p5*ratio <= one .and. Info == 2) Info = 3
//	                    if (Info /= 0) exit main
	//
//	                    ! tests for termination and stringent tolerances.
	//
//	                    if (Nfev >= Maxfev) Info = 5
//	                    if (abs(actred) <= epsmch .and. &
//	                        prered <= epsmch .and. p5*ratio <= one) &
//	                        Info = 6
//	                    if (delta <= epsmch*xnorm) Info = 7
//	                    if (gnorm <= epsmch) Info = 8
//	                    if (Info /= 0) exit main
	//
//	                    if (ratio >= p0001) exit inner
	//
//	                end do inner ! end of the inner loop. repeat if iteration unsuccessful.
	//
//	            end do outer ! end of the outer loop.
	//
//	        end block main
	//
//	        ! termination, either normal or user imposed.
	//
//	        if (iflag < 0) Info = iflag
//	        iflag = 0
//	        if (Nprint > 0) call fcn(m, n, x, Fvec, iflag)
	//
//	    end subroutine lmdif
	// !*****************************************************************************************
	//
	// !*****************************************************************************************
	// !>
	// ! the purpose of lmdif1 is to minimize the sum of the squares of
	// ! m nonlinear functions in n variables by a modification of the
	// ! levenberg-marquardt algorithm. this is done by using the more
	// ! general least-squares solver lmdif. the user must provide a
	// ! subroutine which calculates the functions. the jacobian is
	// ! then calculated by a forward-difference approximation.
	//
//	    subroutine lmdif1(fcn, m, n, x, Fvec, Tol, Info, Iwa, Wa, Lwa)
//	        implicit none
	//
//	        procedure(func2) :: fcn !! the user-supplied subroutine which
//	                                !! calculates the functions.
//	        integer, intent(in) :: m !! a positive integer input variable set to the number
//	                                !! of functions.
//	        integer, intent(in) :: n !! a positive integer input variable set to the number
//	                                !! of variables. n must not exceed m.
//	        integer, intent(out) :: Info !! an integer output variable. if the user has
//	                                    !! terminated execution, info is set to the (negative)
//	                                    !! value of iflag. see description of fcn. otherwise,
//	                                    !! info is set as follows:
//	                                    !!
//	                                    !!  * ***info = 0***  improper input parameters.
//	                                    !!  * ***info = 1***  algorithm estimates that the relative error
//	                                    !!    in the sum of squares is at most tol.
//	                                    !!  * ***info = 2***  algorithm estimates that the relative error
//	                                    !!    between x and the solution is at most tol.
//	                                    !!  * ***info = 3***  conditions for info = 1 and info = 2 both hold.
//	                                    !!  * ***info = 4***  fvec is orthogonal to the columns of the
//	                                    !!    jacobian to machine precision.
//	                                    !!  * ***info = 5***  number of calls to fcn has reached or
//	                                    !!    exceeded 200*(n+1).
//	                                    !!  * ***info = 6***  tol is too small. no further reduction in
//	                                    !!    the sum of squares is possible.
//	                                    !!  * ***info = 7***  tol is too small. no further improvement in
//	                                    !!    the approximate solution x is possible.
//	        integer, intent(in) :: Lwa !! a positive integer input variable not less than
//	                                  !! m*n+5*n+m.
//	        integer, intent(inout) :: Iwa(n) !! an integer work array of length n.
//	        real(wp), intent(in) :: Tol !! a nonnegative input variable. termination occurs
//	                                   !! when the algorithm estimates either that the relative
//	                                   !! error in the sum of squares is at most tol or that
//	                                   !! the relative error between x and the solution is at
//	                                   !! most tol.
//	        real(wp), intent(inout) :: x(n) !! an array of length n. on input x must contain
//	                                       !! an initial estimate of the solution vector. on output x
//	                                       !! contains the final estimate of the solution vector.
//	        real(wp), intent(out) :: Fvec(m) !! an output array of length m which contains
//	                                        !! the functions evaluated at the output x.
//	        real(wp), intent(inout) :: Wa(Lwa) !! a work array of length lwa.
	//
		public static void lmdif1() {
			
		}
//	        integer :: maxfev, mode, mp5n, nfev, nprint
//	        real(wp) :: epsfcn, ftol, gtol, xtol
	//
//	        real(wp), parameter :: factor = 1.0e2_wp
	//
//	        Info = 0
	//
//	        ! check the input parameters for errors.
	//
//	        if (n > 0 .and. m >= n .and. Tol >= zero .and. Lwa >= m*n + 5*n + m) then
	//
//	            ! call lmdif.
	//
//	            maxfev = 200*(n + 1)
//	            ftol = Tol
//	            xtol = Tol
//	            gtol = zero
//	            epsfcn = zero
//	            mode = 1
//	            nprint = 0
//	            mp5n = m + 5*n
//	            call lmdif(fcn, m, n, x, Fvec, ftol, xtol, gtol, maxfev, epsfcn, Wa(1), &
//	                       mode, factor, nprint, Info, nfev, Wa(mp5n + 1), m, Iwa, &
//	                       Wa(n + 1), Wa(2*n + 1), Wa(3*n + 1), Wa(4*n + 1), Wa(5*n + 1))
//	            if (Info == 8) Info = 4
//	        end if
	//
//	    end subroutine lmdif1
	// !*****************************************************************************************
	//
	// !*****************************************************************************************
	// !>
	// ! given an m by n matrix a, an n by n nonsingular diagonal
	// ! matrix d, an m-vector b, and a positive number delta,
	// ! the problem is to determine a value for the parameter
	// ! par such that if x solves the system
	// !```
	// ! a*x = b , sqrt(par)*d*x = 0 ,
	// !```
	// ! in the least squares sense, and dxnorm is the euclidean
	// ! norm of d*x, then either par is zero and
	// !```
	// ! (dxnorm-delta) <= 0.1*delta ,
	// !```
	// ! or par is positive and
	// !```
	// ! abs(dxnorm-delta) <= 0.1*delta .
	// !```
	// ! this subroutine completes the solution of the problem
	// ! if it is provided with the necessary information from the
	// ! qr factorization, with column pivoting, of a. that is, if
	// ! a*p = q*r, where p is a permutation matrix, q has orthogonal
	// ! columns, and r is an upper triangular matrix with diagonal
	// ! elements of nonincreasing magnitude, then lmpar expects
	// ! the full upper triangle of r, the permutation matrix p,
	// ! and the first n components of (q transpose)*b. on output
	// ! lmpar also provides an upper triangular matrix s such that
	// !```
	// ! t t t
	// ! p *(a *a + par*d*d)*p = s *s .
	// !```
	// ! s is employed within lmpar and may be of separate interest.
	// !
	// ! only a few iterations are generally needed for convergence
	// ! of the algorithm. if, however, the limit of 10 iterations
	// ! is reached, then the output par will contain the best
	// ! value obtained so far.
	//
//	    subroutine lmpar(n, r, Ldr, Ipvt, Diag, Qtb, Delta, Par, x, Sdiag, Wa1, Wa2)
//	        implicit none
	//
//	        integer, intent(in) :: n !! a positive integer input variable set to the order of r.
//	        integer, intent(in) :: Ldr !! a positive integer input variable not less than n
//	                                  !! which specifies the leading dimension of the array r.
//	        integer, intent(in) :: Ipvt(n) !! an integer input array of length n which defines the
//	                                      !! permutation matrix p such that a*p = q*r. column j of p
//	                                      !! is column ipvt(j) of the identity matrix.
//	        real(wp) :: Delta !! a positive input variable which specifies an upper
//	                          !! bound on the euclidean norm of d*x.
//	        real(wp), intent(inout) :: Par !! a nonnegative variable. on input par contains an
//	                                      !! initial estimate of the levenberg-marquardt parameter.
//	                                      !! on output par contains the final estimate.
//	        real(wp), intent(inout) :: r(Ldr, n) !! an n by n array. on input the full upper triangle
//	                                            !! must contain the full upper triangle of the matrix r.
//	                                            !! on output the full upper triangle is unaltered, and the
//	                                            !! strict lower triangle contains the strict upper triangle
//	                                            !! (transposed) of the upper triangular matrix s.
//	        real(wp), intent(in) :: Diag(n) !! an input array of length n which must contain the
//	                                       !! diagonal elements of the matrix d.
//	        real(wp), intent(in) :: Qtb(n) !! an input array of length n which must contain the first
//	                                      !! n elements of the vector (q transpose)*b.
//	        real(wp), intent(out) :: x(n) !! an output array of length n which contains the least
//	                                     !! squares solution of the system a*x = b, sqrt(par)*d*x = 0,
//	                                     !! for the output par.
//	        real(wp), intent(out) :: Sdiag(n) !! an output array of length n which contains the
//	                                         !! diagonal elements of the upper triangular matrix s.
//	        real(wp), intent(inout) :: Wa1(n) !! work array of length n.
//	        real(wp), intent(inout) :: Wa2(n) !! work array of length n.
	//
		public static void lmpar() {
			
		}
//	        integer :: i, iter, j, jm1, jp1, k, l, nsing
//	        real(wp) :: dxnorm, fp, gnorm, parc, parl, paru, sum, temp
	//
//	        real(wp), parameter :: p1 = 1.0e-1_wp
//	        real(wp), parameter :: p001 = 1.0e-3_wp
//	        real(wp), parameter :: dwarf = dpmpar(2) !! the smallest positive magnitude
	//
//	        ! compute and store in x the gauss-newton direction. if the
//	        ! jacobian is rank-deficient, obtain a least squares solution.
	//
//	        nsing = n
//	        do j = 1, n
//	            Wa1(j) = Qtb(j)
//	            if (r(j, j) == zero .and. nsing == n) nsing = j - 1
//	            if (nsing < n) Wa1(j) = zero
//	        end do
//	        if (nsing >= 1) then
//	            do k = 1, nsing
//	                j = nsing - k + 1
//	                Wa1(j) = Wa1(j)/r(j, j)
//	                temp = Wa1(j)
//	                jm1 = j - 1
//	                if (jm1 >= 1) then
//	                    do i = 1, jm1
//	                        Wa1(i) = Wa1(i) - r(i, j)*temp
//	                    end do
//	                end if
//	            end do
//	        end if
//	        do j = 1, n
//	            l = Ipvt(j)
//	            x(l) = Wa1(j)
//	        end do
	//
//	        ! initialize the iteration counter.
//	        ! evaluate the function at the origin, and test
//	        ! for acceptance of the gauss-newton direction.
	//
//	        iter = 0
//	        do j = 1, n
//	            Wa2(j) = Diag(j)*x(j)
//	        end do
//	        dxnorm = enorm(n, Wa2)
//	        fp = dxnorm - Delta
//	        if (fp <= p1*Delta) then
//	            ! termination.
//	            if (iter == 0) Par = zero
//	        else
	//
//	            ! if the jacobian is not rank deficient, the newton
//	            ! step provides a lower bound, parl, for the zero of
//	            ! the function. otherwise set this bound to zero.
	//
//	            parl = zero
//	            if (nsing >= n) then
//	                do j = 1, n
//	                    l = Ipvt(j)
//	                    Wa1(j) = Diag(l)*(Wa2(l)/dxnorm)
//	                end do
//	                do j = 1, n
//	                    sum = zero
//	                    jm1 = j - 1
//	                    if (jm1 >= 1) then
//	                        do i = 1, jm1
//	                            sum = sum + r(i, j)*Wa1(i)
//	                        end do
//	                    end if
//	                    Wa1(j) = (Wa1(j) - sum)/r(j, j)
//	                end do
//	                temp = enorm(n, Wa1)
//	                parl = ((fp/Delta)/temp)/temp
//	            end if
	//
//	            ! calculate an upper bound, paru, for the zero of the function.
	//
//	            do j = 1, n
//	                sum = zero
//	                do i = 1, j
//	                    sum = sum + r(i, j)*Qtb(i)
//	                end do
//	                l = Ipvt(j)
//	                Wa1(j) = sum/Diag(l)
//	            end do
//	            gnorm = enorm(n, Wa1)
//	            paru = gnorm/Delta
//	            if (paru == zero) paru = dwarf/min(Delta, p1)
	//
//	            ! if the input par lies outside of the interval (parl,paru),
//	            ! set par to the closer endpoint.
	//
//	            Par = max(Par, parl)
//	            Par = min(Par, paru)
//	            if (Par == zero) Par = gnorm/dxnorm
	//
//	            ! beginning of an iteration.
//	            do
	//
//	                iter = iter + 1
	//
//	                ! evaluate the function at the current value of par.
	//
//	                if (Par == zero) Par = max(dwarf, p001*paru)
//	                temp = sqrt(Par)
//	                do j = 1, n
//	                    Wa1(j) = temp*Diag(j)
//	                end do
//	                call qrsolv(n, r, Ldr, Ipvt, Wa1, Qtb, x, Sdiag, Wa2)
//	                do j = 1, n
//	                    Wa2(j) = Diag(j)*x(j)
//	                end do
//	                dxnorm = enorm(n, Wa2)
//	                temp = fp
//	                fp = dxnorm - Delta
	//
//	                ! if the function is small enough, accept the current value
//	                ! of par. also test for the exceptional cases where parl
//	                ! is zero or the number of iterations has reached 10.
	//
//	                if (abs(fp) <= p1*Delta .or. parl == zero .and. fp <= temp .and. &
//	                    temp < zero .or. iter == 10) then
//	                    if (iter == 0) Par = zero
//	                    exit
//	                else
	//
//	                    ! compute the newton correction.
	//
//	                    do j = 1, n
//	                        l = Ipvt(j)
//	                        Wa1(j) = Diag(l)*(Wa2(l)/dxnorm)
//	                    end do
//	                    do j = 1, n
//	                        Wa1(j) = Wa1(j)/Sdiag(j)
//	                        temp = Wa1(j)
//	                        jp1 = j + 1
//	                        if (n >= jp1) then
//	                            do i = jp1, n
//	                                Wa1(i) = Wa1(i) - r(i, j)*temp
//	                            end do
//	                        end if
//	                    end do
//	                    temp = enorm(n, Wa1)
//	                    parc = ((fp/Delta)/temp)/temp
	//
//	                    ! depending on the sign of the function, update parl or paru.
	//
//	                    if (fp > zero) parl = max(parl, Par)
//	                    if (fp < zero) paru = min(paru, Par)
	//
//	                    ! compute an improved estimate for par.
	//
//	                    Par = max(parl, Par + parc)
	//
//	                end if
	//
//	            end do ! end of an iteration.
	//
//	        end if
	//
//	    end subroutine lmpar
	// !*****************************************************************************************
	//
	// !*****************************************************************************************
	// !>
	// ! the purpose of lmstr is to minimize the sum of the squares of
	// ! m nonlinear functions in n variables by a modification of
	// ! the levenberg-marquardt algorithm which uses minimal storage.
	// ! the user must provide a subroutine which calculates the
	// ! functions and the rows of the jacobian.
	//
//	    subroutine lmstr(fcn, m, n, x, Fvec, Fjac, Ldfjac, Ftol, Xtol, Gtol, Maxfev, &
//	                     Diag, Mode, Factor, Nprint, Info, Nfev, Njev, Ipvt, Qtf, &
//	                     Wa1, Wa2, Wa3, Wa4)
//	        implicit none
	//
//	        procedure(fcn_lmstr) :: fcn !! user-supplied subroutine which
//	                                    !! calculates the functions and the rows of the jacobian.
//	        integer, intent(in) :: m !! a positive integer input variable set to the number
//	                                !! of functions.
//	        integer, intent(in) :: n !! a positive integer input variable set to the number
//	                                !! of variables. n must not exceed m.
//	        integer, intent(in) :: Ldfjac !! a positive integer input variable not less than n
//	                                     !! which specifies the leading dimension of the array fjac.
//	        integer, intent(in) :: Maxfev !! a positive integer input variable. termination
//	                                     !! occurs when the number of calls to fcn with iflag = 1
//	                                     !! has reached maxfev.
//	        integer, intent(in) :: Mode !! an integer input variable. if mode = 1, the
//	                                   !! variables will be scaled internally. if mode = 2,
//	                                   !! the scaling is specified by the input diag. other
//	                                   !! values of mode are equivalent to mode = 1.
//	        integer, intent(in) :: Nprint !! an integer input variable that enables controlled
//	                                     !! printing of iterates if it is positive. in this case,
//	                                     !! fcn is called with iflag = 0 at the beginning of the first
//	                                     !! iteration and every nprint iterations thereafter and
//	                                     !! immediately prior to return, with x and fvec available
//	                                     !! for printing. if nprint is not positive, no special calls
//	                                     !! of fcn with iflag = 0 are made.
//	        integer, intent(out) :: Info !! an integer output variable. if the user has
//	                                    !! terminated execution, info is set to the (negative)
//	                                    !! value of iflag. see description of fcn. otherwise,
//	                                    !! info is set as follows:
//	                                    !!
//	                                    !!  * ***info = 0***  improper input parameters.
//	                                    !!  * ***info = 1***  both actual and predicted relative reductions
//	                                    !!    in the sum of squares are at most ftol.
//	                                    !!  * ***info = 2***  relative error between two consecutive iterates
//	                                    !!    is at most xtol.
//	                                    !!  * ***info = 3***  conditions for info = 1 and info = 2 both hold.
//	                                    !!  * ***info = 4***  the cosine of the angle between fvec and any
//	                                    !!    column of the jacobian is at most gtol in
//	                                    !!    absolute value.
//	                                    !!  * ***info = 5***  number of calls to fcn with iflag = 1 has
//	                                    !!    reached maxfev.
//	                                    !!  * ***info = 6***  ftol is too small. no further reduction in
//	                                    !!    the sum of squares is possible.
//	                                    !!  * ***info = 7***  xtol is too small. no further improvement in
//	                                    !!    the approximate solution x is possible.
//	                                    !!  * ***info = 8***  gtol is too small. fvec is orthogonal to the
//	                                    !!    columns of the jacobian to machine precision.
//	        integer, intent(out) :: Nfev !! an integer output variable set to the number of
//	                                    !! calls to fcn with iflag = 1.
//	        integer, intent(out) :: Njev !! an integer output variable set to the number of
//	                                    !! calls to fcn with iflag = 2.
//	        integer, intent(out) :: Ipvt(n) !! an integer output array of length n. ipvt
//	                                       !! defines a permutation matrix p such that jac*p = q*r,
//	                                       !! where jac is the final calculated jacobian, q is
//	                                       !! orthogonal (not stored), and r is upper triangular.
//	                                       !! column j of p is column ipvt(j) of the identity matrix.
//	        real(wp), intent(in) :: Ftol !! a nonnegative input variable. termination
//	                                    !! occurs when both the actual and predicted relative
//	                                    !! reductions in the sum of squares are at most ftol.
//	                                    !! therefore, ftol measures the relative error desired
//	                                    !! in the sum of squares.
//	        real(wp), intent(in) :: Xtol !! a nonnegative input variable. termination
//	                                    !! occurs when the relative error between two consecutive
//	                                    !! iterates is at most xtol. therefore, xtol measures the
//	                                    !! relative error desired in the approximate solution.
//	        real(wp), intent(in) :: Gtol !! a nonnegative input variable. termination
//	                                    !! occurs when the cosine of the angle between fvec and
//	                                    !! any column of the jacobian is at most gtol in absolute
//	                                    !! value. therefore, gtol measures the orthogonality
//	                                    !! desired between the function vector and the columns
//	                                    !! of the jacobian.
//	        real(wp), intent(in) :: Factor !! a positive input variable used in determining the
//	                                      !! initial step bound. this bound is set to the product of
//	                                      !! factor and the euclidean norm of diag*x if nonzero, or else
//	                                      !! to factor itself. in most cases factor should lie in the
//	                                      !! interval (.1,100.). 100. is a generally recommended value.
//	        real(wp), intent(inout) :: x(n) !! an array of length n. on input x must contain
//	                                       !! an initial estimate of the solution vector. on output x
//	                                       !! contains the final estimate of the solution vector.
//	        real(wp), intent(out) :: Fvec(m) !! an output array of length m which contains
//	                                        !! the functions evaluated at the output x.
//	        real(wp), intent(out) :: Fjac(Ldfjac, n) !! an output n by n array. the upper triangle of fjac
//	                                                !! contains an upper triangular matrix r such that
//	                                                !!```
//	                                                !!        t     t           t
//	                                                !!       p *(jac *jac)*p = r *r,
//	                                                !!```
//	                                                !! where p is a permutation matrix and jac is the final
//	                                                !! calculated jacobian. column j of p is column ipvt(j)
//	                                                !! (see below) of the identity matrix. the lower triangular
//	                                                !! part of fjac contains information generated during
//	                                                !! the computation of r.
//	        real(wp), intent(inout) :: Diag(n) !! an array of length n. if mode = 1 (see
//	                                          !! below), diag is internally set. if mode = 2, diag
//	                                          !! must contain positive entries that serve as
//	                                          !! multiplicative scale factors for the variables.
//	        real(wp), intent(out) :: Qtf(n) !! an output array of length n which contains
//	                                       !! the first n elements of the vector (q transpose)*fvec.
//	        real(wp), intent(inout) :: Wa1(n) !! work array of length n.
//	        real(wp), intent(inout) :: Wa2(n) !! work array of length n.
//	        real(wp), intent(inout) :: Wa3(n) !! work array of length n.
//	        real(wp), intent(inout) :: Wa4(m) !! work array of length m.
	//
		public static void lmstr() {
			
		}
//	        integer :: i, iflag, iter, j, l
//	        real(wp) :: actred, delta, dirder, fnorm, &
//	                    fnorm1, gnorm, par, pnorm, prered, &
//	                    ratio, sum, temp, temp1, temp2, xnorm
//	        logical :: sing
	//
//	        real(wp), parameter :: p1 = 1.0e-1_wp
//	        real(wp), parameter :: p5 = 5.0e-1_wp
//	        real(wp), parameter :: p25 = 2.5e-1_wp
//	        real(wp), parameter :: p75 = 7.5e-1_wp
//	        real(wp), parameter :: p0001 = 1.0e-4_wp
	//
//	        Info = 0
//	        iflag = 0
//	        Nfev = 0
//	        Njev = 0
	//
//	        main : block
	//
//	            ! check the input parameters for errors.
	//
//	            if (n <= 0 .or. m < n .or. Ldfjac < n .or. Ftol < zero .or. &
//	                Xtol < zero .or. Gtol < zero .or. Maxfev <= 0 .or. Factor <= zero) &
//	                exit main
//	            if (Mode == 2) then
//	                do j = 1, n
//	                    if (Diag(j) <= zero) exit main
//	                end do
//	            end if
	//
//	            ! evaluate the function at the starting point
//	            ! and calculate its norm.
	//
//	            iflag = 1
//	            call fcn(m, n, x, Fvec, Wa3, iflag)
//	            Nfev = 1
//	            if (iflag < 0) exit main
//	            fnorm = enorm(m, Fvec)
	//
//	            ! initialize levenberg-marquardt parameter and iteration counter.
	//
//	            par = zero
//	            iter = 1
	//
//	            ! beginning of the outer loop.
	//
//	            outer : do
	//
//	                ! if requested, call fcn to enable printing of iterates.
	//
//	                if (Nprint > 0) then
//	                    iflag = 0
//	                    if (mod(iter - 1, Nprint) == 0) call fcn(m, n, x, Fvec, Wa3, iflag)
//	                    if (iflag < 0) exit main
//	                end if
	//
//	                ! compute the qr factorization of the jacobian matrix
//	                ! calculated one row at a time, while simultaneously
//	                ! forming (q transpose)*fvec and storing the first
//	                ! n components in qtf.
	//
//	                do j = 1, n
//	                    Qtf(j) = zero
//	                    do i = 1, n
//	                        Fjac(i, j) = zero
//	                    end do
//	                end do
//	                iflag = 2
//	                do i = 1, m
//	                    call fcn(m, n, x, Fvec, Wa3, iflag)
//	                    if (iflag < 0) exit main
//	                    temp = Fvec(i)
//	                    call rwupdt(n, Fjac, Ldfjac, Wa3, Qtf, temp, Wa1, Wa2)
//	                    iflag = iflag + 1
//	                end do
//	                Njev = Njev + 1
	//
//	                ! if the jacobian is rank deficient, call qrfac to
//	                ! reorder its columns and update the components of qtf.
	//
//	                sing = .false.
//	                do j = 1, n
//	                    if (Fjac(j, j) == zero) sing = .true.
//	                    Ipvt(j) = j
//	                    Wa2(j) = enorm(j, Fjac(1, j))
//	                end do
//	                if (sing) then
//	                    call qrfac(n, n, Fjac, Ldfjac, .true., Ipvt, n, Wa1, Wa2, Wa3)
//	                    do j = 1, n
//	                        if (Fjac(j, j) /= zero) then
//	                            sum = zero
//	                            do i = j, n
//	                                sum = sum + Fjac(i, j)*Qtf(i)
//	                            end do
//	                            temp = -sum/Fjac(j, j)
//	                            do i = j, n
//	                                Qtf(i) = Qtf(i) + Fjac(i, j)*temp
//	                            end do
//	                        end if
//	                        Fjac(j, j) = Wa1(j)
//	                    end do
//	                end if
	//
//	                ! on the first iteration and if mode is 1, scale according
//	                ! to the norms of the columns of the initial jacobian.
	//
//	                if (iter == 1) then
//	                    if (Mode /= 2) then
//	                        do j = 1, n
//	                            Diag(j) = Wa2(j)
//	                            if (Wa2(j) == zero) Diag(j) = one
//	                        end do
//	                    end if
	//
//	                    ! on the first iteration, calculate the norm of the scaled x
//	                    ! and initialize the step bound delta.
	//
//	                    do j = 1, n
//	                        Wa3(j) = Diag(j)*x(j)
//	                    end do
//	                    xnorm = enorm(n, Wa3)
//	                    delta = Factor*xnorm
//	                    if (delta == zero) delta = Factor
//	                end if
	//
//	                ! compute the norm of the scaled gradient.
	//
//	                gnorm = zero
//	                if (fnorm /= zero) then
//	                    do j = 1, n
//	                        l = Ipvt(j)
//	                        if (Wa2(l) /= zero) then
//	                            sum = zero
//	                            do i = 1, j
//	                                sum = sum + Fjac(i, j)*(Qtf(i)/fnorm)
//	                            end do
//	                            gnorm = max(gnorm, abs(sum/Wa2(l)))
//	                        end if
//	                    end do
//	                end if
	//
//	                ! test for convergence of the gradient norm.
	//
//	                if (gnorm <= Gtol) Info = 4
//	                if (Info /= 0) exit main
	//
//	                ! rescale if necessary.
	//
//	                if (Mode /= 2) then
//	                    do j = 1, n
//	                        Diag(j) = max(Diag(j), Wa2(j))
//	                    end do
//	                end if
	//
//	                ! beginning of the inner loop.
	//
//	                inner : do
	//
//	                    ! determine the levenberg-marquardt parameter.
	//
//	                    call lmpar(n, Fjac, Ldfjac, Ipvt, Diag, Qtf, delta, par, Wa1, Wa2, Wa3, Wa4)
	//
//	                    ! store the direction p and x + p. calculate the norm of p.
	//
//	                    do j = 1, n
//	                        Wa1(j) = -Wa1(j)
//	                        Wa2(j) = x(j) + Wa1(j)
//	                        Wa3(j) = Diag(j)*Wa1(j)
//	                    end do
//	                    pnorm = enorm(n, Wa3)
	//
//	                    ! on the first iteration, adjust the initial step bound.
	//
//	                    if (iter == 1) delta = min(delta, pnorm)
	//
//	                    ! evaluate the function at x + p and calculate its norm.
	//
//	                    iflag = 1
//	                    call fcn(m, n, Wa2, Wa4, Wa3, iflag)
//	                    Nfev = Nfev + 1
//	                    if (iflag < 0) exit main
	//
//	                    fnorm1 = enorm(m, Wa4)
	//
//	                    ! compute the scaled actual reduction.
	//
//	                    actred = -one
//	                    if (p1*fnorm1 < fnorm) actred = one - (fnorm1/fnorm)**2
	//
//	                    ! compute the scaled predicted reduction and
//	                    ! the scaled directional derivative.
	//
//	                    do j = 1, n
//	                        Wa3(j) = zero
//	                        l = Ipvt(j)
//	                        temp = Wa1(l)
//	                        do i = 1, j
//	                            Wa3(i) = Wa3(i) + Fjac(i, j)*temp
//	                        end do
//	                    end do
//	                    temp1 = enorm(n, Wa3)/fnorm
//	                    temp2 = (sqrt(par)*pnorm)/fnorm
//	                    prered = temp1**2 + temp2**2/p5
//	                    dirder = -(temp1**2 + temp2**2)
	//
//	                    ! compute the ratio of the actual to the predicted
//	                    ! reduction.
	//
//	                    ratio = zero
//	                    if (prered /= zero) ratio = actred/prered
	//
//	                    ! update the step bound.
	//
//	                    if (ratio <= p25) then
//	                        if (actred >= zero) temp = p5
//	                        if (actred < zero) temp = p5*dirder/(dirder + p5*actred)
//	                        if (p1*fnorm1 >= fnorm .or. temp < p1) temp = p1
//	                        delta = temp*min(delta, pnorm/p1)
//	                        par = par/temp
//	                    elseif (par == zero .or. ratio >= p75) then
//	                        delta = pnorm/p5
//	                        par = p5*par
//	                    end if
	//
//	                    ! test for successful iteration.
	//
//	                    if (ratio >= p0001) then
	//
//	                        ! successful iteration. update x, fvec, and their norms.
	//
//	                        do j = 1, n
//	                            x(j) = Wa2(j)
//	                            Wa2(j) = Diag(j)*x(j)
//	                        end do
//	                        do i = 1, m
//	                            Fvec(i) = Wa4(i)
//	                        end do
//	                        xnorm = enorm(n, Wa2)
//	                        fnorm = fnorm1
//	                        iter = iter + 1
//	                    end if
	//
//	                    ! tests for convergence.
	//
//	                    if (abs(actred) <= Ftol .and. prered <= Ftol .and. &
//	                        p5*ratio <= one) Info = 1
//	                    if (delta <= Xtol*xnorm) Info = 2
//	                    if (abs(actred) <= Ftol .and. prered <= Ftol .and. &
//	                        p5*ratio <= one .and. Info == 2) Info = 3
//	                    if (Info /= 0) exit main
	//
//	                    ! tests for termination and stringent tolerances.
	//
//	                    if (Nfev >= Maxfev) Info = 5
//	                    if (abs(actred) <= epsmch .and. prered <= epsmch .and. &
//	                        p5*ratio <= one) Info = 6
//	                    if (delta <= epsmch*xnorm) Info = 7
//	                    if (gnorm <= epsmch) Info = 8
//	                    if (Info /= 0) exit main
	//
//	                    if (ratio >= p0001) exit inner
	//
//	                end do inner ! end of the inner loop. repeat if iteration unsuccessful.
	//
//	            end do outer ! end of the outer loop.
	//
//	        end block main
	//
//	        ! termination, either normal or user imposed.
	//
//	        if (iflag < 0) Info = iflag
//	        iflag = 0
//	        if (Nprint > 0) call fcn(m, n, x, Fvec, Wa3, iflag)
	//
//	    end subroutine lmstr
	// !*****************************************************************************************
	//
	// !*****************************************************************************************
	// !>
	// ! the purpose of lmstr1 is to minimize the sum of the squares of
	// ! m nonlinear functions in n variables by a modification of
	// ! the levenberg-marquardt algorithm which uses minimal storage.
	// ! this is done by using the more general least-squares solver
	// ! lmstr. the user must provide a subroutine which calculates
	// ! the functions and the rows of the jacobian.
	//
//	    subroutine lmstr1(fcn, m, n, x, Fvec, Fjac, Ldfjac, Tol, Info, Ipvt, Wa, Lwa)
//	        implicit none
	//
//	        procedure(fcn_lmstr) :: fcn !! user-supplied subroutine which
//	                                    !! calculates the functions and the rows of the jacobian.
//	        integer, intent(in) :: m !! a positive integer input variable set to the number
//	                                !! of functions.
//	        integer, intent(in) :: n !! a positive integer input variable set to the number
//	                                !! of variables. n must not exceed m.
//	        integer, intent(in) :: Ldfjac !! a positive integer input variable not less than n
//	                                     !! which specifies the leading dimension of the array fjac.
//	        integer, intent(out) :: Info !! an integer output variable. if the user has
//	                                    !! terminated execution, info is set to the (negative)
//	                                    !! value of iflag. see description of fcn. otherwise,
//	                                    !! info is set as follows:
//	                                    !!
//	                                    !!  * ***info = 0***  improper input parameters.
//	                                    !!  * ***info = 1***  algorithm estimates that the relative error
//	                                    !!           in the sum of squares is at most tol.
//	                                    !!  * ***info = 2***  algorithm estimates that the relative error
//	                                    !!           between x and the solution is at most tol.
//	                                    !!  * ***info = 3***  conditions for info = 1 and info = 2 both hold.
//	                                    !!  * ***info = 4***  fvec is orthogonal to the columns of the
//	                                    !!           jacobian to machine precision.
//	                                    !!  * ***info = 5***  number of calls to fcn with iflag = 1 has
//	                                    !!           reached 100*(n+1).
//	                                    !!  * ***info = 6***  tol is too small. no further reduction in
//	                                    !!           the sum of squares is possible.
//	                                    !!  * ***info = 7***  tol is too small. no further improvement in
//	                                    !!           the approximate solution x is possible.
//	        integer, intent(in) :: Lwa !! a positive integer input variable not less than 5*n+m.
//	        integer, intent(out) :: Ipvt(n) !! an integer output array of length n. ipvt
//	                                       !! defines a permutation matrix p such that jac*p = q*r,
//	                                       !! where jac is the final calculated jacobian, q is
//	                                       !! orthogonal (not stored), and r is upper triangular.
//	                                       !! column j of p is column ipvt(j) of the identity matrix.
//	        real(wp), intent(in) :: Tol !! a nonnegative input variable. termination occurs
//	                                   !! when the algorithm estimates either that the relative
//	                                   !! error in the sum of squares is at most tol or that
//	                                   !! the relative error between x and the solution is at
//	                                   !! most tol.
//	        real(wp), intent(inout) :: x(n) !! an array of length n. on input x must contain
//	                                       !! an initial estimate of the solution vector. on output x
//	                                       !! contains the final estimate of the solution vector.
//	        real(wp), intent(out) :: Fvec(m) !! an output array of length m which contains
//	                                        !! the functions evaluated at the output x.
//	        real(wp), intent(out) :: Fjac(Ldfjac, n) !! an output n by n array. the upper triangle of fjac
//	                                                !! contains an upper triangular matrix r such that
//	                                                !!```
//	                                                !!        t     t           t
//	                                                !!       p *(jac *jac)*p = r *r,
//	                                                !!```
//	                                                !! where p is a permutation matrix and jac is the final
//	                                                !! calculated jacobian. column j of p is column ipvt(j)
//	                                                !! (see below) of the identity matrix. the lower triangular
//	                                                !! part of fjac contains information generated during
//	                                                !! the computation of r.
//	        real(wp), intent(inout) :: Wa(Lwa) !! a work array of length lwa.
	//
		public static void lmstr1() {
			
		}
//	        integer :: maxfev, mode, nfev, njev, nprint
//	        real(wp) :: ftol, gtol, xtol
	//
//	        real(wp), parameter :: factor = 1.0e2_wp
	//
//	        Info = 0
	//
//	        ! check the input parameters for errors.
	//
//	        if (n > 0 .and. m >= n .and. Ldfjac >= n .and. Tol >= zero .and. &
//	            Lwa >= 5*n + m) then
	//
//	            ! call lmstr.
	//
//	            maxfev = 100*(n + 1)
//	            ftol = Tol
//	            xtol = Tol
//	            gtol = zero
//	            mode = 1
//	            nprint = 0
//	            call lmstr(fcn, m, n, x, Fvec, Fjac, Ldfjac, ftol, xtol, gtol, maxfev, &
//	                       Wa(1), mode, factor, nprint, Info, nfev, njev, Ipvt, Wa(n + 1), &
//	                       Wa(2*n + 1), Wa(3*n + 1), Wa(4*n + 1), Wa(5*n + 1))
//	            if (Info == 8) Info = 4
//	        end if
	//
//	    end subroutine lmstr1
	// !*****************************************************************************************

	/**
	 * <p>
	 * this subroutine proceeds from the computed qr factorization ofan m by n
	 * matrix a to accumulate the m by m orthogonal matrix q from its factored form.
	 * 
	 * @param m   int -> A positive integer variable set to the number of rows of a
	 *            and the order of q.
	 *            <p>
	 * @param n   int -> A positive integer variable set to the number of columns of
	 *            a.
	 *            <p>
	 * @param ldq int -> A positive integer variable not less than m which specifies
	 *            the leading dimension of the array q.
	 *            <p>
	 * @param q   double[ldq][m] -> An m by m array. On input the full lower
	 *            trapezoid in the first min(m,n) columns of q contains the factored
	 *            form. On output q has been accumulated into a square matrix.
	 *            <p>
	 * @param wa  double[m] -> A working array of length m.
	 * 
	 * @return Nothing is returned, rather the objects themselves are modified.
	 */
	
	public static void qform(int m, int n, double[][] q, int ldq, double[] wa) {
		
		int i, j, jm1, k, l, minmn, np1; //temporary variables 
		double sum, temp; //temporary variables 
//	    zero out upper triangle of q in the first min(m,n) columns.
		if (m < n) {
			minmn = m;
		} else {
			minmn = n;
		}
		if (minmn >= 2) {
			for (j = 1; j < minmn; j++) {
				jm1 = j;
				for (i = 0; i < jm1; i++) {
					q[i][j] = 0;
				}
			}
		}
//	    initialize remaining columns to those of the identity matrix.
		np1 = n + 1;
		if (m >= np1) { //this portion has not been tested to ensure it works : USE CAUTION:
			for (j = np1-1; j < m; j++) {
				for (i = 0; i < m; i++) {
					q[i][j] = 0;
				}
				q[j][j] = 1;
			}
		}
//	    accumulate q from its factored form.
		for (l = 0; l < minmn; l++) {
			k = minmn - (l+1);
			for (i = k; i < m; i++) {
				wa[i] = q[i][k];
				q[i][k] = 0;
			}
			q[k][k] = 1;
			if (wa[k] != 0) {
				for (j = k; j < m; j++) {
					sum = 0;
					for (i = k; i < m; i++) {
						sum = sum + q[i][j] * wa[i];
					}
					temp = sum / wa[k ];
					for (i = k; i < m; i++) {
						q[i][j] = q[i][j] - temp * wa[i];

					}
				}
			}
		}
	}
	
/**
 * <p>
 * This subroutine uses householder transformations with column pivoting
 * (optional) to compute a qr factorization of the m by n matrix a. That is,
 * qrfac determines an orthogonal matrix q, a permutation matrix p, and an upper
 * trapezoidal matrix r with diagonal elements of nonincreasing magnitude, such
 * that a*p = q*r. The householder transformation for column k, k =
 * 1,2,...,min(m,n), is of the form {t i - (1/u(k))*u*u} where u has zeros in
 * the first k-1 positions. The form of this transformation and the method of
 * pivoting first appeared in the corresponding linpack subroutine.
 * 
 * 
 * 
 * @param m      int -> A positive integer variable set to the number of rows of
 *               a.
 *               <p>
 * @param n      int -> A positive integer variable set to the number of columns
 *               of a.
 *               <p>
 * @param lda    int -> A positive integer variable not less than m which
 *               specifies the leading dimension of the array a.
 *               <p>
 * @param lipvt  int -> A positive integer variable. If pivot is false, then
 *               lipvt may be as small as 1. If pivot is true, then lipvt must
 *               be at least n.
 *               <p>
 * @param ipvt   int[lipvt] -> An integer array of length lipvt. ipvt defines
 *               the permutation matrix p such that a*p = q*r. Column j of p is
 *               column ipvt(j) of the identity matrix. If pivot is false, ipvt
 *               is not referenced.
 *               <p>
 * @param pivot  boolean -> A logical variable. If pivot is set true, then
 *               column pivoting is enforced. If pivot is set false, then no
 *               column pivoting is done.
 *               <p>
 * @param a      double[lda][n] -> A m by n array. On input a contains the
 *               matrix for which the qr factorization is to be computed. On
 *               output the strict upper trapezoidal part of a contains the
 *               strict upper trapezoidal part of r, and the lower trapezoidal
 *               part of a contains a factored form of q (the non-trivial
 *               elements of the u vectors described above).
 *               <p>
 * @param rDiag  double[n] -> An array of length n which contains the diagonal
 *               elements of r.
 *               <p>
 * @param acNorm double[n] -> An array of length n which contains the norms of
 *               the corresponding columns of the input matrix a. If this
 *               information is not needed, then acnorm can coincide with rdiag.
 *               <p>
 * @param wa     double[n] -> A working array of length n. If pivot is false,
 *               then wa can coincide with rdiag.
 *               <p>
 * @return Nothing is returned, rather the objects themselves are modified.
 * 
 * 
 */				
	

	public static void qrfac(int m, int n, double[][] a, int lda, boolean pivot, int[] ipvt, int lipvt,
			double[] rdiag, double[] acnorm, double[] wa) {
				
		int i, j, jp1, k, kmax, minmn, count; //Temporary variables
		double ajnorm, sum, temp; //Temporary variables

		double[] subArray = new double[n]; //Used for sending sub arrays to enorm for evaluation
		double p05 = 0.05; // 0.05, thats it. 
//	    Compute the initial column norms and initialize several arrays.
		for (j = 0; j < n; j++) {
			subArray = new double[n];
			count = 0;
			for (int z = 0; z < n; z++) {
				subArray[count] = a[z][j]; // Sub Array of a column
				count++;
			}
			acnorm[j] = enorm(m, subArray);
			rdiag[j] = acnorm[j];
			wa[j] = rdiag[j];
			if (pivot) {
				ipvt[j] = j;
			}
		}
//		Reduce a to r with householder transformations.
		if (m < n) {
			minmn = m;
		} else {
			minmn = n;
		}
		for (j = 0; j < minmn; j++) {
			if (pivot) {
				kmax = j;
				for (k = j; k < n; k++) {
					if (rdiag[k] > rdiag[kmax]) {
						kmax = k;
					}
				}
				if (kmax != j) {
					for (i = 0; i < m; i++) {
						temp = a[i][j];
						a[i][j] = a[i][kmax];
						a[i][kmax] = temp;
					}
					rdiag[kmax] = rdiag[j];
					wa[kmax] = wa[j];
					k = ipvt[j];
					ipvt[j] = ipvt[kmax];
					ipvt[kmax] = k;
				}
			}
//	        Compute the householder transformation to reduce the
//	        j-th column of a to a multiple of the j-th unit vector.
			subArray = new double[n];
			count = 0;
			for (int z = j; z < n; z++) {
				subArray[count] = a[z][j]; //As the algorithm moves diagonally down the matrix,
				count++;				//a smaller and smaller subarray is sent to enorm
			}
			ajnorm = enorm(m - j, subArray);
			if (ajnorm != 0) {
				if (a[j][j] < 0) {
					ajnorm = -ajnorm;
				}
				for (i = j; i < m; i++) {
					a[i][j] = a[i][j] / ajnorm;
				}
				a[j][j] = a[j][j] + 1.0;
//				Apply the transformation to the remaining columns
//				and update the norms.
				jp1 = j + 1;
				if (n > jp1) {
					for (k = jp1; k < n; k++) {
						sum = 0;
						for (i = j; i < m; i++) {
							sum = sum + a[i][j] * a[i][k];
						}
						temp = sum / a[j][j];
						for (i = j; i < m ; i++) {
							a[i][k] = a[i][k] - temp * a[i][j];
						}
						if (!(!pivot || rdiag[k] == 0)) {
							temp = a[j][k] / rdiag[k];
							if (0 < 1 - (temp * temp)) {
								rdiag[k] = rdiag[k] * Math.sqrt(1 - (temp * temp));
							} else {
								rdiag[k] = 0;
							}
							if (p05 * ((rdiag[k] / wa[k]) * (rdiag[k] / wa[k])) <= epsmch) {
								subArray = new double[n];
								count = 0;
								for (int z = jp1; z < n; z++) {
									subArray[count] = a[z][k]; 
									count++;
								}
								ajnorm = enorm(m - j, subArray);
								rdiag[k] = ajnorm;
								wa[k] = rdiag[k];
							}
						}
					}
				}
			}
			rdiag[j] = -ajnorm;
		}
	}

	// !*****************************************************************************************
	// !>
	// ! given an m by n matrix a, an n by n diagonal matrix d,
	// ! and an m-vector b, the problem is to determine an x which
	// ! solves the system
	// !```
	// ! a*x = b , d*x = 0 ,
	// !```
	// ! in the least squares sense.
	// !
	// ! this subroutine completes the solution of the problem
	// ! if it is provided with the necessary information from the
	// ! qr factorization, with column pivoting, of a. that is, if
	// ! a*p = q*r, where p is a permutation matrix, q has orthogonal
	// ! columns, and r is an upper triangular matrix with diagonal
	// ! elements of nonincreasing magnitude, then qrsolv expects
	// ! the full upper triangle of r, the permutation matrix p,
	// ! and the first n components of (q transpose)*b. the system
	// ! a*x = b, d*x = 0, is then equivalent to
	// !```
	// ! t t
	// ! r*z = q *b , p *d*p*z = 0 ,
	// !```
	// ! where x = p*z. if this system does not have full rank,
	// ! then a least squares solution is obtained. on output qrsolv
	// ! also provides an upper triangular matrix s such that
	// !```
	// ! t t t
	// ! p *(a *a + d*d)*p = s *s .
	// !```
	// ! s is computed within qrsolv and may be of separate interest.
	//
//	    subroutine qrsolv(n, r, Ldr, Ipvt, Diag, Qtb, x, Sdiag, Wa)
//	        implicit none
	//
//	        integer, intent(in) :: n !! a positive integer input variable set to the order of r.
//	        integer, intent(in) :: Ldr !! a positive integer input variable not less than n
//	                                  !! which specifies the leading dimension of the array r.
//	        integer, intent(in) :: Ipvt(n) !! an integer input array of length n which defines the
//	                                      !! permutation matrix p such that a*p = q*r. column j of p
//	                                      !! is column ipvt(j) of the identity matrix.
//	        real(wp), intent(inout) :: r(Ldr, n) !! an n by n array. on input the full upper triangle
//	                                            !! must contain the full upper triangle of the matrix r.
//	                                            !! on output the full upper triangle is unaltered, and the
//	                                            !! strict lower triangle contains the strict upper triangle
//	                                            !! (transposed) of the upper triangular matrix s.
//	        real(wp), intent(in) :: Diag(n) !! an input array of length n which must contain the
//	                                       !! diagonal elements of the matrix d.
//	        real(wp), intent(in) :: Qtb(n) !! an input array of length n which must contain the first
//	                                      !! n elements of the vector (q transpose)*b.
//	        real(wp), intent(out) :: x(n) !! an output array of length n which contains the least
//	                                     !! squares solution of the system a*x = b, d*x = 0.
//	        real(wp), intent(out) :: Sdiag(n) !! an output array of length n which contains the
//	                                         !! diagonal elements of the upper triangular matrix s.
//	        real(wp), intent(inout) :: Wa(n) !! a work array of length n.
	//
	public static void qrsolv() {
		
	}
//	        integer :: i, j, jp1, k, kp1, l, nsing
//	        real(wp) :: cos, cotan, qtbpj, sin, sum, tan, temp
	//
//	        real(wp), parameter :: p5 = 5.0e-1_wp
//	        real(wp), parameter :: p25 = 2.5e-1_wp
	//
//	        ! copy r and (q transpose)*b to preserve input and initialize s.
//	        ! in particular, save the diagonal elements of r in x.
	//
//	        do j = 1, n
//	            do i = j, n
//	                r(i, j) = r(j, i)
//	            end do
//	            x(j) = r(j, j)
//	            Wa(j) = Qtb(j)
//	        end do
	//
//	        ! eliminate the diagonal matrix d using a givens rotation.
	//
//	        do j = 1, n
	//
//	            ! prepare the row of d to be eliminated, locating the
//	            ! diagonal element using p from the qr factorization.
	//
//	            l = Ipvt(j)
//	            if (Diag(l) /= zero) then
//	                do k = j, n
//	                    Sdiag(k) = zero
//	                end do
//	                Sdiag(j) = Diag(l)
	//
//	                ! the transformations to eliminate the row of d
//	                ! modify only a single element of (q transpose)*b
//	                ! beyond the first n, which is initially zero.
	//
//	                qtbpj = zero
//	                do k = j, n
	//
//	                    ! determine a givens rotation which eliminates the
//	                    ! appropriate element in the current row of d.
	//
//	                    if (Sdiag(k) /= zero) then
//	                        if (abs(r(k, k)) >= abs(Sdiag(k))) then
//	                            tan = Sdiag(k)/r(k, k)
//	                            cos = p5/sqrt(p25 + p25*tan**2)
//	                            sin = cos*tan
//	                        else
//	                            cotan = r(k, k)/Sdiag(k)
//	                            sin = p5/sqrt(p25 + p25*cotan**2)
//	                            cos = sin*cotan
//	                        end if
	//
//	                        ! compute the modified diagonal element of r and
//	                        ! the modified element of ((q transpose)*b,0).
	//
//	                        r(k, k) = cos*r(k, k) + sin*Sdiag(k)
//	                        temp = cos*Wa(k) + sin*qtbpj
//	                        qtbpj = -sin*Wa(k) + cos*qtbpj
//	                        Wa(k) = temp
	//
//	                        ! accumulate the tranformation in the row of s.
	//
//	                        kp1 = k + 1
//	                        if (n >= kp1) then
//	                            do i = kp1, n
//	                                temp = cos*r(i, k) + sin*Sdiag(i)
//	                                Sdiag(i) = -sin*r(i, k) + cos*Sdiag(i)
//	                                r(i, k) = temp
//	                            end do
//	                        end if
//	                    end if
//	                end do
//	            end if
	//
//	            ! store the diagonal element of s and restore
//	            ! the corresponding diagonal element of r.
	//
//	            Sdiag(j) = r(j, j)
//	            r(j, j) = x(j)
//	        end do
	//
//	        ! solve the triangular system for z. if the system is
//	        ! singular, then obtain a least squares solution.
	//
//	        nsing = n
//	        do j = 1, n
//	            if (Sdiag(j) == zero .and. nsing == n) nsing = j - 1
//	            if (nsing < n) Wa(j) = zero
//	        end do
//	        if (nsing >= 1) then
//	            do k = 1, nsing
//	                j = nsing - k + 1
//	                sum = zero
//	                jp1 = j + 1
//	                if (nsing >= jp1) then
//	                    do i = jp1, nsing
//	                        sum = sum + r(i, j)*Wa(i)
//	                    end do
//	                end if
//	                Wa(j) = (Wa(j) - sum)/Sdiag(j)
//	            end do
//	        end if
	//
//	        ! permute the components of z back to components of x.
	//
//	        do j = 1, n
//	            l = Ipvt(j)
//	            x(l) = Wa(j)
//	        end do
	//
//	    end subroutine qrsolv
	// !*****************************************************************************************
	
	/**
	 * 
	 * Given an m by n matrix a, this subroutine computes a*q whereq is the product
	 * of 2*(n - 1) transformations { gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1) } and
	 * gv(i), gw(i) are givens rotations in the (i,n) plane which eliminate elements
	 * in the i-th and n-th planes, respectively. q itself is not given, rather the
	 * information to recover the gv, gw rotations is supplied.
	 * 
	 * @param m   int -> A positive integer variable set to the number of rows of a.
	 *            <p>
	 * @param n   int -> A positive integer variable set to the number of columns of
	 *            a.
	 *            <p>
	 * @param lda int -> A positive integer variable not less than m which specifies
	 *            the leading dimension of the array a.
	 *            <p>
	 * @param a   double[lda][n] -> A m by n array. On input a must contain the
	 *            matrix to be postmultiplied by the orthogonal matrix q described
	 *            above. On output a*q has replaced a.
	 *            <p>
	 * @param v   double[n] -> An array of length n. v(i) must contain the
	 *            information necessary to recover the givens rotation gv(i)
	 *            described above.
	 *            <p>
	 * @param w   double[n] -> An array of length n. w(i) must contain the
	 *            information necessary to recover the givens rotation gw(i)
	 *            described above.
	 *            <p>
	 * @return A double[][] array
	 */
	
	public static double[][] r1mpyq(int m, int n, double[][] a, int lda, double[] v, double[] w) {
		
		int i, j, nmj, nm1; //Temporary variables
		double cos = 0, sin = 0, temp; //Temporary variables

//	    Apply the first set of givens rotations to a.
		nm1 = n;
		if (nm1 >= 0) {
			for (nmj = 0; nmj < nm1-1; nmj++) {
				j = n - (nmj+2);
				if (Math.abs(v[j]) > 1) {
					cos = 1 / v[j];
					sin = Math.sqrt(1 - (cos * cos));
				} else {
					sin = v[j];
					cos = Math.sqrt(1 - (sin * sin));
				}
				for (i = 0; i < m; i++) {
					temp = cos * a[i][j] - sin * a[i][n-1];
					a[i][n-1] = sin * a[i][j] + cos * a[i][n-1];
					a[i][j] = temp;
				}
			}
//			Apply the second set of givens rotations to a.
			for (j = 0; j < nm1-1; j++) {
				if (Math.abs(w[j]) > 1) {
					cos = 1 / w[j];
				}
				if (Math.abs(w[j]) > 1) {
					sin = Math.sqrt(1 - (cos * cos));
				}
				if (Math.abs(w[j]) <= 1) {
					sin = w[j];
				}
				if (Math.abs(w[j]) <= 1) {
					cos = Math.sqrt(1 - (sin * sin));
				}
				for (i = 0; i < m; i++) {
					temp = cos * a[i][j] + sin * a[i][n-1];
					a[i][n-1] = -sin * a[i][j] + cos * a[i][n-1];
					a[i][j] = temp;
				}
			}
		}
		return a;
	}

	/**
	 * <p>
	 * Given an m by n lower trapezoidal matrix s, an m-vector u, and an n-vector v,
	 * the problem is to determine an orthogonal matrix q such that { t (s + u*v )*q
	 * } is again lower trapezoidal.
	 * <p>
	 * This subroutine determines q as the product of 2*(n - 1) transformations {
	 * gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1) } where gv(i), gw(i) are givens rotations
	 * in the (i,n) plane which eliminate elements in the i-th and n-th planes,
	 * respectively. q itself is not accumulated, rather the information to recover
	 * the gv, gw rotations is returned.
	 * 
	 * @param m    int -> A positive integer variable set to the number of rows of
	 *             s.
	 *             <p>
	 * @param n    int -> A positive integer variable set to the number of columns
	 *             of s. n must not exceed m.
	 *             <p>
	 * @param ls   int -> A positive integer variable not less than (n*(2*m-n+1))/2.
	 *             <p>
	 * @param sing boolean -> A logical variable. sing is set true if any of the
	 *             diagonal elements of the output s are zero. Otherwise sing is set
	 *             false.
	 *             <p>
	 * @param s    double[ls] -> An array of length ls. On input s must contain the
	 *             lower trapezoidal matrix s stored by columns. On output s
	 *             contains the lower trapezoidal matrix produced as described
	 *             above.
	 *             <p>
	 * @param u    double[m] -> An array of length m which must contain the vector
	 *             u.
	 *             <p>
	 * @param v    double[n] -> An array of length n. On input v must contain the
	 *             vector v. On output v(i) contains the information necessary to
	 *             recover the givens rotation gv(i) described above.
	 *             <p>
	 * @param w    double[m] -> An array of length m. w(i) contains information
	 *             necessary to recover the givens rotation gw(i) described above.
	 *             <p>
	 * @return Nothing is returned, rather the objects themselves are modified.
	 */
	
	public static void r1updt(int m, int n, double[] s, int Ls, double[] u, double[] v, double[] w,
			boolean sing) {
		
		//Being left as is with the array starting at 1 notation fixes.
		
		int i, j, jj, l, nmj, nm1; //Temporary variables
		double cos, cotan, sin, tan, tau, temp; //Temporary variables
		 
		double p5 = 0.5;
		double p25 = 0.25;
		double giant = 1.7976931348623157 * Math.pow(10, 308); //Largest magnitude possible
		
//	    Initialize the diagonal element pointer.
		jj = (n * (2 * m - n + 1)) / 2 - (m - n);
//	    Move the nontrivial part of the last column of s into w.
		l = jj;
		for (i = n; i < m + 1; i++) {
			w[i - 1] = s[l - 1];
			l = l + 1;
		}
//	    Rotate the vector v into a multiple of the n-th unit vector
//	    in such a way that a spike is introduced into w.
		nm1 = n - 1;
		if (nm1 >= 1) {
			for (nmj = 1; nmj < nm1 + 1; nmj++) {
				j = n - nmj;
				jj = jj - (m - j + 1);
				w[j - 1] = 0;
				if (v[j - 1] != 0) {
//	        	    Determine a givens rotation which eliminates the
//	        	    j-th element of v.
					if (Math.abs(v[n - 1]) >= Math.abs(v[j - 1])) {
						tan = v[j - 1] / v[n - 1];
						cos = p5 / Math.sqrt(p25 + p25 * (tan * tan));
						sin = cos * tan;
						tau = sin;
					} else {
						cotan = v[n - 1] / v[j - 1];
						sin = p5 / Math.sqrt(p25 + p25 * (cotan * cotan));
						cos = sin * cotan;
						tau = 1;
						if (Math.abs(cos) * giant > 1) {
							tau = 1 / cos;
						}
					}
//	        	    Apply the transformation to v and store the information
//	        	    necessary to recover the givens rotation.
					v[n - 1] = sin * v[j - 1] + cos * v[n - 1];
					v[j - 1] = tau;
//					Apply the transformation to s and extend the spike in w.
					l = jj;
					for (i = j; i < m + 1; i++) {
						temp = cos * s[l - 1] - sin * w[i - 1];
						w[i - 1] = sin * s[l - 1] + cos * w[i - 1];
						s[l - 1] = temp;
						l = l + 1;
					}
				}
			}
		}
//	    Add the spike from the rank 1 update to w.
		for (i = 0; i < m; i++) {
			w[i] = w[i] + v[n - 1] * u[i];
		}
//	    Eliminate the spike.
		sing = false;
		if (nm1 >= 1) {
			for (j = 1; j < nm1 + 1; j++) {
				if (w[j - 1] != 0) {
//					Determine a givens rotation which eliminates the
//					j-th element of the spike.
					if (Math.abs(s[jj - 1]) >= Math.abs(w[j - 1])) {
						tan = w[j - 1] / s[jj - 1];
						cos = p5 / Math.sqrt(p25 + p25 * (tan * tan));
						sin = cos * tan;
						tau = sin;
					} else {
						cotan = s[jj - 1] / w[j - 1];
						sin = p5 / Math.sqrt(p25 + p25 * (cotan * cotan));
						cos = sin * cotan;
						tau = 1;
						if (Math.abs(cos) * giant > 1)
							tau = 1 / cos;
					}
//					Apply the transformation to s and reduce the spike in w.
					l = jj;
					for (i = j; i < m + 1; i++) {

						temp = cos * s[l - 1] + sin * w[i - 1];
						w[i - 1] = -sin * s[l - 1] + cos * w[i - 1];
						s[l - 1] = temp;
						l = l + 1;
					}
//					Store the information necessary to recover the
//					givens rotation.
					w[j - 1] = tau;
				}
//				Test for zero diagonal elements in the output s.
				if (s[jj - 1] == 0) {
					sing = true;
				}
				jj = jj + (m - j + 1);
			}
		}
//		Move w back into the last column of the output s.
		l = jj;
		for (i = n; i < m + 1; i++) {
			s[l - 1] = w[i - 1];
			l = l + 1;
		}
		if (s[jj - 1] == 0) {
			sing = true;
		}
	}

	public static double[] rwupdt(int n, double[][] r, int ldr, double[] w, double[] b, double alpha, double[] cos,
			double[] sin) {
		// !*****************************************************************************************
		// !>
		// ! given an n by n upper triangular matrix r, this subroutine
		// ! computes the qr decomposition of the matrix formed when a row
		// ! is added to r. if the row is specified by the vector w, then
		// ! rwupdt determines an orthogonal matrix q such that when the
		// ! n+1 by n matrix composed of r augmented by w is premultiplied
		// ! by (q transpose), the resulting matrix is upper trapezoidal.
		// ! the matrix (q transpose) is the product of n transformations
		// !```
		// ! g(n)*g(n-1)* ... *g(1)
		// !```
		// ! where g(i) is a givens rotation in the (i,n+1) plane which
		// ! eliminates elements in the (n+1)-st plane. rwupdt also
		// ! computes the product (q transpose)*c where c is the
		// ! (n+1)-vector (b,alpha). q itself is not accumulated, rather
		// ! the information to recover the g rotations is supplied.
		//
//	    subroutine rwupdt(n, r, Ldr, w, b, Alpha, Cos, Sin)
//	        implicit none
		//
//	        integer, intent(in) :: n !! a positive integer input variable set to the order of r.
//	        integer, intent(in) :: Ldr !! a positive integer input variable not less than n
//	                                  !! which specifies the leading dimension of the array r.
//	        real(wp), intent(inout) :: Alpha !! a variable. on input alpha must contain the
//	                                        !! (n+1)-st element of the vector c. on output alpha contains
//	                                        !! the (n+1)-st element of the vector (q transpose)*c.
//	        real(wp), intent(inout) :: r(Ldr, n) !! an n by n array. on input the upper triangular part of
//	                                            !! r must contain the matrix to be updated. on output r
//	                                            !! contains the updated triangular matrix.
//	        real(wp), intent(in) :: w(n) !! an input array of length n which must contain the row
//	                                    !! vector to be added to r.
//	        real(wp), intent(inout) :: b(n) !! an array of length n. on input b must contain the
//	                                       !! first n elements of the vector c. on output b contains
//	                                       !! the first n elements of the vector (q transpose)*c.
//	        real(wp), intent(out) :: Cos(n) !! an output array of length n which contains the
//	                                       !! cosines of the transforming givens rotations.
//	        real(wp), intent(out) :: Sin(n) !! an output array of length n which contains the
//	                                       !! sines of the transforming givens rotations.
		//
//	        integer :: i, j, jm1
//	        real(wp) :: cotan, rowj, tan, temp
//				
		return w;

//	        real(wp), parameter :: p5 = 5.0e-1_wp
//	        real(wp), parameter :: p25 = 2.5e-1_wp
		//
//	        do j = 1, n
//	            rowj = w(j)
//	            jm1 = j - 1
		//
//	            ! apply the previous transformations to
//	            ! r(i,j), i=1,2,...,j-1, and to w(j).
		//
//	            if (jm1 >= 1) then
//	                do i = 1, jm1
//	                    temp = Cos(i)*r(i, j) + Sin(i)*rowj
//	                    rowj = -Sin(i)*r(i, j) + Cos(i)*rowj
//	                    r(i, j) = temp
//	                end do
//	            end if
		//
//	            ! determine a givens rotation which eliminates w(j).
		//
//	            Cos(j) = one
//	            Sin(j) = zero
//	            if (rowj /= zero) then
//	                if (abs(r(j, j)) >= abs(rowj)) then
//	                    tan = rowj/r(j, j)
//	                    Cos(j) = p5/sqrt(p25 + p25*tan**2)
//	                    Sin(j) = Cos(j)*tan
//	                else
//	                    cotan = r(j, j)/rowj
//	                    Sin(j) = p5/sqrt(p25 + p25*cotan**2)
//	                    Cos(j) = Sin(j)*cotan
//	                end if
		//
//	                ! apply the current transformation to r(j,j), b(j), and alpha.
		//
//	                r(j, j) = Cos(j)*r(j, j) + Sin(j)*rowj
//	                temp = Cos(j)*b(j) + Sin(j)*Alpha
//	                Alpha = -Sin(j)*b(j) + Cos(j)*Alpha
//	                b(j) = temp
//	            end if
//	        end do
		//
//	    end subroutine rwupdt
		// !*****************************************************************************************
		//
	}
	// !*****************************************************************************************
	// end module minpack_module
	// !*****************************************************************************************
}
