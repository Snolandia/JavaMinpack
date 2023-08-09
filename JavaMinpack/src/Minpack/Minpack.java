package Minpack;

public class Minpack {

//	Converted from Fortran, then modified slightly. Work in progress to finish
//	the rest of the functions as currently, only hybrd1/hybrd is fully functional

//**********Disclaimer from original Minpack: **********************\\

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
//	4. waRRANTY DISCLAIMER. THE SOFTwaRE IS SUPPLIED "AS IS"
//	WITHOUT waRRANTY OF ANY KIND. THE COPYRIGHT HOLDER, THE
//	UNITED STATES, THE UNITED STATES DEPARTMENT OF ENERGY, AND
//	THEIR EMPLOYEES: (1) DISCLAIM ANY waRRANTIES, EXPRESS OR
//	IMPLIED, INCLUDING BUT NOT LIMITED TO ANY IMPLIED waRRANTIES
//	OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE
//	OR NON-INFRINGEMENT, (2) DO NOT ASSUME ANY LEGAL LIABILITY
//	OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR
//	USEFULNESS OF THE SOFTwaRE, (3) DO NOT REPRESENT THAT USE OF
//	THE SOFTwaRE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS, (4)
//	DO NOT waRRANT THAT THE SOFTwaRE WILL FUNCTION
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
//	EVEN IF ANY OF SAID PARTIES HAS BEEN waRNED OF THE
//	POSSIBILITY OF SUCH LOSS OR DAMAGES.

	static double n;
	// static double epsmch = Math.ulp(n); //removing as too much precision causes issues
	static double epsmch = 0.0000000000000002220446049250313; // Using this instead due to above statement
	static String[] infoDict = { "info 0 : Input parameters are incorrect.",
			"info 1 : Relative error between two consecutive iterates is at most `xTol`.",
			"info 2 : Iterations has exceded maxFev, the max allowed iterations",
			"info 3 : xTol is too small. No further improvement in the approximate solution `x` is possible.",
			"info 4 : Iteration is not making good progress, measured by the improvement from the last five jacobian evaluations.",
			"info 5 : Tteration is not making good progress, as measured by the improvement from the last ten iterations." };

	/**
	 * <p>
	 * NOT TESTED YET This subroutine checks the gradients of m nonlinear functions
	 * in n variables, evaluated at a point x, for consistency with the functions
	 * themselves.
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
	 * @param ldfJac int -> A positive integer, not less than m, which specifies the
	 *               leading dimension of the array fJac. (i.e., The number of
	 *               functions/equations.ldfJac can just be set to m)
	 *               <p>
	 * @param mode   int -> An integer input set to 1 on the first call and 2 on the
	 *               second. The user must call chkder twice,first with mode = 1 and
	 *               then with mode = 2. mode = 1 -> x must contain the point of
	 *               evaluation. xp is set to a neighboring point. mode = 2 -> fVec
	 *               must contain the functions and the rows of fJac must contain
	 *               the gradients of the respective functions each evaluated at x,
	 *               and fVecp must contain the functions evaluated at xp. Err
	 *               contains measures of correctness of the respective gradients.
	 *               <p>
	 * @param x      double[n] -> Input array (i.e. The variables)
	 *               <p>
	 * @param fVec   double[m] -> An array of length m. When mode = 2, fVec must
	 *               contain the functions evaluated at x.
	 *               <p>
	 * @param fJac   double[ldfJac][n] -> an m by n array. When mode = 2, the rows
	 *               of fJac must contain the gradients of the respective functions
	 *               evaluated at x.
	 *               <p>
	 * @param xp     double[n] -> An array of length n. When mode = 1, xp is set to
	 *               a neighboring point of x.
	 *               <p>
	 * @param fVecp  double[m] -> An array of length m. When mode = 2, fVecp must
	 *               contain the functions evaluated at xp.
	 *               <p>
	 * @param err    double[m] -> An array of length m. When mode = 2, err contains
	 *               measures of correctness of the respective gradients. If there
	 *               is no severe loss of significance, then if err[i] is 1.0 the
	 *               i-th gradient is correct, while if err[i] is 0.0 the i-th
	 *               gradient is incorrect. For values of err between 0.0 and 1.0,
	 *               the categorization is less certain. In general, a value of
	 *               err[i] greater than 0.5 indicates that the i-th gradient is
	 *               probably correct, while a value of err[i] less than 0.5
	 *               indicates that the i-th gradient is probably incorrect.
	 * 
	 * @return Nothing is returned, rather the objects themselves are modified.
	 */
	public static void chkder(int m, int n, double[] x, double[] fVec, double[][] fJac, int ldfJac, double[] xp,
			double[] fVecp, int mode, double[] err) {

		// NOT TESTED YET

		int i, j; // Only used for loops, not important.
		double temp; // Just used as a temp variable
		double eps = Math.sqrt(epsmch); // sqrt of machine epsilon
		double factor = 100.0; // factor, or percent of change that happens.
		double epsf = factor * epsmch; // factor of epsmch, or percent of change that will happen.
		double epslog = Math.log10(epsf); // log 10 of epsf

		if (mode == 2) {
			err = new double[m];
			for (j = 0; j <= n; j++) {
				temp = Math.abs(x[j]);
				if (temp == 0) {
					temp = 1.0;
				}
				for (i = 0; i <= m; i++) {
					err[i] = err[i] + temp * fJac[i][j];
				}
			}
			for (i = 0; i <= m; i++) {
				temp = 1.0;
				if (fVec[i] != 0 && fVecp[i] != 0 && Math.abs(fVecp[i] - fVec[i]) >= epsf * Math.abs(fVec[i])) {
					temp = eps * Math.abs((fVecp[i] - fVec[i]) / eps - err[i])
							/ (Math.abs(fVec[i]) + Math.abs(fVecp[i]));
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
		double alpha, bnorm, gNorm, qnorm, sgNorm, sum, temp; // temporary doubles
//	    First, calculate the gauss-newton direction.
		jj = (n * (n + 1)) / 2;
		for (k = 0; k < n; k++) {
			j = n - (k + 1);
			jp1 = j + 1;
			jj = jj - (k + 1);
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
			gNorm = enorm(n, wa1);
			sgNorm = 0;
			alpha = delta / qnorm;
			if (gNorm != 0) {
//			Calculate the point along the scaled gradient
//			at which the quadratic is minimized.
				for (j = 0; j < n; j++) {
					wa1[j] = (wa1[j] / gNorm) / diag[j];
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
				sgNorm = (gNorm / temp) / temp;
//				Test whether the scaled gradient direction is acceptable.
				alpha = 0;
				if (sgNorm < delta) {
//				The scaled gradient direction is not acceptable.
//				Finally, calculate the point along the dogleg
//				at which the quadratic is minimized.
					bnorm = enorm(n, qtb);
					temp = (bnorm / gNorm) * (bnorm / qnorm) * (sgNorm / delta);
					temp = temp - (delta / qnorm) * ((sgNorm / delta) * (sgNorm / delta))
							+ Math.sqrt((temp - (delta / qnorm)) * (temp - (delta / qnorm)))
							+ (1 - ((delta / qnorm) * (delta / qnorm))) * (1 - ((sgNorm / delta) * (sgNorm / delta)));
					alpha = ((delta / qnorm) * (1 - ((sgNorm / delta) * (sgNorm / delta)))) / temp;
				}
			}
//			Form appropriate convex combination of the gauss-newton
//			direction and the scaled gradient direction.
			if (sgNorm < delta) {
				temp = (1 - alpha) * sgNorm;
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
	 * The euclidean norm is computed by accumulating the sum of  squares in three
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
	 * @param ldfJac int -> A positive integer variable not less than n which
	 *               specifies the leading dimension of the array fJac. I.e. Number
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
	 * @param fVec   double[n] -> An array of length n which must contain the
	 *               functions evaluated at x.
	 *               <p>
	 * @param fJac   double[ldfJac][n] -> A n by n array which contains the
	 *               approximation to the jacobian matrix evaluated at x.
	 *               <p>
	 * @param wa1    double[n] -> A working array of length n.
	 *               <p>
	 * @param wa2    double[n] -> A working array of length n. If ml + mu + 1 is at
	 *               least n, then the jacobian is considered dense, and wa2 is not
	 *               referenced.
	 * @return Nothing is returned, rather the objects themselves are modified.
	 */
	public static void fdjac1(SystemOfEquations fcn, int n, double[] x, double[] fVec, double[][] fJac, int ldfJac,
			int iFlag, int ml, int mu, double epsfcn, double[] wa1, double[] wa2) {

		int i, j, k, msum; // Temporary variables
		double eps, h, temp; // Temporary variables

		// double epsmch = 0.0000000000000002220446049250313; // Removing this and using
		// the global causes fdjac1 to return too small of a gradient change

		if (epsfcn > epsmch) {
			eps = Math.sqrt(epsfcn);
		} else {
			eps = Math.sqrt(epsmch);
		}
		msum = ml + mu + 1;
		if (msum < n) {
//			Computation of banded approximate jacobian.
			for (k = 0; k < msum + 1; k++) {
				for (j = k; j < n; j = +msum) {
					wa2[j] = x[j];
					h = eps * Math.abs(wa2[j]);
					if (h == 0) {
						h = eps;
					}
					x[j] = wa2[j] + h;
				}
				wa1 = fcn.evaluate(x); // This part needs to be changed later
				for (j = k; k < n; k = +msum) {
					x[j] = wa2[j];
					h = eps * Math.abs(wa2[j]);
					if (h == 0) {
						h = eps;
					}
					for (i = 0; i < n; i++) {
						fJac[i][j] = 0;
						if (i >= j - mu && i <= j + ml) {
							fJac[i][j] = (wa1[i] - fVec[i]) / h;
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
					fJac[i][j] = (wa1[i] - fVec[i]) / h;
				}
			}
		}
	}

	/**
	 * 
	 * <p>
	 * NOT TESTED YET This subroutine computes a forward-difference approximation to
	 * the m by n jacobian matrix associated with a specified problem of m functions
	 * in n variables.
	 * 
	 * procedure(func2) :: fcn  the user-supplied subroutine which calculates the
	 * functions.
	 * <p>
	 * 
	 * @param m      int -> A positive integer variable set to the number of
	 *               functions.
	 *               <p>
	 * @param n      int -> A positive integer variable set to the number of
	 *               variables. n must not exceed m.
	 *               <p>
	 * @param ldfJac int -> A positive integer variable not less than m which
	 *               specifies the leading dimension of the array fJac.
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
	 * @param fVec   double[m] -> An array of length m which must contain the
	 *               functions evaluated at x.
	 *               <p>
	 * @param fJac   double[ldfJac, n] -> A m by n array which contains the
	 *               approximation to the jacobian matrix evaluated at x.
	 *               <p>
	 * @param wa     double[m] -> A working array of length m.
	 *               <p>
	 * @return Nothing is returned, rather the objects themselves are modified.
	 */
	public static void fdjac2(SystemOfEquations fcn, int m, int n, double[] x, double[] fVec, double[][] fJac,
			int ldfJac, int iFlag, double epsfcn, double[] wa) {

		// NOT TESTED YET

		int i, j; // used for loops
		double eps, h, temp; // temporary variables

		if (epsfcn > epsmch) {
			eps = Math.sqrt(epsfcn);
		} else {
			eps = Math.sqrt(epsmch);
		}
		for (j = 0; j < n; j++) {
			temp = x[j];
			h = eps * Math.abs(temp);
			if (h == 0) {
				h = eps;
			}
			x[j] = temp + h;
			fcn.evaluate(x);
			if (iFlag < 0) {
				throw new IllegalArgumentException(
						"Class method fdjac2, at call to fcn, has returned an iFlag less than 0");
			}
			x[j] = temp;
			for (i = 0; i < m; i++) {
				fJac[i][j] = (wa[i] - fVec[i]) / h;
			}
		}
	}

	/**
	 * <p>
	 * The purpose of hybrd is to find a zero of a system of n nonlinear functions
	 * in n variables by a modification of the powell hybrid method. the user must
	 * provide a subroutine which calculates the functions. The jacobian is then
	 * calculated by a forward-difference approximation.
	 * 
	 * procedure(func) :: fcn  user-supplied subroutine which calculates the
	 * functions
	 * 
	 * @param n      int -> A positive integer variable set to the number of
	 *               functions and variables.
	 *               <p>
	 * @param maxFev int -> A positive integer variable. Termination occurs when the
	 *               number of calls to `fcn` is at least `maxFev` by the end of an
	 *               iteration.
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
	 * @param nPrint int -> An integer variable that enables controlled printing of
	 *               iterates if it is positive. In this case, `fcn` is called with
	 *               `iFlag = 0` at the beginning of the first iteration and every
	 *               `nPrint` iterations thereafter and immediately prior to return,
	 *               with `x` and `fVec` available for printing. If `nPrint` is not
	 *               positive, no special calls of `fcn` with `iFlag = 0` are made.
	 *               <p>
	 * @param info   int -> An integer variable. If the user has terminated
	 *               execution, `info` is set to the (negative) value of `iFlag`.
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
	 * @param ldfJac int -> A positive integer variable not less than `n` which
	 *               specifies the leading dimension of the array `fJac`.
	 *               <p>
	 * @param lr     int -> A positive integer input variable not less than
	 *               `(n*(n+1))/2`.
	 *               <p>
	 * @param xTol   double -> A nonnegative variable. Termination occurs when the
	 *               relative error between two consecutive iterates is at most
	 *               `xTol`.
	 *               <p>
	 * @param epsfcn double -> A variable used in determining a suitable step length
	 *               for the forward-difference approximation. This approximation
	 *               assumes that the relative errors in the functions are of the
	 *               order of `epsfcn`. If `epsfcn` is less than the machine
	 *               precision, it is assumed that the relative errors in the
	 *               functions are of the order of the machine precision.
	 *               <p>
	 * @param factor double -> A positive variable used in determining the initial
	 *               step bound. This bound is set to the product of `factor` and
	 *               the euclidean norm of `diag*x` if nonzero, or else to `factor`
	 *               itself. In most cases factor should lie in the interval
	 *               (.1,100.). -> (100.0) is a generally recommended value.
	 *               <p>
	 * @param x      double[n] -> Array of length n. On input `x` must contain an
	 *               initial estimate of the solution vector. On output `x` contains
	 *               the final estimate of the solution vector.
	 *               <p>
	 * @param fVec   double[n] -> An output array of length `n` which contains the
	 *               functions evaluated at the output `x`.
	 *               <p>
	 * @param diag   double[n] -> An array of length `n`. If `mode = 1` (see below),
	 *               `diag` is internally set. If `mode = 2`, `diag` must contain
	 *               positive entries that serve as multiplicative scale factors for
	 *               the variables.
	 *               <p>
	 * @param fJac   double[ldfJac][n] -> Array which contains the orthogonal matrix
	 *               `q` produced by the QR factorization of the final approximate
	 *               jacobian.
	 *               <p>
	 * @param r      double[lr] -> An array which contains the upper triangular
	 *               matrix produced by the QR factorization of the final
	 *               approximate jacobian, stored row-wise.
	 *               <p>
	 * @param qtf    double[n] -> An output array of length `n` which contains the
	 *               vector `(q transpose)*fVec`.
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
	public static double[] hybrd(SystemOfEquations fcn, int n, double[] x, double[] fVec, double xTol, int maxFev,
			int ml, int mu, double epsfcn, double[] diag, int mode, double factor, int nPrint, int info, int nFev,
			double[][] fJac, int ldfJac, double[] r, int lr, double[] qtf, double[] wa1, double[] wa2, double[] wa3,
			double[] wa4) {

		int i, iFlag, iter, j, jm1, l, msum, ncFail, ncSuc, nSlow1, nSlow2; // Temporary Variables
		int[] iwa = new int[0]; // Temporary Variables
		boolean jeval, sing; // Temporary Variables
		double actred, delta = 0, fNorm, fNorm1, pNorm, prered, ratio, sum, temp, xNorm = 0; // Temporary Variables
		double[][] qtfHold = new double[1][n]; // Used to create a sub-array to pass

		double p1 = 0.1;
		double p5 = 0.5;
		double p001 = 0.001;
		double p0001 = 0.0001;

		info = 0;
		iFlag = 0;
		nFev = 0;
		boolean main = true;
		boolean inner = true;
		boolean outer = true;


		do { // main loop
//	        Check the input parameters for errors.
			if (n <= 0 || xTol < 0 || maxFev <= 0 || ml < 0 || mu < 0 || factor <= 0 || ldfJac < n
					|| lr < (n * (n + 1)) / 2) {
				throw new IllegalArgumentException("Input parameters do not match as expected");
			}
			if (mode == 2) {
				for (int a = 0; a < n; a++) {
					if (diag[a] <= 0) {
						throw new IllegalArgumentException("diagonal Error at Hybrd");
					}
				}
			}
//	        Evaluate the function at the starting point
//	        and calculate its norm.
			iFlag = 1;
			fVec = fcn.evaluate(x);
			nFev = 1;
			if (iFlag < 0) {
				throw new IllegalArgumentException("iFlag was returned as negative");
			}
			fNorm = enorm(n, fVec);
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

			do { // outer loop
				jeval = true;
//	            Calculate the jacobian matrix.
				iFlag = 2;
				fdjac1(fcn, n, x, fVec, fJac, ldfJac, iFlag, ml, mu, epsfcn, wa1, wa2);
				nFev = nFev + msum;
//	            Compute the qr factorization of the jacobian.
				qrfac(n, n, fJac, ldfJac, false, iwa, 1, wa1, wa2, wa3);
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
//	            Form (q transpose)*fVec and store in qtf.
				for (i = 0; i < n; i++) {
					qtf[i] = fVec[i];
				}
				for (j = 0; j < n; j++) {
					if (fJac[j][j] != 0) {
						sum = 0;
						for (i = j; i < n; i++) {
							sum = sum + fJac[i][j] * qtf[i];
						}
						temp = -sum / fJac[j][j];
						for (i = j; i < n; i++) {
							qtf[i] = qtf[i] + fJac[i][j] * temp;
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
							r[l] = fJac[i][j];
							l = l + n - 1 - i;
						}
					}
					r[l] = wa1[j];
					if (wa1[j] == 0) {
						sing = true;
					}
				}
//				Accumulate the orthogonal factor in fJac.
				qform(n, n, fJac, ldfJac, wa1);
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
//							Add a method to do the printing
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
					nFev = nFev + 1;
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
//						Successful iteration. update x, fVec, and their norms.
						for (j = 0; j < n; j++) {
							x[j] = wa2[j];
							wa2[j] = diag[j] * x[j];
							fVec[j] = wa4[j];
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
					if (nFev >= maxFev) {
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
						return x;
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
							sum = sum + fJac[i][j] * wa4[i];
						}
						wa2[j] = (sum - wa3[j]) / pNorm;
						wa1[j] = diag[j] * ((diag[j] * wa1[j]) / pNorm);
						if (ratio >= p0001) {
							qtf[j] = sum;
						}
					}
//					Compute the qr factorization of the updated jacobian.
					r1updt(n, n, r, lr, wa1, wa2, wa3, sing);
					r1mpyq(n, n, fJac, ldfJac, wa2, wa3);
					for (int z = 0; z < n; z++) {
						qtfHold[0][z] = qtf[z]; // Used to pass a smaller array to r1mpyq
					}
					qtfHold = r1mpyq(1, n, qtfHold, 1, wa2, wa3);
					for (int z = 0; z < n; z++) {
						qtf[z] = qtfHold[0][z]; // Used to pass r1mpyq results back to qtf
					}
					jeval = false;
				} while (inner); // End of the inner loop.
			} while (outer); // End of the outer loop.
		} while (main); // End of main loop
		return x;
	}

	/**
	 * <p>
	 * The purpose of hybrd is to find a zero of a system of n nonlinear functions
	 * in n variables by a modification of the powell hybrid method. the user must
	 * provide a subroutine which calculates the functions. The jacobian is then
	 * calculated by a forward-difference approximation.
	 * 
	 * procedure(func) :: fcn  user-supplied subroutine which calculates the
	 * functions
	 * @param tol   double -> A nonnegative variable. Termination occurs when the
	 *               relative error between two consecutive iterates is at most
	 *               `xTol`.
	 *               <p>
	 * @param x      double[n] -> Array of length n. On input `x` must contain an
	 *               initial estimate of the solution vector. On output `x` contains
	 *               the final estimate of the solution vector.
	 *               <p> 
	 * @return Returns a double[] that contains the solved values, x
	 * 
	 */	
	public static double[] hybrd1(SystemOfEquations fcn, double[] x, double tol) {
		
//	    Check the input parameters for errors.
				
		if (fcn.size() != x.length) {
			System.out.println("not square, exiting now");
			System.exit(0);
		}
		if (tol <= 0) {
			tol = 1.490116119384766E-008;
		}
		int info = 0;
		int n = x.length;
		int lwa = (n * (3 * n + 13)) / 2;
		if (n > 0 & tol >= 0 & lwa >= (n * (3 * n + 13)) / 2) {
//			Call hybrd.
			int maxFev = 200 * (n + 1);
			double xTol = tol;
			int ml = n;
			int mu = n;
			double epsfcn = 0;
			int mode = 2;
			int nPrint = 0;
			int lr = (n * (n + 1)) / 2;
			int nFev = 0;
			double[] fVec = new double[n];
			double[] diag = new double[n];
			for (int j = 0; j < n; j++) {
				diag[j] = 1;
			}
			double factor = 100.0;
			double[][] fJac = new double[n][n];
			int ldfJac = n;
			double[] r = new double[lr];
			double[] qtf = new double[n];
			double[] wa1 = new double[n];
			double[] wa2 = new double[n];
			double[] wa3 = new double[n];
			double[] wa4 = new double[n];

			x = hybrd(fcn, n, x, fVec, xTol, maxFev, ml, mu, epsfcn, diag, mode, factor, nPrint, info, nFev, fJac,
					ldfJac, r, lr, qtf, wa1, wa2, wa3, wa4);

			if (info == 5) {
				info = 4;
			}
		}
		return x;
	}

	/**
	 * The purpose of hybrj is to find a zero of a system of n nonlinear
	 * functions in n variables by a modification of the powell hybrid method. The
	 * user must provide a subroutine which calculates the functions and the
	 * jacobian.
	 * 
	 * subroutine hybrj(fcn, n, x, fVec, fJac, ldfJac, xTol, maxFev, diag, mode, &
	 * factor, nPrint, info, nFev, nJev, r, Lr, qtf, wa1, wa2, & wa3, wa4)
	 * 
	 * implicit none
	 * 
	 * procedure(fcn_hybrj) :: fcn  the user-supplied subroutine which 
	 * calculates the functions and the jacobian 
	 * integer, intent(in) :: n  a
	 * positive integer input variable set to the number  of functions and
	 * variables. 
	 * integer, intent(in) :: ldfJac  a positive integer input variable
	 * not less than n  which specifies the leading dimension of the array fJac.
	 * integer, intent(in) :: maxFev  a positive integer input variable.
	 * termination  occurs when the number of calls to fcn with iFlag = 1  has
	 * reached maxFev.
	 *  integer, intent(in) :: mode  an integer input variable. if
	 * mode = 1, the  variables will be scaled internally. if mode = 2,  the
	 * scaling is specified by the input diag. other  values of mode are
	 * equivalent to mode = 1.
	 *  integer, intent(in) :: nPrint  an integer input
	 * variable that enables controlled  printing of iterates if it is positive.
	 * in this case,  fcn is called with iFlag = 0 at the beginning of the first
	 *  iteration and every nPrint iterations thereafter and  immediately prior
	 * to return, with x and fVec available  for printing. fVec and fJac should
	 * not be altered.  if nPrint is not positive, no special calls of fcn  with
	 * iFlag = 0 are made. integer, intent(out) :: info  an integer output
	 * variable. if the user has  terminated execution, info is set to the
	 * (negative)  value of iFlag. see description of fcn. otherwise,  info is
	 * set as follows:   * ***info = 0*** improper input parameters.  *
	 * ***info = 1*** relative error between two consecutive iterates  is at most
	 * xTol.  * ***info = 2*** number of calls to fcn with iFlag = 1 has 
	 * reached maxFev.  * ***info = 3*** xTol is too small. no further improvement
	 * in  the approximate solution x is possible.  * ***info = 4*** iteration
	 * is not making good progress, as  measured by the improvement from the last
	 *  five jacobian evaluations.  * ***info = 5*** iteration is not making
	 * good progress, as  measured by the improvement from the last  ten
	 * iterations. integer, intent(out) :: nFev  an integer output variable set to
	 * the number of  calls to fcn with iFlag = 1. integer, intent(out) :: nJev 
	 * an integer output variable set to the number of  calls to fcn with iFlag =
	 * 2. integer, intent(in) :: Lr  a positive integer input variable not less
	 * than  (n*(n+1))/2. real(wp), intent(in) :: xTol  a nonnegative input
	 * variable. termination  occurs when the relative error between two
	 * consecutive  iterates is at most xTol. real(wp), intent(in) :: factor  a
	 * positive input variable used in determining the  initial step bound. this
	 * bound is set to the product of  factor and the euclidean norm of diag*x if
	 * nonzero, or else  to factor itself. in most cases factor should lie in the
	 *  interval (.1,100.). 100. is a generally recommended value. real(wp),
	 * intent(inout) :: x(n)  an array of length n. on input x must contain  an
	 * initial estimate of the solution vector. on output x  contains the final
	 * estimate of the solution vector. real(wp), intent(out) :: fVec(n)  an
	 * output array of length n which contains  the functions evaluated at the
	 * output x. real(wp), intent(out) :: fJac(ldfJac, n)  an output n by n array
	 * which contains the  orthogonal matrix q produced by the qr factorization 
	 * of the final approximate jacobian. real(wp), intent(inout) :: diag(n)  an
	 * array of length n. if mode = 1 (see  below), diag is internally set. if
	 * mode = 2, diag  must contain positive entries that serve as 
	 * multiplicative scale factors for the variables. real(wp), intent(out) ::
	 * r(Lr)  an output array of length lr which contains the  upper triangular
	 * matrix produced by the qr factorization  of the final approximate jacobian,
	 * stored rowwise. real(wp), intent(out) :: qtf(n)  an output array of length
	 * n which contains  the vector (q transpose)*fVec. real(wp), intent(inout) ::
	 * wa1(n)  work array of length n. real(wp), intent(inout) :: wa2(n)  work
	 * array of length n. real(wp), intent(inout) :: wa3(n)  work array of length
	 * n. real(wp), intent(inout) :: wa4(n)  work array of length n.
	 */

	public static double[] hybrj(SystemOfEquations fcn, int n, double[] x, double[] fVec, double[][] fJac, int ldfJac,
			double xTol, int maxFev, double[] diag, int mode, double factor, int nPrint, int info, int nFev, int nJev,
			double[] r, int lr, double[] qtf, double[] wa1, double[] wa2, double[] wa3, double[] wa4) {

		int i, iFlag, iter, j, jm1, l, ncfail, ncsuc, nslow1, nslow2 = 0;
		int[] iwa = new int[1]; // need to verify 1 is the correct value
		boolean jeval, sing;
		double actred = 0, delta = 0, fNorm = 0, fNorm1 = 0, pnorm = 0, prered = 0, ratio = 0, sum = 0, temp = 0,
				xnorm = 0;
		double[][] qtfHold = new double[1][n]; // Used to create a sub-array to pass

		double p1 = 1.0e-1;
		double p5 = 5.0e-1;
		double p001 = 1.0e-3;
		double p0001 = 1.0e-4;
		info = 0;
		iFlag = 0;
		nFev = 0;
		nJev = 0;

		do { // main : block
//	            Check the input parameters for errors.
			if (n <= 0 || ldfJac < n || xTol < 0 || maxFev <= 0 || factor <= 0 || lr < (n * (n + 1)) / 2) {
				throw new IllegalArgumentException("Input parameters do not match as expected");
				// Add throw error
			}
			if (mode == 2) {
				for (j = 0; j < n; j++) {
					if (diag[j] <= 0) {
						throw new IllegalArgumentException("diagonal Error at Hybrd");
						// Add throw error
					}
				}
			}
//	            Evaluate the function at the starting point
//	            and calculate its norm.
			iFlag = 1;
			fcn.evaluate(x); // to update at future date
			nFev = 1;
			if (iFlag < 0) {
				throw new IllegalArgumentException("iFlag was returned as negative");
			}
			fNorm = enorm(n, fVec);
//	            Initialize iteration counter and monitors.
			iter = 1;
			ncsuc = 0;
			ncfail = 0;
			nslow1 = 0;
			nslow2 = 0;
//	            Beginning of the outer loop.
			do { // outer : do
				jeval = true;
//	                Calculate the jacobian matrix.
				iFlag = 2;
				fVec = fcn.evaluate(x); // to update at future date
				fJac = fcn.evaluateJacobian(x); // to update at future date
				nJev = nJev + 1;
				if (iFlag < 0) {
					// throw error
				}
//	                Compute the qr factorization of the jacobian.
				qrfac(n, n, fJac, ldfJac, false, iwa, 1, wa1, wa2, wa3); // maybe to do?
//	                On the first iteration and if mode is 1, scale according
//	                to the norms of the columns of the initial jacobian.
				if (iter == 1) {
					if (mode != 2) {
						for (j = 0; j < n; j++) {
							diag[j] = wa2[j];
							if (wa2[j] == 0) {
								diag[j] = 1;
							}
						}
					}
//	                    On the first iteration, calculate the norm of the scaled x
//	                    and initialize the step bound delta.
					for (j = 0; j < n; j++) {
						wa3[j] = diag[j] * x[j];
					}
					xnorm = enorm(n, wa3);
					delta = factor * xnorm;
					if (delta == 0) {
						delta = factor;
					}
				}
//	                Form (q transpose)*fVec and store in qtf.
				for (i = 0; i < n; i++) {
					qtf[i] = fVec[i];
				}
				for (j = 0; j < n; j++) {
					if (fJac[j][j] != 0) {
						sum = 0;
						for (i = j; i < n; i++) {
							sum = sum + fJac[i][j] * qtf[i];
						}
						temp = -sum / fJac[j][j];
						for (i = j; i < n; i++) {
							qtf[i] = qtf[i] + fJac[i][j] * temp;
						}
					}
				}
//	                Copy the triangular factor of the qr factorization into r.
				sing = false;
				for (j = 0; j < n; j++) {
					l = j;
					jm1 = j - 1;
					if (jm1 >= 0) {
						for (i = 0; i < jm1; i++) {
							r[l] = fJac[i][j];
							l = l + n - i;
						}
					}
					r[l] = wa1[j];
					if (wa1[j] == 0) {
						sing = true;
					}
				}
//	                Accumulate the orthogonal factor in fJac.
				qform(n, n, fJac, ldfJac, wa1);
//	                Rescale if necessary.
				if (mode != 2) {
					for (j = 0; j < n; j++) {
						if (diag[j] > wa2[j]) {
							diag[j] = diag[j];
						} else {
							diag[j] = wa2[j];
						}
					}
				}
//	                Beginning of the inner loop.
//	                inner : do
				do {
//	                    If requested, call fcn to enable printing of iterates.
					if (nPrint > 0) {
						iFlag = 0;
						if ((iter - 1) % nPrint == 0) {
							fcn.evaluate(x); // ToDo fix
						}
						if (iFlag < 0) {
							// Add Error Handling
						}
					}
//	                    Determine the direction p.
					dogleg(n, r, lr, diag, qtf, delta, wa1, wa2, wa3);
//	                    Store the direction p and x + p. calculate the norm of p.
					for (j = 0; j < n; j++) {
						wa1[j] = -wa1[j];
						wa2[j] = x[j] + wa1[j];
						wa3[j] = diag[j] * wa1[j];
					}
					pnorm = enorm(n, wa3);
//	                    On the first iteration, adjust the initial step bound.
					if (iter == 1) {
						if (pnorm < delta) {
							delta = pnorm;
						}
//	                    Evaluate the function at x + p and calculate its norm.
						iFlag = 1;
						wa4 = fcn.evaluate(wa2);// (n, wa2, wa4, fJac, ldfJac, iFlag); //I think this is right
						nFev = nFev + 1;
						if (iFlag < 0) {
							// add error handling
						}
						fNorm1 = enorm(n, wa4);
//	                    Compute the scaled actual reduction.
						actred = -1;
						if (fNorm1 < fNorm) {
							actred = 1 - (fNorm1 / fNorm) * (fNorm1 / fNorm);
						}
//	                    Compute the scaled predicted reduction.
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
							prered = 1 - (temp / fNorm) * (temp / fNorm);
						}
//	                    Compute the ratio of the actual to the predicted
//	                    reduction.
						ratio = 0;
						if (prered > 0) {
							ratio = actred / prered;
						}
//	                    Update the step bound.
						if (ratio >= p1) {
							ncfail = 0;
							ncsuc = ncsuc + 1;
							if (ratio >= p5 || ncsuc > 1) {
								if (pnorm / p5 > delta) {
									delta = pnorm / p5;
								}
							}
							if (Math.abs(ratio - 1) <= p1)
								delta = pnorm / p5;
						}
					} else {
						ncsuc = 0;
						ncfail = ncfail + 1;
						delta = p5 * delta;
					}
//	                    Test for successful iteration.
					if (ratio >= p0001) {
						// Successful iteration. update x, fVec, and their norms.
						for (j = 1; j < n; j++) {
							x[j] = wa2[j];
							wa2[j] = diag[j] * x[j];
							fVec[j] = wa4[j];
						}
						xnorm = enorm(n, wa2);
						fNorm = fNorm1;
						iter = iter + 1;
					}
//	                    Determine the progress of the iteration.
					nslow1 = nslow1 + 1;
					if (actred >= p001) {
						nslow1 = 0;
					}
					if (jeval) {
						nslow2 = nslow2 + 1;
					}
					if (actred >= p1) {
						nslow2 = 0;
					}
//	                    Test for convergence.
					if (delta <= xTol * xnorm || fNorm == 0) {
						info = 1;
					}
					if (info != 0) {
						return x;
						// add error handling
					}
//	                    Tests for termination and stringent tolerances.
					if (nFev >= maxFev) {
						info = 2;
					}
					if (p1 * Math.max(p1 * delta, pnorm) <= epsmch * xnorm) {
						info = 3;
					}
					if (nslow2 == 5) {
						info = 4;
					}
					if (nslow1 == 10) {
						info = 5;
					}
					if (info != 0) {
						return x;
						// add error handling
					}
//	                    Criterion for recalculating jacobian.
					if (ncfail == 2) {
						break; // Cycles outer loop
					}
//	                    Calculate the rank one modification to the jacobian
//	                    and update qtf if necessary.
					for (j = 0; j < n; j++) {
						sum = 0;
						for (i = 0; i < n; i++) {
							sum = sum + fJac[i][j] * wa4[i];
						}
						wa2[j] = (sum - wa3[j]) / pnorm;
						wa1[j] = diag[j] * ((diag[j] * wa1[j]) / pnorm);
						if (ratio >= p0001) {
							qtf[j] = sum;
						}
					}
//	                    Compute the qr factorization of the updated jacobian.
					r1updt(n, n, r, lr, wa1, wa2, wa3, sing);
					r1mpyq(n, n, fJac, ldfJac, wa2, wa3);
					for (int z = 0; z < n; z++) {
						qtfHold[0][z] = qtf[z]; // Used to pass a smaller array to r1mpyq
					}
					qtfHold = r1mpyq(1, n, qtfHold, 1, wa2, wa3);
					for (int z = 0; z < n; z++) {
						qtf[z] = qtfHold[0][z]; // Used to pass r1mpyq results back to qtf
					}

					jeval = false;
				} while (true);// End of the inner loop.
			} while (true); // End of the outer loop.
		} while (true); // End block main
	}

	/**
	 * //
	 * *****************************************************************************************
	 * // //
	 * *****************************************************************************************
	 * // > //  the purpose of hybrj1 is to find a zero of a system of //  n
	 * nonlinear functions in n variables by a modification //  of the powell
	 * hybrid method. this is done by using the //  more general nonlinear equation
	 * solver hybrj. the user //  must provide a subroutine which calculates the
	 * functions //  and the jacobian. // // subroutine hybrj1(fcn, n, x, fVec,
	 * fJac, ldfJac, Tol, info, wa, Lwa) // // implicit none // //
	 * procedure(fcn_hybrj) :: fcn  the user-supplied subroutine which // 
	 * calculates the functions and the jacobian // integer, intent(in) :: n  a
	 * positive integer input variable set to the number //  of functions and
	 * variables. // integer, intent(in) :: ldfJac  a positive integer input
	 * variable not less than n //  which specifies the leading dimension of the
	 * array fJac. // integer, intent(out) :: info  an integer output variable. if
	 * the user has //  terminated execution, info is set to the (negative) // 
	 * value of iFlag. see description of fcn. otherwise, //  info is set as
	 * follows: //  //  * ***info = 0*** improper input parameters. //  *
	 * ***info = 1*** algorithm estimates that the relative error //  between x
	 * and the solution is at most tol. //  * ***info = 2*** number of calls to
	 * fcn with iFlag = 1 has //  reached 100*(n+1). //  * ***info = 3*** tol is
	 * too small. no further improvement in //  the approximate solution x is
	 * possible. //  * ***info = 4*** iteration is not making good progress. //
	 * integer, intent(in) :: Lwa  a positive integer input variable not less than
	 * //  (n*(n+13))/2. // real(wp), intent(in) :: Tol  a nonnegative input
	 * variable. termination occurs //  when the algorithm estimates that the
	 * relative error //  between x and the solution is at most tol. // real(wp),
	 * intent(inout) :: x(n)  an array of length n. on input x must contain // 
	 * an initial estimate of the solution vector. on output x //  contains the
	 * final estimate of the solution vector. // real(wp), intent(out) :: fVec(n) 
	 * an output array of length n which contains //  the functions evaluated at
	 * the output x. // real(wp), intent(out) :: fJac(ldfJac, n)  an output n by n
	 * array which contains the //  orthogonal matrix q produced by the qr
	 * factorization //  of the final approximate jacobian. // real(wp),
	 * intent(inout) :: wa(Lwa)  a work array of length lwa.
	 */
	//
	public static void hyberj1(SystemOfEquations fcn, int n, double[] x, double[] fVec, double[][] fJac, int ldfJac,
			double tol, int info, double[] wa, int lwa) {

//			subroutine hybrj1(fcn, n, x, fVec, fJac, ldfJac, Tol, info, wa, Lwa)
		int j, lr, maxFev, mode, nFev = 0, nJev = 0, nPrint;
		double xTol;

		double factor = 100.0;
		double[] wa1 = new double[n];
		double[] wa2 = new double[n];
		double[] wa3 = new double[n];
		double[] wa4 = new double[n];
		double[] wa5 = new double[n];
		double[] r = new double[n];

		info = 0;

//	        Check the input parameters for errors.

		if (n > 0 && ldfJac >= n && tol >= 0 && lwa >= (n * (n + 13)) / 2) {
//	            Call hybrj.
			maxFev = 100 * (n + 1);
			xTol = tol;
			mode = 2;
			for (j = 0; j < n; j++) {
				wa[j] = 1;
			}
			nPrint = 0;
			lr = (n * (n + 1)) / 2;
			hybrj(fcn, n, x, fVec, fJac, ldfJac, xTol, maxFev, wa, mode, factor, nPrint, info, nFev, nJev, r, lr, wa1,
					wa2, wa3, wa4, wa5);
			if (info == 5) {
				info = 4;
			}
		}
	}

	/**
	 * //
	 * *****************************************************************************************
	 * // > //  the purpose of lmder is to minimize the sum of the squares of // 
	 * m nonlinear functions in n variables by a modification of //  the
	 * levenberg-marquardt algorithm. the user must provide a //  subroutine which
	 * calculates the functions and the jacobian. // // subroutine lmder(fcn, m, n,
	 * x, fVec, fJac, ldfJac, fTol, xTol, gTol, maxFev, & // diag, mode, factor,
	 * nPrint, info, nFev, nJev, ipvt, qtf, & // wa1, wa2, wa3, wa4) // // implicit
	 * none // // procedure(fcn_lmder) :: fcn  the user-supplied subroutine which
	 * //  calculates the functions and the jacobian // integer, intent(in) :: m
	 *  a positive integer input variable set to the number //  of functions. //
	 * integer, intent(in) :: n  a positive integer input variable set to the
	 * number //  of variables. n must not exceed m. // integer, intent(in) ::
	 * ldfJac  a positive integer input variable not less than m //  which
	 * specifies the leading dimension of the array fJac. // integer, intent(in) ::
	 * maxFev  a positive integer input variable. termination //  occurs when
	 * the number of calls to fcn with iFlag = 1 //  has reached maxFev. //
	 * integer, intent(in) :: mode  an integer input variable. if mode = 1, the //
	 *  variables will be scaled internally. if mode = 2, //  the scaling is
	 * specified by the input diag. other //  values of mode are equivalent to
	 * mode = 1. // integer, intent(in) :: nPrint  an integer input variable that
	 * enables controlled //  printing of iterates if it is positive. in this
	 * case, //  fcn is called with iFlag = 0 at the beginning of the first // 
	 * iteration and every nPrint iterations thereafter and //  immediately prior
	 * to return, with x, fVec, and fJac //  available for printing. fVec and fJac
	 * should not be //  altered. if nPrint is not positive, no special calls //
	 *  of fcn with iFlag = 0 are made. // integer, intent(out) :: info  an
	 * integer output variable. if the user has //  terminated execution, info is
	 * set to the (negative) //  value of iFlag. see description of fcn.
	 * otherwise, //  info is set as follows: //  //  * ***info = 0***
	 * improper input parameters. //  * ***info = 1*** both actual and predicted
	 * relative reductions //  in the sum of squares are at most fTol. //  *
	 * ***info = 2*** relative error between two consecutive iterates //  is at
	 * most xTol. //  * ***info = 3*** conditions for info = 1 and info = 2 both
	 * hold. //  * ***info = 4*** the cosine of the angle between fVec and any //
	 *  column of the jacobian is at most gTol in //  absolute value. //  *
	 * ***info = 5*** number of calls to fcn with iFlag = 1 has //  reached
	 * maxFev. //  * ***info = 6*** fTol is too small. no further reduction in //
	 *  the sum of squares is possible. //  * ***info = 7*** xTol is too small.
	 * no further improvement in //  the approximate solution x is possible. // 
	 * * ***info = 8*** gTol is too small. fVec is orthogonal to the //  columns
	 * of the jacobian to machine precision. // integer, intent(out) :: nFev  an
	 * integer output variable set to the number of //  calls to fcn with iFlag =
	 * 1. // integer, intent(out) :: nJev  an integer output variable set to the
	 * number of //  calls to fcn with iFlag = 2. // integer, intent(out) ::
	 * ipvt(n)  an integer output array of length n. ipvt //  defines a
	 * permutation matrix p such that jac*p = q*r, //  where jac is the final
	 * calculated jacobian, q is //  orthogonal (not stored), and r is upper
	 * triangular //  with diagonal elements of nonincreasing magnitude. // 
	 * column j of p is column ipvt[j] of the identity matrix. // real(wp),
	 * intent(in) :: fTol  a nonnegative input variable. termination //  occurs
	 * when both the actual and predicted relative //  reductions in the sum of
	 * squares are at most fTol. //  therefore, fTol measures the relative error
	 * desired //  in the sum of squares. // real(wp), intent(in) :: xTol  a
	 * nonnegative input variable. termination //  occurs when the relative error
	 * between two consecutive //  iterates is at most xTol. therefore, xTol
	 * measures the //  relative error desired in the approximate solution. //
	 * real(wp), intent(in) :: gTol  a nonnegative input variable. termination //
	 *  occurs when the cosine of the angle between fVec and //  any column of
	 * the jacobian is at most gTol in absolute //  value. therefore, gTol
	 * measures the orthogonality //  desired between the function vector and the
	 * columns //  of the jacobian. // real(wp), intent(in) :: factor  a
	 * positive input variable used in determining the //  initial step bound.
	 * this bound is set to the product of //  factor and the euclidean norm of
	 * diag*x if nonzero, or else //  to factor itself. in most cases factor
	 * should lie in the //  interval (.1,100.).100. is a generally recommended
	 * value. // real(wp), intent(inout) :: x(n)  an array of length n. on input x
	 * must contain //  an initial estimate of the solution vector. on output x //
	 *  contains the final estimate of the solution vector. // real(wp),
	 * intent(out) :: fVec(m)  an output array of length m which contains // 
	 * the functions evaluated at the output x. // real(wp), intent(out) ::
	 * fJac(ldfJac, n)  an output m by n array. the upper n by n submatrix // 
	 * of fJac contains an upper triangular matrix r with //  diagonal elements of
	 * nonincreasing magnitude such that // ``` //  t t t //  p *(jac *jac)*p
	 * = r *r, // ``` //  where p is a permutation matrix and jac is the final
	 * //  calculated jacobian. column j of p is column ipvt[j] //  (see below)
	 * of the identity matrix. the lower trapezoidal //  part of fJac contains
	 * information generated during //  the computation of r. // real(wp),
	 * intent(inout) :: diag(n)  an array of length n. if mode = 1 (see // 
	 * below), diag is internally set. if mode = 2, diag //  must contain positive
	 * entries that serve as //  multiplicative scale factors for the variables.
	 * // real(wp), intent(out) :: qtf(n)  an output array of length n which
	 * contains //  the first n elements of the vector (q transpose)*fVec. //
	 * real(wp), intent(inout) :: wa1(n)  work array of length n. // real(wp),
	 * intent(inout) :: wa2(n)  work array of length n. // real(wp), intent(inout)
	 * :: wa3(n)  work array of length n. // real(wp), intent(inout) :: wa4(m) 
	 * work array of length m.
	 */
	public static void lmder(SystemOfEquations fcn, int m, int n, double[] x, double[] fVec, double[][] fJac,
			int ldfJac, double fTol, double xTol, double gTol, int maxFev, double[] diag, int mode, double factor,
			int nPrint, int info, int nFev, int nJev, int[] ipvt, double[] qtf, double[] wa1, double[] wa2,
			double[] wa3, double[] wa4) {

//			subroutine lmder(fcn, m, n, x, fVec, fJac, ldfJac, fTol, xTol, gTol, maxFev, &
//                    diag, mode, factor, nPrint, info, nFev, nJev, ipvt, qtf, &
//                    wa1, wa2, wa3, wa4)

		int i, iFlag, iter, j, l;
		double actred, delta = 0, dirder, fNorm, fNorm1, gNorm, par, pnorm, prered, ratio, sum, temp = 0, temp1, temp2,
				xnorm = 0;

		double p1 = 1.0e-1;
		double p5 = 5.0e-1;
		double p25 = 2.5e-1;
		double p75 = 7.5e-1;
		double p0001 = 1.0e-4;

		info = 0;
		iFlag = 0;
		nFev = 0;
		nJev = 0;

		do { // main : block

//	            Check the input parameters for errors.

			if (n > 0 && m >= n && ldfJac >= m && fTol >= 0 && xTol >= 0 && gTol >= 0 && maxFev > 0 && factor > 0) {
				if (mode == 2) {
					for (j = 0; j < n; j++) {
						if (diag[j] <= 0) {
							// throw error exit main
						}
					}
				}
			} else {
				// throw error exit main
			}

//	            Evaluate the function at the starting point
//	            and calculate its norm.

			iFlag = 1;
			fVec = fcn.evaluate(x); // (m, n, x, fVec, fJac, ldfJac, iFlag);
			nFev = 1;
			if (iFlag < 0) {
				// throw error exit main
			}
			fNorm = enorm(m, fVec);

//	            Initialize levenberg-marquardt parameter and iteration counter.

			par = 0;
			iter = 1;

//	            Beginning of the outer loop.

			do { // outer : do

//	                Calculate the jacobian matrix.

				iFlag = 2;
				fVec = fcn.evaluate(x); // fcn(m, n, x, fVec, fJac, ldfJac, iFlag);
				nJev = nJev + 1;
				if (iFlag < 0) {
					// throw error exit main
				}

//	                If requested, call fcn to enable printing of iterates.

				if (nPrint > 0) {
					iFlag = 0;
					if ((iter - 1 % nPrint) == 0) {
						// fcn(m, n, x, fVec, fJac, ldfJac, iFlag);
					}
					if (iFlag < 0) {
						// throw error exit main
					}
				}

//	                Compute the qr factorization of the jacobian.

				qrfac(m, n, fJac, ldfJac, true, ipvt, n, wa1, wa2, wa3);

//	                On the first iteration and if mode is 1, scale according
//	                to the norms of the columns of the initial jacobian.

				if (iter == 1) {
					if (mode != 2) {
						for (j = 0; j < n; j++) {
							diag[j] = wa2[j];
							if (wa2[j] == 0) {
								diag[j] = 1;
							}
						}
					}

//	                    On the first iteration, calculate the norm of the scaled x
//	                    and initialize the step bound delta.

					for (j = 0; j < n; j++) {
						wa3[j] = diag[j] * x[j];
					}
					xnorm = enorm(n, wa3);
					delta = factor * xnorm;
					if (delta == 0) {
						delta = factor;
					}
				}

//	                Form (q transpose)*fVec and store the first n components in
//	                qtf.

				for (i = 0; i < m; i++) {
					wa4[i] = fVec[i];
				}
				for (j = 0; j < n; j++) {
					if (fJac[j][j] != 0) {
						sum = 0;
						for (i = j; i < m; i++) {
							sum = sum + fJac[i][j] * wa4[i];
						}
						temp = -sum / fJac[j][j];
						for (i = j; i < m; i++) {
							wa4[i] = wa4[i] + fJac[i][j] * temp;
						}
					}
					fJac[j][j] = wa1[j];
					qtf[j] = wa4[j];
				}

//	                Compute the norm of the scaled gradient.

				gNorm = 0;
				if (fNorm != 0) {
					for (j = 0; j < n; j++) {
						l = ipvt[j];
						if (wa2[l] != 0) {
							sum = 0;
							for (i = 0; i < j; i++) {
								sum = sum + fJac[i][j] * (qtf[i] / fNorm);
							}
							gNorm = Math.max(gNorm, Math.abs(sum / wa2[l]));
						}
					}
				}

//	                Test for convergence of the gradient norm.

				if (gNorm <= gTol) {
					info = 4;
				}
				if (info != 0) {
					// throw error exit main
				}

//	                Rescale if necessary.

				if (mode != 2) {
					for (j = 0; j < n; j++) {
						diag[j] = Math.max(diag[j], wa2[j]);
					}
				}

//	                Beginning of the inner loop.
				do { // inner : do

//	                    Determine the levenberg-marquardt parameter.

					lmpar(n, fJac, ldfJac, ipvt, diag, qtf, delta, par, wa1, wa2, wa3, wa4);

//	                    Store the direction p and x + p. calculate the norm of p.

					for (j = 0; j < n; j++) {
						wa1[j] = -wa1[j];
						wa2[j] = x[j] + wa1[j];
						wa3[j] = diag[j] * wa1[j];
					}
					pnorm = enorm(n, wa3);

//	                    On the first iteration, adjust the initial step bound.

					if (iter == 1) {
						delta = Math.min(delta, pnorm);
					}

//	                    Evaluate the function at x + p and calculate its norm.

					iFlag = 1;
					wa4 = fcn.evaluate(wa2); // fcn(m, n, wa2, wa4, fJac, ldfJac, iFlag);
					nFev = nFev + 1;
					if (iFlag < 0) {
						// throw error exit main
					}
					fNorm1 = enorm(m, wa4);

//	                    Compute the scaled actual reduction.

					actred = -1;
					if (p1 * fNorm1 < fNorm) {
						actred = 1 - (fNorm1 / fNorm) * (fNorm1 / fNorm);
					}

//	                    Compute the scaled predicted reduction and
//	                    the scaled directional derivative.

					for (j = 0; j < n; j++) {
						wa3[j] = 0;
						l = ipvt[j];
						temp = wa1[l];
						for (i = 0; i < j; i++) {
							wa3[i] = wa3[i] + fJac[i][j] * temp;
						}
					}
					temp1 = enorm(n, wa3) / fNorm;
					temp2 = (Math.sqrt(par) * pnorm) / fNorm;
					prered = temp1 * temp1 + temp2 * temp2 / p5;
					dirder = -(temp1 * temp1 + temp2 * temp2);

//	                    Compute the ratio of the actual to the predicted
//	                    reduction.

					ratio = 0;
					if (prered != 0) {
						ratio = actred / prered;
					}

//	                    Update the step bound.

					if (ratio <= p25) {
						if (actred >= 0) {
							temp = p5;
						}
						if (actred < 0) {
							temp = p5 * dirder / (dirder + p5 * actred);
						}
						if (p1 * fNorm1 >= fNorm || temp < p1) {
							temp = p1;
						}
						delta = temp * Math.min(delta, pnorm / p1);
						par = par / temp;
					} else if (par == 0 || ratio >= p75) {
						delta = pnorm / p5;
						par = p5 * par;
					}

//	                    Test for successful iteration.

					if (ratio >= p0001) {
//	                        Successful iteration. update x, fVec, and their norms.
						for (j = 0; j < n; j++) {
							x[j] = wa2[j];
							wa2[j] = diag[j] * x[j];
						}
						for (i = 0; i < m; i++) {
							fVec[i] = wa4[i];
						}
						xnorm = enorm(n, wa2);
						fNorm = fNorm1;
						iter = iter + 1;
					}

//	                    Tests for convergence.
					if (Math.abs(actred) <= fTol && prered <= fTol && p5 * ratio <= 1) {
						info = 1;
					}
					if (delta <= xTol * xnorm) {
						info = 2;
					}
					if (Math.abs(actred) <= fTol && prered <= fTol && p5 * ratio <= 1 && info == 2) {
						info = 3;
					}
					if (info != 0) {
						// throw error exit main
					}

//	                    Tests for termination and stringent tolerances.
					if (nFev >= maxFev) {
						info = 5;
					}
					if (Math.abs(actred) <= epsmch && prered <= epsmch && p5 * ratio <= 1) {
						info = 6;
					}
					if (delta <= epsmch * xnorm) {
						info = 7;
					}
					if (gNorm <= epsmch) {
						info = 8;
					}
					if (info != 0) {
						// throw error exit main
					}

					if (ratio >= p0001) {
						break; // exit inner
					}

				} while (true); // inner // End of the inner loop. repeat if iteration unsuccessful.

			} while (true); // outer // End of the outer loop

		} while (true); // end block main

//	        Termination, either normal or user imposed.
	}

	/**
	 * //
	 * *****************************************************************************************
	 * // > //  the purpose of lmder1 is to minimize the sum of the squares of //
	 *  m nonlinear functions in n variables by a modification of the // 
	 * levenberg-marquardt algorithm. this is done by using the more //  general
	 * least-squares solver lmder. the user must provide a //  subroutine which
	 * calculates the functions and the jacobian. // // subroutine lmder1(fcn, m, n,
	 * x, fVec, fJac, ldfJac, Tol, info, ipvt, wa, Lwa) // implicit none // //
	 * procedure(fcn_lmder) :: fcn  user-supplied subroutine which // 
	 * calculates the functions and the jacobian. // integer, intent(in) :: m  a
	 * positive integer input variable set to the number //  of functions. //
	 * integer, intent(in) :: n  a positive integer input variable set to the
	 * number //  of variables. n must not exceed m. // integer, intent(in) ::
	 * ldfJac  a positive integer input variable not less than m //  which
	 * specifies the leading dimension of the array fJac. // integer, intent(out) ::
	 * info  an integer output variable. if the user has //  terminated
	 * execution, info is set to the (negative) //  value of iFlag. see
	 * description of fcn. otherwise, //  info is set as follows. //  //  *
	 * ***info = 0*** improper input parameters. //  * ***info = 1*** algorithm
	 * estimates that the relative error //  in the sum of squares is at most tol.
	 * //  * ***info = 2*** algorithm estimates that the relative error // 
	 * between x and the solution is at most tol. //  * ***info = 3*** conditions
	 * for info = 1 and info = 2 both hold. //  * ***info = 4*** fVec is
	 * orthogonal to the columns of the //  jacobian to machine precision. //  *
	 * ***info = 5*** number of calls to fcn with iFlag = 1 has //  reached
	 * 100*(n+1). //  * ***info = 6*** tol is too small. no further reduction in
	 * //  the sum of squares is possible. //  * ***info = 7*** tol is too
	 * small. no further improvement in //  the approximate solution x is
	 * possible. // integer, intent(in) :: Lwa  a positive integer input variable
	 * not less than 5*n+m. // integer, intent(out) :: ipvt(n)  an integer output
	 * array of length n. ipvt //  defines a permutation matrix p such that jac*p
	 * = q*r, //  where jac is the final calculated jacobian, q is // 
	 * orthogonal (not stored), and r is upper triangular //  with diagonal
	 * elements of nonincreasing magnitude. //  column j of p is column ipvt[j] of
	 * the identity matrix. // real(wp), intent(in) :: Tol  a nonnegative input
	 * variable. termination occurs //  when the algorithm estimates either that
	 * the relative //  error in the sum of squares is at most tol or that // 
	 * the relative error between x and the solution is at //  most tol. //
	 * real(wp), intent(inout) :: x(n)  an array of length n. on input x must
	 * contain //  an initial estimate of the solution vector. on output x // 
	 * contains the final estimate of the solution vector. // real(wp), intent(out)
	 * :: fVec(m)  an output array of length m which contains //  the functions
	 * evaluated at the output x. // real(wp), intent(out) :: fJac(ldfJac, n)  an
	 * output m by n array. the upper n by n submatrix //  of fJac contains an
	 * upper triangular matrix r with //  diagonal elements of nonincreasing
	 * magnitude such that // ``` //  t t t //  p *(jac *jac)*p = r *r, //
	 * ``` //  where p is a permutation matrix and jac is the final // 
	 * calculated jacobian. column j of p is column ipvt[j] //  (see below) of the
	 * identity matrix. the lower trapezoidal //  part of fJac contains
	 * information generated during //  the computation of r. // real(wp),
	 * intent(inout) :: wa(Lwa)  a work array of length lwa. // // integer ::
	 * maxFev, mode, nFev, nJev, nPrint // real(wp) :: fTol, gTol, xTol // //
	 * real(wp), parameter :: factor = 100.0_wp
	 */

	public static void lmder1(SystemOfEquations fcn, int m, int n, double[] x, double[] fVec, double[][] fJac,
			int ldfJac, double tol, int info, int[] ipvt, double[] wa, int lwa) {
//			lmder1(fcn, m, n, x, fVec, fJac, ldfJac, Tol, info, ipvt, wa, Lwa)

		int maxFev, mode, nFev = 0, nJev = 0, nPrint;
		double fTol, gTol, xTol;

		double factor = 100.0;

		info = 0;
		double[] wa1 = new double[n];
		double[] wa2 = new double[n];
		double[] wa3 = new double[n];
		double[] wa4 = new double[n];
		double[] wa5 = new double[n];

//	        Check the input parameters for errors.

		if (n > 0 && m >= n && ldfJac >= m && tol >= 0 && lwa >= 5 * n + m) {
//	            Call lmder.
			maxFev = 100 * (n + 1);
			fTol = tol;
			xTol = tol;
			gTol = 0;
			mode = 1;
			nPrint = 0;
			lmder(fcn, m, n, x, fVec, fJac, ldfJac, fTol, xTol, gTol, maxFev, wa, mode, factor, nPrint, info, nFev,
					nJev, ipvt, wa1, wa2, wa3, wa4, wa5);
			if (info == 8) {
				info = 4;
			}
		}

	}

	/**
	 * //
	 * *****************************************************************************************
	 * // > //  the purpose of lmdif is to minimize the sum of the squares of // 
	 * m nonlinear functions in n variables by a modification of //  the
	 * levenberg-marquardt algorithm. the user must provide a //  subroutine which
	 * calculates the functions. the jacobian is //  then calculated by a
	 * forward-difference approximation. // // subroutine lmdif(fcn, m, n, x, fVec,
	 * fTol, xTol, gTol, maxFev, Epsfcn, diag, & // mode, factor, nPrint, info,
	 * nFev, fJac, ldfJac, ipvt, & // qtf, wa1, wa2, wa3, wa4) // implicit none //
	 * // procedure(func2) :: fcn  the user-supplied subroutine which // 
	 * calculates the functions. // integer, intent(in) :: m  a positive integer
	 * input variable set to the number //  of functions. // integer, intent(in)
	 * :: n  a positive integer input variable set to the number //  of
	 * variables. n must not exceed m. // integer, intent(in) :: maxFev  a
	 * positive integer input variable. termination //  occurs when the number of
	 * calls to fcn is at least //  maxFev by the end of an iteration. // integer,
	 * intent(in) :: mode  an integer input variable. if mode = 1, the // 
	 * variables will be scaled internally. if mode = 2, //  the scaling is
	 * specified by the input diag. other //  values of mode are equivalent to
	 * mode = 1. // integer, intent(in) :: nPrint  an integer input variable that
	 * enables controlled //  printing of iterates if it is positive. in this
	 * case, //  fcn is called with iFlag = 0 at the beginning of the first // 
	 * iteration and every nPrint iterations thereafter and //  immediately prior
	 * to return, with x and fVec available //  for printing. if nPrint is not
	 * positive, no special calls //  of fcn with iFlag = 0 are made. // integer,
	 * intent(out) :: info  an integer output variable. if the user has // 
	 * terminated execution, info is set to the (negative) //  value of iFlag. see
	 * description of fcn. otherwise, //  info is set as follows: //  //  *
	 * ***info = 0*** improper input parameters. //  * ***info = 1*** both actual
	 * and predicted relative reductions //  in the sum of squares are at most
	 * fTol. //  * ***info = 2*** relative error between two consecutive iterates
	 * //  is at most xTol. //  * ***info = 3*** conditions for info = 1 and
	 * info = 2 both hold. //  * ***info = 4*** the cosine of the angle between
	 * fVec and any //  column of the jacobian is at most gTol in //  absolute
	 * value. //  * ***info = 5*** number of calls to fcn has reached or // 
	 * exceeded maxFev. //  * ***info = 6*** fTol is too small. no further
	 * reduction in //  the sum of squares is possible. //  * ***info = 7***
	 * xTol is too small. no further improvement in //  the approximate solution x
	 * is possible. //  * ***info = 8*** gTol is too small. fVec is orthogonal to
	 * the //  columns of the jacobian to machine precision. // integer,
	 * intent(out) :: nFev  an integer output variable set to the number of // 
	 * calls to fcn. // integer, intent(in) :: ldfJac  a positive integer input
	 * variable not less than m //  which specifies the leading dimension of the
	 * array fJac. // integer, intent(out) :: ipvt(n)  an integer output array of
	 * length n. ipvt //  defines a permutation matrix p such that jac*p = q*r, //
	 *  where jac is the final calculated jacobian, q is //  orthogonal (not
	 * stored), and r is upper triangular //  with diagonal elements of
	 * nonincreasing magnitude. //  column j of p is column ipvt[j] of the
	 * identity matrix. // real(wp), intent(in) :: fTol  a nonnegative input
	 * variable. termination //  occurs when both the actual and predicted
	 * relative //  reductions in the sum of squares are at most fTol. // 
	 * therefore, fTol measures the relative error desired //  in the sum of
	 * squares. // real(wp), intent(in) :: xTol  a nonnegative input variable.
	 * termination //  occurs when the relative error between two consecutive //
	 *  iterates is at most xTol. therefore, xTol measures the //  relative
	 * error desired in the approximate solution. // real(wp), intent(in) :: gTol 
	 * a nonnegative input variable. termination //  occurs when the cosine of the
	 * angle between fVec and //  any column of the jacobian is at most gTol in
	 * absolute //  value. therefore, gTol measures the orthogonality // 
	 * desired between the function vector and the columns //  of the jacobian. //
	 * real(wp), intent(in) :: Epsfcn  an input variable used in determining a
	 * suitable //  step length for the forward-difference approximation. this //
	 *  approximation assumes that the relative errors in the //  functions are
	 * of the order of epsfcn. if epsfcn is less //  than the machine precision,
	 * it is assumed that the relative //  errors in the functions are of the
	 * order of the machine //  precision. // real(wp), intent(in) :: factor  a
	 * positive input variable used in determining the //  initial step bound.
	 * this bound is set to the product of //  factor and the euclidean norm of
	 * diag*x if nonzero, or else //  to factor itself. in most cases factor
	 * should lie in the //  interval (.1,100.). 100. is a generally recommended
	 * value. // real(wp), intent(inout) :: x(n)  an array of length n. on input x
	 * must contain //  an initial estimate of the solution vector. on output x //
	 *  contains the final estimate of the solution vector. // real(wp),
	 * intent(out) :: fVec(m)  an output array of length m which contains // 
	 * the functions evaluated at the output x. // real(wp), intent(inout) ::
	 * diag(n)  an array of length n. if mode = 1 (see //  below), diag is
	 * internally set. if mode = 2, diag //  must contain positive entries that
	 * serve as //  multiplicative scale factors for the variables. // real(wp),
	 * intent(out) :: fJac(ldfJac, n)  an output m by n array. the upper n by n
	 * submatrix //  of fJac contains an upper triangular matrix r with // 
	 * diagonal elements of nonincreasing magnitude such that // ``` //  t t t
	 * //  p *(jac *jac)*p = r *r, // ``` //  where p is a permutation matrix
	 * and jac is the final //  calculated jacobian. column j of p is column
	 * ipvt[j] //  (see below) of the identity matrix. the lower trapezoidal // 
	 * part of fJac contains information generated during //  the computation of
	 * r. // real(wp), intent(out) :: qtf(n)  an output array of length n which
	 * contains //  the first n elements of the vector (q transpose)*fVec. //
	 * real(wp), intent(inout) :: wa1(n)  work array of length n. // real(wp),
	 * intent(inout) :: wa2(n)  work array of length n. // real(wp), intent(inout)
	 * :: wa3(n)  work array of length n. // real(wp), intent(inout) :: wa4(m) 
	 * work array of length n. // // integer :: i, iFlag, iter, j, l // real(wp) ::
	 * actred, delta, dirder, fNorm, & // fNorm1, gNorm, par, pnorm, prered, & //
	 * ratio, sum, temp, temp1, temp2, xnorm // // real(wp), parameter :: p1 =
	 * 1.0e-1_wp // real(wp), parameter :: p5 = 5.0e-1_wp // real(wp), parameter ::
	 * p25 = 2.5e-1_wp // real(wp), parameter :: p75 = 7.5e-1_wp // real(wp),
	 * parameter :: p0001 = 1.0e-4_wp
	 */
	public static void lmdif(SystemOfEquations fcn, int m, int n, double[] x, double[] fVec, double fTol, double xTol,
			double gTol, int maxFev, double epsfcn, double[] diag, int mode, double factor, int nPrint, int info,
			int nFev, double[][] fJac, int ldfJac, int[] ipvt, double[] qtf, double[] wa1, double[] wa2, double[] wa3,
			double[] wa4) {
		// subroutine lmdif(fcn, m, n, x, fVec, fTol, xTol, gTol, maxFev, Epsfcn, diag,
		// &
//	                     mode, factor, nPrint, info, nFev, fJac, ldfJac, ipvt, &
//	                     qtf, wa1, wa2, wa3, wa4)
		int i, iFlag, iter, j, l;
		double actred, delta = 0, dirder, fNorm, fNorm1, gNorm, par, pnorm, prered, ratio, sum, temp = 0, temp1, temp2,
				xnorm = 0;

		double p1 = 1.0e-1;
		double p5 = 5.0e-1;
		double p25 = 2.5e-1;
		double p75 = 7.5e-1;
		double p0001 = 1.0e-4;

		info = 0;
		iFlag = 0;
		nFev = 0;
		do { // main : block
//	            Check the input parameters for errors.
			if (n > 0 && m >= n && ldfJac >= m && fTol >= 0 && xTol >= 0 && gTol >= 0 && maxFev > 0 && factor > 0) {
				if (mode == 2) {
					for (j = 0; j < n; j++) {
						if (diag[j] <= 0) {
							// throw error exit main
						}
					}
				}
				// throw error exit main
			}
//	            Evaluate the function at the starting point
//	            and calculate its norm.
			iFlag = 1;
			fVec = fcn.evaluate(x);// fcn(m, n, x, fVec, iFlag);
			nFev = 1;
			if (iFlag < 0) {
				// throw error exit main
			}
			fNorm = enorm(m, fVec);
//	            Initialize levenberg-marquardt parameter and iteration counter.
			par = 0;
			iter = 1;
//	            Beginning of the outer loop.
			do { // outer : do
//	                Calculate the jacobian matrix.
				iFlag = 2;
				fdjac2(fcn, m, n, x, fVec, fJac, ldfJac, iFlag, epsfcn, wa4);
				nFev = nFev + n;
				if (iFlag < 0) {
					// throw error exit main
				}
//	                If requested, call fcn to enable printing of iterates.
				if (nPrint > 0) {
					iFlag = 0;
					if (((iter - 1) % nPrint) == 0) {
						// call fcn(m, n, x, fVec, iFlag);
					}
					if (iFlag < 0) {
						// throw error // throw error exit main
					}
				}
//	                Compute the qr factorization of the jacobian.
				qrfac(m, n, fJac, ldfJac, true, ipvt, n, wa1, wa2, wa3);
//	                On the first iteration and if mode is 1, scale according
//	                to the norms of the columns of the initial jacobian.
				if (iter == 1) {
					if (mode != 2) {
						for (j = 0; j < n; j++) {
							diag[j] = wa2[j];
							if (wa2[j] == 0) {
								diag[j] = 1;
							}
						}
					}
//	                    On the first iteration, calculate the norm of the scaled x
//	                    and initialize the step bound delta.
					for (j = 0; j < n; j++) {
						wa3[j] = diag[j] * x[j];
					}
					xnorm = enorm(n, wa3);
					delta = factor * xnorm;
					if (delta == 0) {
						delta = factor;
					}
				}
//	                Form (q transpose)*fVec and store the first n components in
//	                qtf.
				for (i = 0; i < m; i++) {
					wa4[i] = fVec[i];
				}
				for (j = 0; j < n; j++) {
					if (fJac[j][j] != 0) {
						sum = 0;
						for (i = j; i < m; i++) {
							sum = sum + fJac[i][j] * wa4[i];
						}
						temp = -sum / fJac[j][j];
						for (i = j; i < m; i++) {
							wa4[i] = wa4[i] + fJac[i][j] * temp;
						}
					}
					fJac[j][j] = wa1[j];
					qtf[j] = wa4[j];
				}
//	                Compute the norm of the scaled gradient.
				gNorm = 0;
				if (fNorm != 0) {
					for (j = 0; j < n; j++) {
						l = ipvt[j];
						if (wa2[l] != 0) {
							sum = 0;
							for (i = 0; i < j; i++) {
								sum = sum + fJac[i][j] * (qtf[i] / fNorm);
							}
							gNorm = Math.max(gNorm, Math.abs(sum / wa2[l]));
						}
					}
				}
//	                Test for convergence of the gradient norm.
				if (gNorm <= gTol) {
					info = 4;
				}
				if (info != 0) {
					// throw error exit main
				}
//	                Rescale if necessary.
				if (mode != 2) {
					for (j = 0; j < n; j++) {
						diag[j] = Math.max(diag[j], wa2[j]);
					}
				}
//	                Beginning of the inner loop.
				do { // inner : do
//	                    Determine the levenberg-marquardt parameter.
					lmpar(n, fJac, ldfJac, ipvt, diag, qtf, delta, par, wa1, wa2, wa3, wa4);
//	                    Store the direction p and x + p. calculate the norm of p.
					for (j = 0; j < n; j++) {
						wa1[j] = -wa1[j];
						wa2[j] = x[j] + wa1[j];
						wa3[j] = diag[j] * wa1[j];
					}
					pnorm = enorm(n, wa3);
//	                    On the first iteration, adjust the initial step bound.
					if (iter == 1) {
						delta = Math.min(delta, pnorm);
					}
//	                    Evaluate the function at x + p and calculate its norm.
					iFlag = 1;
					wa4 = fcn.evaluate(wa2); // call fcn(m, n, wa2, wa4, iFlag);
					nFev = nFev + 1;
					if (iFlag < 0) {
						// throw error exit main
					}
					fNorm1 = enorm(m, wa4);
//	                    Compute the scaled actual reduction.
					actred = -1;
					if (p1 * fNorm1 < fNorm)
						actred = 1 - (fNorm1 / fNorm) * (fNorm1 / fNorm);
//	                    Compute the scaled predicted reduction and
//	                    the scaled directional derivative.
					for (j = 0; j < n; j++) {
						wa3[j] = 0;
						l = ipvt[j];
						temp = wa1[l];
						for (i = 0; i < j; i++) {
							wa3[i] = wa3[i] + fJac[i][j] * temp;
						}
					}
					temp1 = enorm(n, wa3) / fNorm;
					temp2 = (Math.sqrt(par) * pnorm) / fNorm;
					prered = temp1 * temp1 + temp2 * temp2 / p5;
					dirder = -(temp1 * temp1 + temp2 * temp2);
//	                    Compute the ratio of the actual to the predicted
//	                    reduction.
					ratio = 0;
					if (prered != 0)
						ratio = actred / prered;
//	                    Update the step bound.
					if (ratio <= p25) {
						if (actred >= 0) {
							temp = p5;
						}
						if (actred < 0) {
							temp = p5 * dirder / (dirder + p5 * actred);
						}
						if (p1 * fNorm1 >= fNorm || temp < p1) {
							temp = p1;
						}
						if (delta < pnorm / p1) {
							delta = temp * delta;
						} else {
							delta = temp * pnorm / p1;
						}
						par = par / temp;
					} else if (par == 0 || ratio >= p75) {
						delta = pnorm / p5;
						par = p5 * par;
					}
//	                    Test for successful iteration.
					if (ratio >= p0001) {
//	                        Successful iteration. update x, fVec, and their norms.
						for (j = 0; j < n; j++) {
							x[j] = wa2[j];
							wa2[j] = diag[j] * x[j];
						}
						for (i = 0; i < m; i++) {
							fVec[i] = wa4[i];
						}
						xnorm = enorm(n, wa2);
						fNorm = fNorm1;
						iter = iter + 1;
					}
//	                    Tests for convergence.

					if (Math.abs(actred) <= fTol && prered <= fTol && p5 * ratio <= 1) {
						info = 1;
					}
					if (delta <= xTol * xnorm) {
						info = 2;
					}
					if (Math.abs(actred) <= fTol && prered <= fTol && p5 * ratio <= 1 && info == 2) {
						info = 3;
					}
					if (info != 0) // throw error exit main
//	                    Tests for termination and stringent tolerances.
						if (nFev >= maxFev) {
							info = 5;
						}
					if (Math.abs(actred) <= epsmch && prered <= epsmch && p5 * ratio <= 1) {
						info = 6;
					}
					if (delta <= epsmch * xnorm) {
						info = 7;
					}
					if (gNorm <= epsmch) {
						info = 8;
					}
					if (info != 0) {
						// throw error exit main
					}
					if (ratio >= p0001) {
						break; // exit inner
					}
				} while (true); // } inner  end of the inner loop. repeat if iteration unsuccessful.
			} while (true); // } outer  end of the outer loop.
		} while (true); // end block main
//	         termination, either normal or user imposed.
//	        if (iFlag < 0) info = iFlag
//	        iFlag = 0
//	        if (nPrint > 0) call fcn(m, n, x, fVec, iFlag)
		//
	}

	/**
	 * // //
	 * *****************************************************************************************
	 * // > //  the purpose of lmdif1 is to minimize the sum of the squares of //
	 *  m nonlinear functions in n variables by a modification of the // 
	 * levenberg-marquardt algorithm. this is done by using the more //  general
	 * least-squares solver lmdif. the user must provide a //  subroutine which
	 * calculates the functions. the jacobian is //  then calculated by a
	 * forward-difference approximation. // // // // procedure(func2) :: fcn  the
	 * user-supplied subroutine which //  calculates the functions. // integer,
	 * intent(in) :: m  a positive integer input variable set to the number // 
	 * of functions. // integer, intent(in) :: n  a positive integer input
	 * variable set to the number //  of variables. n must not exceed m. //
	 * integer, intent(out) :: info  an integer output variable. if the user has
	 * //  terminated execution, info is set to the (negative) //  value of
	 * iFlag. see description of fcn. otherwise, //  info is set as follows: // 
	 * //  * ***info = 0*** improper input parameters. //  * ***info = 1***
	 * algorithm estimates that the relative error //  in the sum of squares is at
	 * most tol. //  * ***info = 2*** algorithm estimates that the relative error
	 * //  between x and the solution is at most tol. //  * ***info = 3***
	 * conditions for info = 1 and info = 2 both hold. //  * ***info = 4*** fVec
	 * is orthogonal to the columns of the //  jacobian to machine precision. //
	 *  * ***info = 5*** number of calls to fcn has reached or //  exceeded
	 * 200*(n+1). //  * ***info = 6*** tol is too small. no further reduction in
	 * //  the sum of squares is possible. //  * ***info = 7*** tol is too
	 * small. no further improvement in //  the approximate solution x is
	 * possible. // integer, intent(in) :: Lwa  a positive integer input variable
	 * not less than //  m*n+5*n+m. // integer, intent(inout) :: Iwa(n)  an
	 * integer work array of length n. // real(wp), intent(in) :: Tol  a
	 * nonnegative input variable. termination occurs //  when the algorithm
	 * estimates either that the relative //  error in the sum of squares is at
	 * most tol or that //  the relative error between x and the solution is at //
	 *  most tol. // real(wp), intent(inout) :: x(n)  an array of length n. on
	 * input x must contain //  an initial estimate of the solution vector. on
	 * output x //  contains the final estimate of the solution vector. //
	 * real(wp), intent(out) :: fVec(m)  an output array of length m which
	 * contains //  the functions evaluated at the output x. // real(wp),
	 * intent(inout) :: wa(Lwa)  a work array of length lwa.
	 */
	public static void lmdif1(SystemOfEquations fcn, int m, int n, double[] x, double[] fVec, double tol, int info,
			int[] iwa, double[] wa, int lwa) {
		// subroutine lmdif1(fcn, m, n, x, fVec, Tol, info, Iwa, wa, Lwa)
//	        implicit none

		int maxFev, mode, nFev, nPrint;
		double epsfcn, fTol, gTol, xTol;

		double factor = 1.0e2;
		double[] wa1 = new double[n];
		double[] wa2 = new double[n];
		double[] wa3 = new double[n];
		double[] wa4 = new double[n];
		double[] wa5 = new double[n];
		double[][] fjac = new double[n][m];

		info = 0;

//	        Check the input parameters for errors.

		if (n > 0 && m >= n && tol >= 0 && lwa >= m * n + 5 * n + m) {

//	            Call lmdif.

			maxFev = 200 * (n + 1);
			nFev = 0;
			fTol = tol;
			xTol = tol;
			gTol = 0;
			epsfcn = 0;
			mode = 1;
			nPrint = 0;
			lmdif(fcn, m, n, x, fVec, fTol, xTol, gTol, maxFev, epsfcn, wa, mode, factor, nPrint, info, nFev, fjac, m,
					iwa, wa1, wa2, wa3, wa4, wa5);
			if (info == 8) {
				info = 4;
			}
		}

	}

	/**
	 * //
	 * *****************************************************************************************
	 * // > //  given an m by n matrix a, an n by n nonSingular diagonal // 
	 * matrix d, an m-vector b, and a positive number delta, //  the problem is to
	 * determine a value for the parameter //  par such that if x solves the system
	 * // ``` //  a*x = b , sqrt(par)*d*x = 0 , // ``` //  in the least squares
	 * sense, and dxnorm is the euclidean //  norm of d*x, then either par is zero
	 * and // ``` //  (dxnorm-delta) <= 0.1*delta , // ``` //  or par is
	 * positive and // ``` //  abs(dxnorm-delta) <= 0.1*delta . // ``` //  this
	 * subroutine completes the solution of the problem //  if it is provided with
	 * the necessary information from the //  qr factorization, with column
	 * pivoting, of a. that is, if //  a*p = q*r, where p is a permutation matrix,
	 * q has orthogonal //  columns, and r is an upper triangular matrix with
	 * diagonal //  elements of nonincreasing magnitude, then lmpar expects // 
	 * the full upper triangle of r, the permutation matrix p, //  and the first n
	 * components of (q transpose)*b. on output //  lmpar also provides an upper
	 * triangular matrix s such that // ``` //  t t t //  p *(a *a + par*d*d)*p =
	 * s *s . // ``` //  s is employed within lmpar and may be of separate
	 * interest. //  //  only a few iterations are generally needed for
	 * convergence //  of the algorithm. if, however, the limit of 10 iterations //
	 *  is reached, then the output par will contain the best //  value obtained
	 * so far. // // subroutine lmpar(n, r, Ldr, ipvt, diag, Qtb, Delta, Par, x,
	 * Sdiag, wa1, wa2) // implicit none // // integer, intent(in) :: n  a
	 * positive integer input variable set to the order of r. // integer, intent(in)
	 * :: Ldr  a positive integer input variable not less than n //  which
	 * specifies the leading dimension of the array r. // integer, intent(in) ::
	 * ipvt(n)  an integer input array of length n which defines the // 
	 * permutation matrix p such that a*p = q*r. column j of p //  is column
	 * ipvt[j] of the identity matrix. // real(wp) :: Delta  a positive input
	 * variable which specifies an upper //  bound on the euclidean norm of d*x.
	 * // real(wp), intent(inout) :: Par  a nonnegative variable. on input par
	 * contains an //  initial estimate of the levenberg-marquardt parameter. //
	 *  on output par contains the final estimate. // real(wp), intent(inout) ::
	 * r(Ldr, n)  an n by n array. on input the full upper triangle //  must
	 * contain the full upper triangle of the matrix r. //  on output the full
	 * upper triangle is unaltered, and the //  strict lower triangle contains the
	 * strict upper triangle //  (transposed) of the upper triangular matrix s. //
	 * real(wp), intent(in) :: diag(n)  an input array of length n which must
	 * contain the //  diagonal elements of the matrix d. // real(wp), intent(in)
	 * :: Qtb(n)  an input array of length n which must contain the first //  n
	 * elements of the vector (q transpose)*b. // real(wp), intent(out) :: x(n) 
	 * an output array of length n which contains the least //  squares solution
	 * of the system a*x = b, sqrt(par)*d*x = 0, //  for the output par. //
	 * real(wp), intent(out) :: Sdiag(n)  an output array of length n which
	 * contains the //  diagonal elements of the upper triangular matrix s. //
	 * real(wp), intent(inout) :: wa1(n)  work array of length n. // real(wp),
	 * intent(inout) :: wa2(n)  work array of length n. //
	 */
	public static void lmpar(int n, double[][] r, int ldr, int[] ipvt, double[] diag, double[] qtb, double delta,
			double par, double[] x, double[] sdiag, double[] wa1, double[] wa2) {
		// lmpar(n, r, Ldr, ipvt, diag, Qtb, Delta, Par, x, Sdiag, wa1, wa2)

		int i, iter, j, jm1, jp1, k, l, nSing;
		double dxnorm, fp, gNorm, parc, parl, paru, sum, temp;

		double p1 = 1.0e-1;
		double p001 = 1.0e-3;
		double dwarf = 1.0e-20; // This is mostly a guess, may need adjusting //dpmpar[2]; // the smallest
								// positive magnitude

//	        Compute and store in x the gauss-newton direction. if the
//	        jacobian is rank-deficient, obtain a least squares solution.

		nSing = n;
		for (j = 0; j < n; j++) {
			wa1[j] = qtb[j];
			if (r[j][j] == 0 && nSing == n) {
				nSing = j - 1;// might need fixing
			}
			if (nSing < n) {
				wa1[j] = 0;
			}
		}
		if (nSing >= 1) {
			for (k = 0; k < nSing; k++) {
				j = nSing - k + 1; // might need fixing
				wa1[j] = wa1[j] / r[j][j];
				temp = wa1[j];
				jm1 = j - 1;
				if (jm1 >= 0) {
					for (i = 0; i < jm1; i++) {
						wa1[i] = wa1[i] - r[i][j] * temp;
					}
				}
			}
		}
		for (j = 0; j < n; j++) {
			l = ipvt[j];
			x[l] = wa1[j];
		}

//	        Initialize the iteration counter.
//	        evaluate the function at the origin, and test
//	        for acceptance of the gauss-newton direction.

		iter = 0;
		for (j = 0; j < n; j++) {
			wa2[j] = diag[j] * x[j];
		}
		dxnorm = enorm(n, wa2);
		fp = dxnorm - delta;
		if (fp <= p1 * delta) {
//	            Termination.
			if (iter == 0) {
				par = 0;
			}
		} else {

//	            If the jacobian is not rank deficient, the newton
//	            step provides a lower bound, parl, for the zero of
//	            the function. otherwise set this bound to zero.

			parl = 0;
			if (nSing >= n) {
				for (j = 0; j < n; j++) {
					l = ipvt[j];
					wa1[j] = diag[l] * (wa2[l] / dxnorm);
				}
				for (j = 0; j < n; j++) {
					sum = 0;
					jm1 = j - 1;
					if (jm1 >= 0) {
						for (i = 0; i < jm1; i++) {
							sum = sum + r[i][j] * wa1[i];
						}
					}
					wa1[j] = (wa1[j] - sum) / r[j][j];
				}
				temp = enorm(n, wa1);
				parl = ((fp / delta) / temp) / temp;
			}

//	            Calculate an upper bound, paru, for the zero of the function.

			for (j = 0; j < n; j++) {
				sum = 0;
				for (i = 0; i < j; i++) {
					sum = sum + r[i][j] * qtb[i];
				}
				l = ipvt[j];
				wa1[j] = sum / diag[l];
			}
			gNorm = enorm(n, wa1);
			paru = gNorm / delta;
			if (paru == 0) {
				if (delta > p1) {
					paru = dwarf / p1;
				} else {
					paru = dwarf / delta;
				}
			}
//	            If the input par lies outside of the interval (parl,paru),
//	            set par to the closer endpoint.
			if (parl > par) {
				par = parl;
			}
			if (par > paru) {
				par = paru;
			}
			if (par == 0) {
				par = gNorm / dxnorm;
			}

//	            Beginning of an iteration.
			do {

				iter = iter + 1;

//	                Evaluate the function at the current value of par.

				if (par == 0) {
					if (dwarf > p001 * paru) {
						par = dwarf;
					} else {
						par = p001 * paru;
					}
				}
				temp = Math.sqrt(par);
				for (j = 0; j < n; j++) {
					wa1[j] = temp * diag[j];
				}
				qrsolv(n, r, ldr, ipvt, wa1, qtb, x, sdiag, wa2);
				for (j = 0; j < n; j++) {
					wa2[j] = diag[j] * x[j];
				}
				dxnorm = enorm(n, wa2);
				temp = fp;
				fp = dxnorm - delta;

//	                If the function is small enough, accept the current value
//	                of par. also test for the exceptional cases where parl
//	                is zero or the number of iterations has reached 10.

				if (Math.abs(fp) <= p1 * delta || parl == 0 && fp <= temp && temp < 0 || iter == 10) {
					if (iter == 0) {
						par = 0;
					}
					// exit maybe a break or continue?
				} else {

//	                    Compute the Newton correction.

					for (j = 0; j < n; j++) {
						l = ipvt[j];
						wa1[j] = diag[l] * (wa2[l] / dxnorm);
					}
					for (j = 0; j < n; j++) {
						wa1[j] = wa1[j] / sdiag[j];
						temp = wa1[j];
						jp1 = j + 1;
						if (n >= jp1) {
							for (i = jp1; i < n; i++) {
								wa1[i] = wa1[i] - r[i][j] * temp;
							}
						}
					}
					temp = enorm(n, wa1);
					parc = ((fp / delta) / temp) / temp;

//	                     depending on the sign of the function, update parl or paru.

					if (fp > 0) {
						if (par > parl) {
							parl = par;
						}
					}
					if (fp < 0) {
						if (par < paru) {
							paru = par;
						}
					}
//	                    Compute an improved estimate for par.
					if (parl > par + parc) {
						par = parl;
					} else {
						par = par + parc;
					}

				}

			} while (true); // end of an iteration.

		}
		while (true)
			;
	}// while(true);

	/**
	 * //
	 * *****************************************************************************************
	 * // > //  the purpose of lmstr is to minimize the sum of the squares of // 
	 * m nonlinear functions in n variables by a modification of //  the
	 * levenberg-marquardt algorithm which uses minimal storage. //  the user must
	 * provide a subroutine which calculates the //  functions and the rows of the
	 * jacobian. // // // implicit none // // procedure(fcn_lmstr) :: fcn 
	 * user-supplied subroutine which //  calculates the functions and the rows of
	 * the jacobian. // integer, intent(in) :: m  a positive integer input
	 * variable set to the number //  of functions. // integer, intent(in) :: n 
	 * a positive integer input variable set to the number //  of variables. n
	 * must not exceed m. // integer, intent(in) :: ldfJac  a positive integer
	 * input variable not less than n //  which specifies the leading dimension of
	 * the array fJac. // integer, intent(in) :: maxFev  a positive integer input
	 * variable. termination //  occurs when the number of calls to fcn with iFlag
	 * = 1 //  has reached maxFev. // integer, intent(in) :: mode  an integer
	 * input variable. if mode = 1, the //  variables will be scaled internally.
	 * if mode = 2, //  the scaling is specified by the input diag. other // 
	 * values of mode are equivalent to mode = 1. // integer, intent(in) :: nPrint
	 *  an integer input variable that enables controlled //  printing of
	 * iterates if it is positive. in this case, //  fcn is called with iFlag = 0
	 * at the beginning of the first //  iteration and every nPrint iterations
	 * thereafter and //  immediately prior to return, with x and fVec available
	 * //  for printing. if nPrint is not positive, no special calls //  of fcn
	 * with iFlag = 0 are made. // integer, intent(out) :: info  an integer output
	 * variable. if the user has //  terminated execution, info is set to the
	 * (negative) //  value of iFlag. see description of fcn. otherwise, // 
	 * info is set as follows: //  //  * ***info = 0*** improper input
	 * parameters. //  * ***info = 1*** both actual and predicted relative
	 * reductions //  in the sum of squares are at most fTol. //  * ***info =
	 * 2*** relative error between two consecutive iterates //  is at most xTol.
	 * //  * ***info = 3*** conditions for info = 1 and info = 2 both hold. // 
	 * * ***info = 4*** the cosine of the angle between fVec and any //  column of
	 * the jacobian is at most gTol in //  absolute value. //  * ***info = 5***
	 * number of calls to fcn with iFlag = 1 has //  reached maxFev. //  *
	 * ***info = 6*** fTol is too small. no further reduction in //  the sum of
	 * squares is possible. //  * ***info = 7*** xTol is too small. no further
	 * improvement in //  the approximate solution x is possible. //  * ***info
	 * = 8*** gTol is too small. fVec is orthogonal to the //  columns of the
	 * jacobian to machine precision. // integer, intent(out) :: nFev  an integer
	 * output variable set to the number of //  calls to fcn with iFlag = 1. //
	 * integer, intent(out) :: nJev  an integer output variable set to the number
	 * of //  calls to fcn with iFlag = 2. // integer, intent(out) :: ipvt(n) 
	 * an integer output array of length n. ipvt //  defines a permutation matrix
	 * p such that jac*p = q*r, //  where jac is the final calculated jacobian, q
	 * is //  orthogonal (not stored), and r is upper triangular. //  column j
	 * of p is column ipvt[j] of the identity matrix. // real(wp), intent(in) ::
	 * fTol  a nonnegative input variable. termination //  occurs when both the
	 * actual and predicted relative //  reductions in the sum of squares are at
	 * most fTol. //  therefore, fTol measures the relative error desired //  in
	 * the sum of squares. // real(wp), intent(in) :: xTol  a nonnegative input
	 * variable. termination //  occurs when the relative error between two
	 * consecutive //  iterates is at most xTol. therefore, xTol measures the //
	 *  relative error desired in the approximate solution. // real(wp),
	 * intent(in) :: gTol  a nonnegative input variable. termination //  occurs
	 * when the cosine of the angle between fVec and //  any column of the
	 * jacobian is at most gTol in absolute //  value. therefore, gTol measures
	 * the orthogonality //  desired between the function vector and the columns
	 * //  of the jacobian. // real(wp), intent(in) :: factor  a positive input
	 * variable used in determining the //  initial step bound. this bound is set
	 * to the product of //  factor and the euclidean norm of diag*x if nonzero,
	 * or else //  to factor itself. in most cases factor should lie in the // 
	 * interval (.1,100.). 100. is a generally recommended value. // real(wp),
	 * intent(inout) :: x(n)  an array of length n. on input x must contain // 
	 * an initial estimate of the solution vector. on output x //  contains the
	 * final estimate of the solution vector. // real(wp), intent(out) :: fVec(m) 
	 * an output array of length m which contains //  the functions evaluated at
	 * the output x. // real(wp), intent(out) :: fJac(ldfJac, n)  an output n by n
	 * array. the upper triangle of fJac //  contains an upper triangular matrix r
	 * such that // ``` //  t t t //  p *(jac *jac)*p = r *r, // ``` // 
	 * where p is a permutation matrix and jac is the final //  calculated
	 * jacobian. column j of p is column ipvt[j] //  (see below) of the identity
	 * matrix. the lower triangular //  part of fJac contains information
	 * generated during //  the computation of r. // real(wp), intent(inout) ::
	 * diag(n)  an array of length n. if mode = 1 (see //  below), diag is
	 * internally set. if mode = 2, diag //  must contain positive entries that
	 * serve as //  multiplicative scale factors for the variables. // real(wp),
	 * intent(out) :: qtf(n)  an output array of length n which contains //  the
	 * first n elements of the vector (q transpose)*fVec. // real(wp), intent(inout)
	 * :: wa1(n)  work array of length n. // real(wp), intent(inout) :: wa2(n) 
	 * work array of length n. // real(wp), intent(inout) :: wa3(n)  work array of
	 * length n. // real(wp), intent(inout) :: wa4(m)  work array of length m. //
	 */
	public static void lmstr(SystemOfEquations fcn, int m, int n, double[] x, double[] fVec, double[][] fJac,
			int ldfJac, double fTol, double xTol, double gTol, int maxFev, double[] diag, int mode, double factor,
			int nPrint, int info, int nFev, int nJev, int[] ipvt, double[] qtf, double[] wa1, double[] wa2,
			double[] wa3, double[] wa4) {
		// subroutine lmstr(fcn, m, n, x, fVec, fJac, ldfJac, fTol, xTol, gTol, maxFev,
		// &
//	                     diag, mode, factor, nPrint, info, nFev, nJev, ipvt, qtf, &
//	                     wa1, wa2, wa3, wa4)

		int i, iFlag, iter, j, l;
		double actred, delta = 0, dirder, fNorm, fNorm1, gNorm, par, pnorm, prered, ratio, sum, temp = 0, temp1, temp2,
				xnorm = 0;
		boolean sing;

		double p1 = 1.0e-1;
		double p5 = 5.0e-1;
		double p25 = 2.5e-1;
		double p75 = 7.5e-1;
		double p0001 = 1.0e-4;

		info = 0;
		iFlag = 0;
		nFev = 0;
		nJev = 0;

		do { // main : block

//	            Check the input parameters for errors.

			if (n <= 0 || m < n || ldfJac < n || fTol < 0 || xTol < 0 || gTol < 0 || maxFev <= 0 || factor <= 0) {
				// thorw.exception();
			}
			if (mode == 2) {
				for (j = 0; j < n; j++) {
					if (diag[j] <= 0) {
						// throw error
					}
				}
			}

//	            Evaluate the function at the starting point
//	            and calculate its norm.

			iFlag = 1;
			fVec = fcn.evaluate(x); // fcn(m, n, x, fVec, wa3, iFlag)
			nFev = 1;
			if (iFlag < 0) {
				// throw error exit main
			}
			fNorm = enorm(m, fVec);

//	            Initialize levenberg-marquardt parameter and iteration counter.

			par = 0;
			iter = 1;

//	            Beginning of the outer loop.

			do { // outer : do

//	                If requested, call fcn to enable printing of iterates.

				if (nPrint > 0) {
					iFlag = 0;
					if (((iter - 1) % nPrint) == 0) {
						// fcn(m, n, x, fVec, wa3, iFlag)
					}
					if (iFlag < 0) {
						// throw error exit main
					}
				}

//	                Compute the qr factorization of the jacobian matrix
//	                calculated one row at a time, while simultaneously
//	                forming (q transpose)*fVec and storing the first
//	                n components in qtf.
				for (j = 0; j < n; j++) {
					qtf[j] = 0;
					for (i = 0; i < n; i++) {
						fJac[i][j] = 0;
					}
				}
				iFlag = 2;
				for (i = 0; i < m; i++) {
					fVec = fcn.evaluate(x); // fcn(m, n, x, fVec, wa3, iFlag)
					if (iFlag < 0) {
						// throw error exit main
					}
					temp = fVec[i];
					rwupdt(n, fJac, ldfJac, wa3, qtf, temp, wa1, wa2);
					iFlag = iFlag + 1;
				}
				nJev = nJev + 1;
				//
//	                If the jacobian is rank deficient, call qrfac to
//	                reorder its columns and update the components of qtf.

				sing = false;
				for (j = 0; j < n; j++) {
					if (fJac[j][j] == 0) {
						sing = true;
					}
					ipvt[j] = j;
					double[] fJacTemp = new double[j];
					for (i = 0; i < j; i++) {
						fJacTemp[i] = fJac[0][i];
					}
					wa2[j] = enorm(j, fJacTemp); // should be correct? Wa2(j) = enorm(j, Fjac(1, j)) originally this,
													// with the loop starting at 1
				}
				if (sing) {
					qrfac(n, n, fJac, ldfJac, true, ipvt, n, wa1, wa2, wa3);
					for (j = 0; j < n; j++) {
						if (fJac[j][j] != 0) {
							sum = 0;
							for (i = j; i < n; i++) {
								sum = sum + fJac[i][j] * qtf[i];
							}
							temp = -sum / fJac[j][j];
							for (i = j; i < n; i++) {
								qtf[i] = qtf[i] + fJac[i][j] * temp;
							}
						}
						fJac[j][j] = wa1[j];
					}
				}

//	                On the first iteration and if mode is 1, scale according
//	                to the norms of the columns of the initial jacobian.

				if (iter == 1) {
					if (mode != 2) {
						for (j = 0; j < n; j++) {
							diag[j] = wa2[j];
							if (wa2[j] == 0) {
								diag[j] = 1;
							}
						}
					}

//	                    On the first iteration, calculate the norm of the scaled x
//	                    and initialize the step bound delta.

					for (j = 0; j < n; j++) {
						wa3[j] = diag[j] * x[j];
					}
					xnorm = enorm(n, wa3);
					delta = factor * xnorm;
					if (delta == 0) {
						delta = factor;
					}
				}

//	                Compute the norm of the scaled gradient.

				gNorm = 0;
				if (fNorm != 0) {
					for (j = 0; j < n; j++) {
						l = ipvt[j];
						if (wa2[l] != 0) {
							sum = 0;
							for (i = 0; i < j; i++) {
								sum = sum + fJac[i][j] * (qtf[i] / fNorm);
							}
							if (Math.abs(sum / wa2[l]) > gNorm) {
								gNorm = Math.abs(sum / wa2[l]);
							}
						}
					}
				}

//	                Test for convergence of the gradient norm.

				if (gNorm <= gTol) {
					info = 4;
				}
				if (info != 0) {
					// throw error exit
				}

//	                Rescale if necessary.

				if (mode != 2) {
					for (j = 0; j < n; j++) {
						if (wa2[j] > diag[j]) {
							diag[j] = wa2[j];
						}
					}
				}

//	                Beginning of the inner loop.

				do { // inner : do

//	                    Determine the levenberg-marquardt parameter.
					//
					lmpar(n, fJac, ldfJac, ipvt, diag, qtf, delta, par, wa1, wa2, wa3, wa4);
					//
//	                    Store the direction p and x + p. calculate the norm of p.
					//
					for (j = 0; j < n; j++) {
						wa1[j] = -wa1[j];
						wa2[j] = x[j] + wa1[j];
						wa3[j] = diag[j] * wa1[j];
					}
					pnorm = enorm(n, wa3);

//	                    On the first iteration, adjust the initial step bound.

					if (iter == 1) {
						if (delta > pnorm) {
							delta = pnorm;
						}
					}
//	                    Evaluate the function at x + p and calculate its norm.
					//
					iFlag = 1;
					wa4 = fcn.evaluate(x); // call fcn(m, n, wa2, wa4, wa3, iFlag)
					nFev = nFev + 1;
					if (iFlag < 0) {
						/// exit main
					}
					fNorm1 = enorm(m, wa4);

//	                    Compute the scaled actual reduction.

					actred = -1;
					if (p1 * fNorm1 < fNorm)
						actred = 1 - (fNorm1 / fNorm) * (fNorm1 / fNorm);

//	                    Compute the scaled predicted reduction and
//	                    the scaled directional derivative.

					for (j = 0; j < n; j++) {
						wa3[j] = 0;
						l = ipvt[j];
						temp = wa1[l];
						for (i = 0; i < j; i++) {
							wa3[i] = wa3[i] + fJac[i][j] * temp;
						}
					}
					temp1 = enorm(n, wa3) / fNorm;
					temp2 = (Math.sqrt(par) * pnorm) / fNorm;
					prered = temp1 * temp1 + temp2 * temp2 / p5;
					dirder = -(temp1 * temp1 + temp2 * temp2);

//	                    Compute the ratio of the actual to the predicted
//	                    reduction.

					ratio = 0;
					if (prered != 0) {
						ratio = actred / prered;
					}
//	                    Update the step bound.

					if (ratio <= p25) {
						if (actred >= 0) {
							temp = p5;
						}
						if (actred < 0) {
							temp = p5 * dirder / (dirder + p5 * actred);
						}
						if (p1 * fNorm1 >= fNorm || temp < p1) {
							temp = p1;
						}
						if (delta > pnorm / p1) {
							delta = temp * pnorm / p1;
						} else {
							delta = delta * temp;
						}
						par = par / temp;
					} else if (par == 0 || ratio >= p75) {
						delta = pnorm / p5;
						par = p5 * par;
					}

//	                    Test for successful iteration.

					if (ratio >= p0001) {

//	                        Successful iteration. update x, fVec, and their norms.

						for (j = 0; j < n; j++) {
							x[j] = wa2[j];
							wa2[j] = diag[j] * x[j];
						}
						for (i = 0; i < m; i++) {
							fVec[i] = wa4[i];
						}
						xnorm = enorm(n, wa2);
						fNorm = fNorm1;
						iter = iter + 1;
					}

//	                    Tests for convergence.

					if (Math.abs(actred) <= fTol && prered <= fTol && p5 * ratio <= 1) {
						info = 1;
					}
					if (delta <= xTol * xnorm) {
						info = 2;
					}
					if (Math.abs(actred) <= fTol && prered <= fTol && p5 * ratio <= 1 && info == 2) {
						info = 3;
					}
					if (info != 0) {
						// exit main
					}
					// Tests for termination and stringent tolerances.

					if (nFev >= maxFev) {
						info = 5;
					}
					if (Math.abs(actred) <= epsmch && prered <= epsmch && p5 * ratio <= 1) {
						info = 6;
					}
					if (delta <= epsmch * xnorm) {
						info = 7;
					}
					if (gNorm <= epsmch) {
						info = 8;
					}
					if (info != 0) {
						// exit main
					}
					if (ratio >= p0001) {
						// exit inner
					}

				} while (true); // } inner  end of the inner loop. repeat if iteration unsuccessful.
			} while (true); // //} outer  end of the outer loop.
		} while (true); // //end block main
//	         termination, either normal or user imposed.
		//
//	        if (iFlag < 0) info = iFlag
//	        iFlag = 0
//	        if (nPrint > 0) call fcn(m, n, x, fVec, wa3, iFlag)
		//
	}

	/**
	 * //
	 * *****************************************************************************************
	 * // > //  the purpose of lmstr1 is to minimize the sum of the squares of //
	 *  m nonlinear functions in n variables by a modification of //  the
	 * levenberg-marquardt algorithm which uses minimal storage. //  this is done
	 * by using the more general least-squares solver //  lmstr. the user must
	 * provide a subroutine which calculates //  the functions and the rows of the
	 * jacobian. // // subroutine lmstr1(fcn, m, n, x, fVec, fJac, ldfJac, Tol,
	 * info, ipvt, wa, Lwa) // implicit none // // procedure(fcn_lmstr) :: fcn 
	 * user-supplied subroutine which //  calculates the functions and the rows of
	 * the jacobian. // integer, intent(in) :: m  a positive integer input
	 * variable set to the number //  of functions. // integer, intent(in) :: n 
	 * a positive integer input variable set to the number //  of variables. n
	 * must not exceed m. // integer, intent(in) :: ldfJac  a positive integer
	 * input variable not less than n //  which specifies the leading dimension of
	 * the array fJac. // integer, intent(out) :: info  an integer output
	 * variable. if the user has //  terminated execution, info is set to the
	 * (negative) //  value of iFlag. see description of fcn. otherwise, // 
	 * info is set as follows: //  //  * ***info = 0*** improper input
	 * parameters. //  * ***info = 1*** algorithm estimates that the relative
	 * error //  in the sum of squares is at most tol. //  * ***info = 2***
	 * algorithm estimates that the relative error //  between x and the solution
	 * is at most tol. //  * ***info = 3*** conditions for info = 1 and info = 2
	 * both hold. //  * ***info = 4*** fVec is orthogonal to the columns of the //
	 *  jacobian to machine precision. //  * ***info = 5*** number of calls to
	 * fcn with iFlag = 1 has //  reached 100*(n+1). //  * ***info = 6*** tol is
	 * too small. no further reduction in //  the sum of squares is possible. //
	 *  * ***info = 7*** tol is too small. no further improvement in //  the
	 * approximate solution x is possible. // integer, intent(in) :: Lwa  a
	 * positive integer input variable not less than 5*n+m. // integer, intent(out)
	 * :: ipvt(n)  an integer output array of length n. ipvt //  defines a
	 * permutation matrix p such that jac*p = q*r, //  where jac is the final
	 * calculated jacobian, q is //  orthogonal (not stored), and r is upper
	 * triangular. //  column j of p is column ipvt[j] of the identity matrix. //
	 * real(wp), intent(in) :: Tol  a nonnegative input variable. termination
	 * occurs //  when the algorithm estimates either that the relative // 
	 * error in the sum of squares is at most tol or that //  the relative error
	 * between x and the solution is at //  most tol. // real(wp), intent(inout)
	 * :: x(n)  an array of length n. on input x must contain //  an initial
	 * estimate of the solution vector. on output x //  contains the final
	 * estimate of the solution vector. // real(wp), intent(out) :: fVec(m)  an
	 * output array of length m which contains //  the functions evaluated at the
	 * output x. // real(wp), intent(out) :: fJac(ldfJac, n)  an output n by n
	 * array. the upper triangle of fJac //  contains an upper triangular matrix r
	 * such that // ``` //  t t t //  p *(jac *jac)*p = r *r, // ``` // 
	 * where p is a permutation matrix and jac is the final //  calculated
	 * jacobian. column j of p is column ipvt[j] //  (see below) of the identity
	 * matrix. the lower triangular //  part of fJac contains information
	 * generated during //  the computation of r. // real(wp), intent(inout) ::
	 * wa(Lwa)  a work array of length lwa. //
	 */
	public static void lmstr1(SystemOfEquations fcn, int m, int n, double[] x, double[] fVec, double[][] fJac,
			int ldfJac, double tol, int info, int[] ipvt, double[] wa, int lwa) {

		int maxFev, mode, nFev = 0, nJev = 0, nPrint;
		double fTol, gTol, xTol;
		double[] wa1 = new double[n];// real(wp), intent(inout) :: wa1(n)  work array
		double[] wa2 = new double[n];// real(wp), intent(inout) :: wa2(n)  work array
		double[] wa3 = new double[n];// real(wp), intent(inout) :: wa3(n)  work array
		double[] wa4 = new double[n];// real(wp), intent(inout) :: wa4(n)  work array
		double[] wa5 = new double[n];// real(wp), intent(inout) :: wa4(n)  work array

		double factor = 1.0e2;

		info = 0;

		// check the input parameters for errors.

		if (n > 0 & m >= n & ldfJac >= n & tol >= 0 & lwa >= (5 * n + m)) {
			// call lmstr.
			maxFev = 100 * (n + 1);
			fTol = tol;
			xTol = tol;
			gTol = 0;
			mode = 1;
			nPrint = 0;
			lmstr(fcn, m, n, x, fVec, fJac, ldfJac, fTol, xTol, gTol, maxFev, wa, mode, factor, nPrint, info, nFev,
					nJev, ipvt, wa1, wa2, wa3, wa4, wa5);
			if (info == 8) {
				info = 4;
			}
		}
	}

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

		int i, j, jm1, k, l, minmn, np1; // temporary variables
		double sum, temp; // temporary variables
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
		if (m >= np1) { // this portion has not been tested to ensure it works : USE CAUTION:
			for (j = np1 - 1; j < m; j++) {
				for (i = 0; i < m; i++) {
					q[i][j] = 0;
				}
				q[j][j] = 1;
			}
		}
//	    accumulate q from its factored form.
		for (l = 0; l < minmn; l++) {
			k = minmn - (l + 1);
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
					temp = sum / wa[k];
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
	 *               column ipvt[j] of the identity matrix. If pivot is false, ipvt
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
	 * @param rdiag  double[n] -> An array of length n which contains the diagonal
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
	public static void qrfac(int m, int n, double[][] a, int lda, boolean pivot, int[] ipvt, int lipvt, double[] rdiag,
			double[] acnorm, double[] wa) {

		int i, j, jp1, k, kmax, minmn, count; // Temporary variables
		double ajnorm, sum, temp; // Temporary variables

		double[] subArray = new double[n]; // Used for sending sub arrays to enorm for evaluation
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
				subArray[count] = a[z][j]; // As the algorithm moves diagonally down the matrix,
				count++; // a smaller and smaller subarray is sent to enorm
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
						for (i = j; i < m; i++) {
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

	/*
	 * ****************************************************************************
	 * ************* >  given an m by n matrix a, an n by n diagonal matrix d, 
	 * and an m-vector b, the problem is to determine an x which  solves the system
	 * ```  a*x = b , d*x = 0 , ```  in the least squares sense.   this
	 * subroutine completes the solution of the problem  if it is provided with the
	 * necessary information from the  qr factorization, with column pivoting, of
	 * a. that is, if  a*p = q*r, where p is a permutation matrix, q has orthogonal
	 *  columns, and r is an upper triangular matrix with diagonal  elements of
	 * nonincreasing magnitude, then qrsolv expects  the full upper triangle of r,
	 * the permutation matrix p,  and the first n components of (q transpose)*b.
	 * the system  a*x = b, d*x = 0, is then equivalent to ```  t t  r*z = q *b
	 * , p *d*p*z = 0 , ```  where x = p*z. if this system does not have full
	 * rank,  then a least squares solution is obtained. on output qrsolv  also
	 * provides an upper triangular matrix s such that ```  t t t  p *(a *a +
	 * d*d)*p = s *s . ```  s is computed within qrsolv and may be of separate
	 * interest.
	 * 
	 * subroutine qrsolv(n, r, Ldr, ipvt, diag, Qtb, x, Sdiag, wa) implicit none
	 * 
	 * integer, intent(in) :: n  a positive integer input variable set to the
	 * order of r. integer, intent(in) :: Ldr  a positive integer input variable
	 * not less than n  which specifies the leading dimension of the array r.
	 * integer, intent(in) :: ipvt(n)  an integer input array of length n which
	 * defines the  permutation matrix p such that a*p = q*r. column j of p  is
	 * column ipvt[j] of the identity matrix. real(wp), intent(inout) :: r(Ldr, n)
	 *  an n by n array. on input the full upper triangle  must contain the full
	 * upper triangle of the matrix r.  on output the full upper triangle is
	 * unaltered, and the  strict lower triangle contains the strict upper
	 * triangle  (transposed) of the upper triangular matrix s. real(wp),
	 * intent(in) :: diag(n)  an input array of length n which must contain the 
	 * diagonal elements of the matrix d. real(wp), intent(in) :: Qtb(n)  an input
	 * array of length n which must contain the first  n elements of the vector (q
	 * transpose)*b. real(wp), intent(out) :: x(n)  an output array of length n
	 * which contains the least  squares solution of the system a*x = b, d*x = 0.
	 * real(wp), intent(out) :: Sdiag(n)  an output array of length n which
	 * contains the  diagonal elements of the upper triangular matrix s. real(wp),
	 * intent(inout) :: wa(n)  a work array of length n.
	 */
	public static void qrsolv(int n, double[][] r, int ldr, int[] ipvt, double[] diag, double[] qtb, double[] x,
			double[] sdiag, double[] wa) {

		// !!!!!!!UNTESTED!!!!!!!!!!!!

		int i, j, jp1, k, kp1, l, nSing;
		double cos, cotan, qtbpj, sin, sum, tan, temp;

		double p5 = 5.0e-1;
		double p25 = 2.5e-1;

		// Copy r and (q transpose)*b to preserve input and initialize s.
		// in particular, save the diagonal elements of r in x.

		for (j = 0; j < n; j++) {
			for (i = j; i < n; i++) {
				r[i][j] = r[j][i];
			}
			x[j] = r[j][j];
			wa[j] = qtb[j];
		}

		//  eliminate the diagonal matrix d using a givens rotation.

		for (j = 0; j < n; j++) {

			// prepare the row of d to be eliminated, locating the
			// diagonal element using p from the qr factorization.

			l = ipvt[j];
			if (diag[l] != 0) {
				for (k = j; k < n; k++) {
					sdiag[k] = 0;
				}
				sdiag[j] = diag[l];

				// the transformations to eliminate the row of d
				// modify only a single element of (q transpose)*b
				// beyond the first n, which is initially zero.

				qtbpj = 0;
				for (k = j; k < n; k++) {

					// determine a givens rotation which eliminates the
					// appropriate element in the current row of d.

					if (sdiag[k] != 0) {
						if (Math.abs(r[k][k]) >= Math.abs(sdiag[k])) {
							tan = sdiag[k] / r[k][k];
							cos = p5 / Math.sqrt(p25 + p25 * tan * tan);
							sin = cos * tan;
						} else {
							cotan = r[k][k] / sdiag[k];
							sin = p5 / Math.sqrt(p25 + p25 * cotan * cotan);
							cos = sin * cotan;
						}

						// compute the modified diagonal element of r and
						// the modified element of ((q transpose)*b,0).

						r[k][k] = cos * r[k][k] + sin * sdiag[k];
						temp = cos * wa[k] + sin * qtbpj;
						qtbpj = -sin * wa[k] + cos * qtbpj;
						wa[k] = temp;

						// accumulate the tranformation in the row of s.

						kp1 = k + 1;
						if (n >= kp1) {
							for (i = kp1; i < n; i++) {
								temp = cos * r[i][k] + sin * sdiag[i];
								sdiag[i] = -sin * r[i][k] + cos * sdiag[i];
								r[i][k] = temp;
							}
						}
					}
				}
			}

			// store the diagonal element of s and restore
			// the corresponding diagonal element of r.

			sdiag[j] = r[j][j];
			r[j][j] = x[j];
		}

		// solve the triangular system for z. if the system is
		// singular, then obtain a least squares solution.

		nSing = n;
		for (j = 0; j < n; j++) {
			if (sdiag[j] == 0 & nSing == n) {
				nSing = j - 1;
			}
			if (nSing < n) {
				wa[j] = 0;
			}
		}
		if (nSing >= 1) { // might need to be changed o >= 0
			for (k = 0; k < nSing; k++) {
				j = nSing - k + 1;
				sum = 0;
				jp1 = j + 1;
				if (nSing >= jp1) {
					for (i = jp1; i < nSing; i++) {
						sum = sum + r[i][j] * wa[i];
					}
				}
				wa[j] = (wa[j] - sum) / sdiag[j];
			}
		}

		// permute the components of z back to components of x.

		for (j = 0; j < n; j++) {
			l = ipvt[j];
			x[l] = wa[j];
		}
	}

	/**
	 * 
	 * Given an m by n matrix a, this subroutine computes a*q whereq is the product
	 * of 2*(n - 1) transformations { gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1) } and
	 * gv[i], gw[i] are givens rotations in the (i,n) plane which eliminate elements
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
	 * @param v   double[n] -> An array of length n. v[i] must contain the
	 *            information necessary to recover the givens rotation gv[i]
	 *            described above.
	 *            <p>
	 * @param w   double[n] -> An array of length n. w[i] must contain the
	 *            information necessary to recover the givens rotation gw[i]
	 *            described above.
	 *            <p>
	 * @return A double[][] array
	 */
	public static double[][] r1mpyq(int m, int n, double[][] a, int lda, double[] v, double[] w) {

		int i, j, nmj, nm1; // Temporary variables
		double cos = 0, sin = 0, temp; // Temporary variables

//	    Apply the first set of givens rotations to a.
		nm1 = n;
		if (nm1 >= 0) {
			for (nmj = 0; nmj < nm1 - 1; nmj++) {
				j = n - (nmj + 2);
				if (Math.abs(v[j]) > 1) {
					cos = 1 / v[j];
					sin = Math.sqrt(1 - (cos * cos));
				} else {
					sin = v[j];
					cos = Math.sqrt(1 - (sin * sin));
				}
				for (i = 0; i < m; i++) {
					temp = cos * a[i][j] - sin * a[i][n - 1];
					a[i][n - 1] = sin * a[i][j] + cos * a[i][n - 1];
					a[i][j] = temp;
				}
			}
//			Apply the second set of givens rotations to a.
			for (j = 0; j < nm1 - 1; j++) {
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
					temp = cos * a[i][j] + sin * a[i][n - 1];
					a[i][n - 1] = -sin * a[i][j] + cos * a[i][n - 1];
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
	 * gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1) } where gv[i], gw[i] are givens rotations
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
	 *             vector v. On output v[i] contains the information necessary to
	 *             recover the givens rotation gv[i] described above.
	 *             <p>
	 * @param w    double[m] -> An array of length m. w[i] contains information
	 *             necessary to recover the givens rotation gw[i] described above.
	 *             <p>
	 * @return Nothing is returned, rather the objects themselves are modified.
	 */
	public static void r1updt(int m, int n, double[] s, int Ls, double[] u, double[] v, double[] w, boolean sing) {

		// Being left as is with the array starting at 1 notation fixes.

		int i, j, jj, l, nmj, nm1; // Temporary variables
		double cos, cotan, sin, tan, tau, temp; // Temporary variables

		double p5 = 0.5;
		double p25 = 0.25;
		double giant = 1.7976931348623157 * Math.pow(10, 308); // Largest magnitude possible

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

	/*
	 * ****************************************************************************
	 * ************* >  given an n by n upper triangular matrix r, this subroutine
	 *  computes the qr decomposition of the matrix formed when a row  is added to
	 * r. if the row is specified by the vector w, then  rwupdt determines an
	 * orthogonal matrix q such that when the  n+1 by n matrix composed of r
	 * augmented by w is premultiplied  by (q transpose), the resulting matrix is
	 * upper trapezoidal.  the matrix (q transpose) is the product of n
	 * transformations ```  g(n)*g(n-1)* ... *g(1) ```  where g[i] is a givens
	 * rotation in the (i,n+1) plane which  eliminates elements in the (n+1)-st
	 * plane. rwupdt also  computes the product (q transpose)*c where c is the 
	 * (n+1)-vector (b,alpha). q itself is not accumulated, rather  the information
	 * to recover the g rotations is supplied.
	 * 
	 * subroutine rwupdt(n, r, Ldr, w, b, Alpha, Cos, Sin) implicit none
	 * 
	 * integer, intent(in) :: n  a positive integer input variable set to the
	 * order of r. integer, intent(in) :: Ldr  a positive integer input variable
	 * not less than n  which specifies the leading dimension of the array r.
	 * real(wp), intent(inout) :: Alpha  a variable. on input alpha must contain
	 * the  (n+1)-st element of the vector c. on output alpha contains  the
	 * (n+1)-st element of the vector (q transpose)*c. real(wp), intent(inout) ::
	 * r(Ldr, n)  an n by n array. on input the upper triangular part of  r must
	 * contain the matrix to be updated. on output r  contains the updated
	 * triangular matrix. real(wp), intent(in) :: w(n)  an input array of length n
	 * which must contain the row  vector to be added to r. real(wp),
	 * intent(inout) :: b(n)  an array of length n. on input b must contain the 
	 * first n elements of the vector c. on output b contains  the first n
	 * elements of the vector (q transpose)*c. real(wp), intent(out) :: Cos(n)  an
	 * output array of length n which contains the  cosines of the transforming
	 * givens rotations. real(wp), intent(out) :: Sin(n)  an output array of
	 * length n which contains the  sines of the transforming givens rotations.
	 * 
	 * 
	 */
	public static void rwupdt(int n, double[][] r, int ldr, double[] w, double[] b, double alpha, double[] cos,

//!!!!!!!!!!!!!! CURRENTLY UNTESTED

			double[] sin) {

		int i, j, jm1;
		double cotan, rowj, tan, temp;

		double p5 = 5.0e-1;
		double p25 = 2.5e-1;

		for (j = 0; j < n; j++) {
		}
		rowj = w[j];
		jm1 = j - 1;

		// Apply the previous transformations to
		// r(i,j), i=1,2,...,j-1, and to w[j].

		if (jm1 >= 0) {
			for (i = 0; i < jm1; i++) {
				temp = cos[i] * r[i][j] + sin[i] * rowj;
				rowj = -sin[i] * r[i][j] + cos[i] * rowj;
				r[i][j] = temp;
			}
		}

		// Determine a givens rotation which eliminates w[j].

		cos[j] = 1;
		sin[j] = 0;
		if (rowj != 0) {
			if (Math.abs(r[j][j]) >= Math.abs(rowj)) {
				tan = rowj / r[j][j];
				cos[j] = p5 / Math.sqrt(p25 + p25 * tan * tan);
				sin[j] = cos[j] * tan;
			} else {
				cotan = r[j][j] / rowj;
				sin[j] = p5 / Math.sqrt(p25 + p25 * cotan * cotan);
				cos[j] = sin[j] * cotan;
			}

			// apply the current transformation to r(j,j), b[j], and alpha.

			r[j][j] = cos[j] * r[j][j] + sin[j] * rowj;
			temp = cos[j] * b[j] + sin[j] * alpha;
			alpha = -sin[j] * b[j] + cos[j] * alpha;
			b[j] = temp;
		}
	}
}
