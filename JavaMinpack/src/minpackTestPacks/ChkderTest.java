package minpackTestPacks;

import java.util.Arrays;

import Minpack.Minpack;

public class ChkderTest {
//	!*****************************************************************************************
//	!>
//	!  This program tests the ability of [[chkder]] to detect
//	!  inconsistencies between functions and their first derivatives.
//	!  fourteen test function vectors and jacobians are used. eleven of
//	!  the tests are false(f), i.e. there are inconsistencies between
//	!  the function vectors and the corresponding jacobians. three of
//	!  the tests are true(t), i.e. there are no inconsistencies. the
//	!  driver reads in data, calls chkder and prints out information
//	!  required by and received from chkder.
//
//	program test_chkder
		static int ncases = 14;
		static int[] nprobs  = {1,2,3,4,5,6,7,8,9,10,11,12,13,14};
		static int[] ns      = {2,4,2,4,3,9,7,10,10,10,10,10,10,10};
	//
		static int i, ldfJac, lnp = 0, mode, n, nprob, icase;
		static double cp;
		static int[] na = new int[ncases], np = new int[ncases];
		static double[] errmax = new double[ncases], errmin = new double[ncases];
		static double[] diff, err, fVec1, fVec2, x1, x2;
		static double[][] fJac;
	//
		static boolean[] a = {false, false, false, true, false, false, false,true, false, false, false, false, true, false};
		static double one = 1.0;
		static double tol = 1.490116119384766E-008; //double tol = Math.sqrt(dpmpar(1)); //!! abstol for matching previously generated solutions
		static double solution_reltol = 1.0e-4; //!! reltol for matching previously generated solutions
	//
		static int[] info_original = {1}; //! not used here
	
		public static void chkderTest() {
//
//	    use minpack_module, only: wp, dpmpar, chkder
//	    use iso_fortran_env, only: nwrite => output_unit
//
//	    implicit none
//
//	    ! originally from file23
	    
//
	    cp = 1.23e-1;
//
	    for( icase = 0;icase< ncases+1;icase++) {
	    	
	        if (icase == ncases) {
	            System.out.println ("SUMMARY OF " + lnp + " TESTS OF CHKDER");
	            System.out.println (" NPROB   N    STATUS     ERRMIN         ERRMAX");
	            for(i = 0;i< lnp;i++) {
	                System.out.println ("" + np[i] + "         " + na[i] + "      " + a[i] + "      " + errmin[i] + "         " + errmax[i]);
	            }
	            break; //stop
	        }else {
	        	
	            nprob = nprobs[icase];
	            n = ns[icase];
	            ldfJac = n;

	            System.out.println(">>>   " + nprobs[icase]);
	            diff = new double[n];
            	err= new double[n];
            	fJac= new double[n][n];
            	fVec1= new double[n];
            	fVec2= new double[n];
            	x1= new double[n];
            	x2= new double[n];
	            
	           

	            initpt(n, x1, nprob, one);
	            for(i = 0;i< n;i++) {
	                x1[i] = x1[i] + cp;
	                cp = -cp;
	            }
	            System.out.println (" PROBLEM : " + nprob + " WITH DIMENSION"+ n + " IS  " + a[nprob-1]);
	            mode = 1;
	            Minpack.chkder(n, n, x1, fVec1, fJac, ldfJac, x2, fVec2, mode, err);
	            mode = 2;
	            vecfcn(n, x1, fVec1, nprob);
	            errjac(n, x1, fJac, ldfJac, nprob);
	            vecfcn(n, x2, fVec2, nprob);
	            Minpack.chkder(n, n, x1, fVec1, fJac, ldfJac, x2, fVec2, mode, err);
//	            
	            errmin[nprob-1] = err[0];
	            errmax[nprob-1] = err[0];
	            for( i = 0;i< n;i++) {
	                diff[i] = fVec2[i] - fVec1[i];
	                if (errmin[nprob-1] > err[i]) {
	                	errmin[nprob-1] = err[i];
	                }
	                if (errmax[nprob-1] < err[i]) {
	                	errmax[nprob-1] = err[i];
	                }
	            }
	            np[nprob-1] = nprob;
	            lnp = nprob;
	            na[nprob-1] = n;
	            System.out.println (" FIRST FUNCTION VECTOR " + fVec1[i-1]  + n);
	            System.out.println (" FUNCTION DIFFERENCE VECTOR " + diff[i-1] + n);
	            System.out.println (" ERROR VECTOR " + err[i-1] +  n);
	            compareSolutions(nprob, diff, solution_reltol, tol);

	        }
	    	}
//
		}
//	    contains
//	!*****************************************************************************************
//
//	!*****************************************************************************************
//	!>
//	!  Compare with previously generated solutions.
	public static void compareSolutions(int ic, double[] x, double relTol, double absTol){
//	    subroutine compare_solutions(ic, x, reltol, abstol)
//
//	    implicit none
//
//	    integer,intent(in) :: ic !! problem number (index is `solution` vector)
		
//	    real(wp),dimension(:),intent(in) :: x !! computed `x` vector from the method
	
//	    real(wp),intent(in) :: reltol !! relative tolerance for `x` to pass
		
//	    real(wp),intent(in) :: abstol !! absolute tolerance for `x` to pass
		
//
//	    real(wp),dimension(size(x)) :: diff, absdiff, reldiff
		double[] diff = new double[x.length], absDiff = new double[x.length], relDiff = new double[x.length], icF;
//
//	    if (info_original(ic)<5) {//    Ignore any where the original minpack failed
			icF = solution(ic);
			for(int i =0;i<x.length;i++) {
				diff[i] = icF[i] - x[i];
				absDiff[i] = Math.abs(icF[i] - x[i]);
				if(absDiff[i]>absTol) {
					relDiff[i] = absDiff[i];
					if(icF[i] != 0) {
						relDiff[i] = absDiff[i] / Math.abs(icF[i]);
					}
					if(relDiff[i] > relTol) {
						System.out.println("(A)"+ "Failed case");
		                System.out.println("Expected x: " + Arrays.toString(icF));
		                System.out.println("Computed x: " + Arrays.toString(x));
		                System.out.println("absdiff: " + absDiff[i]);
		                System.out.println( "reldiff: " + relDiff[i]);
					}
				}
			}
//	        diff = solution(ic) - x;
//	        absDiff = Math.abs(diff);
//	        if ((absDiff>absTol)) { // if (any(absdiff>abstol)) { // ! first do an absolute diff
//	            ! also do a rel diff if the abs diff fails (also protect for divide by zero)
//	            relDiff = absDiff;
//	            where (solution(ic) /= 0.0;) reldiff = absdiff / abs(solution(ic))
//	            if ((relDiff > relTol)) { //if any reldiff > reltol then print
//	                System.out.println(nwrite,"(A)") "Failed case"
//	                System.out.println(nwrite, "(//5x, a//(5x, 5d15.7))") "Expected x: ", solution(ic)
//	                System.out.println(nwrite, "(/5x, a//(5x, 5d15.7))")  "Computed x: ", x
//	                System.out.println(nwrite, "(/5x, a//(5x, 5d15.7))")  "absdiff: ", absdiff
//	                System.out.println(nwrite, "(/5x, a//(5x, 5d15.7))")  "reldiff: ", reldiff
	                //error stop ! test failed
//	}
//	}
//	}
	}
//	    end subroutine compare_solutions
//	!*****************************************************************************************
//
//	!*****************************************************************************************
//	!>
//	!  Replaced statement function in original code.
		public static double dfloat(int i) {
//	    pure elemental function dfloat[i] result(f)
//	        implicit none
//	        integer, intent(in) :: i
//	        double f
//	        f = real(i, wp)
			double f = i;
			return f;
		}
//	    end function dfloat
//	!*****************************************************************************************
//
//	!*****************************************************************************************
//	!>
//	!  Get expected `diff` vectors for each case.
		public static double[] solution(int nprob) {
//	    pure function solution(nprob) result(x)
//
//	        implicit none
//
//	        integer,intent(in) :: nprob
//	        real(wp),dimension(:),allocatable :: x
			double[] x;
//
			switch(nprob) {
//	        select case (nprob)
		        case(1): x = new double[]  {-0.1604855093262358E-07,0.4763689633868751E-06};
		        return x;
		        case(2): x = new double[]  {0.2138763655068487E-06,-0.2512328678427878E-07,-0.3578105500778861E-07,
		                      0.4754114257821129E-06};
		        return x;
		        case(3): x = new double[]  {0.3214806338291964E-04,-0.7057517459330143E-08};
		        return x;
		        case(4): x = new double[]  {0.2322078535144101E-03,0.5335169362297165E-04,0.2089913541567512E-03,
		                      0.4808346034224087E-04};
		        return x;
		        case(5): x = new double[]  {0.1832842144722235E-07,-0.1319622122686326E-06,0.1832842821958280E-08};
		        return x;
		        case(6): x = new double[]  {0.4515008482641747E-06,0.8252608125758343E-06,0.1047075926408070E-05,
		                      0.1220878203866960E-05,0.1363612746274612E-05,0.1485299534920159E-05,
		                      0.1592014982065848E-05,0.1687697022134671E-05,0.1775024010441939E-05};
		        return x;
		        case(7): x = new double[]  {0.1542483058641908E-07,0.2060287401794980E-07,0.2376945576476608E-07,
		                      0.5349154558187408E-07,0.9704076181504817E-07,0.1576633008593120E-06,
		                      0.2112184365188341E-06};
		        return x;
		        case(8): x = new double[]  {0.8012354335562577E-07,0.8378922800034161E-07,0.8012354335562577E-07,
		                      0.8378922800034161E-07,0.8012354335562577E-07,0.8378922800034161E-07,
		                      0.8012354335562577E-07,0.8378922800034161E-07,0.8012354335562577E-07,
		                      0.1065043608861060E-09};
		        return x;
		        case(9): x = new double[]  {0.5774599409757997E-08,-0.7078711450336783E-08,0.7631400400498478E-08,
		                      -0.7053519379685014E-08,0.7658129463905539E-08,-0.7038501670386665E-08,
		                      0.7685261871337445E-08,-0.7047089523037897E-08,0.6495039228671118E-08,
		                      -0.2818530409065545E-08};
		        return x;
		        case(10):x = new double[]  {0.3294705064327275E-08,0.8148107188965525E-09,0.5413627796047038E-08,
		                      0.2381044444943470E-08,0.6401980501280491E-08,0.2764786941056308E-08,
		                      0.6166095051218790E-08,0.1882141179021524E-08,0.4645276802106579E-08,
		                      0.9133731548871538E-09};
		        return x;
		        case(11):x = new double[]  {0.3284536753689338E-08,0.1864164600462459E-08,0.3268772807984988E-08,
		                      0.3333951337225471E-08,0.3253008529213730E-08,0.4803738740122299E-08,
		                      0.3237243362264053E-08,0.6273523034394657E-08,0.3221478195314376E-08,
		                      0.7743309993202274E-08};
		        return x;
		        case(12):x = new double[]  {0.2249654149636626E-02,0.4499298869632185E-02,0.6748936313670129E-02,
		                      0.8998581033665687E-02,0.1124821836128831E-01,0.1349786319769919E-01,
		                      0.1574750046711415E-01,0.1799714541994035E-01,0.2024678350426257E-01,
		                      0.2249642740935087E-01};
		        return x;
		        case(13):x = new double[]  {0.9923452193305593E-07,0.3484660204833290E-07,0.8616620350565540E-07,
		                      0.3484660204833290E-07,0.8616620350565540E-07,0.3484660204833290E-07,
		                      0.8616620350565540E-07,0.3484660204833290E-07,0.8616620350565540E-07,
		                      0.6831460996892247E-07};
		        return x;
		        case(14):x = new double[]  {0.3598775801805232E-06,0.2186061109910042E-06,0.3905816612359558E-06,
		                      0.2493101911582585E-06,0.4212857422913885E-06,0.2800142722136911E-06,
		                      0.4311392522993174E-06,0.2800142722136911E-06,0.4311392522993174E-06,
		                      0.2591637029425442E-06};
		        return x;
		        default:
		        	x = new double[] {10};
		            return x; //error stop "invalid case"
			}
}
//	    end function solution
//	!*****************************************************************************************
//
//	!*****************************************************************************************
//	!>
//	!  This subroutine is derived from [[vecjac]] which defines the
//	!  jacobian matrices of fourteen test functions. the problem
//	!  dimensions are as described in the prologue comments of vecfcn.
//	!  various errors are deliberately introduced to provide a test
//	!  for chkder.
		public static void errjac(int n, double[] x, double[][] fJac, int ldfJac, int nProb) {
//	    subroutine errjac(n, x, fJac, LdfJac, Nprob)
//	        implicit none
//
//	        integer, intent(in) :: n !! a positive integer variable.
//	        integer, intent(in) :: LdfJac !! a positive integer variable not less than n
//	                                      !! which specifies the leading dimension of the array fJac.
//	        integer, intent(in) :: Nprob !! a positive integer variable which defines the
//	                                     !! number of the problem. nprob must not exceed 14.
//	        real(wp), intent(in) :: x(n) !! an array of length n.
//	        real(wp), intent(out) :: fJac(LdfJac, n) !! an n by n array. on output fJac contains the
//	                                                 !! jacobian matrix, with various errors deliberately
//	                                                 !! introduced, of the nprob function evaluated at x.
//
	        double zero = 0.0e0;
	        double one = 1.0e0;
	        double two = 2.0e0;
	        double three = 3.0e0;
	        double four = 4.0e0;
	        double five = 5.0e0;
	        double six = 6.0e0;
	        double eight = 8.0e0;
	        double ten = 1.0e1;
	        double fiftn = 1.5e1;
	        double twenty = 2.0e1;
	        double hundrd = 1.0e2;
	        double c1 = 1.0e4;
	        double c3 = 2.0e2;
	        double c4 = 2.02e1;
	        double c5 = 1.98e1;
	        double c6 = 1.8e2;
	        double c9 = 2.9e1;
	        double tpi = eight*Math.atan(one);
//
	        int i, j, k, k1, k2, ml, mu;
	        double h, prod, sum, sum1, sum2, temp, temp1, temp2,
	                    temp3, temp4, ti, tj, tk;
//
//	        fJac = new double[n][n]; // fJac(1:n,1:n) = zero
	        for(i = 0;i<n;i++) {
	        	for(j = 0;i<n;i++) {
	    	    	fJac[i][j]=0;
	    	    }
		    }
//
//	        Jacobian routine selector.
//
	        switch(nProb) {
//	        select case (nprob)
	        case (2):
//	            Powell singular function with sign reversal affecting element
//	            ! (3,3).
	            for (k = 0;k< 4;k++){
	                for( j = 0;j< 4;j++){
	                    fJac[k][j] = 0;
	                }
	            }
	            fJac[0][0] = one;
	            fJac[0][1] = ten;
	            fJac[1][2] = Math.sqrt(five);
	            fJac[1][3] = -fJac[1][2];
	            fJac[2][1] = two*(x[1] - two*x[2]);
	            fJac[2][2] = two*fJac[2][1];
	            fJac[3][0] = two*Math.sqrt(ten)*(x[0] - x[3]);
	            fJac[3][3] = -fJac[3][0];
	            break;
	        
	        case (3):
//	            ! powell badly scaled function with the sign of the jacobian
//	            ! reversed.
	            fJac[0][0] = -c1*x[1];
	            fJac[0][1] = -c1*x[0];
	            fJac[1][0] = Math.exp(-x[0]);
	            fJac[1][1] = Math.exp(-x[1]);
	            break;
	        case (4):
	            // wood function without error.
	            for (k = 0;k< 4;k++){
	                for( j = 0;j< 4;j++){
	                    fJac[k][j] = zero;
	                }
	            }
	            temp1 = x[1] - three*x[0]*x[0];
	            temp2 = x[3] - three*x[2]*x[2];
	            fJac[0][0] = -c3*temp1 + one;
	            fJac[0][1] = -c3*x[0];
	            fJac[1][0] = -two*c3*x[0];
	            fJac[1][1] = c3 + c4;
	            fJac[1][3] = c5;
	            fJac[2][2] = -c6*temp2 + one;
	            fJac[2][3] = -c6*x[2];
	            fJac[3][1] = c5;
	            fJac[3][2] = -two*c6*x[2];
	            fJac[3][3] = c6 + c4;
	            		break;
	        case (5):
	            // helical valley function with multiplicative error affecting
	            // elements (2,1) and (2,2).
	            temp = x[0]*x[0] + x[1]*x[1];
	            temp1 = tpi*temp;
	            temp2 = Math.sqrt(temp);
	            fJac[0][0] = hundrd*x[1]/temp1;
	            fJac[0][1] = -hundrd*x[0]/temp1;
	            fJac[0][2] = ten;
	            fJac[1][0] = five*x[0]/temp2;
	            fJac[1][1] = five*x[1]/temp2;
	            fJac[1][2] = zero;
	            fJac[2][0] = zero;
	            fJac[2][1] = zero;
	            fJac[2][2] = one;
	            		break;
	        case (6):
	            // watson function with sign reversals affecting the computation of
	            // temp1.
	            for(k=0;k<n;k++){
	                for(j=k;j<n;j++){
	                    fJac[k][j] = zero;
	                }
	            }
	            for(i = 0;i< 29;i++) {
	                ti = dfloat(i+1)/c9;
	                sum1 = zero;
	                temp = one;
	                for(j = 1;j< n;j++) {
	                    sum1 = sum1 + dfloat(j)*temp*x[j];
	                    temp = ti*temp;
	                }
	                sum2 = zero;
	                temp = one;
	                for(j=0;j<n;j++){
	                    sum2 = sum2 + temp*x[j];
	                    temp = ti*temp;
	                }
	                temp1 = two*(sum1 + sum2*sum2 + one);
	                temp2 = two*sum2;
	                temp = ti*ti;
	                tk = one;
	                for(k=0;k<n;k++){
	                    tj = tk;
	                    for(j=k;j<n;j++){
	                        fJac[k][j] = fJac[k][j] + tj*((dfloat(k)/ti - temp2)*(dfloat(j)/ti - temp2) - temp1);
	                        tj = ti*tj;
	                    }
	                    tk = temp*tk;
	                }
	            }
	            fJac[0][0] = fJac[0][0] + six*x[0]*x[0] - two*x[1] + three;
	            fJac[0][1] = fJac[0][1] - two*x[0];
	            fJac[1][1] = fJac[1][1] + one;
	            for(k=0;k<n;k++){
	                for(j=k;j<n;j++){
	                    fJac[j][k] = fJac[k][j];
	                }
	            }
	            break;
	        case (7):
	            // chebyquad function with jacobian twice correct size.
	            tk = one/dfloat(n);
	            for(j=0;j<n;j++){
	                temp1 = one;
	                temp2 = two*x[j] - one;
	                temp = two*temp2;
	                temp3 = zero;
	                temp4 = two;
	                for(k=0;k<n;k++){
	                    fJac[k][j] = two*tk*temp4;
	                    ti = four*temp2 + temp*temp4 - temp3;
	                    temp3 = temp4;
	                    temp4 = ti;
	                    ti = temp*temp2 - temp1;
	                    temp1 = temp2;
	                    temp2 = ti;
	                }
	            }
	        break;
	        case (8):
	            // brown almost-linear function without error.
	            prod = one;
	            for(j=0;j<n;j++){
	                prod = x[j]*prod;
	                for(k=0;k<n;k++){
	                    fJac[k][j] = one;
	                }
	                fJac[j][j] = two;
	            }
	            for(j=0;j<n;j++){
	                temp = x[j];
	                if (temp == zero) {
	                    temp = one;
	                    prod = one;
	                    for(k=0;k<n;k++){
	                        if (k != j) {
	                        	prod = x[k]*prod;
	                        }
	                    }
	                }
	                fJac[n-1][j] = prod/temp;
	            }
	            break;
	        case (9):
	            // discrete boundary value function with multiplicative error
	            // affecting the jacobian diagonal.
	            h = one/dfloat(n + 1);
	            for(k=0;k<n;k++){
	                temp = three*(x[k] + dfloat(k+1)*h + one)*(x[k] + dfloat(k+1)*h + one);
	                for(j=0;j<n;j++){
	                    fJac[k][j] = zero;
	                }
	                fJac[k][k] = four + temp*h*h;
	                if (k != 0) {
	                	fJac[k][k - 1] = -one;
	                }
	                if (k != n-1) {
	                	fJac[k][k + 1] = -one;
	                }
	            }
	        break;
	        case (10):
	            // discrete integral equation function with sign error affecting
	            // the jacobian diagonal.
	            h = one/dfloat(n + 1);
	            for(k=0;k<n;k++){
	                tk = dfloat(k+1)*h;
	                for(j=0;j<n;j++){
	                    tj = dfloat(j+1)*h;
	                    temp = three*(x[j] + tj + one)*(x[j] + tj + one);
	                    fJac[k][j] = h*Math.min(tj*(one - tk), tk*(one - tj))*temp/two;
	                }
	                fJac[k][k] = fJac[k][k] - one;
	            }
	        break;
	        case (11):
	            // trigonometric function with sign errors affecting the
	            // offdiagonal elements of the jacobian.
	            for(j=0;j<n;j++){
	                temp = Math.sin(x[j]);
	                for(k=0;k<n;k++){
	                    fJac[k][j] = -temp;
	                }
	                fJac[j][j] = dfloat(j + 2)*temp - Math.cos(x[j]);
	            }
	        break;
	        case (12):
//	            // variably dimensioned function with operation error affecting
//	            // the upper triangular elements of the jacobian.
	            sum = zero;
	            for(j=0;j<n;j++){
	                sum = sum + dfloat(j+1)*(x[j] - one);
	            }
	            temp = one + six*sum*sum;
	            for(k=0;k<n;k++){
	                for(j=k;j<n;j++){
	                    fJac[k][j] = dfloat((k+1)*(j+1))/temp;
	                    fJac[j][k] = fJac[k][j];
	                }
	                fJac[k][k] = fJac[k][k] + one;
	            }
	                break;
	        case (13):
//	            // broyden tridiagonal function without error.
	            for(k=0;k<n;k++){
	                for(j=0;j<n;j++){
	                    fJac[k][j] = zero;
	                }
	                fJac[k][k] = three - four*x[k];
	                if (k != 0) {
	                	fJac[k][k - 1] = -one;
	                }
	                if (k != n-1) {
	                	fJac[k][k + 1] = -two;
	                }
	            }
	        break;
	        case (14):
//	            // broyden banded function with sign error affecting the jacobian
//	            // diagonal.
	            ml = 5;
	            mu = 2;
	            for(k=0;k<n;k++){
	                for(j=0;j<n;j++){
	                    fJac[k][j] = zero;
	                }
	                k1 = Math.max(0, k - ml);
	                k2 = Math.min(k + mu, n);
	                for(j = k1;j< k2;j++) {
	                    if (j != k) {
	                    	fJac[k][j] = -(one + two*x[j]);
	                    }
	                }
	                fJac[k][k] = two - fiftn*x[k]*x[k];
	             }
	             break;
	        default:
//	            // rosenbrock function with sign reversal affecting element (1,1).
	            fJac[0][0] = one;
	            fJac[0][1] = zero;
	            fJac[1][0] = -twenty*x[0];
	            fJac[1][1] = ten;
	            break;
	        
	        }
	}
//	    end subroutine errjac
//	!*****************************************************************************************
//
//	!*****************************************************************************************
//	!>
//	!  This subroutine specifies the standard starting points for
//	!  the functions defined by subroutine vecfcn. the subroutine
//	!  returns in x a multiple (factor) of the standard starting
//	!  point. for the sixth function the standard starting point is
//	!  zero, so in this case, if factor is not unity, then the
//	!  subroutine returns the vector  x(j) = factor, j=1,...,n.
		public static double[] initpt(int n, double[] x, int nProb, double factor) {
//	    subroutine initpt(n, x, Nprob, Factor)
//	        implicit none
//
//	        integer,intent(in) :: n !! a positive integer input variable.
//	        real(wp),intent(out) :: x(n) !! an output array of length n which contains the standard
//	                                     !! starting point for problem nprob multiplied by factor.
//	        integer,intent(in) :: Nprob !! a positive integer input variable which defines the
//	                                    !! number of the problem. nprob must not exceed 14.
//	        real(wp),intent(in) :: Factor !! an input variable which specifies the multiple of
//	                                      !! the standard starting point. if factor is unity, no
//	                                      !! multiplication is performed.
//
	        double zero = 0.0;
	        double half = 5.0e-1;
	        double one = 1.0;
	        double three = 3.0;
	        double c1 = 1.2;

	        int j;
	        double h, tj;
//
//	        x(1:n) = zero
//	        x = new double[n];
	        for(int i = 0;i<n;i++) {
		    	x[i] = 0;
		    }
//
//	        // selection of initial point.
//
	        switch(nProb) {
	        case (2):
	            // powell singular function.
	            x[0] = three;
	            x[1] = -one;
	            x[2] = zero;
	            x[3] = one;
	        break;
	        case (3):
	            // powell badly scaled function.
	            x[0] = zero;
	            x[1] = one;
	        break;
	        case (4):
	            // wood function.
	            x[0] = -three;
	            x[1] = -one;
	            x[2] = -three;
	            x[3] = -one;
	        break;
	        case (5):
	            // helical valley function.
	            x[0] = -one;
	            x[1] = zero;
	            x[2] = zero;
	        break;
	        case (6):
	            // watson function.
	            for(j=0;j<n;j++){
	                x[j] = zero;
	            }
	        break;
	        case (7):
	            // chebyquad function.
	            h = one/dfloat(n + 1);
	            for(j=0;j<n;j++){
	                x[j] = dfloat(j+1)*h;
	            }
	        break;
	        case (8):
	            // brown almost-linear function.
	            for(j=0;j<n;j++){
	                x[j] = half;
	            }
	        break;
	        case (9):
	        case (10):
	            // discrete boundary value and integral equation functions.
	            h = one/dfloat(n + 1);
	            for(j=0;j<n;j++){
	                tj = dfloat(j+1)*h;
	                x[j] = tj*(tj - one);
	            }
	        break;
	        case (11):
	            // trigonometric function.
	            h = one/dfloat(n);
	            for(j=0;j<n;j++){
	                x[j] = h;
	            }
	        break;
	        case (12):
	            // variably dimensioned function.
	            h = one/dfloat(n);
	            for(j=0;j<n;j++){
	                x[j] = one - dfloat(j+1)*h;
	            }
	        break;
	        case (13):
	        case (14):
	            // broyden tridiagonal and banded functions.
	            for(j=0;j<n;j++){
	                x[j] = -one;
	            }
	        break;
	        default:
	            // rosenbrock function.
	            x[0] = -c1;
	            x[1] = one;
	        break;
//	        end select
	        }
//
//	        // Compute multiple of initial point.
//
	        if (factor != one) {
	            if (nProb == 6) {
	                for(j=0;j<n;j++){
	                    x[j] = factor;
	                }
	            }else {
	                for(j=0;j<n;j++){
	                    x[j] = factor*x[j];
	                }
	            }
	        }
	    return x;
}
//	    end subroutine initpt
//	!*****************************************************************************************
//
//	!*****************************************************************************************
//	!>
//	!  This subroutine defines fourteen test functions. the first
//	!  five test functions are of dimensions 2,4,2,4,3, respectively,
//	!  while the remaining test functions are of variable dimension
//	!  n for any n greater than or equal to 1 (problem 6 is an
//	!  exception to this, since it does not allow n = 1).
public static void vecfcn(int n, double[] x, double[] fVec, int nProb) {
//	    subroutine vecfcn(n, x, fVec, Nprob)
//	        implicit none
//
//	        integer,intent(in) :: n !! a positive integer input variable.
//	        integer,intent(in) :: Nprob !! a positive integer input variable which defines the
//	                                    !! number of the problem. nprob must not exceed 14.
//	        real(wp),intent(in) :: x(n) !!  an input array of length n.
//	        real(wp),intent(out) :: fVec(n) !! an output array of length n which contains the nprob
//	                                        !! function vector evaluated at x.
//
	        double zero = 0.0;
	        double one = 1.0;
	        double two = 2.0;
	        double three = 3.0;
	        double five = 5.0;
	        double eight = 8.0;
	        double ten = 10.0;
	        double c1 = 1.0e4;
	        double c2 = 1.0001;
	        double c3 = 2.0e2;
	        double c4 = 2.02e1;
	        double c5 = 1.98e1;
	        double c6 = 1.8e2;
	        double c7 = 2.5e-1;
	        double c8 = 5.0e-1;
	        double c9 = 2.9e1;
	        double tpi = eight*Math.atan(one);

	        int i, iev, j, k, k1, k2, kp1, ml, mu;
	        double h, prod, sum, sum1, sum2, temp, temp1, temp2, ti, tj = 0, tk;
//
//	        fVec(1:n) = zero
//	        fVec = new double[n];
	        for(i = 0;i<n;i++) {
		    	fVec[i]=0;
		    }
//
//	        // problem selector.
//
//	        select case (nprob)
	        switch(nProb) {
	        case (2):
	            // powell singular function.
	            fVec[0] = x[0] + ten*x[1];
	            fVec[1] = Math.sqrt(five)*(x[2] - x[3]);
	            fVec[2] = (x[1] - two*x[2])*(x[1] - two*x[2]);
	            fVec[3] = Math.sqrt(ten)*(x[0] - x[3])*(x[0] - x[3]);
	            break;
	        case (3):
	            // powell badly scaled function.
	            fVec[0] = c1*x[0]*x[1] - one;
	            fVec[1] = Math.exp(-x[0]) + Math.exp(-x[1]) - c2;
	            break;
	        case (4):
	            // wood function.
	            temp1 = x[1] - x[0]*x[0];
	            temp2 = x[3] - x[2]*x[2];
	            fVec[0] = -c3*x[0]*temp1 - (one - x[0]);
	            fVec[1] = c3*temp1 + c4*(x[1] - one) + c5*(x[3] - one);
	            fVec[2] = -c6*x[2]*temp2 - (one - x[2]);
	            fVec[3] = c6*temp2 + c4*(x[3] - one) + c5*(x[1] - one);
	            break;
	        case (5):
	            // helical valley function.
	            temp1 = Math.signum(x[1]) * Math.abs(c7);
	            if (x[0] > zero) temp1 = Math.atan(x[1]/x[0])/tpi;
	            if (x[0] < zero) temp1 = Math.atan(x[1]/x[0])/tpi + c8;
	            temp2 = Math.sqrt(x[0]*x[0] + x[1]*x[1]);
	            fVec[0] = ten*(x[2] - ten*temp1);
	            fVec[1] = ten*(temp2 - one);
	            fVec[2] = x[2];
	            break;
	        case (6):
	            // watson function.
	            for(k=0;k<n;k++){
	                fVec[k] = zero;
	            }
	            for( i = 0;i< 29;i++) {
	                ti = dfloat(i+1)/c9;
	                sum1 = zero;
	                temp = one;
	                for( j = 1;j< n;j++) {
	                    sum1 = sum1 + dfloat(j)*temp*x[j];
	                    temp = ti*temp;
	                }
	                sum2 = zero;
	                temp = one;
	                for(j=0;j<n;j++){
	                    sum2 = sum2 + temp*x[j];
	                    temp = ti*temp;
	                }
	                temp1 = sum1 - sum2*sum2 - one;
	                temp2 = two*ti*sum2;
	                temp = one/ti;
	                for(k=0;k<n;k++){
	                    fVec[k] = fVec[k] + temp*(dfloat(k) - temp2)*temp1;
	                    temp = ti*temp;
	                }
	            }
	            temp = x[1] - x[0]*x[0] - one;
	            fVec[0] = fVec[0] + x[0]*(one - two*temp);
	            fVec[1] = fVec[1] + temp;
	            break;
	        case (7):
	            // chebyquad function.
	            for(k=0;k<n;k++){
	                fVec[k] = zero;
	            }
	            for(j=0;j<n;j++){
	                temp1 = one;
	                temp2 = two*x[j] - one;
	                temp = two*temp2;
	                for( i = 0;i< n;i++) {
	                    fVec[i] = fVec[i] + temp2;
	                    ti = temp*temp2 - temp1;
	                    temp1 = temp2;
	                    temp2 = ti;
	                }
	            }
	            tk = one/dfloat(n);
	            iev = -1;
	            for(k=0;k<n;k++){
	                fVec[k] = tk*fVec[k];
	                if (iev > 0) {
	                	fVec[k] = fVec[k] + one/(dfloat(k+1)*dfloat(k+1) - one);
	                }
	                iev = -iev;
	            }
	            break;
	        case (8):
	            // brown almost-linear function.
	            sum = -dfloat(n + 1);
	            prod = one;
	            for(j=0;j<n;j++){
	                sum = sum + x[j];
	                prod = x[j]*prod;
	            }
	            for(k=0;k<n;k++){
	                fVec[k] = x[k] + sum;
	            }
	            fVec[n-1] = prod - one;
	            break;
	        case (9):
	            // discrete boundary value function.
	            h = one/dfloat(n + 1);
	            for(k=0;k<n;k++){
	                temp = (x[k] + dfloat(k+1)*h + one)*(x[k] + dfloat(k+1)*h + one)*(x[k] + dfloat(k+1)*h + one);
	                temp1 = zero;
	                if (k != 0) {
	                	temp1 = x[k - 1];
	                }
	                temp2 = zero;
	                if (k != n-1) {
	                	temp2 = x[k + 1];
	                }
	                fVec[k] = two*x[k] - temp1 - temp2 + temp*h*h/two;
	            }
	            break;
	        case (10):
	            // discrete integral equation function.
	            h = one/dfloat(n + 1);
	            for(k=0;k<n;k++){
	                tk = dfloat(k+1)*h;
	                sum1 = zero;
	                for( j = 0;j< k+1;j++){
	                    tj = dfloat(j+1)*h;
	                    temp = (x[j] + tj + one)*(x[j] + tj + one)*(x[j] + tj + one);
	                    sum1 = sum1 + tj*temp;
	                }
	                sum2 = zero;
	                kp1 = k + 1;
	                if (n > kp1) {
	                    for(j = kp1;j< n;j++) {
	                        tj = dfloat(j+1)*h;
	                        temp = (x[j] + tj + one)*(x[j] + tj + one)*(x[j] + tj + one);
	                        sum2 = sum2 + (one - tj)*temp;
	                    }
	                }
	                fVec[k] = x[k] + h*((one - tk)*sum1 + tk*sum2)/two;
	            }
	            break;
	        case (11):
	            // trigonometric function.
	            sum = zero;
	            for(j=0;j<n;j++){
	                fVec[j] = Math.cos(x[j]);
	                sum = sum + fVec[j];
	            }
	            for(k=0;k<n;k++){
	                fVec[k] = dfloat(n + k+1) - Math.sin(x[k]) - sum - dfloat(k+1)*fVec[k];
	            }
	            break;
	        case (12):
	            // variably dimensioned function.
	            sum = zero;
	            for(j=0;j<n;j++){
	                sum = sum + dfloat(j+1)*(x[j] - one);
	            }
	            temp = sum*(one + two*sum*sum);
	            for(k=0;k<n;k++){
	                fVec[k] = x[k] - one + dfloat(k+1)*temp;
	            }
	            break;
	        case (13):
	            // broyden tridiagonal function.
	            for(k=0;k<n;k++){
	                temp = (three - two*x[k])*x[k];
	                temp1 = zero;
	                if (k != 0){
	        			temp1 = x[k - 1];
	        }
	                temp2 = zero;
	                if (k != n-1) temp2 = x[k + 1];
	                fVec[k] = temp - temp1 - two*temp2 + one;
	            }
	        break;
	        case (14):
	            // broyden banded function.
	            ml = 5;
	            mu = 2;
	            for(k=0;k<n;k++){
	                k1 = Math.max(0, k - ml);
	                k2 = Math.min(k + mu, n);
	                temp = zero;
	                for(j = k1;j<k2;j++){
	                    if (j != k) {
	        				temp = temp + x[j]*(one + x[j]);
	        }
	                }
	                fVec[k] = x[k]*(two + five*x[k]*x[k]) + one - temp;
	            }
	            break;
	        default:
	            // rosenbrock function.
	            fVec[0] = one - x[0];
	            fVec[1] = ten*(x[1] - x[0]*x[0]);
	            break;
//	        end select
	        }
//
//	    end subroutine vecfcn
//
//	end program test_chkder
}
}
