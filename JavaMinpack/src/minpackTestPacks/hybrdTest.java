package minpackTestPacks;

import java.util.Arrays;

import Minpack.Minpack;
import Minpack.SystemOfEquations;
import Minpack.systemInterface;

//
public class hybrdTest {
//	!*****************************************************************************************
//	!>
//	// This program tests codes for the solution of n nonlinear
//	// equations in n variables. it consists of a driver and an
//	// interface // subroutine fcn. the driver reads in data, calls the
//	// nonlinear equation solver, and finally prints out information
//	// on the performance of the solver. this is only a sample driver,
//	// many other drivers are possible. the interface // subroutine fcn
//	// is necessary to take into account the forms of calling
//	// sequences used by the function // subroutines in the various
//	// nonlinear equation solvers.
//
//	program test_hybrd
//
//	    use Math.minpack_module, only: wp, dpmpar, enorm, hybrd1
//	    use iso_fortran_env, only: output_unit
//
//	    implicit none
//
//	    //originally from file21
	static int ncases = 22;
	static int[] nProbs = new int[] { 1, 2, 3, 4, 5, 6, 6, 7, 7, 7, 7, 7, 8, 8, 8, 9, 10, 10, 11, 12, 13, 14 };
	static int[] ns = new int[] { 2, 4, 2, 4, 3, 6, 9, 5, 6, 7, 8, 9, 10, 30, 40, 10, 1, 10, 10, 10, 10, 10 };
	static int[] ntriess = new int[] { 3, 3, 2, 3, 3, 2, 2, 3, 3, 3, 1, 1, 3, 1, 1, 3, 3, 3, 3, 3, 3, 3 };

	static int[] info_original = new int[] { 1, 1, 1, 4, 4, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			1, 4, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

	static int i, ic, info, k, lwa, n, nFev, nProb, ntries, icase;
	static int[] na = new int[55], nf = new int[55], np = new int[55], nx = new int[55];
	static double[] fnm = new double[55];
	static double factor, fnorm1, fnorm2;
	static double[] fVec, wa, x;
	static SystemOfEquations fcn = new SystemOfEquations();

//	    int nwrite = output_unit // logical output unit
	static double one = 1.0;
	static double ten = 10.0;
	static double tol = 1.490116119384766E-008; // sqrt dmpar(1)
	static double solution_reltol = 1.0e-4; // reltol for matching previously generated solutions

	public static void testHybrd() {

		System.out.println(tol);
		systemInterface.func f = (i) -> ffcn(i);

		ic = 0;
		for (icase = 0; icase < ncases + 1; icase++) {
			if(icase == 1) {
				break;
			}

			if (icase == ncases) {
//				System.out.println("(A,I3,A/)" + "1SUMMARY OF " + ic + " CALLS TO HYBRD1");
//				System.out.println("(A/)" + " nProb   N    nFev  INFO  FINAL L2 NORM");
				for (i = 0; i < ic; i++) {
//					System.out.println("(i4, i6, i7, i6, 1x, d15.7)" + np[i] + " " + na[i] + " " + nf[i] + " " + nx[i]
//							+ " " + fnm[i]);
				}
				break;
			} else {
				nProb = nProbs[icase];
				n = ns[icase];
				fcn.useFunctionMethod(f,n);
				lwa = (n * (3 * n + 13)) / 2;
				ntries = ntriess[icase];

				fVec = new double[n];
				wa = new double[lwa];
				x = new double[n];

				factor = one;
				for (k = 0; k < ntries; k++) {
					ic = ic + 1;
					initpt(n, x, nProb, factor);
					vecfcn(n, x, fVec, nProb);
					fnorm1 = Minpack.enorm(n, fVec);
	//				System.out.println(
	//						"(////5x,A,I5,5X,A,I5,5X//)" + " PROBLEM" + " " + nProb + " " + " DIMENSION" + " " + n);
					nFev = 0;
					System.out.println("-------");
					System.out.println(Arrays.toString(x));
					x = Minpack.hybrd1(fcn, x, tol); // hybrd1(fcn, n, x, fVec, tol, info, wa, lwa);
					fnorm2 = Minpack.enorm(n, fVec);
					np[ic-1] = nProb;
					na[ic-1] = n;
					nf[ic-1] = nFev;
					nx[ic-1] = info;
					fnm[ic-1] = fnorm2;
	//				System.out.println("(5X,A,D15.7//5X,A,D15.7//5X,A,I10//5X,A,18X,I10//5X,A//,*(5X,5D15.7/))"
	//						+ " INITIAL L2 NORM OF THE RESIDUALS" + " " + fnorm1 + " " + " FINAL L2 NORM OF THE RESIDUALS  "
	//						+ " " + fnorm2 + " " + " NUMBER OF FUNCTION EVALUATIONS  " + " " + nFev + " "
	//						+ " EXIT PARAMETER" + " " + info + " " + " FINAL APPROXIMATE SOLUTION" + " " + x[0]);
					factor = ten * factor;
					compareSolutions(ic, x, solution_reltol, tol);
			}
		}
	}
	}
//
//	    contains
//	!*****************************************************************************************
//
//	!*****************************************************************************************
//	!>
//	// Compare with previously generated solutions.
//
//	    // subroutine compare_solutions(ic, x, reltol, abstol)
//
	public static void compareSolutions(int ic, double[] x, double relTol, double absTol) {

		double[] diff = new double[x.length], absDiff = new double[x.length], relDiff = new double[x.length], icF;
		boolean failed = false;
		if (info_original[ic-1] < 5) {// Ignore any where the original Math.minpack failed
			icF = solution(ic);
			for (int i = 0; i < x.length; i++) {
				diff[i] = icF[i] - x[i];
				absDiff[i] = Math.abs(icF[i] - x[i]);
				if (absDiff[i] > absTol) {
					relDiff[i] = absDiff[i];
					if (icF[i] != 0) {
						relDiff[i] = absDiff[i] / Math.abs(icF[i]);
					}
					if (relDiff[i] > relTol) {
						failed = true;
					}
				}
			}
			if(failed) {
				System.out.println("(A)" + "Failed case : " + nProb + "/" + ic);
				System.out.println("Math.expected x: " + Arrays.toString(icF));
				System.out.println("Computed x: " + Arrays.toString(x));
				System.out.println("absdiff: " + Arrays.toString(absDiff));
				System.out.println("reldiff: " + Arrays.toString(relDiff));
				System.out.println("nfev : " + nFev);
			}
		}
	}

//	    implicit none
//
//	    integer,intent(in) :: ic // problem number (index is `solution` vector)
//	    real(wp),dimension(:),intent(in) :: x // computed `x` vector from the method
//	    real(wp),intent(in) :: reltol // relative tolerance for `x` to pass
//	    real(wp),intent(in) :: abstol // absolute tolerance for `x` to pass
//
//	    real(wp),dimension(size(x)) :: diff, absdiff, reldiff
//
//	    if (info_original(ic)<5) {    //ignore any where the original Math.minpack failed
//	        diff = solution(ic) - x
//	        absdiff = abs(diff)
//	        if (any(absdiff>abstol)) { //first do an absolute diff
//	            //also do a rel diff if the abs diff fails (also protect for divide by zero)
//	            reldiff = absdiff
//	            where (solution(ic) != 0.0) reldiff = absdiff / abs(solution(ic))
//	            if (any(reldiff > reltol)) {
//	                System.out.println("(A)") "Failed case"
//	                System.out.println( "(//5x + " " + a//(5x + " " + 5d15.7))") "Math.expected x: " + " " + solution(ic)
//	                System.out.println( "(/5x + " " + a//(5x + " " + 5d15.7))")  "Computed x: " + " " + x
//	                System.out.println( "(/5x + " " + a//(5x + " " + 5d15.7))")  "absdiff: " + " " + absdiff
//	                System.out.println( "(/5x + " " + a//(5x + " " + 5d15.7))")  "reldiff: " + " " + reldiff
//	                error stop //test failed
//	            }
//	        }
//	    }
//
//	    } compare_solutions
//	!*****************************************************************************************
//
//	!*****************************************************************************************
//	!>
//	// The calling sequence of fcn should be identical to the
//	// calling sequence of the function // subroutine in the nonlinear
//	// equation solver. fcn should only the testing function
//	// // subroutine vecfcn with the appropriate value of problem
//	// number (nProb).
//
//	    // subroutine fcn(n, x, fVec, Iflag)
//	        implicit none
//
//	        integer, intent(in) :: n // the number of variables.
//	        real(wp), intent(in) :: x(n) // independent variable vector
//	        real(wp), intent(out) :: fVec[n] // value of function at `x`
//	        integer, intent(inout) :: iflag // set to <0 to terMath.minate execution
	public static double[] ffcn(double[] x) {
		vecfcn(n, x, fVec, nProb);
		nFev = nFev + 1;
		return fVec;
	}

//	    } fcn
//	!*****************************************************************************************
//
//	!*****************************************************************************************
//	!>
//	// Replaced statement function in original code.
//
	public static double dfloat(int i) {
		double f = i;
		return f;
	}

//	    pure elemental function dfloat(i) result(f)
//	        implicit none
//	        integer, intent(in) :: i
//	        double f
//	        f = real(i, wp)
//	    end function dfloat
//	!*****************************************************************************************
//
//	!*****************************************************************************************
//	!>
//	// Get Math.expected `x` vectors for each case.
//
//	    pure function solution(nProb) result(x)
	public static double[] solution(int nProb) {
//	        implicit none
//
//	        integer,intent(in) :: nProb
		double[] x;
//
		switch (nProb) {
		case (1):
			x = new double[] { 0.1000000000000000E+01, 0.1000000000000000E+01 };
			return x;
		case (2):
			x = new double[] { 0.1000000000000000E+01, 0.1000000000000000E+01 };
			return x;
		case (3):
			x = new double[] { 0.1000000000000000E+01, 0.1000000000000054E+01 };
			return x;
		case (4):
			x = new double[] { -0.1717232736993598E-17, 0.1717232736993598E-18, 0.4885791716489701E-17,
					0.4885791716489701E-17 };
			return x;
		case (5):
			x = new double[] { 0.1137338840805565E-17, -0.1137338840805565E-18, 0.1509231876962185E-17,
					0.1509231876962185E-17 };
			return x;
		case (6):
			x = new double[] { 0.2071578632741476E-17, -0.2071578632741476E-18, 0.3365608460018520E-17,
					0.3365608460018520E-17 };
			return x;
		case (7):
			x = new double[] { 0.1098159327798559E-04, 0.9106146740037904E+01 };
			return x;
		case (8):
			x = new double[] { 0.1098159288127784E-04, 0.9106146743611500E+01 };
			return x;
		case (9):
			x = new double[] { -0.9679740249513952E+00, 0.9471391408446626E+00, -0.9695163103174519E+00,
					0.9512476657649955E+00 };
			return x;
		case (10):
			x = new double[] { -0.9679740249362032E+00, 0.9471391408151544E+00, -0.9695163103329795E+00,
					0.9512476657950136E+00 };
			return x;
		case (11):
			x = new double[] { -0.9679740247487649E+00, 0.9471391404515958E+00, -0.9695163105216826E+00,
					0.9512476661606586E+00 };
			return x;
		case (12):
			x = new double[] { 0.1000000000000010E+01, -0.1612103012704913E-13, -0.8125674485239315E-34 };
			return x;
		case (13):
			x = new double[] { 0.1000000000004274E+01, 0.1388633807362216E-10, 0.0000000000000000E+00 };
			return x;
		case (14):
			x = new double[] { 0.9999999999999840E+00, 0.1238719384388127E-14, 0.0000000000000000E+00 };
			return x;
		case (15):
			x = new double[] { -0.1572508640134011E-01, 0.1012434869369118E+01, -0.2329916259567960E+00,
					0.1260430087800365E+01, -0.1513728922723441E+01, 0.9929964324318560E+00 };
			return x;
		case (16):
			x = new double[] { -0.1572508639053690E-01, 0.1012434869369907E+01, -0.2329916259623205E+00,
					0.1260430087870332E+01, -0.1513728922830754E+01, 0.9929964324984680E+00 };
			return x;
		case (17):
			x = new double[] { -0.1530703652147214E-04, 0.9997897039319488E+00, 0.1476396369354227E-01,
					0.1463423282995085E+00, 0.1000821103004075E+01, -0.2617731140517609E+01, 0.4104403164477184E+01,
					-0.3143612278555425E+01, 0.1052626408009917E+01 };
			return x;
		case (18):
			x = new double[] { -0.1530703713913476E-04, 0.9997897039319459E+00, 0.1476396369267279E-01,
					0.1463423283016577E+00, 0.1000821102994390E+01, -0.2617731140495362E+01, 0.4104403164446781E+01,
					-0.3143612278533631E+01, 0.1052626408003180E+01 };
			return x;
		case (19):
			x = new double[] { 0.8375125649983552E-01, 0.3127292952224503E+00, 0.5000000000008663E+00,
					0.6872707047760241E+00, 0.9162487435008237E+00 };
			return x;
		case (20):
			x = new double[] { 0.6872707047770207E+00, 0.9162487435004409E+00, 0.4999999999996554E+00,
					0.8375125649942240E-01, 0.3127292952234606E+00 };
			return x;
		case (21):
			x = new double[] { 0.5000000001403082E+00, 0.6872707047172968E+00, 0.8375125651152422E-01,
					0.3127292951388482E+00, 0.9162487434920225E+00 };
			return x;
		case (22):
			x = new double[] { 0.6687659094768218E-01, 0.3666822990106460E+00, 0.2887406733471217E+00,
					0.7112593271119373E+00, 0.6333177005294922E+00, 0.9331234090531204E+00 };
			return x;
		case (23):
			x = new double[] { 0.9331234090539151E+00, 0.2887406731191550E+00, 0.6687659094608109E-01,
					0.7112593268807791E+00, 0.3666822992420550E+00, 0.6333177007580146E+00 };
			return x;
		case (24):
			x = new double[] { 0.3666822992460977E+00, 0.6333177007631426E+00, 0.7112593268766961E+00,
					0.6687659094559376E-01, 0.9331234090527033E+00, 0.2887406731157666E+00 };
			return x;
		case (25):
			x = new double[] { 0.5806914977912385E-01, 0.2351716115667845E+00, 0.3380440957994889E+00,
					0.4999999991484015E+00, 0.6619559063505691E+00, 0.7648283868041619E+00, 0.9419308505514702E+00 };
			return x;
		case (26):
			x = new double[] { 0.3380440947502161E+00, 0.4999999997939603E+00, 0.2351716123793344E+00,
					0.7648283876237569E+00, 0.9419308503662596E+00, 0.5806914962097134E-01, 0.6619559054655015E+00 };
			return x;
		case (27):
			x = new double[] { -0.3267366079625573E+02, -0.2996209843926172E+02, -0.8587775264169514E+02,
					0.2222113097994968E+02, 0.5957249137089175E+02, -0.1038025653217158E+01, 0.8600842862942351E+02 };
			return x;
		case (28):
			x = new double[] { 0.4985640222318974E-01, 0.1986351285003365E+00, 0.2698288337443381E+00,
					0.4992723176748156E+00, 0.5007277255753518E+00, 0.7301712224978171E+00, 0.8013649159179719E+00,
					0.9501436601762751E+00 };
			return x;
		case (29):
			x = new double[] { 0.4420534615691015E-01, 0.1994906721904692E+00, 0.2356191086780681E+00,
					0.4160469078466623E+00, 0.4999999996232831E+00, 0.5839530926716184E+00, 0.7643808911925417E+00,
					0.8005093278089829E+00, 0.9557946538314642E+00 };
			return x;
		case (30):
			x = new double[] { 0.1000000000000008E+01, 0.1000000000000008E+01, 0.1000000000000008E+01,
					0.1000000000000008E+01, 0.1000000000000008E+01, 0.1000000000000008E+01, 0.1000000000000008E+01,
					0.1000000000000008E+01, 0.1000000000000008E+01, 0.9999999999999193E+00 };
			return x;
		case (31):
			x = new double[] { 0.1000000000000002E+01, 0.1000000000000002E+01, 0.1000000000000002E+01,
					0.1000000000000002E+01, 0.1000000000000002E+01, 0.1000000000000002E+01, 0.1000000000000002E+01,
					0.1000000000000002E+01, 0.1000000000000002E+01, 0.9999999999999805E+00 };
			return x;
		case (32):
			x = new double[] { 0.1000000000000126E+01, 0.1000000000000126E+01, 0.1000000000000126E+01,
					0.1000000000000126E+01, 0.1000000000000126E+01, 0.1000000000000126E+01, 0.1000000000000126E+01,
					0.1000000000000126E+01, 0.1000000000000126E+01, 0.9999999999987337E+00 };
			return x;
		case (33):
			x = new double[] { 0.1000000000000116E+01, 0.1000000000000116E+01, 0.1000000000000116E+01,
					0.1000000000000116E+01, 0.1000000000000116E+01, 0.1000000000000116E+01, 0.1000000000000116E+01,
					0.1000000000000116E+01, 0.1000000000000116E+01, 0.1000000000000116E+01, 0.1000000000000116E+01,
					0.1000000000000116E+01, 0.1000000000000116E+01, 0.1000000000000116E+01, 0.1000000000000116E+01,
					0.1000000000000116E+01, 0.1000000000000116E+01, 0.1000000000000116E+01, 0.1000000000000116E+01,
					0.1000000000000116E+01, 0.1000000000000116E+01, 0.1000000000000116E+01, 0.1000000000000116E+01,
					0.1000000000000116E+01, 0.1000000000000116E+01, 0.1000000000000116E+01, 0.1000000000000116E+01,
					0.1000000000000116E+01, 0.1000000000000116E+01, 0.9999999999965015E+00 };
			return x;
		case (34):
			x = new double[] { 0.9999999999999820E+00, 0.9999999999999820E+00, 0.9999999999999820E+00,
					0.9999999999999820E+00, 0.9999999999999820E+00, 0.9999999999999820E+00, 0.9999999999999820E+00,
					0.9999999999999820E+00, 0.9999999999999820E+00, 0.9999999999999820E+00, 0.9999999999999820E+00,
					0.9999999999999820E+00, 0.9999999999999820E+00, 0.9999999999999820E+00, 0.9999999999999820E+00,
					0.9999999999999820E+00, 0.9999999999999820E+00, 0.9999999999999820E+00, 0.9999999999999820E+00,
					0.9999999999999820E+00, 0.9999999999999820E+00, 0.9999999999999820E+00, 0.9999999999999820E+00,
					0.9999999999999820E+00, 0.9999999999999820E+00, 0.9999999999999820E+00, 0.9999999999999820E+00,
					0.9999999999999820E+00, 0.9999999999999820E+00, 0.9999999999999820E+00, 0.9999999999999820E+00,
					0.9999999999999820E+00, 0.9999999999999820E+00, 0.9999999999999820E+00, 0.9999999999999820E+00,
					0.9999999999999820E+00, 0.9999999999999820E+00, 0.9999999999999820E+00, 0.9999999999999820E+00,
					0.1000000000000726E+01 };
			return x;
		case (35):
			x = new double[] { -0.4316498251876486E-01, -0.8157715653538729E-01, -0.1144857143805310E+00,
					-0.1409735768625996E+00, -0.1599086961819857E+00, -0.1698772023127759E+00, -0.1690899837812081E+00,
					-0.1552495352218312E+00, -0.1253558916789345E+00, -0.7541653368589182E-01 };
			return x;
		case (36):
			x = new double[] { -0.4316498251881876E-01, -0.8157715653546944E-01, -0.1144857143805966E+00,
					-0.1409735768626190E+00, -0.1599086961819499E+00, -0.1698772023126901E+00, -0.1690899837811062E+00,
					-0.1552495352217907E+00, -0.1253558916789970E+00, -0.7541653368596339E-01 };
			return x;
		case (37):
			x = new double[] { -0.4316498254524553E-01, -0.8157715658128553E-01, -0.1144857144024714E+00,
					-0.1409735768723229E+00, -0.1599086963003020E+00, -0.1698772022538506E+00, -0.1690899837944776E+00,
					-0.1552495352060321E+00, -0.1253558916432041E+00, -0.7541653366609002E-01 };
			return x;
		case (38):
			x = new double[] { -0.1528138835625800E+00 };
			return x;
		case (39):
			x = new double[] { -0.1528138835625801E+00 };
			return x;
		case (40):
			x = new double[] { -0.1528138835625801E+00 };
			return x;
		case (41):
			x = new double[] { -0.4316498251876487E-01, -0.8157715653538730E-01, -0.1144857143805310E+00,
					-0.1409735768625996E+00, -0.1599086961819857E+00, -0.1698772023127759E+00, -0.1690899837812081E+00,
					-0.1552495352218312E+00, -0.1253558916789344E+00, -0.7541653368589175E-01 };
			return x;
		case (42):
			x = new double[] { -0.4316498251881876E-01, -0.8157715653546946E-01, -0.1144857143805966E+00,
					-0.1409735768626190E+00, -0.1599086961819498E+00, -0.1698772023126901E+00, -0.1690899837811062E+00,
					-0.1552495352217907E+00, -0.1253558916789970E+00, -0.7541653368596334E-01 };
			return x;
		case (43):
			x = new double[] { -0.4316498251876519E-01, -0.8157715653538752E-01, -0.1144857143805303E+00,
					-0.1409735768625980E+00, -0.1599086961819844E+00, -0.1698772023127748E+00, -0.1690899837812073E+00,
					-0.1552495352218307E+00, -0.1253558916789341E+00, -0.7541653368589159E-01 };
			return x;
		case (44):
			x = new double[] { 0.5526115715943968E-01, 0.5695713693095779E-01, 0.5889020336090119E-01,
					0.6113556183519356E-01, 0.6377799292665667E-01, 0.6700414432471043E-01, 0.2079417258113421E+00,
					0.1642681193175131E+00, 0.8643704095817571E-01, 0.9133212907808361E-01 };
			return x;
		case (45):
			x = new double[] { 0.3439628896235289E-01, 0.3503231575416022E-01, 0.3571919583574593E-01,
					0.3646522422001942E-01, 0.3728091174083566E-01, 0.3817986258974846E-01, 0.3918014109819012E-01,
					0.4030650261419996E-01, 0.1797201916815169E+00, 0.1562408814749922E+00 };
			return x;
		case (46):
			x = new double[] { 0.1888395221036672E+02, 0.2516777354434988E+02, 0.1888527511739392E+02,
					0.1888602114554983E+02, 0.1888683683452867E+02, 0.1888773578345333E+02, 0.1888873606083097E+02,
					0.1888986242379029E+02, 0.1902927611240455E+02, 0.1900579680335179E+02 };
			return x;
		case (47):
			x = new double[] { 0.9999999999999992E+00, 0.9999999999999984E+00, 0.9999999999999977E+00,
					0.9999999999999969E+00, 0.9999999999999961E+00, 0.9999999999999953E+00, 0.9999999999999946E+00,
					0.9999999999999938E+00, 0.9999999999999930E+00, 0.9999999999999922E+00 };
			return x;
		case (48):
			x = new double[] { 0.1000000000000000E+01, 0.1000000000000000E+01, 0.1000000000000000E+01,
					0.1000000000000000E+01, 0.1000000000000001E+01, 0.1000000000000001E+01, 0.1000000000000001E+01,
					0.1000000000000001E+01, 0.1000000000000001E+01, 0.1000000000000001E+01 };
			return x;
		case (49):
			x = new double[] { 0.1000000000000000E+01, 0.1000000000000000E+01, 0.1000000000000000E+01,
					0.1000000000000000E+01, 0.1000000000000000E+01, 0.1000000000000000E+01, 0.1000000000000000E+01,
					0.1000000000000000E+01, 0.1000000000000000E+01, 0.1000000000000000E+01 };
			return x;
		case (50):
			x = new double[] { -0.5707221307212121E+00, -0.6818069509055232E+00, -0.7022100775689857E+00,
					-0.7055106309936168E+00, -0.7049061557572888E+00, -0.7014966060124587E+00, -0.6918893211477919E+00,
					-0.6657965141985400E+00, -0.5960351099566767E+00, -0.4164122574358191E+00 };
			return x;
		case (51):
			x = new double[] { -0.5707221320932939E+00, -0.6818069495100820E+00, -0.7022100764111258E+00,
					-0.7055106298493696E+00, -0.7049061556844529E+00, -0.7014966070294095E+00, -0.6918893223739674E+00,
					-0.6657965143753547E+00, -0.5960351092038981E+00, -0.4164122574142932E+00 };
			return x;
		case (52):
			x = new double[] { -0.5707221320171143E+00, -0.6818069499829604E+00, -0.7022100760171542E+00,
					-0.7055106298955310E+00, -0.7049061557301967E+00, -0.7014966070327222E+00, -0.6918893223590803E+00,
					-0.6657965144072677E+00, -0.5960351090088830E+00, -0.4164122575177334E+00 };
			return x;
		case (53):
			x = new double[] { -0.4283028636053099E+00, -0.4765964242962535E+00, -0.5196524638125549E+00,
					-0.5580993246169652E+00, -0.5925061569509362E+00, -0.6245036821428087E+00, -0.6232394714478015E+00,
					-0.6213938418388717E+00, -0.6204535966122983E+00, -0.5864692707477792E+00 };
			return x;
		case (54):
			x = new double[] { -0.4283028634881692E+00, -0.4765964236396713E+00, -0.5196524642776766E+00,
					-0.5580993248351936E+00, -0.5925061568131795E+00, -0.6245036817962692E+00, -0.6232394720687789E+00,
					-0.6213938417874499E+00, -0.6204535965224117E+00, -0.5864692707287930E+00 };
			return x;
		case (55):
			x = new double[] { -0.4283028635608067E+00, -0.4765964243232715E+00, -0.5196524637037395E+00,
					-0.5580993248328234E+00, -0.5925061568292707E+00, -0.6245036822076749E+00, -0.6232394714256790E+00,
					-0.6213938418143937E+00, -0.6204535966527651E+00, -0.5864692707189498E+00 };
			return x;

		default:
			x = new double[] { 10 };
			return x;
		}
	}

//	    end function solution
//	!*****************************************************************************************
//
//	!*****************************************************************************************
//	!>
//	// This // subroutine defines fourteen test functions. the first
//	// five test functions are of dimensions 2,4,2,4,3, respectively,
//	// while the remaining test functions are of variable dimension
//	// n for any n greater than or equal to 1 (problem 6 is an
//	// exception to this, Math.since it does not allow n = 1).
	public static void vecfcn(int n, double[] x, double[] fVec, int nProb) {
//	    // subroutine vecfcn(n, x, fVec, nProb)
//	        implicit none
//
//	        integer, intent(in) :: n // a positive integer input variable.
//	        integer, intent(in) :: nProb // a positive integer input variable which defines the
//	                                     // number of the problem. nProb must not exceed 14.
//	        real(wp), intent(in) :: x(n) // an input array of length n.
//	        real(wp), intent(out) :: fVec[n] // an output array of length n which contains the nProb
//	                                         // function vector evaluated at x.
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

		int i, iev, ivar, j, k, k1, k2, kp1, ml, mu;
		double h, prod, sum, sum1, sum2, temp, temp1, temp2, ti, tj, tk, tpi;

		// fVec(1:n) = zero
		for (i = 0; i < n; i++) {
			fVec[i] = 0;
		}
		// PROBLEM SELECTOR.
		switch (nProb) {
		case (2):
			// POWELL Math.sinGULAR FUNCTION.
			fVec[0] = x[0] + ten * x[1];
			fVec[1] = Math.sqrt(five) * (x[2] - x[3]);
			fVec[2] = (x[1] - two * x[2]) * (x[1] - two * x[2]);
			fVec[3] = Math.sqrt(ten) * (x[0] - x[3]) * (x[0] - x[3]);
			break;
		case (3):
			// POWELL BADLY SCALED FUNCTION.
			fVec[0] = c1 * x[0] * x[1] - one;
			fVec[1] = Math.exp(-x[0]) + Math.exp(-x[1]) - c2;
			break;
		case (4):
			// WOOD FUNCTION.
			temp1 = x[1] - x[0] * x[0];
			temp2 = x[3] - x[2] * x[2];
			fVec[0] = -c3 * x[0] * temp1 - (one - x[0]);
			fVec[1] = c3 * temp1 + c4 * (x[1] - one) + c5 * (x[3] - one);
			fVec[2] = -c6 * x[2] * temp2 - (one - x[2]);
			fVec[3] = c6 * temp2 + c4 * (x[3] - one) + c5 * (x[1] - one);
			break;
		case (5):
			// HELICAL VALLEY FUNCTION.
			tpi = eight * Math.atan(one);
			if (x[0] > zero) {
				temp1 = Math.atan(x[1] / x[0]) / tpi;
			} else if (x[0] < zero) {
				temp1 = Math.atan(x[1] / x[0]) / tpi + c8;
			} else {
				temp1 = Math.signum(x[1]) * Math.abs(c7); // (c7, x[1]); //does this ever happen?
			}
			temp2 = Math.sqrt(x[0] * x[0] + x[1] * x[1]);
			fVec[0] = ten * (x[2] - ten * temp1);
			fVec[1] = ten * (temp2 - one);
			fVec[2] = x[2];
			break;
		case (6):
			// WATSON FUNCTION.
			for (k = 0; k < n; k++) {
				fVec[k] = zero;
			}
			for (i = 0; i < 29; i++) {
				ti = dfloat(i + 1) / c9;
				sum1 = zero;
				temp = one;
				for (j = 1; j < n; j++) {
					sum1 = sum1 + dfloat(j) * temp * x[j];
					temp = ti * temp;
				}
				sum2 = zero;
				temp = one;
				for (j = 0; j < n; j++) {
					sum2 = sum2 + temp * x[j];
					temp = ti * temp;
				}
				temp1 = sum1 - sum2 * sum2 - one;
				temp2 = two * ti * sum2;
				temp = one / ti;
				for (k = 0; k < n; k++) {
					fVec[k] = fVec[k] + temp * (dfloat(k) - temp2) * temp1;
					temp = ti * temp;
				}
			}
			temp = x[1] - x[0] * x[0] - one;
			fVec[0] = fVec[0] + x[0] * (one - two * temp);
			fVec[1] = fVec[1] + temp;
			break;
		case (7):
			// CHEBYQUAD FUNCTION.
			for (k = 0; k < n; k++) {
				fVec[k] = zero;
			}
			for (j = 0; j < n; j++) {
				temp1 = one;
				temp2 = two * x[j] - one;
				temp = two * temp2;
				for (i = 0; i < n; i++) {
					fVec[i] = fVec[i] + temp2;
					ti = temp * temp2 - temp1;
					temp1 = temp2;
					temp2 = ti;
				}
			}
			tk = one / dfloat(n);
			iev = -1;
			for (k = 0; k < n; k++) {
				fVec[k] = tk * fVec[k];
				if (iev > 0)
					fVec[k] = fVec[k] + one / (dfloat(k + 1) * dfloat(k + 1) - one);
				iev = -iev;
			}
			break;
		case (8):
			// BROWN ALMOST-LINEAR FUNCTION.
			sum = -dfloat(n + 1);
			prod = one;
			for (j = 0; j < n; j++) {
				sum = sum + x[j];
				prod = x[j] * prod;
			}
			for (k = 0; k < n; k++) {
				fVec[k] = x[k] + sum;
			}
			fVec[n-1] = prod - one;
			break;
		case (9):
			// DISCRETE BOUNDARY VALUE FUNCTION.
			h = one / dfloat(n + 1);
			for (k = 0; k < n; k++) {
				temp = (x[k] + dfloat(k + 1) * h + one) * (x[k] + dfloat(k + 1) * h + one)
						* (x[k] + dfloat(k + 1) * h + one);
				temp1 = zero;
				if (k != 0) { 
					temp1 = x[k - 1];
				}
				temp2 = zero;
				if (k != n-1) {
					temp2 = x[k + 1];
				}
				fVec[k] = two * x[k] - temp1 - temp2 + temp * h * h / two;
			}
			break;
		case (10):
			// DISCRETE INTEGRAL EQUATION FUNCTION.
			h = one / dfloat(n + 1);
			for (k = 0; k < n; k++) {
				tk = dfloat(k + 1) * h;
				sum1 = zero;
				for (j = 0; j < k; j++) {
					tj = dfloat(j + 1) * h;
					temp = (x[j] + tj + one) * (x[j] + tj + one) * (x[j] + tj + one);
					sum1 = sum1 + tj * temp;
				}
				sum2 = zero;
				kp1 = k + 1;
				if (n >= kp1) {
					for (j = kp1; j < n; j++) {
						tj = dfloat(j + 1) * h;
						temp = (x[j] + tj + one) * (x[j] + tj + one) * (x[j] + tj + one);
						sum2 = sum2 + (one - tj) * temp;
					}
				}
				fVec[k] = x[k] + h * ((one - tk) * sum1 + tk * sum2) / two;
			}
			break;
		case (11):
			// TRIGONOMETRIC FUNCTION.
			sum = zero;
			for (j = 0; j < n; j++) {
				fVec[j] = Math.cos(x[j]);
				sum = sum + fVec[j];
			}
			for (k = 0; k < n; k++) {
				fVec[k] = dfloat(n + k + 1) - Math.sin(x[k]) - sum - dfloat(k + 1) * fVec[k];
			}
			break;
		case (12):
			// VARIABLY DIMENSIONED FUNCTION.
			sum = zero;
			for (j = 0; j < n; j++) {
				sum = sum + dfloat(j + 1) * (x[j] - one);
			}
			temp = sum * (one + two * sum * sum);
			for (k = 0; k < n; k++) {
				fVec[k] = x[k] - one + dfloat(k + 1) * temp;
			}
			break;
		case (13):
			// BROYDEN TRIDIAGONAL FUNCTION.
			for (k = 0; k < n; k++) {
				temp = (three - two * x[k]) * x[k];
				temp1 = zero;
				if (k != 0) {
					temp1 = x[k - 1];
				}
				temp2 = zero;
				if (k != n-1) {
					temp2 = x[k + 1];
				}
				fVec[k] = temp - temp1 - two * temp2 + one;
			}
			break;
		case (14):
			// BROYDEN BANDED FUNCTION.
			ml = 5;
			mu = 2;
			for (k = 0; k < n; k++) {
				k1 = Math.max(0, k - ml);
				k2 = Math.min(k + mu, n);
				temp = zero;
				for (j = k1; j < k2; j++) {
					if (j != k) {
						temp = temp + x[j] * (one + x[j]);
					}
				}
				fVec[k] = x[k] * (two + five * x[k] * x[k]) + one - temp;
			}
			break;
		default:
			// ROSENBROCK FUNCTION.
			fVec[0] = one - x[0];
			fVec[1] = ten * (x[1] - x[0] * x[0]);
			break;
		}
	} // vecfcn
//	!*****************************************************************************************
//
//	!*****************************************************************************************
//	!>
//	// This // subroutine specifies the standard starting points for
//	// the functions defined by // subroutine vecfcn. the // subroutine
//	// returns in x a multiple (factor) of the standard starting
//	// point. for the sixth function the standard starting point is
//	// zero, so in this case, if factor is not unity, { the
//	// // subroutine returns the vector  x[j] = factor, j=1,...,n.

	public static void initpt(int n, double[] x, int nProb, double factor) {
		//// subroutine initpt(n, x, nProb, Factor)
		// implicit none

		// integer, intent(in) :: n // a positive integer input variable.
		// real(wp), intent(out) :: x(n) // an output array of length n which contains
		// the standard
		// starting point for problem nProb multiplied by factor.
		// integer, intent(in) :: nProb // a positive integer input variable which
		// defines the
		// number of the problem. nProb must not exceed 14.
		// real(wp), intent(in) :: Factor // an input variable which specifies the
		// multiple of
		// the standard starting point. if factor is unity, no
		// multiplication is performed.

		int ivar, j;
		double h, tj;

		double zero = 0.0;
		double half = 0.5;
		double one = 1.0;
		double three = 3.0;
		double c1 = 1.2;

		// x(1:n) = zero
		for (i = 0; i < n; i++) {
			x[i] = 0;
		}

		// selection of initial point.

		switch (nProb) {
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
			for (j = 0; j < n; j++) {
				x[j] = zero;
			}
		break;
		case (7):
			// chebyquad function.
			h = one / dfloat(n + 1);
			for (j = 0; j < n; j++) {
				x[j] = dfloat(j + 1) * h;
			}
			break;
		case (8):
			// brown almost-linear function.
			for (j = 0; j < n; j++) {
				x[j] = half;
			}
		break;
		case (9):
		case (10):
			// discrete boundary value and integral equation functions.
			h = one / dfloat(n + 1);
			for (j = 0; j < n; j++) {
				tj = dfloat(j + 1) * h;
				x[j] = tj * (tj - one);
			}
			break;
		case (11):
			// trigonometric function.
			h = one / dfloat(n);
			for (j = 0; j < n; j++) {
				x[j] = h;
			}
			break;
		case (12):
			// variably dimensioned function.
			h = one / dfloat(n);
			for (j = 0; j < n; j++) {
				x[j] = one - dfloat(j + 1) * h;
			}
			break;
		case (13):
		case (14):
			// broyden tridiagonal and banded functions.
			for (j = 0; j < n; j++) {
				x[j] = -one;
			}
		break;
		default:
			// rosenbrock function.
			x[0] = -c1;
			x[1] = one;
			break;
		}

		// compute multiple of initial point.

		if (factor != one) {
			if (nProb == 6) {
				for (j = 0; j < n; j++) {
					x[j] = factor;
				}
			} else {
				for (j = 0; j < n; j++) {
					x[j] = factor * x[j];
				}
			}
		}

	}// initpt
//	!*****************************************************************************************
//
//	!*****************************************************************************************
//	end program test_hybrd
//	!*****************************************************************************************
}
