package minpackTestPacks;

import java.util.Arrays;

import Minpack.Minpack;
import Minpack.SystemOfEquations;
import Minpack.systemInterface;

//
public class hybrjTest {
//	!*****************************************************************************************
//	!>
//	// This program tests codes for the solution of n nonlinear
//	// equations in n variables. it consists of a driver and an
//	// interface subroutine fcn. the driver reads in data, calls the
//	// nonlinear equation solver, and finally prints out information
//	// on the performance of the solver. this is only a sample driver,
//	// many other drivers are possible. the interface subroutine fcn
//	// is necessary to take into account the forms of calling
//	// sequences used by the function and jacobian subroutines in
//	// the various nonlinear equation solvers.
//
//	program test_hybrj
//
//	    use Math.minpack_module, only: wp, dpmpar, enorm, hybrj1
//	    use iso_fortran_env, only: nwrite => output_unit
//
//	    implicit none
//
//	    //originally from file21
	static int ncases = 22;
	static int[] nProbs = new int[] { 1, 2, 3, 4, 5, 6, 6, 7, 7, 7, 7, 7, 8, 8, 8, 9, 10, 10, 11, 12, 13, 14 };
	static int[] ns = new int[] { 2, 4, 2, 4, 3, 6, 9, 5, 6, 7, 8, 9, 10, 30, 40, 10, 1, 10, 10, 10, 10, 10 };
	static int[] ntriess = new int[] { 3, 3, 2, 3, 3, 2, 2, 3, 3, 3, 1, 1, 3, 1, 1, 3, 3, 3, 3, 3, 3, 3 };

	static int[] info_original = { 1, 1, 1, 4, 4, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 1, 1, 1, 1, 1, 1, 1, 1, 4, 4,
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

	static int i, ic, info, k, n, nFev, nJev, nProb, ntries, icase, lwa, ldfJac;
	static int[] na = new int[55], nf = new int[55], nj = new int[55], np = new int[55], nx = new int[55];
	static double[] fnm = new double[55];
	static double factor, fnorm1, fnorm2;
	static double[][] fJac;
	static double[] fVec, wa, x;

	static double one = 1.0;
	static double ten = 10.0;
	static double tol = 1.490116119384766E-008; // sqrt dmpar(1)
	static double solution_reltol = 1.0e-4; // reltol for matching previously generated solutions
	static SystemOfEquations fcn = new SystemOfEquations();

//
	public static void testHybrj() {
		ic = 0;
		systemInterface.func f = (i) -> ffcn(i);
		
		systemInterface.jacob j = (i) -> jfcn(i);
		

		for (icase = 0; icase < ncases + 1; icase++) {
			if (icase == ncases + 1) {
				System.out.println("(A,I3,A/)" + "1SUMMARY OF " + ic + " CALLS TO HYBRJ1");
				System.out.println("(A/)" + " nProb   N    nFev   nJev  INFO  FINAL L2 NORM");
				for (i = 0; i < ic; i++) {
					System.out.println("(I4,I6,2I7,I6,1X,D15.7)" + np[i] + " " + na[i] + " " + nf[i] + " " + nj[i] + " "
							+ nx[i] + " " + fnm[i]);
				}
				break;
			} else {
				nProb = nProbs[icase];
				n = ns[icase];
				fcn.useFunctionMethod(f,n);
				fcn.useJacobianMethod(j,n,n);
				ldfJac = n;
				lwa = (n * (n + 13)) / 2;
				ntries = ntriess[icase];

				fJac = new double[n][n];
				fVec = new double[n];
				wa = new double[lwa];
				x = new double[n];

				factor = one;
				for (k = 0; k < ntries; k++) {
					ic = ic + 1;
					initpt(n, x, nProb, factor);
					vecfcn(n, x, fVec, nProb);
					fnorm1 = Minpack.enorm(n, fVec);
					System.out.println(
							"(////5X,A,I5,5X,A,I5,5X//)" + " PROBLEM" + " " + nProb + " " + " DIMENSION" + " " + n);
					nFev = 0;
					nJev = 0;
					Minpack.hybrj1(fcn, n, x, fVec, fJac, ldfJac, tol, info, wa, lwa);
					fnorm2 = Minpack.enorm(n, fVec);
					np[ic-1] = nProb;
					na[ic-1] = n;
					nf[ic-1] = nFev;
					nj[ic-1] = nJev;
					nx[ic-1] = info;
					fnm[ic-1] = fnorm2;
					System.out.println("(5X,A,D15.7//5X,A,D15.7//5X,A,I10//5X,A,I10//5X,A,18X,I10//5X,A//*(5X,5D15.7/))"
							+ " INITIAL L2 NORM OF THE RESIDUALS" + " " + fnorm1 + " "
							+ " FINAL L2 NORM OF THE RESIDUALS  " + " " + fnorm2 + " "
							+ " NUMBER OF FUNCTION EVALUATIONS  " + " " + nFev + " "
							+ " NUMBER OF JACOBIAN EVALUATIONS  " + " " + nJev + " " + " EXIT PARAMETER" + " " + info
							+ " " + " FINAL APPROXIMATE SOLUTION" + " " + x[0]);
					factor = ten * factor;
					compareSolutions(ic, x, solution_reltol, tol);
					;

				}
			}
		}
	}

//	    contains
//	!*****************************************************************************************
//
//	!*****************************************************************************************
//	!>
//	// Compare with previously generated solutions.
//	
	public static void compareSolutions(int ic, double[] x, double relTol, double absTol) {

		double[] diff = new double[x.length], absDiff = new double[x.length], relDiff = new double[x.length], icF;
		if (info_original[ic-1] < 5) {// Ignore any where the original minpack failed
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
						System.out.println("(A)" + "Failed case");
						System.out.println("Math.expected x: " + Arrays.toString(icF));
						System.out.println("Computed x: " + Arrays.toString(x));
						System.out.println("absdiff: " + absDiff[i]);
						System.out.println("reldiff: " + relDiff[i]);
					}
				}
			}
		}
	}

//	   
//	!*****************************************************************************************
//
//	!*****************************************************************************************
//	!>
//	// The calling sequence of fcn should be identical to the
//	// calling sequence of the function subroutine in the nonlinear
//	// equation solver. fcn should only the testing function
//	// and jacobian subroutines vecfcn and vecjac with the
//	// appropriate value of problem number (nProb).
	public static double[] ffcn(double[] x) {
		double[] fVec = new double[n];
		double[][] fJac = new double[ldfJac][n];
		vecfcn(n, x, fVec, nProb);
		nFev = nFev + 1;
		return fVec;

	}

	public static double[][] jfcn(double[] x) {
		double[] fVec = new double[n];
		double[][] fJac = new double[ldfJac][n];
		vecjac(n, x, fJac, ldfJac, nProb);
		nJev = nJev + 1;
		return fJac;

	}

//	    subroutine fcn(n, x, fVec, fJac, ldfJac, iFlag)
//	        implicit none
//
//	        int n;
//	        int ldfJac;
//	        int iFlag;
//	        double[] fVec = new double[n];
//	        double[][] fJac = new double[ldfJac][n];
//
//	        switch (iFlag)
//	        case (1)
//	            vecfcn(n, x, fVec, nProb)
//	            nFev = nFev + 1
//	        case(2)
//	            vecjac(n, x, fJac, ldfJac, nProb)
//	            nJev = nJev + 1
//	        default:
//	            error stop "invalid iFlag value"
//	        } 
//
//	    }  // fcn
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
//	        int nProb
//	        real(wp),dimension(:),allocatable :: x
		double[] x;

		switch (nProb) {
		case (1):
			x = new double[] { 0.1000000000000000E+01, 0.1000000000000000E+01 };
			return x;
		case (2):
			x = new double[] { 0.1000000000000000E+01, 0.9999999999999971E+00 };
			return x;
		case (3):
			x = new double[] { 0.1000000000000000E+01, 0.1000000000000000E+01 };
			return x;
		case (4):
			x = new double[] { 0.2695066833053901E-17, -0.2695066833053901E-18, 0.3749757049013966E-17,
					0.3749757049013966E-17 };
			return x;
		case (5):
			x = new double[] { 0.2961047510523875E-17, -0.2961047510523875E-18, 0.3392094461844849E-17,
					0.3392094461844849E-17 };
			return x;
		case (6):
			x = new double[] { -0.1125674354180310E-16, 0.1125674354180310E-17, -0.6482313483583572E-17,
					-0.6482313483583572E-17 };
			return x;
		case (7):
			x = new double[] { 0.1098159327798296E-04, 0.9106146740038449E+01 };
			return x;
		case (8):
			x = new double[] { 0.1098159288132353E-04, 0.9106146743612113E+01 };
			return x;
		case (9):
			x = new double[] { -0.9679740249512818E+00, 0.9471391408444433E+00, -0.9695163103175638E+00,
					0.9512476657652131E+00 };
			return x;
		case (10):
			x = new double[] { -0.9679740250937678E+00, 0.9471391411205137E+00, -0.9695163101751340E+00,
					0.9512476654887939E+00 };
			return x;
		case (11):
			x = new double[] { -0.9679740249230577E+00, 0.9471391407896745E+00, -0.9695163103461089E+00,
					0.9512476658204914E+00 };
			return x;
		case (12):
			x = new double[] { 0.1000000000000010E+01, -0.1612103288815217E-13, 0.0000000000000000E+00 };
			return x;
		case (13):
			x = new double[] { 0.1000000000004274E+01, 0.1388617146084722E-10, 0.0000000000000000E+00 };
			return x;
		case (14):
			x = new double[] { 0.1000000000050375E+01, 0.3234657151609015E-10, 0.0000000000000000E+00 };
			return x;
		case (15):
			x = new double[] { -0.1572508640131874E-01, 0.1012434869369120E+01, -0.2329916259568086E+00,
					0.1260430087800509E+01, -0.1513728922723659E+01, 0.9929964324319899E+00 };
			return x;
		case (16):
			x = new double[] { -0.1572508640142940E-01, 0.1012434869369113E+01, -0.2329916259567570E+00,
					0.1260430087799810E+01, -0.1513728922722583E+01, 0.9929964324313200E+00 };
			return x;
		case (17):
			x = new double[] { -0.1530703902928690E-04, 0.9997897039319343E+00, 0.1476396368990648E-01,
					0.1463423283093984E+00, 0.1000821102959179E+01, -0.2617731140413794E+01, 0.4104403164336241E+01,
					-0.3143612278455328E+01, 0.1052626407979455E+01 };
			return x;
		case (18):
			x = new double[] { 0.3119911862752667E+00, 0.1089342426782224E+01, 0.1379293238883866E+01,
					-0.1175426477217227E+02, 0.6265967591560398E+02, -0.1636004875532691E+03, 0.2326239209287591E+03,
					-0.1697823090157466E+03, 0.5076761702794597E+02 };
			return x;
		case (19):
			x = new double[] { 0.8375125649943561E-01, 0.3127292952232932E+00, 0.5000000000000018E+00,
					0.6872707047767043E+00, 0.9162487435005652E+00 };
			return x;
		case (20):
			x = new double[] { 0.8375125650172599E-01, 0.4999999999984706E+00, 0.3127292952187522E+00,
					0.9162487434957841E+00, 0.6872707047852671E+00 };
			return x;
		case (21):
			x = new double[] { 0.6872707048164862E+00, 0.4999999999708427E+00, 0.8375125649245271E-01,
					0.3127292952359536E+00, 0.9162487434842648E+00 };
			return x;
		case (22):
			x = new double[] { 0.6687659094732229E-01, 0.3666822992416477E+00, 0.2887406731168656E+00,
					0.7112593268831344E+00, 0.6333177007583523E+00, 0.9331234090526778E+00 };
			return x;
		case (23):
			x = new double[] { 0.9331234090548961E+00, 0.3666822992072173E+00, 0.6687659093602476E-01,
					0.7112593268461027E+00, 0.2887406731541705E+00, 0.6333177008015886E+00 };
			return x;
		case (24):
			x = new double[] { 0.6687659094738703E-01, 0.7112593269091415E+00, 0.2887406731888937E+00,
					0.9331234090619587E+00, 0.6333177007462377E+00, 0.3666822991463812E+00 };
			return x;
		case (25):
			x = new double[] { 0.5806914960717641E-01, 0.2351716124057688E+00, 0.3380440947078516E+00,
					0.4999999999993442E+00, 0.6619559052930899E+00, 0.7648283875938304E+00, 0.9419308503929387E+00 };
			return x;
		case (26):
			x = new double[] { 0.3380440947234102E+00, 0.2351716123658910E+00, 0.9419308503765472E+00,
					0.7648283876531939E+00, 0.6619559052401212E+00, 0.5806914961691421E-01, 0.5000000000239224E+00 };
			return x;
		case (27):
			x = new double[] { -0.4649136147100040E+02, -0.1039521377302678E+02, -0.7517577291666804E+01,
					0.1217934366521229E+02, 0.4297678475778002E+02, 0.2373113109255284E+02, 0.4306103920017470E+02 };
			return x;
		case (28):
			x = new double[] { 0.4985639998615101E-01, 0.1986351265136680E+00, 0.2698288232492614E+00,
					0.4992722954902931E+00, 0.5007277045097069E+00, 0.7301711767507386E+00, 0.8013648734863320E+00,
					0.9501436000138489E+00 };
			return x;
		case (29):
			x = new double[] { 0.4420534627964380E-01, 0.1994906715249946E+00, 0.2356191097165776E+00,
					0.4160469066798081E+00, 0.5000000008711383E+00, 0.5839530917260407E+00, 0.7643808911354805E+00,
					0.8005093282727184E+00, 0.9557946537935981E+00 };
			return x;
		case (30):
			x = new double[] { 0.1000000000000006E+01, 0.1000000000000006E+01, 0.1000000000000006E+01,
					0.1000000000000006E+01, 0.1000000000000006E+01, 0.1000000000000006E+01, 0.1000000000000006E+01,
					0.1000000000000006E+01, 0.1000000000000006E+01, 0.9999999999999419E+00 };
			return x;
		case (31):
			x = new double[] { 0.1000000000000000E+01, 0.1000000000000000E+01, 0.1000000000000000E+01,
					0.1000000000000000E+01, 0.1000000000000000E+01, 0.1000000000000000E+01, 0.1000000000000000E+01,
					0.1000000000000000E+01, 0.1000000000000000E+01, 0.9999999999999953E+00 };
			return x;
		case (32):
			x = new double[] { 0.9794303033498583E+00, 0.9794303033498583E+00, 0.9794303033498583E+00,
					0.9794303033498583E+00, 0.9794303033498583E+00, 0.9794303033498583E+00, 0.9794303033498583E+00,
					0.9794303033498583E+00, 0.9794303033498583E+00, 0.1205696966501417E+01 };
			return x;
		case (33):
			x = new double[] { 0.9999999999999868E+00, 0.9999999999999868E+00, 0.9999999999999868E+00,
					0.9999999999999868E+00, 0.9999999999999868E+00, 0.9999999999999868E+00, 0.9999999999999868E+00,
					0.9999999999999868E+00, 0.9999999999999868E+00, 0.9999999999999868E+00, 0.9999999999999868E+00,
					0.9999999999999868E+00, 0.9999999999999868E+00, 0.9999999999999868E+00, 0.9999999999999868E+00,
					0.9999999999999868E+00, 0.9999999999999868E+00, 0.9999999999999868E+00, 0.9999999999999868E+00,
					0.9999999999999868E+00, 0.9999999999999868E+00, 0.9999999999999868E+00, 0.9999999999999868E+00,
					0.9999999999999868E+00, 0.9999999999999868E+00, 0.9999999999999868E+00, 0.9999999999999868E+00,
					0.9999999999999868E+00, 0.9999999999999868E+00, 0.1000000000000383E+01 };
			return x;
		case (34):
			x = new double[] { 0.1000000000000006E+01, 0.1000000000000006E+01, 0.1000000000000006E+01,
					0.1000000000000006E+01, 0.1000000000000006E+01, 0.1000000000000006E+01, 0.1000000000000006E+01,
					0.1000000000000006E+01, 0.1000000000000006E+01, 0.1000000000000006E+01, 0.1000000000000006E+01,
					0.1000000000000006E+01, 0.1000000000000006E+01, 0.1000000000000006E+01, 0.1000000000000006E+01,
					0.1000000000000006E+01, 0.1000000000000006E+01, 0.1000000000000006E+01, 0.1000000000000006E+01,
					0.1000000000000006E+01, 0.1000000000000006E+01, 0.1000000000000006E+01, 0.1000000000000006E+01,
					0.1000000000000006E+01, 0.1000000000000006E+01, 0.1000000000000006E+01, 0.1000000000000006E+01,
					0.1000000000000006E+01, 0.1000000000000006E+01, 0.1000000000000006E+01, 0.1000000000000006E+01,
					0.1000000000000006E+01, 0.1000000000000006E+01, 0.1000000000000006E+01, 0.1000000000000006E+01,
					0.1000000000000006E+01, 0.1000000000000006E+01, 0.1000000000000006E+01, 0.1000000000000006E+01,
					0.9999999999997405E+00 };
			return x;
		case (35):
			x = new double[] { -0.4316498251876486E-01, -0.8157715653538729E-01, -0.1144857143805310E+00,
					-0.1409735768625996E+00, -0.1599086961819857E+00, -0.1698772023127759E+00, -0.1690899837812081E+00,
					-0.1552495352218312E+00, -0.1253558916789345E+00, -0.7541653368589182E-01 };
			return x;
		case (36):
			x = new double[] { -0.4316498251881878E-01, -0.8157715653546950E-01, -0.1144857143805966E+00,
					-0.1409735768626191E+00, -0.1599086961819499E+00, -0.1698772023126901E+00, -0.1690899837811062E+00,
					-0.1552495352217907E+00, -0.1253558916789970E+00, -0.7541653368596339E-01 };
			return x;
		case (37):
			x = new double[] { -0.4316498254522300E-01, -0.8157715658124411E-01, -0.1144857144024344E+00,
					-0.1409735768722910E+00, -0.1599086963002226E+00, -0.1698772022538783E+00, -0.1690899837944877E+00,
					-0.1552495352060589E+00, -0.1253558916432355E+00, -0.7541653366610668E-01 };
			return x;
		case (38):
			x = new double[] { -0.1528138835625800E+00 };
			return x;
		case (39):
			x = new double[] { -0.1528138835625801E+00 };
			return x;
		case (40):
			x = new double[] { -0.1528138835625800E+00 };
			return x;
		case (41):
			x = new double[] { -0.4316498251876487E-01, -0.8157715653538729E-01, -0.1144857143805310E+00,
					-0.1409735768625996E+00, -0.1599086961819857E+00, -0.1698772023127759E+00, -0.1690899837812080E+00,
					-0.1552495352218311E+00, -0.1253558916789344E+00, -0.7541653368589175E-01 };
			return x;
		case (42):
			x = new double[] { -0.4316498251881876E-01, -0.8157715653546944E-01, -0.1144857143805966E+00,
					-0.1409735768626190E+00, -0.1599086961819498E+00, -0.1698772023126901E+00, -0.1690899837811062E+00,
					-0.1552495352217907E+00, -0.1253558916789970E+00, -0.7541653368596334E-01 };
			return x;
		case (43):
			x = new double[] { -0.4316498251876519E-01, -0.8157715653538752E-01, -0.1144857143805303E+00,
					-0.1409735768625981E+00, -0.1599086961819844E+00, -0.1698772023127748E+00, -0.1690899837812073E+00,
					-0.1552495352218307E+00, -0.1253558916789341E+00, -0.7541653368589157E-01 };
			return x;
		case (44):
			x = new double[] { 0.5526154410956152E-01, 0.5695755464278050E-01, 0.5889066688359926E-01,
					0.6113606294453851E-01, 0.6377855838223215E-01, 0.6700482468993153E-01, 0.2079410281099167E+00,
					0.1642671068222599E+00, 0.8643947917821507E-01, 0.9133506311839429E-01 };
			return x;
		case (45):
			x = new double[] { 0.3439628896235700E-01, 0.3503231575416286E-01, 0.3571919583574922E-01,
					0.3646522422002401E-01, 0.3728091174083832E-01, 0.3817986258974627E-01, 0.3918014109818273E-01,
					0.4030650261421058E-01, 0.1797201916815176E+00, 0.1562408814749914E+00 };
			return x;
		case (46):
			x = new double[] { 0.1888395221037194E+02, 0.2516777354435942E+02, 0.1888527511766734E+02,
					0.1888602114534231E+02, 0.1888683683396955E+02, 0.1888773578345518E+02, 0.1888873606137402E+02,
					0.1888986242382218E+02, 0.1902927611206524E+02, 0.1900579680367850E+02 };
			return x;
		case (47):
			x = new double[] { 0.9999999999999993E+00, 0.9999999999999986E+00, 0.9999999999999979E+00,
					0.9999999999999972E+00, 0.9999999999999964E+00, 0.9999999999999958E+00, 0.9999999999999950E+00,
					0.9999999999999943E+00, 0.9999999999999937E+00, 0.9999999999999929E+00 };
			return x;
		case (48):
			x = new double[] { 0.1000000000000017E+01, 0.1000000000000033E+01, 0.1000000000000050E+01,
					0.1000000000000066E+01, 0.1000000000000083E+01, 0.1000000000000100E+01, 0.1000000000000116E+01,
					0.1000000000000133E+01, 0.1000000000000150E+01, 0.1000000000000166E+01 };
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
			x = new double[] { -0.5707221320993724E+00, -0.6818069494976288E+00, -0.7022100764193087E+00,
					-0.7055106298465235E+00, -0.7049061556842444E+00, -0.7014966070281702E+00, -0.6918893223746194E+00,
					-0.6657965143749420E+00, -0.5960351092084393E+00, -0.4164122574088334E+00 };
			return x;
		case (52):
			x = new double[] { -0.5707221320171143E+00, -0.6818069499829605E+00, -0.7022100760171542E+00,
					-0.7055106298955310E+00, -0.7049061557301967E+00, -0.7014966070327222E+00, -0.6918893223590803E+00,
					-0.6657965144072679E+00, -0.5960351090088830E+00, -0.4164122575177334E+00 };
			return x;
		case (53):
			x = new double[] { -0.4283028636053096E+00, -0.4765964242962532E+00, -0.5196524638125551E+00,
					-0.5580993246169653E+00, -0.5925061569509360E+00, -0.6245036821428090E+00, -0.6232394714478015E+00,
					-0.6213938418388717E+00, -0.6204535966122983E+00, -0.5864692707477790E+00 };
			return x;
		case (54):
			x = new double[] { -0.4283028634881692E+00, -0.4765964236396711E+00, -0.5196524642776768E+00,
					-0.5580993248351936E+00, -0.5925061568131795E+00, -0.6245036817962691E+00, -0.6232394720687791E+00,
					-0.6213938417874499E+00, -0.6204535965224117E+00, -0.5864692707287930E+00 };
			return x;
		case (55):
			x = new double[] { -0.4283028635608067E+00, -0.4765964243232715E+00, -0.5196524637037395E+00,
					-0.5580993248328234E+00, -0.5925061568292707E+00, -0.6245036822076749E+00, -0.6232394714256790E+00,
					-0.6213938418143938E+00, -0.6204535966527650E+00, -0.5864692707189498E+00 };
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
//	// This subroutine defines the jacobian matrices of fourteen
//	// test functions. the problem dimensions are as described
//	// in the prologue comments of vecfcn.
	public static void vecjac(int n, double[] x, double[][] fJac, int ldfJac, int nProb) {
//	    subroutine vecjac(n, x, fJac, ldfJac, nProb)
//	        implicit none
//
//	        int n // a positive integer variable.
//	        int ldfJac // a positive integer variable not less than n
//	                                     // which specifies the leading dimension of the array fJac.
//	        int nProb // a positive integer variable which defines the
//	                                    // number of the problem. nProb must not exceed 14.
//	        double x(n) // an array of length n.
//	        double fJac(ldfJac, n) // an n by n array. on output fJac contains the
//	                                                // jacobian matrix of the nProb function evaluated at x.

		double zero = 0.0;
		double one = 1.0;
		double two = 2.0;
		double three = 3.0;
		double four = 4.0;
		double five = 5.0;
		double six = 6.0;
		double eight = 8.0;
		double ten = 10.0;
		double fiftn = 15.0;
		double twenty = 20.0;
		double hundrd = 100.0;
		double c1 = 1.0e4;
		double c3 = 2.0e2;
		double c4 = 2.02e1;
		double c5 = 1.98e1;
		double c6 = 1.8e2;
		double c9 = 2.9e1;

		int i, j, k, k1, k2, ml, mu;
		double h, prod, sum, sum1, sum2, temp, temp1, temp2, temp3, temp4, ti, tj, tk, tpi;

//	        fJac(1:n,1:n) = zero
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				fJac[i][j] = 0;
			}
		}
		// jacobian routine selector.

		switch (nProb) {
		case (2):
			// powell Math.singular function.
			for (k = 0; k < 4; k++) {
				for (j = 0; j < 4; j++) {
					fJac[k][j] = zero;
				}
			}
			fJac[0][0] = one;
			fJac[0][1] = ten;
			fJac[1][2] = Math.sqrt(five);
			fJac[1][3] = -fJac[1][2];
			fJac[2][1] = two * (x[1] - two * x[2]);
			fJac[2][2] = -two * fJac[2][1];
			fJac[3][0] = two * Math.sqrt(ten) * (x[0] - x[3]);
			fJac[3][3] = -fJac[3][0];
			break;
		case (3):
			// powell badly scaled function.
			fJac[0][0] = c1 * x[1];
			fJac[0][1] = c1 * x[0];
			fJac[1][0] = -Math.exp(-x[0]);
			fJac[1][1] = -Math.exp(-x[1]);
			break;
		case (4):
			// wood function.
			for (k = 0; k < 4; k++) {
				for (j = 0; j < 4; j++) {
					fJac[k][j] = zero;
				}
			}
			temp1 = x[1] - three * x[0] * x[0];
			temp2 = x[3] - three * x[2] * x[2];
			fJac[0][0] = -c3 * temp1 + one;
			fJac[0][1] = -c3 * x[0];
			fJac[1][0] = -two * c3 * x[0];
			fJac[1][1] = c3 + c4;
			fJac[1][3] = c5;
			fJac[2][2] = -c6 * temp2 + one;
			fJac[2][3] = -c6 * x[2];
			fJac[3][1] = c5;
			fJac[3][2] = -two * c6 * x[2];
			fJac[3][3] = c6 + c4;
			break;
		case (5):
			// helical valley function.
			tpi = eight * Math.atan(one);
			temp = x[0] * x[0] + x[1] * x[1];
			temp1 = tpi * temp;
			temp2 = Math.sqrt(temp);
			fJac[0][0] = hundrd * x[1] / temp1;
			fJac[0][1] = -hundrd * x[0] / temp1;
			fJac[0][2] = ten;
			fJac[1][0] = ten * x[0] / temp2;
			fJac[1][1] = ten * x[1] / temp2;
			fJac[1][2] = zero;
			fJac[2][0] = zero;
			fJac[2][1] = zero;
			fJac[2][2] = one;
			break;
		case (6):
			// watson function.
			for (k = 0; k < n; k++) {
				for (j = k; j < n; j++) {
					fJac[k][j] = zero;
				}
			}
			for (i = 0; i < 29; i++) {
				ti = dfloat(i+1) / c9;
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
				temp1 = two * (sum1 - sum2 * sum2 - one);
				temp2 = two * sum2;
				temp = ti * ti;
				tk = one;
				for (k = 0; k < n; k++) {
					tj = tk;
					for (j = k; j < n; j++) {
						fJac[k][j] = fJac[k][j]
								+ tj * ((dfloat(k) / ti - temp2) * (dfloat(j) / ti - temp2) - temp1);
						tj = ti * tj;
					}
					tk = temp * tk;
				}
			}
			fJac[0][0] = fJac[0][0] + six * x[0] * x[0] - two * x[1] + three;
			fJac[0][1] = fJac[0][1] - two * x[0];
			fJac[1][1] = fJac[1][1] + one;
			for (k = 0; k < n; k++) {
				for (j = k; j < n; j++) {
					fJac[j][k] = fJac[k][j];
				}
			}
			break;
		case (7):
			// chebyquad function.
			tk = one / dfloat(n);
			for (j = 0; j < n; j++) {
				temp1 = one;
				temp2 = two * x[j] - one;
				temp = two * temp2;
				temp3 = zero;
				temp4 = two;
				for (k = 0; k < n; k++) {
					fJac[k][j] = tk * temp4;
					ti = four * temp2 + temp * temp4 - temp3;
					temp3 = temp4;
					temp4 = ti;
					ti = temp * temp2 - temp1;
					temp1 = temp2;
					temp2 = ti;
				}
			}
			break;
		case (8):
			// brown almost-linear function.
			prod = one;
			for (j = 0; j < n; j++) {
				prod = x[j] * prod;
				for (k = 0; k < n; k++) {
					fJac[k][j] = one;
				}
				fJac[j][j] = two;
			}
			for (j = 0; j < n; j++) {
				temp = x[j];
				if (temp == zero) {
					temp = one;
					prod = one;
					for (k = 0; k < n; k++) {
						if (k != j)
							prod = x[k] * prod;
					}
				}
				fJac[n][j] = prod / temp;
			}
			break;
		case (9):
			// discrete boundary value function.
			h = one / dfloat(n + 1);
			for (k = 0; k < n; k++) {
				temp = three * (x[k] + dfloat(k+1) * h + one) * (x[k] + dfloat(k+1) * h + one);
				for (j = 0; j < n; j++) {
					fJac[k][j] = zero;
				}
				fJac[k][k] = two + temp * h * h / two;
				if (k != 0)
					fJac[k][k - 1] = -one;
				if (k != n-1)
					fJac[k][k + 1] = -one;
			}
			break;
		case (10):
			// discrete integral equation function.
			h = one / dfloat(n + 1);
			for (k = 0; k < n; k++) {
				tk = dfloat(k+1) * h;
				for (j = 0; j < n; j++) {
					tj = dfloat(j+1) * h;
					temp = three * (x[j] + tj + one) * (x[j] + tj + one);
					fJac[k][j] = h * Math.min(tj * (one - tk), tk * (one - tj)) * temp / two;
				}
				fJac[k][k] = fJac[k][k] + one;
			}
			break;
		case (11):
			// trigonometric function.
			for (j = 0; j < n; j++) {
				temp = Math.sin(x[j]);
				for (k = 0; k < n; k++) {
					fJac[k][j] = temp;
				}
				fJac[j][j] = dfloat(j + 2) * temp - Math.cos(x[j]);
			}
			break;
		case (12):
			// variably dimensioned function.
			sum = zero;
			for (j = 0; j < n; j++) {
				sum = sum + dfloat(j+1) * (x[j] - one);
			}
			temp = one + six * sum * sum;
			for (k = 0; k < n; k++) {
				for (j = k; j < n; j++) {
					fJac[k][j] = dfloat((k+1) * (j+1)) * temp;
					fJac[j][k] = fJac[k][j];
				}
				fJac[k][k] = fJac[k][k] + one;
			}
			break;
		case (13):
			// broyden tridiagonal function.
			for (k = 0; k < n; k++) {
				for (j = 0; j < n; j++) {
					fJac[k][j] = zero;
				}
				fJac[k][k] = three - four * x[k];
				if (k != 0) {
					fJac[k][k - 1] = -one;
				}
				if (k != n-1) {
					fJac[k][k + 1] = -two;
				}
			}
			break;
		case (14):
			// broyden banded function.
			ml = 5;
			mu = 1;
			for (k = 0; k < n; k++) {
				for (j = 0; j < n; j++) {
					fJac[k][j] = zero;
				}
				k1 = Math.max(1, k - ml);
				k2 = Math.min(k + mu, n);
				for (j = k1; j < k2; j++) {
					if (j != k) {
						fJac[k][j] = -(one + two * x[j]);
					}
				}
				fJac[k][k] = two + fiftn * x[k] * x[k];
			}
			break;
		default:
			// rosenbrock function.
			fJac[0][0] = -one;
			fJac[0][1] = zero;
			fJac[1][0] = -twenty * x[0];
			fJac[1][1] = ten;
			break;
		}

	} // vecjac
//	!*****************************************************************************************
//
//	!*****************************************************************************************
//	!>
//	// This subroutine specifies the standard starting points for
//	// the functions defined by subroutine vecfcn. the subroutine
//	// returns in x a multiple (factor) of the standard starting
//	// point. for the sixth function the standard starting point is
//	// zero, so in this case, if factor is not unity, { the
//	// subroutine returns the vector  x[j] = factor, j=1,...,n.

	public static void initpt(int n, double[] x, int nProb, double factor) {
//	    subroutine initpt(n, x, nProb, Factor)
//	        implicit none
//
//	        int n // a positive integer input variable.
//	        int nProb // a positive integer input variable which defines the
//	                                    // number of the problem. nProb must not exceed 14.
//	        double Factor // an input variable which specifies the multiple of
//	                                      // the standard starting point. if factor is unity, no
//	                                      // multiplication is performed.
//	        double x(n) // an output array of length n which contains the standard
//	                                     // starting point for problem nProb multiplied by factor.

		int j;
		double h, tj;

		double zero = 0.0;
		double half = 0.5;
		double one = 1.0;
		double three = 3.0;
		double c1 = 1.2;

//	        x(1:n) = zero
		for (i = 0; i < n; i++) {
			x[i] = 0;
		}

		// selection of initial point.

		switch (nProb) {
		case (2):
			// powell Math.singular function.
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
				x[j] = dfloat(j+1) * h;
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
				tj = dfloat(j+1) * h;
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
				x[j] = one - dfloat(j+1) * h;
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

	} // initpt
//	!*****************************************************************************************
//
//	!*****************************************************************************************
//	!>
//	// This subroutine defines fourteen test functions. the first
//	// five test functions are of dimensions 2,4,2,4,3, respectively,
//	// while the remaining test functions are of variable dimension
//	// n for any n greater than or equal to 1 (problem 6 is an
//	// exception to this, Math.since it does not allow n = 1).

	public static void vecfcn(int n, double[] x, double[] fVec, int nProb) {
//	    subroutine vecfcn(n, x, fVec, nProb)
//	        implicit none

//	        int n // a positive integer input variable.
//	        double[] x(n) // an input array of length n.
//	        double[] fVec[n] // an output array of length n which contains the nProb
//	                                        // function vector evaluated at x.
//	        int nProb // a positive integer input variable which defines the
		// number of the problem. nProb must not exceed 14.

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

		double tpi = eight * Math.atan(one);

		int i, iev, j, k, k1, k2, kp1, ml, mu;
		double h, prod, sum, sum1, sum2, temp, temp1, temp2, ti, tj, tk;

//	        fVec(1:n) = zero;
		for (i = 0; i < n; i++) {
			fVec[i] = 0;
		}
		// problem selector.

		switch (nProb) {
		case (2):
			// powell Math.singular function.
			fVec[0] = x[0] + ten * x[1];
			fVec[1] = Math.sqrt(five) * (x[2] - x[3]);
			fVec[2] = (x[1] - two * x[2]) * (x[1] - two * x[2]);
			fVec[3] = Math.sqrt(ten) * (x[0] - x[3]) * (x[0] - x[3]);
			break;
		case (3):
			// powell badly scaled function.
			fVec[0] = c1 * x[0] * x[1] - one;
			fVec[1] = Math.exp(-x[0]) + Math.exp(-x[1]) - c2;
			break;
		case (4):
			// wood function.
			temp1 = x[1] - x[0] * x[0];
			temp2 = x[3] - x[2] * x[2];
			fVec[0] = -c3 * x[0] * temp1 - (one - x[0]);
			fVec[1] = c3 * temp1 + c4 * (x[1] - one) + c5 * (x[3] - one);
			fVec[2] = -c6 * x[2] * temp2 - (one - x[2]);
			fVec[3] = c6 * temp2 + c4 * (x[3] - one) + c5 * (x[1] - one);
			break;
		case (5):
			// helical valley function.
			temp1 = Math.signum(x[1]) * Math.abs(c7); // sign(c7, x[1]);
			if (x[0] > zero)
				temp1 = Math.atan(x[1] / x[0]) / tpi;
			if (x[0] < zero)
				temp1 = Math.atan(x[1] / x[0]) / tpi + c8;
			temp2 = Math.sqrt(x[0] * x[0] + x[1] * x[1]);
			fVec[0] = ten * (x[2] - ten * temp1);
			fVec[1] = ten * (temp2 - one);
			fVec[2] = x[2];
			break;
		case (6):
			// watson function.
			for (k = 0; k < n; k++) {
				fVec[k] = zero;
			}
			for (i = 0; i < 29; i++) {
				ti = dfloat(i+1) / c9;
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
			// chebyquad function.
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
				if (iev > 0) {
					fVec[k] = fVec[k] + one / (dfloat(k+1) * dfloat(k+1) - one);
				}
				iev = -iev;
			}
			break;
		case (8):
			// brown almost-linear function.
			sum = -dfloat(n + 1);
			prod = one;
			for (j = 0; j < n; j++) {
				sum = sum + x[j];
				prod = x[j] * prod;
			}
			for (k = 0; k < n; k++) {
				fVec[k] = x[k] + sum;
			}
			fVec[n] = prod - one;
			break;
		case (9):
			// discrete boundary value function.
			h = one / dfloat(n + 1);
			for (k = 0; k < n; k++) {
				temp = (x[k] + dfloat(k+1) * h + one) * (x[k] + dfloat(k+1) * h + one) * (x[k] + dfloat(k+1) * h + one);
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
			// discrete integral equation function.
			h = one / dfloat(n + 1);
			for (k = 0; k < n; k++) {
				tk = dfloat(k+1) * h;
				sum1 = zero;
				for (j = 0; j < k; j++) {
					tj = dfloat(j+1) * h;
					temp = (x[j] + tj + one) * (x[j] + tj + one) * (x[j] + tj + one);
					sum1 = sum1 + tj * temp;
				}
				sum2 = zero;
				kp1 = k + 1;
				if (n >= kp1) {
					for (j = kp1; j < n; j++) {
						tj = dfloat(j+1) * h;
						temp = (x[j] + tj + one) * (x[j] + tj + one) * (x[j] + tj + one);
						sum2 = sum2 + (one - tj) * temp;
					}
				}
				fVec[k] = x[k] + h * ((one - tk) * sum1 + tk * sum2) / two;
			}
			break;
		case (11):
			// trigonometric function.
			sum = zero;
			for (j = 0; j < n; j++) {
				fVec[j] = Math.cos(x[j]);
				sum = sum + fVec[j];
			}
			for (k = 0; k < n; k++) {
				fVec[k] = dfloat(n + k+1) - Math.sin(x[k]) - sum - dfloat(k+1) * fVec[k];
			}
			break;
		case (12):
			// variably dimensioned function.
			sum = zero;
			for (j = 0; j < n; j++) {
				sum = sum + dfloat(j+1) * (x[j] - one);
			}
			temp = sum * (one + two * sum * sum);
			for (k = 0; k < n; k++) {
				fVec[k] = x[k] - one + dfloat(k+1) * temp;
			}
			break;
		case (13):
			// broyden tridiagonal function.
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
			// broyden banded function.
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
			// rosenbrock function.
			fVec[0] = one - x[0];
			fVec[1] = ten * (x[1] - x[0] * x[0]);
			break;
		}

	} // vecfcn

//	end program test_hybrj
}
