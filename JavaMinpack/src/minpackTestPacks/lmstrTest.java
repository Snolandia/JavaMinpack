package minpackTestPacks;

import java.util.Arrays;

import Minpack.Minpack;

public class lmstrTest {
//	!*****************************************************************************************
//	!>
//	!  This program tests codes for the least-squares solution of
//	!  m nonlinear equations in n variables. it consists of a driver
//	!  and an interface subroutine fcn. the driver reads in data,
//	!  calls the nonlinear least-squares solver, and finally prints
//	!  out information on the performance of the solver. this is
//	!  only a sample driver, many other drivers are possible. the
//	!  interface subroutine fcn is necessary to take into account the
//	!  forms of calling sequences used by the function and jacobian
//	!  subroutines in the various nonlinear least-squares solvers.
//
	static int ncases = 28;
	static int[] nProbs  = new int[] {1,1,2,2,3,3,4,5,6,7,8,9,10,11,11,11,12,13,14,15,15,15,15,16,16,16,17,18};
	static int[] ns      = new int[] {5,5,5,5,5,5,2,3,4,2,3,4,3,6,9,12,3,2,4,1,8,9,10,10,30,40,5,11};
	static int[] ms      =  new int[] {10,50,10,50,10,50,2,3,4,2,15,11,16,31,31,31,10,10,20,8,8,9,10,10,30,40,33,65};
	static int[] ntriess =  new int[] {1,1,1,1,1,1,3,3,3,3,3,3,2,3,3,3,1,1,3,3,1,1,1,3,1,1,1,1};

	static int[] info_original =  new int[] {2,3,1,1,1,1,4,2,2,2,2,2,4,4,4,1,1,1,1,1,1,1,1,5,3,5,
                                                        1,1,1,3,3,2,3,2,3,2,1,1,1,1,4,1,1,1,2,1,2,2,2,2,2,1,1};
                                                        //original `info` from the original minpack

	static int[] iwa;
	static double[] wa;
	static double[][] fjac;
	static double[] fVec;
	static double[] x;
	static int i, ic, info, k, ldfJac, lwa, m, n, nFev, nJev, nProb, ntries, icase;
	static int[] ma = new int[53] , na =  new int[53] , nf =  new int[53] , nj =  new int[53] , np = new int[53] , nx = new int[53]; 
	static double[] fnm = new double[53];
	static double factor, fnorm1, fnorm2;

	static  double one = 1.0;
	static double ten = 10.0;
	static double tol = Math.sqrt(1.0e-16);//dpmpar(1))
	static double solution_abstol = 1.0e-5;    //abstol for matching previously generated solutions
	static double solution_reltol = 1.0e-4;    //reltol for matching previously generated solutions
	
 
		public static void testLmstr(){
//	    use minpack_module, only: wp, dpmpar, enorm, lmstr1
//	    use iso_fortran_env, only: nwrite => output_unit
//
//	    implicit none
//
//	    //originally from file22
	    

	    ic = 0;
	    for(icase = 0;icase< ncases+1;icase++){

	        if (icase == ncases+1) {
	        	System.out.println("(A,I3,A)" + "1SUMMARY OF " + ic + " CALLS TO LMSTR1");
	        	System.out.println("(A/)" +    " nProb   N    M   nFev  nJev  INFO  FINAL L2 NORM");
	            for(i=0;i<ic;i++){ 
	                System.out.println("(3I5, 3I6, 1X, D15.7)"+ np[i] + "  " +  na[i] + "  " +  ma[i] + "  " +  nf[i] + "  " +  nj[i] + "  " + nx[i] + "  " + fnm[i]);
	            }
	            break; //stop
	        }else{
	            nProb = nProbs[icase];
	            n = ns[icase];
	            m = ms[icase];
	            lwa = 5*n+m;
	            ldfJac = n;
				iwa = new int[n];
				wa = new double[lwa];
				double[][] fJac = new double[n][n];
				fVec = new double[m];
				x = new double[n];

	            ntries = ntriess[icase];
	            factor = one;
	            for(k = 0;k<ntries;k++){
	                ic = ic + 1;
	                initpt(n, x, nProb, factor);
	                ssqfcn(m, n, x, fVec, nProb);
	                fnorm1 = Minpack.enorm(m, fVec);
	                System.out.println( "(////5X,A,I5,5X,A,2I5,5X//)" + " PROBLEM  " + nProb + " DIMENSIONS " + n + " " + m);
	                nFev = 0;
	                nJev = 0;
	                Minpack.lmstr1(fcn, m, n, x, fVec, fjac, ldfJac, tol, info, iwa, wa, lwa);
	                ssqfcn(m, n, x, fVec, nProb);
	                fnorm2 = Minpack.enorm(m, fVec);
	                np[ic] = nProb;
	                na[ic] = n;
	                ma[ic] = m;
	                nf[ic] = nFev;
	                nj[ic] = nJev;
	                nx[ic] = info;
	                fnm[ic] = fnorm2;
	                System.out.println("(5x,a,d15.7//5x,a,d15.7//5x,a,i10//5x,a,i10//5x,a,18x,i10//5x,a//*(5x,5d15.7/))" +
	                        " INITIAL L2 NORM OF THE RESIDUALS"+ fnorm1+
	                        " FINAL L2 NORM OF THE RESIDUALS  "+ fnorm2+ 
	                        " NUMBER OF FUNCTION EVALUATIONS  "+ nFev+ 
	                        " NUMBER OF JACOBIAN EVALUATIONS  "+ nJev+   
	                        " EXIT PARAMETER " + info + 
	                        " FINAL APPROXIMATE SOLUTION " + x[1]);
	                factor = ten*factor;
	                compareSolutions(ic, x, solution_reltol, solution_abstol);
	            }
	        }
	    }
		}
//	    contains
//	
//
			public static void compareSolutions(int ic, double[] x, double relTol, double absTol){

				double[] diff = new double[x.length], absDiff = new double[x.length], relDiff = new double[x.length], icF;
		//
			    if (info_original[ic]<5) {//    Ignore any where the original minpack failed
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
			    }

			}

//	!*****************************************************************************************
//	!>
//	!  The calling sequence of fcn should be identical to the
//	!  calling sequence of the function subroutine in the nonlinear
//	!  least squares solver. if iflag = 1, fcn should only call the
//	!  testing function subroutine ssqfcn. if iflag = i, i >= 2,
//	!  fcn should only call subroutine ssqjac to calculate the
//	!  (i-1)-st row of the jacobian. (the ssqjac subroutine provided
//	!  here for testing purposes calculates the entire jacobian
//	!  matrix and is therefore called only when iflag = 2.) each
//	!  call to ssqfcn or ssqjac should specify the appropriate
//	!  value of problem number (nProb).
	public static void fcn(int m, int n, double[] x, double[] fVec, double[] fJrom, int iFlag){
//	subroutine fcn(m, n, x, fVec, Fjrow, Iflag)
//	    implicit none
//
//	    int m
//	    int n
//	    integer,intent(inout) :: Iflag
//	    real(wp),intent(in) :: x(n)
//	    real(wp),intent(inout) :: fVec[m]
//	    real(wp),intent(inout) :: Fjrow(n)
//
	    int ldfJac = 65;
//
	    int j;
		double[][] temp = new double[ldfJac][40];
//	    real(wp), save :: temp(ldfJac, 40)
//	        !//this array is filled when FCN is called with IFLAG=2.
//	        !//When FCN is called with IFLAG=2,3,..., the argument array
//	        !//FJROW is filled with a row of TEMP. This will work only if
//	        !//TEMP is given the SAVE attribute, which was not done in
//	        !//the original code.
//
//	    switch(iFlag){
//	    case(1):
//	        ssqfcn(m, n, x, fVec, nProb);
//	        nFev = nFev + 1;
//	        break;
//	    case(2):
//	        //populate the temp array
//	        ssqjac(m, n, x, temp, ldfJac, nProb);
//	        nJev = nJev + 1;
//		break;
//	    }
//
//	    //for iflag = 2,3,... get the row from temp
//	    for(j=0;j<n;j++){
//	        fJrow[j] = temp(iFlag - 1, j);
//	    }
}
//	end subroutine fcn
//	!*****************************************************************************************
//
//	!*****************************************************************************************
//	!>
//	!  Replaced statement function in original code.
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
//	!  Get expected `x` vectors for each case.
		public static double[] solution(int nProb){
//	    pure function solution(nProb) result(x)
//
//	        implicit none
//			
			double[] x;
//	        integer,intent(in) :: nProb
//	        real(wp),dimension(:),allocatable :: x
//
	        switch (nProb){
	        case (  1): x = new double[] {-0.1000000000000000E+01,-0.1000000000000000E+01,-0.1000000000000000E+01,
	                         -0.1000000000000000E+01,-0.9999999999999992E+00};
			break;
	        case (  2): x = new double[] {-0.1000000000000001E+01,-0.1000000000000000E+01,-0.9999999999999960E+00,
	                         -0.9999999999999982E+00,-0.1000000000000001E+01};
			break;
	        case (  3): x = new double[] { 0.2561515109810807E+03,-0.1557989602750597E+02, 0.1294492228753740E+03,
	                         -0.7289948013753003E+01,-0.1168073476708643E+03};
			break;
	        case (  4): x = new double[] { 0.9059979768060809E+02, 0.7667581483950413E+02, 0.8463689104480578E+02,
	                          0.3883790741975204E+02,-0.1306368054405490E+03};
			break;
	        case (  5): x = new double[] { 0.1000000000000000E+01, 0.2146513294878052E+03,-0.1024194005964338E+03,
	                         -0.3046699664951840E+02, 0.1000000000000000E+01};
			break;
	        case (  6): x = new double[] { 0.1000000000000000E+01,-0.1453732603711001E+03, 0.1438958372244028E+03,
	                         -0.3522751577398922E+02, 0.1000000000000000E+01};
			break;
	        case (  7): x = new double[] { 0.1000000000000000E+01, 0.1000000000000000E+01};
			break;
	        case (  8): x = new double[] { 0.1000000000000000E+01, 0.1000000000000000E+01};
			break;
	        case (  9): x = new double[] { 0.1000000000000000E+01, 0.1000000000000000E+01};
			break;
	        case ( 10): x = new double[] { 0.1000000000000000E+01,-0.6243302423970055E-17, 0.0000000000000000E+00};
			break;
	        case ( 11): x = new double[] { 0.1000000000000000E+01, 0.6563910805155555E-20, 0.0000000000000000E+00};
			break;
	        case ( 12): x = new double[] { 0.1000000000000000E+01,-0.1972152263052530E-29, 0.0000000000000000E+00};
			break;
	        case ( 13): x = new double[] { 0.1652117596168386E-16,-0.1652117596168386E-17, 0.2643388153869423E-17,
	                          0.2643388153869423E-17};
			break;
	        case ( 14): x = new double[] { 0.2065146995210483E-16,-0.2065146995210483E-17, 0.3304235192336774E-17,
	                          0.3304235192336774E-17};
			break;
	        case ( 15): x = new double[] { 0.1290716872006552E-16,-0.1290716872006552E-17, 0.2065146995210478E-17,
	                          0.2065146995210478E-17};
			break;
	        case ( 16): x = new double[] { 0.1141248446549926E+02,-0.8968279137315114E+00};
			break;
	        case ( 17): x = new double[] { 0.1141300466147463E+02,-0.8967960386859540E+00};
			break;
	        case ( 18): x = new double[] { 0.1141278178578780E+02,-0.8968051074920977E+00};
			break;
	        case ( 19): x = new double[] { 0.8241057657583328E-01, 0.1133036653471502E+01, 0.2343694638941156E+01};
			break;
	        case ( 20): x = new double[] { 0.8406666738183292E+00,-0.1588480332595652E+09,-0.1643786716535350E+09};
			break;
	        case ( 21): x = new double[] { 0.8406666738676454E+00,-0.1589461672055181E+09,-0.1644649068577709E+09};
			break;
	        case ( 22): x = new double[] { 0.1928078104762493E+00, 0.1912626533540697E+00, 0.1230528010469306E+00,
	                          0.1360532211505162E+00};
			break;
	        case ( 23): x = new double[] { 0.7286754737688021E+06,-0.1407588031293926E+02,-0.3297779778420302E+08,
	                         -0.2057159419780570E+08};
			break;
	        case ( 24): x = new double[] { 0.1928081061389676E+00, 0.1912560109262705E+00, 0.1230515496222497E+00,
	                          0.1360501459183326E+00};
			break;
	        case ( 25): x = new double[] { 0.5609636476947878E-02, 0.6181346345409401E+04, 0.3452236345945471E+03};
			break;
	        case ( 26): x = new double[] { 0.8796346510290913E-11, 0.3464217981622328E+05, 0.9148762438218395E+03};
			break;
	        case ( 27): x = new double[] {-0.1572496150837810E-01, 0.1012434882329655E+01,-0.2329917223876732E+00,
	                          0.1260431011028184E+01,-0.1513730313944207E+01, 0.9929972729184212E+00};
			break;
	        case ( 28): x = new double[]  {-0.1572519013866755E-01, 0.1012434858601051E+01,-0.2329915458438306E+00,
	                          0.1260429320891634E+01,-0.1513727767065756E+01, 0.9929957342632842E+00};
			break;
	        case ( 29): x = new double[] {-0.1572470197125892E-01, 0.1012434909256583E+01,-0.2329919227616438E+00,
	                          0.1260432929295550E+01,-0.1513733204527069E+01, 0.9929990192232202E+00};
			break;
	        case ( 30): x = new double[] {-0.1530706441665223E-04, 0.9997897039345965E+00, 0.1476396349111473E-01,
	                          0.1463423301458067E+00, 0.1000821094549045E+01,-0.2617731120707204E+01,
	                          0.4104403139436332E+01,-0.3143612262364274E+01, 0.1052626403788084E+01};
			break;
	        case ( 31): x = new double[] {-0.1530703649598812E-04, 0.9997897039319464E+00, 0.1476396369371614E-01,
	                          0.1463423282979424E+00, 0.1000821103011242E+01,-0.2617731140534388E+01,
	                          0.4104403164498381E+01,-0.3143612278569126E+01, 0.1052626408013491E+01};
			break;
	        case ( 32): x = new double[] {-0.1530703652139991E-04, 0.9997897039319481E+00, 0.1476396369357189E-01,
	                          0.1463423282991793E+00, 0.1000821103005637E+01,-0.2617731140521392E+01,
	                          0.4104403164482094E+01,-0.3143612278558673E+01, 0.1052626408010778E+01};
			break;
	        case ( 33): x = new double[] {-0.6602659327962938E-08, 0.1000001644118327E+01,-0.5639321470774273E-03,
	                          0.3478205400521048E+00,-0.1567315002542552E+00, 0.1052815158301577E+01,
	                         -0.3247271095327381E+01, 0.7288434784001158E+01,-0.1027184809891813E+02,
	                          0.9074113537386525E+01,-0.4541375419278813E+01, 0.1012011879768094E+01};
			break;
	        case ( 34): x = new double[] {-0.6638060464411604E-08, 0.1000001644117862E+01,-0.5639322103635158E-03,
	                          0.3478205405045170E+00,-0.1567315041009225E+00, 0.1052815177233611E+01,
	                         -0.3247271153548033E+01, 0.7288434898122469E+01,-0.1027184824156328E+02,
	                          0.9074113647268060E+01,-0.4541375466778046E+01, 0.1012011888568983E+01};
			break;
	        case ( 35): x = new double[] {-0.6637951126821797E-08, 0.1000001644117862E+01,-0.5639322100134774E-03,
	                          0.3478205405003073E+00,-0.1567315040668739E+00, 0.1052815177081149E+01,
	                         -0.3247271153132584E+01, 0.7288434897408339E+01,-0.1027184824078557E+02,
	                          0.9074113646748899E+01,-0.4541375466584788E+01, 0.1012011888538411E+01};
			break;
	        case ( 36): x = new double[] { 0.9999999999999999E+00, 0.1000000000000000E+02, 0.1000000000000000E+01};
			break;
	        case ( 37): x = new double[] { 0.2578199266368066E+00, 0.2578299767645460E+00};
			break;
	        case ( 38): x = new double[] {-0.1159124627062127E+02, 0.1320248655827936E+02,-0.4035748532873218E+00,
	                          0.2367362096235223E+00};
			break;
	        case ( 39): x = new double[] {-0.1159592742718941E+02, 0.1320418669261399E+02,-0.4034173629293835E+00,
	                          0.2367711433462491E+00};
			break;
	        case ( 40): x = new double[] {-0.1159026120687546E+02, 0.1320206344847849E+02,-0.4036926117131069E+00,
	                          0.2366618932408272E+00};
			break;
	        case ( 41): x = new double[] { 0.5000000000000000E+00};
			break;
	        case ( 42): x = new double[] { 0.9817314924683995E+00};
			break;
	        case ( 43): x = new double[] { 0.9817314852933997E+00};
			break;
	        case ( 44): x = new double[] { 0.4315366485873597E-01, 0.1930916378432680E+00, 0.2663285938126974E+00,
	                          0.4999993346289134E+00, 0.5000006653710866E+00, 0.7336714061873026E+00,
	                          0.8069083621567320E+00, 0.9568463351412640E+00};
	        case ( 45): x = new double[] { 0.4420534613578272E-01, 0.1994906723098810E+00, 0.2356191084710600E+00,
	                          0.4160469078925981E+00, 0.5000000000000000E+00, 0.5839530921074020E+00,
	                          0.7643808915289400E+00, 0.8005093276901191E+00, 0.9557946538642172E+00};
			break;
	        case ( 46): x = new double[] { 0.5962026717535887E-01, 0.1667087838059411E+00, 0.2391710188135120E+00,
	                          0.3988852903460409E+00, 0.3988836678709105E+00, 0.6011163321290895E+00,
	                          0.6011147096539592E+00, 0.7608289811864880E+00, 0.8332912161940590E+00,
	                          0.9403797328246412E+00};
			break;
	        case ( 47): x = new double[] { 0.9794303033498616E+00, 0.9794303033498616E+00, 0.9794303033498616E+00,
	                          0.9794303033498616E+00, 0.9794303033498616E+00, 0.9794303033498616E+00,
	                          0.9794303033498616E+00, 0.9794303033498616E+00, 0.9794303033498616E+00,
	                          0.1205696966501385E+01};
			break;
	        case ( 48): x = new double[] { 0.9794303033498644E+00, 0.9794303033498644E+00, 0.9794303033498644E+00,
	                          0.9794303033498644E+00, 0.9794303033498644E+00, 0.9794303033498644E+00,
	                          0.9794303033498644E+00, 0.9794303033498644E+00, 0.9794303033498644E+00,
	                          0.1205696966501354E+01};
			break;
	        case ( 49): x = new double[] { 0.9794303033498626E+00, 0.9794303033498626E+00, 0.9794303033498626E+00,
	                          0.9794303033498626E+00, 0.9794303033498626E+00, 0.9794303033498626E+00,
	                          0.9794303033498626E+00, 0.9794303033498626E+00, 0.9794303033498626E+00,
	                          0.1205696966501374E+01};
			break;
	        case ( 50): x = new double[] { 0.9977542164428299E+00, 0.9977542164428299E+00, 0.9977542164428299E+00,
	                          0.9977542164428299E+00, 0.9977542164428299E+00, 0.9977542164428299E+00,
	                          0.9977542164428299E+00, 0.9977542164428299E+00, 0.9977542164428299E+00,
	                          0.9977542164428299E+00, 0.9977542164428299E+00, 0.9977542164428299E+00,
	                          0.9977542164428299E+00, 0.9977542164428299E+00, 0.9977542164428299E+00,
	                          0.9977542164428299E+00, 0.9977542164428299E+00, 0.9977542164428299E+00,
	                          0.9977542164428299E+00, 0.9977542164428299E+00, 0.9977542164428299E+00,
	                          0.9977542164428299E+00, 0.9977542164428299E+00, 0.9977542164428299E+00,
	                          0.9977542164428299E+00, 0.9977542164428299E+00, 0.9977542164428299E+00,
	                          0.9977542164428299E+00, 0.9977542164428299E+00, 0.1067373506715082E+01};
			break;
	        case ( 51): x = new double[] { 0.1000000000000028E+01, 0.1000000000000028E+01, 0.1000000000000028E+01,
	                          0.1000000000000028E+01, 0.1000000000000028E+01, 0.1000000000000028E+01,
	                          0.1000000000000028E+01, 0.1000000000000028E+01, 0.1000000000000028E+01,
	                          0.1000000000000028E+01, 0.1000000000000028E+01, 0.1000000000000028E+01,
	                          0.1000000000000028E+01, 0.1000000000000028E+01, 0.1000000000000028E+01,
	                          0.1000000000000028E+01, 0.1000000000000028E+01, 0.1000000000000028E+01,
	                          0.1000000000000028E+01, 0.1000000000000028E+01, 0.1000000000000028E+01,
	                          0.1000000000000028E+01, 0.1000000000000028E+01, 0.1000000000000028E+01,
	                          0.1000000000000028E+01, 0.1000000000000028E+01, 0.1000000000000028E+01,
	                          0.1000000000000028E+01, 0.1000000000000028E+01, 0.1000000000000028E+01,
	                          0.1000000000000028E+01, 0.1000000000000028E+01, 0.1000000000000028E+01,
	                          0.1000000000000028E+01, 0.1000000000000028E+01, 0.1000000000000028E+01,
	                          0.1000000000000028E+01, 0.1000000000000028E+01, 0.1000000000000028E+01,
	                          0.9999999999988959E+00};
			break;
	        case ( 52): x = new double[] { 0.3754100492440249E+00, 0.1935846545431039E+01,-0.1464686767487118E+01,
	                          0.1286753391104384E-01, 0.2212270118130775E-01};
			break;
	        case ( 53): x = new double[] { 0.1309976638100963E+01, 0.4315524807599997E+00, 0.6336612616028594E+00,
	                          0.5994285609916951E+00, 0.7541797682724487E+00, 0.9043000823785183E+00,
	                          0.1365799495210074E+01, 0.4823731997481072E+01, 0.2398684751048711E+01,
	                          0.4568875547914517E+01, 0.5675342062730520E+01};
			break;
	        default:
	            x = new double[] {1};//error stop 'invalid case'
				break;
	        }
			return x;
}
//	end function solution
//	!*****************************************************************************************
//
//	!*****************************************************************************************
//	!>
//	!  This subroutine defines the jacobian matrices of eighteen
//	!  nonlinear least squares problems. the problem dimensions are
//	!  as described in the prologue comments of ssqfcn.
public static void ssqjac(int m, int n, double[] x, double[][] fJac, int ldfJac, int nProb) {
//	subroutine ssqjac(m, n, x, Fjac, ldfJac, nProb)
//	    implicit none
//
//	    integer,intent(in) :: m !! positive integer input variable.
//	    integer,intent(in) :: n !! positive integer input variable. n must not exceed m.
//	    integer,intent(in) :: ldfJac !! a positive integer input variable not less than m
//	                                 !! which specifies the leading dimension of the array fjac.
//	    integer,intent(in) :: nProb !! a positive integer variable which defines the
//	                                !! number of the problem. nProb must not exceed 18.
//	    real(wp),intent(in) :: x(n) !! an input array of length n.
//	    real(wp),intent(out) :: Fjac(ldfJac, n) !! an m by n output array which contains the jacobian
//	                                            !! matrix of the nProb function evaluated at x.
//
	    double zero = 0.0;
	    double one = 1.0;
	    double two = 2.0;
	    double three = 3.0;
	    double four = 4.0;
	    double five = 5.0;
	    double eight = 8.0;
	    double ten = 10.0;
	    double c14 = 14.0;
	    double c20 = 20.0;
	    double c29 = 29.0;
	    double c45 = 45.0;
	    double c100 = 100.0;

	    double[] v = {4.0, 2.0, 1.0, 5.0e-1, 2.5e-1, 1.67e-1,
	                                   1.25e-1, 1.0e-1, 8.33e-2, 7.14e-2, 6.25e-2};

	    int i, ivar, j, k, mm1, nm1;
	    double div, dx, prod, s2, temp, ti, tmp1, tmp2, tmp3, tmp4, tpi;

//	    Fjac(1:m, 1:n) = zero
	    fJac = new double[m][n];
	    for(i = 0;i<m;i++) {
	    	for(j = 0;i<n;i++) {
		    	fJac[i][j] = 0;
		    }
	    }

	    //JACOBIAN ROUTINE SELECTOR.

	    switch(nProb) {
	    case (2):
	        //LINEAR FUNCTION - RANK 1.
	        for(j=0;j<n;j++){
	            for(i=0;i<m;i++){
	                fJac[i][j] = dfloat(i+1)*dfloat(j+1);
	            }
	        }
	    break;
	    case (3):
	        //LINEAR FUNCTION - RANK 1 WITH ZERO COLUMNS AND ROWS.
	        for(j=0;j<n;j++){
	            for(i=0;i<m;i++){
	                fJac[i][j] = zero;
	            }
	        }
	        nm1 = n - 1;
	        mm1 = m - 1;
	        if (nm1 >= 2) {
	            for( j = 1;j< nm1;j++){
	                for( i = 1;i< mm1;i++){
	                    fJac[i][j] = dfloat(i)*dfloat(j+1);
	                }
	            }
	        }
		    break;
	    case (4):
	        //ROSENBROCK FUNCTION.
	        fJac[1][1] = -c20*x[0];
	        fJac[1][2] = ten;
	        fJac[2][1] = -one;
	        fJac[2][2] = zero;
		    break;
	    case (5):
	        //HELICAL VALLEY FUNCTION.
	        tpi = eight*Math.atan(one);
	        temp = x[0]*x[0] + x[1]*x[1];
	        tmp1 = tpi*temp;
	        tmp2 = Math.sqrt(temp);
	        fJac[1][1] = c100*x[1]/tmp1;
	        fJac[1][2] = -c100*x[0]/tmp1;
	        fJac[1][3] = ten;
	        fJac[2][1] = ten*x[0]/tmp2;
	        fJac[2][2] = ten*x[1]/tmp2;
	        fJac[2][3] = zero;
	        fJac[3][1] = zero;
	        fJac[3][2] = zero;
	        fJac[3][3] = one;
		    break;
	    case (6):
	        //POWELL SINGULAR FUNCTION.
	        for(j = 0;j< 4;j++){
	            for(i=0;i< 4;i++){
	                fJac[i][j] = zero;
	            }
	        }
	        fJac[1][1] = one;
	        fJac[1][2] = ten;
	        fJac[2][3] = Math.sqrt(five);
	        fJac[2][4] = -fJac[2][3];
	        fJac[3][2] = two*(x[1] - two*x[2]);
	        fJac[3][3] = -two*fJac[3][2];
	        fJac[4][1] = two*Math.sqrt(ten)*(x[0] - x[3]);
	        fJac[4][4] = -fJac[4][1];
		    break;
	    case (7):
	        //FREUDENSTEIN AND ROTH FUNCTION.
	        fJac[1][1] = one;
	        fJac[1][2] = x[1]*(ten - three*x[1]) - two;
	        fJac[2][1] = one;
	        fJac[2][2] = x[1]*(two + three*x[1]) - c14;
		    break;
	    case (8):
	        //BARD FUNCTION.
	        for(i=0;i<15;i++){
	            tmp1 = dfloat(i+1);
	            tmp2 = dfloat(16 - (i+1));
	            tmp3 = tmp1;
	            if (i > 8) {tmp3 = tmp2;}
	            tmp4 = (x[1]*tmp2 + x[2]*tmp3)*(x[1]*tmp2 + x[2]*tmp3);
	            fJac[i][1] = -one;
	            fJac[i][2] = tmp1*tmp2/tmp4;
	            fJac[i][3] = tmp1*tmp3/tmp4;
	        }
	    break;
	    case (9):
	        //KOWALIK AND OSBORNE FUNCTION.
	        for(i=0;i< 11;i++){ 
	            tmp1 = v[i]*(v[i] + x[1]);
	            tmp2 = v[i]*(v[i] + x[2]) + x[3];
	            fJac[i][1] = -tmp1/tmp2;
	            fJac[i][2] = -v[i]*x[0]/tmp2;
	            fJac[i][3] = fJac[i][1]*fJac[i][2];
	            fJac[i][4] = fJac[i][3]/v[i];
	        }
	    break;
	    case (10):
	        //MEYER FUNCTION.
	        for(i=0;i< 16;i++){
	            temp = five*dfloat(i+1) + c45 + x[2];
	            tmp1 = x[1]/temp;
	            tmp2 = Math.exp(tmp1);
	            fJac[i][1] = tmp2;
	            fJac[i][2] = x[0]*tmp2/temp;
	            fJac[i][3] = -tmp1*fJac[i][2];
	        }
	    break;
	    case (11):
	        //WATSON FUNCTION.
	        for(i=0;i<29;i++){
	            div = dfloat(i+1)/c29;
	            s2 = zero;
	            dx = one;
	            for(j=0;j<n;j++){
	                s2 = s2 + dx*x[j];
	                dx = div*dx;
	            }
	            temp = two*div*s2;
	            dx = one/div;
	            for(j=0;j<n;j++){
	                fJac[i][j] = dx*(dfloat(j) - temp);
	                dx = div*dx;
	            }
	        }
	        for(j=0;j<n;j++){
	            for( i = 29;i<31;i++){
	                fJac[i][j] = zero;
	            }
	        }
	        fJac[30][1] = one;
	        fJac[31][1] = -two*x[0];
	        fJac[31][2] = one;
		    break;
	    case (12):
	        //BOX 3-DIMENSIONAL FUNCTION.
	        for(i=0;i<m;i++){
	            temp = dfloat(i+1);
	            tmp1 = temp/ten;
	            fJac[i][1] = -tmp1*Math.exp(-tmp1*x[0]);
	            fJac[i][2] = tmp1*Math.exp(-tmp1*x[1]);
	            fJac[i][3] = Math.exp(-temp) - Math.exp(-tmp1);
	        }
	    break;
	    case (13):
	        //JENNRICH AND SAMPSON FUNCTION.
	        for(i=0;i<m;i++){
	            temp = dfloat(i+1);
	            fJac[i][1] = -temp*Math.exp(temp*x[0]);
	            fJac[i][2] = -temp*Math.exp(temp*x[1]);
	        }
	    break;
	    case (14):
	        //BROWN AND DENNIS FUNCTION.
	        for(i=0;i<m;i++){
	            temp = dfloat(i+1)/five;
	            ti = Math.sin(temp);
	            tmp1 = x[0] + temp*x[1] - Math.exp(temp);
	            tmp2 = x[2] + ti*x[3] - Math.cos(temp);
	            fJac[i][1] = two*tmp1;
	            fJac[i][2] = temp*fJac[i][1];
	            fJac[i][3] = two*tmp2;
	            fJac[i][4] = ti*fJac[i][3];
	        }
	    break;
	    case (15):
	        //CHEBYQUAD FUNCTION.
	        dx = one/dfloat(n);
	        for(j=0;j<n;j++){
	            tmp1 = one;
	            tmp2 = two*x[j] - one;
	            temp = two*tmp2;
	            tmp3 = zero;
	            tmp4 = two;
	            for(i=0;i<m;i++){
	                fJac[i][j] = dx*tmp4;
	                ti = four*tmp2 + temp*tmp4 - tmp3;
	                tmp3 = tmp4;
	                tmp4 = ti;
	                ti = temp*tmp2 - tmp1;
	                tmp1 = tmp2;
	                tmp2 = ti;
	            }
	        }
		    break;
	    case (16):
	        //BROWN ALMOST-LINEAR FUNCTION.
	        prod = one;
	        for(j=0;j<n;j++){
	            prod = x[j]*prod;
	            for(i=0;i<n;i++){
	                fJac[i][j] = one;
	            }
	            fJac[j][j] = two;
	        }
	        for(j=0;j<n;j++){
	            temp = x[j];
	            if (temp == zero) {
	                temp = one;
	                prod = one;
	                for( k = 0;k<n;k++){
	                    if (k != j) {prod = x[k]*prod;}
	                }
	            }
	            fJac[n][j] = prod/temp;
	        }
		    break;
	    case (17):
	        //OSBORNE 1 FUNCTION.
	        for(i=0;i< 33;i++){ 
	            temp = ten*dfloat(i);
	            tmp1 = Math.exp(-x[3]*temp);
	            tmp2 = Math.exp(-x[4]*temp);
	            fJac[i][1] = -one;
	            fJac[i][2] = -tmp1;
	            fJac[i][3] = -tmp2;
	            fJac[i][4] = temp*x[1]*tmp1;
	            fJac[i][5] = temp*x[2]*tmp2;
	        }
	    break;
	    case (18):
	        //OSBORNE 2 FUNCTION.
	        for(i=0;i< 65;i++){
	            temp = dfloat(i)/ten;
	            tmp1 = Math.exp(-x[4]*temp);
	            tmp2 = Math.exp(-x[5]*(temp - x[8])*x[8]);
	            tmp3 = Math.exp(-x[6]*(temp - x[9])*x[9]);
	            tmp4 = Math.exp(-x[7]*(temp - x[10])*x[10]);
	            fJac[i][1] = -tmp1;
	            fJac[i][2] = -tmp2;
	            fJac[i][3] = -tmp3;
	            fJac[i][4] = -tmp4;
	            fJac[i][5] = temp*x[0]*tmp1;
	            fJac[i][6] = x[1]*(temp - x[8])*x[8]*tmp2;
	            fJac[i][7] = x[2]*(temp - x[9])*x[9]*tmp3;
	            fJac[i][8] = x[3]*(temp - x[10])*x[10]*tmp4;
	            fJac[i][9] = -two*x[1]*x[5]*(temp - x[8])*tmp2;
	            fJac[i][10] = -two*x[2]*x[6]*(temp - x[9])*tmp3;
	            fJac[i][11] = -two*x[3]*x[7]*(temp - x[10])*tmp4;
	        }
	    break;
	    default:
	        //LINEAR FUNCTION - FULL RANK.
	        temp = two/dfloat(m);
	        for(j=0;j<n;j++){
	            for(i=0;i<m;i++){
	                fJac[i][j] = -temp;
	            }
	            fJac[j][j] = fJac[j][j] + one;
	        }
		    break;
	    }
}
//	end subroutine ssqjac
//	!*****************************************************************************************
//
//	!*****************************************************************************************
//	!>
//	!  This subroutine specifies the standard starting points for the
//	!  functions defined by subroutine ssqfcn. the subroutine returns
//	!  in x a multiple (factor) of the standard starting point. for
//	!  the 11th function the standard starting point is zero, so in
//	!  this case, if factor is not unity, then the subroutine returns
//	!  the vector  x[j] = factor, j=1,...,n.
public static void initpt(int n, double[] x, int nProb, double factor) {
//	subroutine initpt(n, x, nProb, factor)
//	    implicit none
//
//	    integer,intent(in) :: n !! a positive integer input variable.
//	    integer,intent(in) :: nProb !! a positive integer input variable which defines the
//	                                !! number of the problem. nProb must not exceed 18.
//	    real(wp),intent(in) :: factor !! an input variable which specifies the multiple of
//	                                  !! the standard starting point. if factor is unity, no
//	                                  !! multiplication is performed.
//	    real(wp),intent(out) :: x(n) !! an output array of length n which contains the standard
//	                                 !! starting point for problem nProb multiplied by factor.
//
	    double zero = 0.0;
	    double half = 0.5;
	    double one = 1.0;
	    double two = 2.0;
	    double three = 3.0;
	    double five = 5.0;
	    double seven = 7.0;
	    double ten = 10.0;
	    double twenty = 20.0;
	    double twntf = 25.0;

	    double c1 = 1.2;
	    double c2 = 2.5e-1;
	    double c3 = 3.9e-1;
	    double c4 = 4.15e-1;
	    double c5 = 2.0e-2;
	    double c6 = 4.0e3;
	    double c7 = 2.5e2;
	    double c8 = 3.0e-1;
	    double c9 = 4.0e-1;
	    double c10 = 1.5;
	    double c11 = 1.0e-2;
	    double c12 = 1.3;
	    double c13 = 6.5e-1;
	    double c14 = 7.0e-1;
	    double c15 = 6.0e-1;
	    double c16 = 4.5;
	    double c17 = 5.5;

	    int ivar, j;
	    double h;

//	    x(1:n) = zero
	    x = new double[n];
	    for(int i = 0;i<n;i++) {
	    	x[i] = 0;
	    }

	    //SELECTION OF INITIAL POINT.

//	    select case (nProb)
		switch(nProb){
	    case (4):
	        //ROSENBROCK FUNCTION.
	        x[0] = -c1;
	        x[1] = one;
	        break;
	    case (5):
	        //HELICAL VALLEY FUNCTION.
	        x[0] = -one;
	        x[1] = zero;
	        x[2] = zero;
	        break;
	    case (6):
	        //POWELL SINGULAR FUNCTION.
	        x[0] = three;
	        x[1] = -one;
	        x[2] = zero;
	        x[3] = one;
	        break;
	    case (7):
	        //FREUDENSTEIN AND ROTH FUNCTION.
	        x[0] = half;
	        x[1] = -two;
	        break;
	    case (8):
	        //BARD FUNCTION.
	        x[0] = one;
	        x[1] = one;
	        x[2] = one;
	        break;
	    case (9):
	        //KOWALIK AND OSBORNE FUNCTION.
	        x[0] = c2;
	        x[1] = c3;
	        x[2] = c4;
	        x[3] = c3;
	        break;
	    case (10):
	        //MEYER FUNCTION.
	        x[0] = c5;
	        x[1] = c6;
	        x[2] = c7;
	        break;
	    case (11):
	        //WATSON FUNCTION.
	        for(j=0;j<n;j++){
	            x[j] = zero;
	        }
        break;
	    case (12):
	        //BOX 3-DIMENSIONAL FUNCTION.
	        x[0] = zero;
	        x[1] = ten;
	        x[2] = twenty;
	        break;
	    case (13):
	        //JENNRICH AND SAMPSON FUNCTION.
	        x[0] = c8;
	        x[1] = c9;
	        break;
	    case (14):
	        //BROWN AND DENNIS FUNCTION.
	        x[0] = twntf;
	        x[1] = five;
	        x[2] = -five;
	        x[3] = -one;
	        break;
	    case (15):
	        //CHEBYQUAD FUNCTION.
	        h = one/dfloat(n + 1);
	        for(j=0;j<n;j++){
	            x[j] = dfloat(j+1)*h;
	        }
	        break;
	    case (16):
	        //BROWN ALMOST-LINEAR FUNCTION.
	        for(j=0;j<n;j++){
	            x[j] = half;
	        }
        break;
	    case (17):
	        //OSBORNE 1 FUNCTION.
	        x[0] = half;
	        x[1] = c10;
	        x[2] = -one;
	        x[3] = c11;
	        x[4] = c5;
	        break;
	    case (18):
	        //OSBORNE 2 FUNCTION.
	        x[0] = c12;
	        x[1] = c13;
	        x[2] = c13;
	        x[3] = c14;
	        x[4] = c15;
	        x[5] = three;
	        x[6] = five;
	        x[7] = seven;
	        x[8] = two;
	        x[9] = c16;
	        x[10] = c17;
	        break;
	    default:
	        //LINEAR FUNCTION - FULL RANK OR RANK 1.
	        for(j=0;j<n;j++){
	            x[j] = one;
	        }
	        break;
	    }

	    //COMPUTE MULTIPLE OF INITIAL POINT.

	    if (factor != one) {
	        if (nProb == 11) {
	            for(j=0;j<n;j++){
	                x[j] = factor;
	            }
	        }else{
	            for(j=0;j<n;j++){
	                x[j] = factor*x[j];
	            }
	        }
	    }
}
//	end subroutine initpt
//	!*****************************************************************************************
//
//	!*****************************************************************************************
//	!>
//	!  This subroutine defines the functions of eighteen nonlinear
//	!  least squares problems. the allowable values of (m,n) for
//	!  functions 1,2 and 3 are variable but with m >= n.
//	!  for functions 4,5,6,7,8,9 and 10 the values of (m,n) are
//	!  (2,2),(3,3),(4,4),(2,2),(15,3),(11,4) and (16,3), respectively.
//	!  function 11 (watson) has m = 31 with n usually 6 or 9.
//	!  however, any n, n = 2,...,31, is permitted.
//	!  functions 12,13 and 14 have n = 3,2 and 4, respectively, but
//	!  allow any m >= n, with the usual choices being 10,10 and 20.
//	!  function 15 (chebyquad) allows m and n variable with m >= n.
//	!  function 16 (brown) allows n variable with m = n.
//	!  for functions 17 and 18, the values of (m,n) are
//	!  (33,5) and (65,11), respectively.
public static void ssqfcn(int m, int n, double[] x, double[] fVec, int nProb) {
//	subroutine ssqfcn(m, n, x, fVec, nProb)
//	    implicit none
//
//	    integer,intent(in) :: m !! a positive integer input variable.
//	    integer,intent(in) :: n !! a positive integer input variable. n must not exceed m.
//	    integer,intent(in) :: nProb !! a positive integer input variable which defines the
//	                                !! number of the problem. nProb must not exceed 18.
//	    real(wp),intent(in) :: x(n) !! an input array of length n.
//	    real(wp),intent(out) :: fVec[m] !! an output array of length m which contains the nProb
//	                                    !! function evaluated at x.
//
	    double zero = 0.0;
	    double zp25 = 2.5e-1;
	    double zp5 = 5.0e-1;
	    double one = 1.0;
	    double two = 2.0;
	    double five = 5.0;
	    double eight = 8.0;
	    double ten = 10.0;
	    double c13 = 13.0;
	    double c14 = 14.0;
	    double c29 = 29.0;
	    double c45 = 45.0;

	    double[] v =  {4.0, 2.0, 1.0, 5.0e-1,
	                                    2.5e-1, 1.67e-1, 1.25e-1, 1.0e-1, 8.33e-2, 7.14e-2,
	                                    6.25e-2};
	    double[] y1 = {1.4e-1, 1.8e-1, 2.2e-1, 2.5e-1, 2.9e-1, 3.2e-1,
	                                    3.5e-1, 3.9e-1, 3.7e-1, 5.8e-1, 7.3e-1, 9.6e-1,
	                                    1.34, 2.1, 4.39};
	    double[] y2 = {1.957e-1, 1.947e-1,
	                                    1.735e-1, 1.6e-1, 8.44e-2, 6.27e-2, 4.56e-2, 3.42e-2,
	                                    3.23e-2, 2.35e-2, 2.46e-2};
	    double[] y3 = {3.478e4, 2.861e4, 2.365e4, 1.963e4,
	                                    1.637e4, 1.372e4, 1.154e4, 9.744e3, 8.261e3, 7.03e3,
	                                    6.005e3, 5.147e3, 4.427e3, 3.82e3, 3.307e3, 2.872e3};
	    double[] y4 = {8.44e-1,
	                                    9.08e-1, 9.32e-1, 9.36e-1, 9.25e-1, 9.08e-1, 8.81e-1,
	                                    8.5e-1, 8.18e-1, 7.84e-1, 7.51e-1, 7.18e-1, 6.85e-1,
	                                    6.58e-1, 6.28e-1, 6.03e-1, 5.8e-1, 5.58e-1, 5.38e-1,
	                                    5.22e-1, 5.06e-1, 4.9e-1, 4.78e-1, 4.67e-1, 4.57e-1,
	                                    4.48e-1, 4.38e-1, 4.31e-1, 4.24e-1, 4.2e-1, 4.14e-1,
	                                    4.11e-1, 4.06e-1};
	    double[] y5 = {1.366,
	                                    1.191, 1.112, 1.013, 9.91e-1, 8.85e-1, 8.31e-1,
	                                    8.47e-1, 7.86e-1, 7.25e-1, 7.46e-1, 6.79e-1, 6.08e-1,
	                                    6.55e-1, 6.16e-1, 6.06e-1, 6.02e-1, 6.26e-1, 6.51e-1,
	                                    7.24e-1, 6.49e-1, 6.49e-1, 6.94e-1, 6.44e-1, 6.24e-1,
	                                    6.61e-1, 6.12e-1, 5.58e-1, 5.33e-1, 4.95e-1, 5.0e-1,
	                                    4.23e-1, 3.95e-1, 3.75e-1, 3.72e-1, 3.91e-1, 3.96e-1,
	                                    4.05e-1, 4.28e-1, 4.29e-1, 5.23e-1, 5.62e-1, 6.07e-1,
	                                    6.53e-1, 6.72e-1, 7.08e-1, 6.33e-1, 6.68e-1, 6.45e-1,
	                                    6.32e-1, 5.91e-1, 5.59e-1, 5.97e-1, 6.25e-1, 7.39e-1,
	                                    7.1e-1, 7.29e-1, 7.2e-1, 6.36e-1, 5.81e-1, 4.28e-1,
	                                    2.92e-1, 1.62e-1, 9.8e-2, 5.4e-2};

	    int i, iev, ivar, j, nm1;
	    double div, dx, prod, sum, s1, s2, temp, ti, tmp1, tmp2, tmp3, tmp4, tpi;

//	    fVec(1:m) = zero;
	    fVec = new double[m];
for(i = 0;i<m;i++) {
	    	fVec[i] = 0;
	    }

	    //FUNCTION ROUTINE SELECTOR.

	    //select case (nProb)
		switch(nProb){		
	    case (2):
	        //LINEAR FUNCTION - RANK 1.
	        sum = zero;
	        for(j=0;j<n;j++){
	            sum = sum + dfloat(j+1)*x[j];
	        }
	        for(i=0;i<m;i++){
	            fVec[i] = dfloat(i+1)*sum - one;
	        }
	        break;
	    case (3):
	        //LINEAR FUNCTION - RANK 1 WITH ZERO COLUMNS AND ROWS.
	        sum = zero;
	        nm1 = n - 1;
	        if (nm1 >= 2) {
	            for(j = 1;j< nm1;j++){
	                sum = sum + dfloat(j+1)*x[j];
	            }
	        }
	        for(i=0;i<m;i++){
	            fVec[i] = dfloat(i)*sum - one;
	        }
	        fVec[m] = -one;
	        break;
	    case (4):
	        //ROSENBROCK FUNCTION.
	        fVec[0] = ten*(x[1] - x[0]*x[0]);
	        fVec[1] = one - x[0];
	    case (5):
	        //HELICAL VALLEY FUNCTION.
	        tpi = eight*Math.atan(one);
	        tmp1 = Math.signum(x[1]) * Math.abs(zp25); //, x[1]);
	        if (x[0] > zero) {tmp1 = Math.atan(x[1]/x[0])/tpi;}
	        if (x[0] < zero) {tmp1 = Math.atan(x[1]/x[0])/tpi + zp5;}
	        tmp2 = Math.sqrt(x[0]*x[0] + x[1]*x[1]);
	        fVec[1] = ten*(x[2] - ten*tmp1);
	        fVec[2] = ten*(tmp2 - one);
	        fVec[3] = x[2];
	        break;
	    case (6):
	        //POWELL SINGULAR FUNCTION.
	        fVec[1] = x[0] + ten*x[1];
	        fVec[2] = Math.sqrt(five)*(x[2] - x[3]);
	        fVec[3] = (x[1] - two*x[2])*(x[1] - two*x[2]);
	        fVec[4] = Math.sqrt(ten)*(x[0] - x[3])*(x[0] - x[3]);
	        break;
	    case (7):
	        //FREUDENSTEIN AND ROTH FUNCTION.
	        fVec[1] = -c13 + x[0] + ((five - x[1])*x[1] - two)*x[1];
	        fVec[2] = -c29 + x[0] + ((one + x[1])*x[1] - c14)*x[1];
	        break;
	    case (8):
	        //BARD FUNCTION.
	        for(i=0;i< 15;i++){ 
	            tmp1 = dfloat(i+1);
	            tmp2 = dfloat(16 - (i + 1));
	            tmp3 = tmp1;
	            if (i > 8) {tmp3 = tmp2;}
	            fVec[i] = y1[i] - (x[0] + tmp1/(x[1]*tmp2 + x[2]*tmp3));
	        }
	    break;
	    case (9):
	        //KOWALIK AND OSBORNE FUNCTION.
	        for(i=0;i<11;i++){
	            tmp1 = v[i]*(v[i] + x[1]);
	            tmp2 = v[i]*(v[i] + x[2]) + x[3];
	            fVec[i] = y2[i] - x[0]*tmp1/tmp2;
	        }
	    break;
	    case (10):
	        //MEYER FUNCTION.
	        for(i=0;i<16;i++){
	            temp = five*dfloat(i+1) + c45 + x[2];
	            tmp1 = x[1]/temp;
	            tmp2 = Math.exp(tmp1);
	            fVec[i] = x[0]*tmp2 - y3[i];
	        }
	    break;
	    case (11):
	        //WATSON FUNCTION.
	        for(i=0;i<29;i++){ 
	            div = dfloat(i+1)/c29;
	            s1 = zero;
	            dx = one;
	            for(j = 1;j< n;j++){
	                s1 = s1 + dfloat(j)*dx*x[j];
	                dx = div*dx;
	            }
	            s2 = zero;
	            dx = one;
	            for(j=0;j<n;j++){
	                s2 = s2 + dx*x[j];
	                dx = div*dx;
	            }
	            fVec[i] = s1 - s2*s2 - one;
	        }
	        fVec[30] = x[0];
	        fVec[31] = x[1] - x[0]*x[0] - one;
	        break;
	    case (12):
	        //BOX 3-DIMENSIONAL FUNCTION.
	        for(i=0;i<m;i++){
	            temp = dfloat(i+1);
	            tmp1 = temp/ten;
	            fVec[i] = Math.exp(-tmp1*x[0]) - Math.exp(-tmp1*x[1]) + (Math.exp(-temp) - Math.exp(-tmp1))*x[2];
	        }
	    break;
	    case (13):
	        //JENNRICH AND SAMPSON FUNCTION.
	        for(i=0;i<m;i++){
	            temp = dfloat(i+1);
	            fVec[i] = two + two*temp - Math.exp(temp*x[0]) - Math.exp(temp*x[1]);
	        }
	    break;
	    case (14):
	        //BROWN AND DENNIS FUNCTION.
	        for(i=0;i<m;i++){
	            temp = dfloat(i+1)/five;
	            tmp1 = x[0] + temp*x[1] - Math.exp(temp);
	            tmp2 = x[2] + Math.sin(temp)*x[3] - Math.cos(temp);
	            fVec[i] = tmp1*tmp1 + tmp2*tmp2;
	        }
	    break;
	    case (15):
	        //CHEBYQUAD FUNCTION.
	        for(i=0;i<m;i++){
	            fVec[i] = zero;
	        }
	        for(j=0;j<n;j++){
	            tmp1 = one;
	            tmp2 = two*x[j] - one;
	            temp = two*tmp2;
	            for(i=0;i<m;i++){
	                fVec[i] = fVec[i] + tmp2;
	                ti = temp*tmp2 - tmp1;
	                tmp1 = tmp2;
	                tmp2 = ti;
	            }
	        }
	        dx = one/dfloat(n);
	        iev = -1;
	        for(i=0;i<m;i++){
	            fVec[i] = dx*fVec[i];
	            if (iev > 0) fVec[i] = fVec[i] + one/(dfloat(i+1)*dfloat(i+1) - one);
	            iev = -iev;
	        }
	        break;
	    case (16):
	        //BROWN ALMOST-LINEAR FUNCTION.
	        sum = -dfloat(n + 1);
	        prod = one;
	        for(j=0;j<n;j++){
	            sum = sum + x[j];
	            prod = x[j]*prod;
	        }
	        for(i=0;i<n;i++){
	            fVec[i] = x[i] + sum;
	        }
	        fVec[n] = prod - one;
	        break;
	    case (17):
	        //OSBORNE 1 FUNCTION.
	        for(i=0;i< 33;i++){
	            temp = ten*dfloat(i);
	            tmp1 = Math.exp(-x[3]*temp);
	            tmp2 = Math.exp(-x[4]*temp);
	            fVec[i] = y4[i] - (x[0] + x[1]*tmp1 + x[2]*tmp2);
	        }
	    break;
	    case (18):
	        //OSBORNE 2 FUNCTION.
	        for(i=0;i<65 ;i++){ 
	            temp = dfloat(i)/ten;
	            tmp1 = Math.exp(-x[4]*temp);
	            tmp2 = Math.exp(-x[5]*(temp - x[8])*(temp - x[8]));
	            tmp3 = Math.exp(-x[6]*(temp - x[9])*(temp - x[9]));
	            tmp4 = Math.exp(-x[7]*(temp - x[10])*(temp - x[10]));
	            fVec[i] = y5[i] - (x[0]*tmp1 + x[1]*tmp2 + x[2]*tmp3 + x[3]*tmp4);
	        }
	    break;
	    default:
	        //LINEAR FUNCTION - FULL RANK.
	        sum = zero;
	        for(j=0;j<n;j++){
	            sum = sum + x[j];
	        }
	        temp = two*sum/dfloat(m) + one;
	        for(i=0;i<m;i++){
	            fVec[i] = -temp;
	            if (i <= n) {fVec[i] = fVec[i] + x[i];}
	        }
	        break;
	    } //end select
}
}
