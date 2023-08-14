package minpackTestPacks;
//
public class lmdifTest {
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
//	program test_lmdif
//
//	    use minpack_module, only: wp, dpmpar, enorm, lmdif1
//	    use iso_fortran_env, only: nwrite => output_unit
//
//	    implicit none
//
//	    ! originally from file22
//	    int ncases = 28
//	    int[] nProbs  = new int[ncases] [1,1,2,2,3,3,4,5,6,7,8,9,10,11,11,11,12,13,14,15,15,15,15,16,16,16,17,18]
//	    int[] ns = new int[ncases]  [5,5,5,5,5,5,2,3,4,2,3,4,3,6,9,12,3,2,4,1,8,9,10,10,30,40,5,11]
//	    int[] ms = new int[ncases]  [10,50,10,50,10,50,2,3,4,2,15,11,16,31,31,31,10,10,20,8,8,9,10,10,30,40,33,65]
//	    int[] nTriess = new int[ncases]  [1,1,1,1,1,1,3,3,3,3,3,3,2,3,3,3,1,1,3,3,1,1,1,3,1,1,1,1]
//
//	    int[] info_original = new int[53]  [1,1,1,1,1,1,2,2,2,2,2,2,5,5,5,1,1,1,1,1,1,1,1,5,1,4,&
//	                                                        1,1,1,1,2,2,2,2,2,2,1,5,1,5,1,1,1,1,2,1,1,2,1,2,2,1,1]
//	                                                        !! original `info` from the original minpack
//
//	    int i, ic, info, k, m, n, NFEv, NJEv, nProb, nTries, iCase, lwa
//	    double factor, fnorm1, fnorm2
//	    int ma(53), na(53), nf(53), nj(53), np(53), nx(53)
//	    double fnm(53)
//	    int[] iwa
//	    double[] fVec, wa, x
//
//	    double one = 1.0;
//	    double ten = 10.0;
//	    double tol = Math.sqrt(dpmpar(1))
//	    real(wp), parameter :: solution_reltol = 1.0e-3; !! reltol for matching previously generated solutions
//
//	    ic = 0
//	    do iCase = 1, ncases+1
//
//	        if (iCase == ncases+1) {
//	            write (nwrite, '(A,I3,A/)') '1SUMMARY OF ', ic, ' CALLS TO LMDIF1'
//	            write (nwrite, '(A/)')      ' nProb   N    M   NFEV  NJEV  INFO  FINAL L2 NORM'
//	            for(i=0;i< ;i++){ic
//	                write (nwrite, '(3I5,3I6,1X,D15.7)') np(i), na(i), ma(i), nf(i), nj(i), nx[i], fnm(i)
//	            }
//	            stop
//	        } else {
//	            nProb = nProbs(iCase)
//	            n = ns(iCase)
//	            m = ms(iCase)
//	            lwa = m*n+5*n+m
//	            nTries = nTriess(iCase)
//	            if (allocated(iwa))  deallocate(iwa);  allocate(iwa(n))
//	            if (allocated(fVec)) deallocate(fVec); allocate(fVec(m))
//	            if (allocated(wa))   deallocate(wa);   allocate(wa(lwa))
//	            if (allocated(x))    deallocate(x);    allocate(x(n))
//	            factor = one
//	            do k = 1, nTries
//	                ic = ic + 1
//	                call initpt(n, x, nProb, factor)
//	                call ssqfcn(m, n, x, fVec, nProb)
//	                fnorm1 = enorm(m, fVec)
//	                write (nwrite, '(////5X,A,I5,5X,A,2I5,5X//)') ' PROBLEM', nProb, ' DIMENSIONS', n, m
//	                NFEv = 0
//	                NJEv = 0
//	                call lmdif1(fcn, m, n, x, fVec, tol, info, iwa, wa, lwa)
//	                call ssqfcn(m, n, x, fVec, nProb)
//	                fnorm2 = enorm(m, fVec)
//	                np[ic] = nProb
//	                na[ic] = n
//	                ma[ic] = m
//	                nf[ic] = NFEv
//	                NJEv = NJEv/n
//	                nj[ic] = NJEv
//	                nx[ic] = info
//	                fnm[ic] = fnorm2
//	                write(nwrite,'(5X,A,D15.7//5X,A,D15.7//5X,A,I10//5X,A,I10//5X,A,18X,I10//5X,A//*(5X,5D15.7/))') &
//	                             ' INITIAL L2 NORM OF THE RESIDUALS', fnorm1, &
//	                             ' FINAL L2 NORM OF THE RESIDUALS  ', fnorm2, &
//	                             ' NUMBER OF FUNCTION EVALUATIONS  ', NFEv,   &
//	                             ' NUMBER OF JACOBIAN EVALUATIONS  ', NJEv,   &
//	                             ' EXIT PARAMETER', info,                     &
//	                             ' FINAL APPROXIMATE SOLUTION', x(1:n)
//	                factor = ten*factor
//	                call compare_solutions(ic, x, solution_reltol, tol)
//	            }
//	        }
//	    }
//
//	    contains
//	!*****************************************************************************************
//
//	!*****************************************************************************************
//	!>
//	!  Compare with previously generated solutions.
//
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
//
//	    if (info_original[ic]<5) {    ! ignore any where the original minpack failed
//	        diff = solution[ic] - x
//	        absdiff = abs(diff)
//	        if (any(absdiff>abstol)) { ! first do an absolute diff
//	            ! also do a rel diff if the abs diff fails (also protect for divide by zero)
//	            reldiff = absdiff
//	            where (solution[ic] != 0.0;) reldiff = absdiff / abs(solution[ic])
//	            if (any(reldiff > reltol)) {
//	                write(nwrite,'(A)') 'Failed case'
//	                write(nwrite, '(//5x, a//(5x, 5d15.7))') 'Math.expected x: ', solution[ic]
//	                write(nwrite, '(/5x, a//(5x, 5d15.7))')  'Computed x: ', x
//	                write(nwrite, '(/5x, a//(5x, 5d15.7))')  'absdiff: ', absdiff
//	                write(nwrite, '(/5x, a//(5x, 5d15.7))')  'reldiff: ', reldiff
//	                error stop ! test failed
//	            }
//	        }
//	    }
//
//	    end subroutine compare_solutions
//	!*****************************************************************************************
//
//	!*****************************************************************************************
//	!>
//	!  The calling sequence of fcn should be identical to the
//	!  calling sequence of the function subroutine in the nonlinear
//	!  least-squares solver. fcn should only call the testing
//	!  function subroutine ssqfcn with the appropriate value of
//	!  problem number (nProb).
//
//	    subroutine fcn(m, n, x, fVec, Iflag)
//	        implicit none
//
//	        integer,intent(in) :: m
//	        integer,intent(in) :: n
//	        real(wp),intent(in) :: x(n)
//	        real(wp),intent(out) :: fVec(m)
//	        integer,intent(inout) :: Iflag
//
//	        call ssqfcn(m, n, x, fVec, nProb)
//
//	        switch (Iflag)
//	        case(1)
//	            NFEv = NFEv + 1
//	        case(2)
//	            NJEv = NJEv + 1
//	       default:
//	            error stop 'invalid iflag value'
//	        }
//
//	    end subroutine fcn
//	!*****************************************************************************************
//
//	!*****************************************************************************************
//	!>
//	!  Replaced statement function in original code.
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
//	!  Get Math.expected `x` vectors for each case.
//
//	    pure function solution(nProb) result(x)
//
	        public static double[] solution(int nProb) {

	        double[] x;

	        switch (nProb) {
	        case (1): 
	        	x = new double[] {-0.1000000029802320E+01,-0.1000000029802320E+01,-0.1000000029802320E+01,
	                         -0.1000000029802320E+01,-0.9999999552965200E+00};
	        return x;
	        case (  2): x = new double[] {-0.9999998927116482E+00,-0.9999998927116488E+00,-0.9999998927116555E+00,
	                         -0.9999998927116560E+00,-0.9999998927116568E+00};
	        return x;
	        case (  3): x = new double[] {-0.1677968180239693E+03,-0.8339840901198468E+02, 0.2211100430795781E+03,
	                         -0.4119920450599233E+02,-0.3275936360479385E+02};
	        return x;
	        case (  4): x = new double[] {-0.2029999900022674E+02,-0.9649999500113370E+01,-0.1652451975264496E+03,
	                         -0.4324999750056676E+01, 0.1105330585100652E+03};
	        return x;
	        case (  5): x = new double[] { 0.1000000000000000E+01,-0.2103615324224772E+03, 0.3212042081132130E+02,
	                          0.8113456824980642E+02, 0.1000000000000000E+01};
	        return x;
	        case (  6): x = new double[] { 0.1000000000000000E+01, 0.1103865983923618E+03,-0.1465594395269163E+03,
	                          0.5473401240776926E+02, 0.1000000000000000E+01};
	        return x;
	        case (  7): x = new double[] { 0.1000000000000000E+01, 0.1000000000000000E+01};
	        return x;
	        case (  8): x = new double[] { 0.1000000000000000E+01, 0.1000000000000000E+01};
	        return x;
	        case (  9): x = new double[] { 0.1000000000000000E+01, 0.1000000000000000E+01};
	        return x;
	        case ( 10): x = new double[] { 0.1000000000000000E+01, 0.1425619853622940E-16, 0.0000000000000000E+00};
	        return x;
	        case ( 11): x = new double[] { 0.1000000000000000E+01, 0.7098020938940133E-18, 0.0000000000000000E+00};
	        return x;
	        case ( 12): x = new double[] { 0.1000000000000000E+01,-0.2725107730298909E-22, 0.0000000000000000E+00};
	        return x;
	        case ( 13): x = new double[] { 0.5263054387860989E-10,-0.5263054387860989E-11, 0.2465981630582231E-10,
	                          0.2465981630582231E-10};
	        return x;
	        case ( 14): x = new double[] { 0.4417018305110184E-10,-0.4417018305113451E-11, 0.2043852798973299E-10,
	                          0.2043852798973299E-10};
	        case ( 15): x = new double[] { 0.7114477974627450E-10,-0.7114477974627451E-11, 0.3363831937585708E-10,
	                          0.3363831937585708E-10};
	        return x;
	        case ( 16): x = new double[] { 0.1141248445844032E+02,-0.8968279141291072E+00};
	        return x;
	        case ( 17): x = new double[] { 0.1141300460911224E+02,-0.8967960407402507E+00};
	        return x;
	        case ( 18): x = new double[] { 0.1141278171916741E+02,-0.8968051105401963E+00};
	        return x;
	        case ( 19): x = new double[] { 0.8241057720241220E-01, 0.1133036677062726E+01, 0.2343694616119322E+01};
	        return x;
	        case ( 20): x = new double[] { 0.8406666710204237E+00,-0.2370847506164708E+09,-0.2404886744537081E+09};
	        return x;
	        case ( 21): x = new double[] { 0.8406666930558223E+00,-0.6605876218531588E+08,-0.6952091093079293E+08};
	        return x;
	        case ( 22): x = new double[] { 0.1928078051588323E+00, 0.1912627645045558E+00, 0.1230528191311660E+00,
	                          0.1360532725211912E+00};
	        return x;
	        case ( 23): x = new double[] { 0.1829224893098868E+06,-0.1407587255998039E+02,-0.8278551064687881E+07,
	                         -0.5164170812975379E+07};
	        return x;
	        case ( 24): x = new double[] { 0.3731181508014032E-02, 0.6367392058467618E+03, 0.7288615223592633E+01,
	                          0.4922699528439693E+01};
	        return x;
	        case ( 25): x = new double[] { 0.5609634098324512E-02, 0.6181346698953550E+04, 0.3452236464967932E+03};
	        return x;
	        case ( 26): x = new double[] { 0.1382351771873629E+01,-0.3663634184932995E+04,-0.2257365229134962E+02};
	        return x;
	        case ( 27): x = new double[] {-0.1389015072214595E-23, 0.1013638416971579E+01,-0.2443033770679000E+00,
	                          0.1373770737764654E+01,-0.1685639342611184E+01, 0.1098096981457304E+01};
	        return x;
	        case ( 28): x = new double[] {-0.1572518809742218E-01, 0.1012434868019557E+01,-0.2329916478346281E+00,
	                          0.1260429642080260E+01,-0.1513728144449624E+01, 0.9929958881927142E+00};
	        return x;
	        case ( 29): x = new double[] {-0.1572475363171166E-01, 0.1012434901340758E+01,-0.2329918373092776E+00,
	                          0.1260432420000366E+01,-0.1513732491253852E+01, 0.9929986220470630E+00};
	        return x;
	        case ( 30): x = new double[] { 0.1403007020665606E-21, 0.9997896480662041E+00, 0.1478229170675247E-01,
	                          0.1463069883566064E+00, 0.1001011648885689E+01,-0.2618209254498177E+01,
	                          0.4105099755986520E+01,-0.3144132399727149E+01, 0.1052791707674771E+01};
	        return x;
	        case ( 31): x = new double[] {-0.1502306463900362E-04, 0.9997897169117623E+00, 0.1476367815080281E-01,
	                          0.1463488075636044E+00, 0.1000791909805600E+01,-0.2617666211006538E+01, 
	                          0.4104328767857482E+01,-0.3143569813449599E+01, 0.1052617084764119E+01};
	        return x;
	        case ( 32): x = new double[] {-0.1532378591782743E-04, 0.9997896601115883E+00, 0.1476311309182103E-01,
	                          0.1463524931728402E+00, 0.1000777881470954E+01,-0.2617638010535779E+01,
	                          0.4104295790382706E+01,-0.3143549498653720E+01, 0.1052611820110282E+01};
	        return x;
	        case ( 33): x = new double[] {-0.1847463656530098E-22, 0.1000001659949695E+01,-0.5664287022818456E-03,
	                          0.3478754629617173E+00,-0.1572556993994226E+00, 0.1055542916650735E+01,
	                         -0.3255800432321387E+01, 0.7305174796577345E+01,-0.1029263063796092E+02,
	                          0.9089960730505242E+01,-0.4548149719224188E+01, 0.1013254861288710E+01};
	        return x;
	        case ( 34): x = new double[] { 0.9744479990467999E-07, 0.1000001546850994E+01,-0.5599459691289794E-03,
	                          0.3477832466116677E+00,-0.1565911336856186E+00, 0.1052711654130422E+01,
	                         -0.3248141797201425E+01, 0.7291658238611263E+01,-0.1027711647271800E+02,
	                          0.9078801606718631E+01,-0.4543583141978685E+01, 0.1012443953789911E+01};
	        return x;
	        case ( 35): x = new double[] {-0.2091128558991212E-07, 0.1000001538991780E+01,-0.5571965521936125E-03,
	                          0.3477159330232413E+00,-0.1559410543073800E+00, 0.1049320617322875E+01,
	                         -0.3237500657632093E+01, 0.7270621188722325E+01,-0.1025070733588128E+02,
	                          0.9058370428814367E+01,-0.4534697862977898E+01, 0.1010781878687835E+01};
	        return x;
	        case ( 36): x = new double[] { 0.1000000000000000E+01, 0.9999999999999996E+01, 0.9999999999999999E+00};
	        return x;
	        case ( 37): x = new double[] { 0.2578199268103105E+00, 0.2578299761926529E+00};
	        return x;
	        case ( 38): x = new double[] {-0.1157261477453999E+02, 0.1319584188052311E+02,-0.4077307965700902E+00,
	                          0.2333612117062785E+00};
	        return x;
	        case ( 39): x = new double[] {-0.1159528811684483E+02, 0.1320395000564385E+02,-0.4034298473767633E+00,
	                          0.2367736424151024E+00};
	        return x;
	        case ( 40): x = new double[] {-0.1158599992682185E+02, 0.1320046721521291E+02,-0.4030467979128602E+00,
	                          0.2371622176669948E+00};
	        return x;
	        case ( 41): x = new double[] { 0.5000000053094538E+00};
	        return x;
	        case ( 42): x = new double[] { 0.9817314947348010E+00};
	        return x;
	        case ( 43): x = new double[] { 0.9817314875630089E+00};
	        return x;
	        case ( 44): x = new double[] { 0.4315366607966317E-01, 0.1930916404198016E+00, 0.2663285955379775E+00,
	                          0.4999993382086716E+00, 0.5000006675747786E+00, 0.7336714089175319E+00,
	                          0.8069083676132508E+00, 0.9568463407766595E+00};
	        return x;
	        case ( 45): x = new double[] { 0.4420534613578277E-01, 0.1994906723098810E+00, 0.2356191084710600E+00,
	                          0.4160469078925980E+00, 0.5000000000000001E+00, 0.5839530921074019E+00,
	                          0.7643808915289401E+00, 0.8005093276901190E+00, 0.9557946538642172E+00};
	        return x;
	        case ( 46): x = new double[] { 0.5962027126608570E-01, 0.1667087889853926E+00, 0.2391710264712834E+00,
	                          0.3988852939227040E+00, 0.3988836749314369E+00, 0.6011163388953045E+00,
	                          0.6011147164581175E+00, 0.7608289915836429E+00, 0.8332912255167690E+00,
	                          0.9403797400270459E+00};
	        return x;
	        case ( 47): x = new double[] {-0.5469127990081513E-01,-0.5469127990081425E-01,-0.5469127990081792E-01,
	                         -0.5469127990081825E-01,-0.5469127990081792E-01,-0.5469127990081825E-01,
	                         -0.5469127990081792E-01,-0.5469127990081792E-01,-0.5469127990081391E-01,
	                          0.1154691279900732E+02};
	        return x;
	        case ( 48): x = new double[] { 0.9794303033498610E+00, 0.9794303033498610E+00, 0.9794303033498610E+00,
	                          0.9794303033498610E+00, 0.9794303033498610E+00, 0.9794303033498610E+00,
	                          0.9794303033498610E+00, 0.9794303033498610E+00, 0.9794303033498610E+00,
	                          0.1205696966501391E+01};
	        return x;
	        case ( 49): x = new double[] {-0.5465647913616722E-02,-0.5465647913616175E-02,-0.5465647913616722E-02,
	                         -0.5465647913616722E-02,-0.5465647913616330E-02,-0.5465647913616175E-02,
	                         -0.5465647913616450E-02,-0.5465647913616450E-02,-0.5465647913616330E-02,
	                          0.1105465647913852E+02};
	        return x;
	        case ( 50): x = new double[] { 0.1000000000000012E+01, 0.1000000000000012E+01, 0.1000000000000012E+01,
	                          0.1000000000000012E+01, 0.1000000000000012E+01, 0.1000000000000012E+01,
	                          0.1000000000000012E+01, 0.1000000000000012E+01, 0.1000000000000012E+01,
	                          0.1000000000000012E+01, 0.1000000000000012E+01, 0.1000000000000012E+01,
	                          0.1000000000000012E+01, 0.1000000000000012E+01, 0.1000000000000012E+01,
	                          0.1000000000000012E+01, 0.1000000000000012E+01, 0.1000000000000012E+01,
	                          0.1000000000000012E+01, 0.1000000000000012E+01, 0.1000000000000012E+01,
	                          0.1000000000000012E+01, 0.1000000000000012E+01, 0.1000000000000012E+01,
	                          0.1000000000000012E+01, 0.1000000000000012E+01, 0.1000000000000012E+01,
	                          0.1000000000000012E+01, 0.1000000000000012E+01, 0.9999999999996652E+00};
	        return x;
	        case ( 51): x = new double[] { 0.1000000000000012E+01, 0.1000000000000012E+01, 0.1000000000000012E+01,
	                          0.1000000000000012E+01, 0.1000000000000012E+01, 0.1000000000000012E+01,
	                          0.1000000000000012E+01, 0.1000000000000012E+01, 0.1000000000000012E+01,
	                          0.1000000000000012E+01, 0.1000000000000012E+01, 0.1000000000000012E+01,
	                          0.1000000000000012E+01, 0.1000000000000012E+01, 0.1000000000000012E+01,
	                          0.1000000000000012E+01, 0.1000000000000012E+01, 0.1000000000000012E+01,
	                          0.1000000000000012E+01, 0.1000000000000012E+01, 0.1000000000000012E+01,
	                          0.1000000000000012E+01, 0.1000000000000012E+01, 0.1000000000000012E+01,
	                          0.1000000000000012E+01, 0.1000000000000012E+01, 0.1000000000000012E+01,
	                          0.1000000000000012E+01, 0.1000000000000012E+01, 0.1000000000000012E+01,
	                          0.1000000000000012E+01, 0.1000000000000012E+01, 0.1000000000000012E+01,
	                          0.1000000000000012E+01, 0.1000000000000012E+01, 0.1000000000000012E+01,
	                          0.1000000000000012E+01, 0.1000000000000012E+01, 0.1000000000000012E+01,
	                          0.9999999999995192E+00};
	        return x;
	        case ( 52): x = new double[] { 0.3754100561705969E+00, 0.1935847320681608E+01,-0.1464687548211473E+01,
	                          0.1286753549696225E-01, 0.2212269807171213E-01};
	        return x;
	        case ( 53): x = new double[] { 0.1309976637986268E+01, 0.4315524809261695E+00, 0.6336612618356566E+00,
	                          0.5994285609242268E+00, 0.7541797688621015E+00, 0.9043000800785158E+00,
	                          0.1365799494307133E+01, 0.4823732003709051E+01, 0.2398684750507853E+01,
	                          0.4568875548318950E+01, 0.5675342062895008E+01};
	        return x;

	        default:
	            x = new double[] {10};
				return x;
	        }
	     }
//
//	    end function solution
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
//		
	        public static void ssqfcn(int m, int n, double[] x, double[] fVec, int nProb) {
//	    subroutine ssqfcn(m, n, x, fVec, nProb)
//	        implicit none
//
//	        integer,intent(in) :: m !! positive integer input variable.
//	        integer,intent(in) :: n !! positive integer input variable. n must not exceed m.
//	        integer,intent(in) :: nProb !! a positive integer input variable which defines the
//	                                    !! number of the problem. nProb must not exceed 18.
//	        real(wp),intent(in) :: x(n) !! an input array of length n.
//	        real(wp),intent(out) :: fVec(m) !! an output array of length m which contains the nProb
//	                                        !! function evaluated at x.
//
	        double zero  = 0.0;
	        double zp25  = 2.5e-1;
	        double zp5   = 5.0e-1;
	        double one   = 1.0;
	        double two   = 2.0;
	        double five  = 5.0;
	        double eight = 8.0;
	        double ten   = 10.0;
	        double c13   = 13.0;
	        double c14   = 14.0;
	        double c29   = 29.0;
	        double c45   = 45.0;

	        double[] v  = new double[] {4.0e0, 2.0e0, 1.0e0, 5.0e-1, 
	                2.5e-1, 1.67e-1, 1.25e-1, 1.0e-1, 8.33e-2, 7.14e-2, 
	                6.25e-2};
	        double[] y1 = new double[]  {1.4e-1, 1.8e-1, 2.2e-1, 2.5e-1, 
	                2.9e-1, 3.2e-1, 3.5e-1, 3.9e-1, 3.7e-1, 5.8e-1, 
	                7.3e-1, 9.6e-1, 1.34e0, 2.1e0, 4.39e0};
	        double[] y2 = new double[]  {1.957e-1, 1.947e-1, 
	                1.735e-1, 1.6e-1, 8.44e-2, 6.27e-2, 4.56e-2, 3.42e-2,  
	                3.23e-2, 2.35e-2, 2.46e-2};
	        double[] y3 = new double[]  {3.478e4, 2.861e4, 2.365e4, 1.963e4,  
	                1.637e4, 1.372e4, 1.154e4, 9.744e3, 8.261e3, 7.03e3,   
	                6.005e3, 5.147e3, 4.427e3, 3.82e3, 3.307e3, 2.872e3};
	        double[] y4 =  new double[] {8.44e-1,
	                9.08e-1, 9.32e-1, 9.36e-1, 9.25e-1, 9.08e-1, 8.81e-1,  
	                8.5e-1, 8.18e-1, 7.84e-1, 7.51e-1, 7.18e-1, 6.85e-1,   
	                6.58e-1, 6.28e-1, 6.03e-1, 5.8e-1, 5.58e-1, 5.38e-1,   
	                5.22e-1, 5.06e-1, 4.9e-1, 4.78e-1, 4.67e-1, 4.57e-1,   
	                4.48e-1, 4.38e-1, 4.31e-1, 4.24e-1, 4.2e-1, 4.14e-1,   
	                4.11e-1, 4.06e-1};
	        double[] y5 =  new double[] {1.366e0, 
	                1.191e0, 1.112e0, 1.013e0, 9.91e-1, 8.85e-1, 8.31e-1,  
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

//	        fVec(1:m) = zero
	        fVec = new double[m];
	        for(i=0;i<m;i++) {
	        	fVec[i] = 0;
	        }
//	        ! function routine selector.

	        switch (nProb) {
	        case (2):
//	            ! LINEAR FUNCTION - RANK 1.
	            sum = zero;
	            for(j=0;j<n;j++){
	                sum = sum + dfloat(j)*x[j];
	            }
	            for(i=0;i<m;i++){
	                fVec[i] = dfloat(i)*sum - one;
	            }
		        break;
	        case (3):
//	            ! LINEAR FUNCTION - RANK 1 WITH ZERO COLUMNS AND ROWS.
	            sum = zero;
	            nm1 = n - 1;
	            if (nm1 >= 2) {
	               for(j=1;j<nm1;j++){
	                    sum = sum + dfloat(j)*x[j];
	                }
	            }
	            for(i=0;i<m;i++){
	                fVec[i] = dfloat(i - 1)*sum - one;
	            }
	            fVec[m] = -one;
		        break;
	        case (4):
//	            ! ROSENBROCK FUNCTION.
	            fVec[0] = ten*(x[1] - x[0]*x[0]);
	            fVec[1] = one - x[0];
		        break;
	        case (5):
//	            ! HELICAL VALLEY FUNCTION.
	            tpi = eight*Math.atan(one);
	            tmp1 = Math.signum(x[1]) * Math.abs(zp25);
	            if (x[0] > zero) tmp1 = Math.atan(x[1]/x[0])/tpi;
	            if (x[0] < zero) tmp1 = Math.atan(x[1]/x[0])/tpi + zp5;
	            tmp2 = Math.sqrt(x[0]*x[0] + x[1]*x[1]);
	            fVec[0] = ten*(x[2] - ten*tmp1);
	            fVec[1] = ten*(tmp2 - one);
	            fVec[2] = x[2];
		        break;
	        case (6):
//	            ! POWELL Math.sinGULAR FUNCTION.
	            fVec[0] = x[0] + ten*x[1];
	            fVec[1] = Math.sqrt(five)*(x[2] - x[3]);
	            fVec[2] = (x[1] - two*x[2])*(x[1] - two*x[2]);
	            fVec[3] = Math.sqrt(ten)*(x[0] - x[3])*(x[0] - x[3]);
		        break;
	        case (7):
//	            ! FREUDENSTEIN AND ROTH FUNCTION.
	            fVec[0] = -c13 + x[0] + ((five - x[1])*x[1] - two)*x[1];
	            fVec[1] = -c29 + x[0] + ((one + x[1])*x[1] - c14)*x[1];
		        break;
	        case (8):
//	            ! BARD FUNCTION.
	            for(i=0;i< 15;i++){
	                tmp1 = dfloat(i);
	                tmp2 = dfloat(16 - i);
	                tmp3 = tmp1;
	                if (i > 8) {tmp3 = tmp2;}
	                fVec[i] = y1[i] - (x[0] + tmp1/(x[1]*tmp2 + x[2]*tmp3));
	            }
	        break;
	        case (9):
//	            ! KOWALIK AND OSBORNE FUNCTION.
	            for(i=0;i< 11;i++){
	                tmp1 = v[i]*(v[i] + x[1]);
	                tmp2 = v[i]*(v[i] + x[2]) + x[3];
	                fVec[i] = y2[i] - x[0]*tmp1/tmp2;
	            }
	        break;
	        case (10):
//	            ! MEYER FUNCTION.
	            for(i=0;i<16 ;i++){
	                temp = five*dfloat(i) + c45 + x[2];
	                tmp1 = x[1]/temp;
	                tmp2 = Math.exp(tmp1);
	                fVec[i] = x[0]*tmp2 - y3[i];
	            }
	        break;
	        case (11):
//	            ! WATSON FUNCTION.
	            for(i=0;i<29 ;i++){
	                div = dfloat(i)/c29;
	                s1 = zero;
	                dx = one;
	               for(j=1;j<n;j++){
	                    s1 = s1 + dfloat(j - 1)*dx*x[j];
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
//	            ! BOX 3-DIMENSIONAL FUNCTION.
	            for(i=0;i<m;i++){
	                temp = dfloat(i);
	                tmp1 = temp/ten;
	                fVec[i] = Math.exp(-tmp1*x[0]) - Math.exp(-tmp1*x[1]) + (Math.exp(-temp) - Math.exp(-tmp1))*x[2];
	            }
	        break;
	        case (13):
//	            ! JENNRICH AND SAMPSON FUNCTION.
	            for(i=0;i<m;i++){
	                temp = dfloat(i);
	                fVec[i] = two + two*temp - Math.exp(temp*x[0]) - Math.exp(temp*x[1]);
	            }
	        break;
	        case (14):
//	            ! BROWN AND DENNIS FUNCTION.
	            for(i=0;i<m;i++){
	                temp = dfloat(i)/five;
	                tmp1 = x[0] + temp*x[1] - Math.exp(temp);
	                tmp2 = x[2] + Math.sin(temp)*x[3] - Math.cos(temp);
	                fVec[i] = tmp1*tmp1 + tmp2*tmp2;
	            }
	        break;
	        case (15):
//	            ! CHEBYQUAD FUNCTION.
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
	                if (iev > 0) { fVec[i] = fVec[i] + one/(dfloat(i)*dfloat(i) - one);}
	                iev = -iev;
	            }
		        break;
	        case (16):
//	            ! BROWN ALMOST-LINEAR FUNCTION.
	            sum = -dfloat(n + 1);
	            prod = one;
	            for(j=0;j<n;j++){
	                sum = sum + x[j];
	                prod = x[j]*prod;
	            }
	            for(i=0;i< n;i++){
	                fVec[i] = x[i] + sum;
	            }
	            fVec[n] = prod - one;
		        break;
	        case (17):
//	            ! OSBORNE 1 FUNCTION.
	            for(i=0;i< 33;i++){
	                temp = ten*dfloat(i - 1);
	                tmp1 = Math.exp(-x[3]*temp);
	                tmp2 = Math.exp(-x[4]*temp);
	                fVec[i] = y4[i] - (x[0] + x[1]*tmp1 + x[2]*tmp2);
	            }
	        break;
	        case (18):
//	            ! OSBORNE 2 FUNCTION.
	            for(i=0;i< 65;i++){
	                temp = dfloat(i - 1)/ten;
	                tmp1 = Math.exp(-x[4]*temp);
	                tmp2 = Math.exp(-x[5]*(temp - x[8])*(temp - x[8]));
	                tmp3 = Math.exp(-x[6]*(temp - x[9])*(temp - x[9]));
	                tmp4 = Math.exp(-x[7]*(temp - x[10])*(temp - x[10]));
	                fVec[i] = y5[i] - (x[0]*tmp1 + x[1]*tmp2 + x[2]*tmp3 + x[3]*tmp4);
	            }
	        break;
	       default:
//	            ! LINEAR FUNCTION - FULL RANK.
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
	        }

	        }
//	    end subroutine ssqfcn
//	!*****************************************************************************************
//
//	!*****************************************************************************************
//	!>
//	!  This subroutine specifies the standard starting points for the
//	!  functions defined by subroutine ssqfcn. the subroutine returns
//	!  in x a multiple (factor) of the standard starting point. for
//	!  the 11th function the standard starting point is zero, so in
//	!  this case, if factor is not unity, { the subroutine returns
//	!  the vector  x[j] = factor, j=1,...,n.
		public static void initpt(int n,double[] x,int nProb, double factor){
//	    subroutine initpt(n, x, nProb, factor)
//	        implicit none
//
//	        integer,intent(in) :: n !! a positive integer input variable.
//	        integer,intent(in) :: nProb !! an input variable which specifies the multiple of
//	                                    !! the standard starting point. if factor is unity, no
//	                                    !! multiplication is performed.
//	        real(wp),intent(in) :: factor !! an input variable which specifies the multiple of
//	                                      !! the standard starting point. if factor is unity, no
//	                                      !! multiplication is performed.
//	        real(wp),intent(out) :: x(n) !! an output array of length n which contains the standard
//	                                     !! starting point for problem nProb multiplied by factor.

	        double zero   = 0.0;
	        double half   = 0.5;
	        double one    = 1.0;
	        double two    = 2.0;
	        double three  = 3.0;
	        double five   = 5.0;
	        double seven  = 7.0;
	        double ten    = 10.0;
	        double twenty = 20.0;
	        double twntf  = 25.0;

	        double c1  = 1.2;
	        double c2  = 2.5e-1;
	        double c3  = 3.9e-1;
	        double c4  = 4.15e-1;
	        double c5  = 2.0e-2;
	        double c6  = 4.0e3;
	        double c7  = 2.5e2;
	        double c8  = 3.0e-1;
	        double c9  = 4.0e-1;
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

//	        x(1:n) = zero
	        x = new double[n];
	        for(int i=0;i<n;i++) {
	        	x[i] = 0;
	        }

//	        ! selection of initial point.

	        switch (nProb) {
	        case (4):
//	            ! ROSENBROCK FUNCTION.
	            x[0] = -c1;
	            x[1] = one;
				break;
	        case (5):
//	            ! HELICAL VALLEY FUNCTION.
	            x[0] = -one;
	            x[1] = zero;
	            x[2] = zero;
				break;
	        case (6):
//	            ! POWELL Math.sinGULAR FUNCTION.
	            x[0] = three;
	            x[1] = -one;
	            x[2] = zero;
	            x[3] = one;
				break;
	        case (7):
//	            ! FREUDENSTEIN AND ROTH FUNCTION.
	            x[0] = half;
	            x[1] = -two;
				break;
	        case (8):
//	            ! BARD FUNCTION.
	            x[0] = one;
	            x[1] = one;
	            x[2] = one;
				break;
	        case (9):
//	            ! KOWALIK AND OSBORNE FUNCTION.
	            x[0] = c2;
	            x[1] = c3;
	            x[2] = c4;
	            x[3] = c3;
				break;
	        case (10):
//	            ! MEYER FUNCTION.
	            x[0] = c5;
	            x[1] = c6;
	            x[2] = c7;
				break;
	        case (11):
//	            ! WATSON FUNCTION.
	            for(j=0;j<n;j++){
	                x[j] = zero;
	            }
				break;
	        case (12):
//	            ! BOX 3-DIMENSIONAL FUNCTION.
	            x[0] = zero;
	            x[1] = ten;
	            x[2] = twenty;
				break;
	        case (13):
//	            ! JENNRICH AND SAMPSON FUNCTION.
	            x[0] = c8;
	            x[1] = c9;
				break;
	        case (14):
//	            ! BROWN AND DENNIS FUNCTION.
	            x[0] = twntf;
	            x[1] = five;
	            x[2] = -five;
	            x[3] = -one;
				break;
	        case (15):
//	            ! CHEBYQUAD FUNCTION.
	            h = one/dfloat(n + 1);
	            for(j=0;j<n;j++){
	                x[j] = dfloat(j)*h;
	            }
				break;
	        case (16):
//	            ! BROWN ALMOST-LINEAR FUNCTION.
	            for(j=0;j<n;j++){
	                x[j] = half;
	            }
				break;
	        case (17):
//	            ! OSBORNE 1 FUNCTION.
	            x[0] = half;
	            x[1] = c10;
	            x[2] = -one;
	            x[3] = c11;
	            x[4] = c5;
				break;
	        case (18):
//	            ! OSBORNE 2 FUNCTION.
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
//	            ! LINEAR FUNCTION - FULL RANK OR RANK 1.
	            for(j=0;j<n;j++){
	                x[j] = one;
	            }
				break;
	        }

//	        ! compute multiple of initial point.
	        if (factor != one) {
	            if (nProb == 10) {
	                for(j=0;j<n;j++){
	                    x[j] = factor;
	                }
	            } else {
	                for(j=0;j<n;j++){
	                    x[j] = factor*x[j];
	                }
	            }
	        }
		}
//
//	    end subroutine initpt
//	!*****************************************************************************************
//
//	!*****************************************************************************************
//	    end program test_lmdif
////	!*****************************************************************************************
}
