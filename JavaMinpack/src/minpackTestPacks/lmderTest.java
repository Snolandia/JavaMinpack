package minpackTestPacks;
//
public class lmderTest {
//	!*****************************************************************************************
//	!>
//	!  This program tests codes for the least-squares solution of
//	!  M nonlinear equations in n variables. it consists of a driver
//	!  And an interface subroutine fcn. the driver reads in data,
//	!  Calls the nonlinear least-squares solver, and finally prints
//	!  Out information on the performance of the solver. this is
//	!  Only a sample driver, many other drivers are possible. the
//	!  Interface subroutine fcn is necessary to take into account the
//	!  Forms of calling sequences used by the function and jacobian
//	!  Subroutines in the various nonlinear least-squares solvers.
//
//	program test_lmder
//
//	    use minpack_module, only: wp, dpmpar, enorm, lmder1
//	    use iso_fortran_env, only: output_unit
//
//	    implicit none
//
//	    ! originally from file22
//	    integer,parameter :: ncases = 28
//	    integer,dimension(ncases),parameter :: nProbs  = [1,1,2,2,3,3,4,5,6,7,8,9,10,11,11,11,12,13,14,15,15,15,15,16,16,16,17,18]
//	    integer,dimension(ncases),parameter :: ns      = [5,5,5,5,5,5,2,3,4,2,3,4,3,6,9,12,3,2,4,1,8,9,10,10,30,40,5,11]
//	    integer,dimension(ncases),parameter :: ms      = [10,50,10,50,10,50,2,3,4,2,15,11,16,31,31,31,10,10,20,8,8,9,10,10,30,40,33,65]
//	    integer,dimension(ncases),parameter :: ntriess = [1,1,1,1,1,1,3,3,3,3,3,3,2,3,3,3,1,1,3,3,1,1,1,3,1,1,1,1]
//
//	    integer,dimension(53),parameter :: info_original = [3,3,1,1,1,1,4,2,2,2,2,2,4,4,4,1,1,1,1,1,1,
//	                                                        1,1,5,2,5,1,1,1,3,1,3,3,3,2,2,1,1,1,1,4,1,
//	                                                        1,1,2,1,2,2,2,2,2,1,1] !! original `info` from the original minpack
//
//	    int i, ic, info, k, ldfjac, lwa, m, n, NFEv, NJEv, nProb, ntries, icase, iunit
//	   double factor, fnorm1, fnorm2
//	    int ma(53), na(53), nf(53), nj(53), np(53), nx(53)
//	   double fnm(53)
//	    integer,dimension(:),allocatable :: iwa
//	    real(wp),dimension(:),allocatable :: fVec, wa, x
//	    real(wp),dimension(:,:),allocatable :: fjac
//
//	    integer, parameter :: nwrite = output_unit ! logical output unit
//	    double one = 1.0
//	    double ten = 10.0
//	    double tol = Math.sqrt(dpmpar(1))
//	    double solution_reltol = 1.0e-4 !! reltol for matching previously generated solutions
//
//	    ic = 0
//	    do icase = 1, ncases+1
//
//	        if (icase == ncases+1) {
//	            write (nwrite, '(A,I3,A/)') '1SUMMARY OF ', ic, ' CALLS TO LMDER1'
//	            write (nwrite, '(A/)')      ' nProb   N    M   NFEV  NJEV  INFO  FINAL L2 NORM'
//	            for(i=0;i< ;i++){ic
//	                write (nwrite, '(3I5,3I6,1X,D15.7)') np(i), na(i), ma(i), nf(i), nj(i), nx[i], fnm(i)
//	            }
//	            stop
//	        }else{
//
//	            nProb = nProbs(icase)
//	            n = ns(icase)
//	            m = ms(icase)
//	            lwa = 5*n+m
//	            ldfjac = m
//
//	            if (allocated(fjac)) deallocate(fjac); allocate(fJac[m,n))
//	            if (allocated(fVec)) deallocate(fVec); allocate(fVec[m])
//	            if (allocated(wa))   deallocate(wa);   allocate(wa(lwa))
//	            if (allocated(x))    deallocate(x);    allocate(x(n))
//	            if (allocated(iwa))  deallocate(iwa);  allocate(iwa(n))
//
//	            ntries = ntriess(icase)
//
//	            factor = one
//	            for(k=0;k<n;k++){tries
//	                ic = ic + 1
//	                call initpt(n, x, nProb, factor)
//	                call ssqfcn(m, n, x, fVec, nProb)
//	                fnorm1 = enorm(m, fVec)
//	                write (nwrite, '(////5X,A,I5,5X,A,2I5,5X//)') ' PROBLEM', nProb, ' DIMENSIONS', n, m
//	                NFEv = 0
//	                NJEv = 0
//	                call lmder1(fcn, m, n, x, fVec, fjac, ldfjac, tol, info, iwa, wa, lwa)
//	                call ssqfcn(m, n, x, fVec, nProb)
//	                fnorm2 = enorm(m, fVec)
//	                np(ic) = nProb
//	                na(ic) = n
//	                ma(ic) = m
//	                nf(ic) = NFEv
//	                nj(ic) = NJEv
//	                nx(ic) = info
//	                fnm(ic) = fnorm2
//	                write (nwrite, '(5X,A,D15.7//5X,A,D15.7//5X,A,I10//5X,A,I10//5X,A,18X,I10//5X,A//,*(5X,5D15.7/))') 
//	                            ' INITIAL L2 NORM OF THE RESIDUALS', fnorm1, 
//	                            ' FINAL L2 NORM OF THE RESIDUALS  ', fnorm2, 
//	                            ' NUMBER OF FUNCTION EVALUATIONS  ', NFEv, 
//	                            ' NUMBER OF JACOBIAN EVALUATIONS  ', NJEv, 
//	                            ' EXIT PARAMETER', info, 
//	                            ' FINAL APPROXIMATE SOLUTION',x(1:n)
//	                factor = ten*factor
//	                call compare_solutions(ic, x, solution_reltol, tol)
//	            }
//	        }
//	    }
//
//	contains
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
//	    if (info_original(ic)<5) {    ! ignore any where the original minpack failed
//	        diff = solution(ic) - x
//	        absdiff = abs(diff)
//	        if (any(absdiff>abstol)) { ! first do an absolute diff
//	            ! also do a rel diff if the abs diff fails (also protect for divide by zero)
//	            reldiff = absdiff
//	            where (solution(ic) != 0.0) reldiff = absdiff / abs(solution(ic))
//	            if (any(reldiff > reltol)) {
//	                write(nwrite,'(A)') 'Failed case'
//	                write(nwrite, '(//5x, a//(5x, 5d15.7))') 'Math.expected x: ', solution(ic)
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
//	!  function and jacobian subroutines ssqfcn and ssqjac with
//	!  the appropriate value of problem number (nProb).
//
//	    subroutine fcn(m, n, x, fVec, Fjac, Ldfjac, Iflag)
//
//	        implicit none
//
//	        integer, intent(in) :: m !! the number of functions.
//	        integer, intent(in) :: n !! the number of variables.
//	        integer, intent(in) :: ldfjac !! leading dimension of the array fjac.
//	        integer, intent(inout) :: iflag !! if iflag = 1 calculate the functions at x and
//	                                       !! return this vector in fVec. do not alter fjac.
//	                                       !! if iflag = 2 calculate the jacobian at x and
//	                                       !! return this matrix in fjac. do not alter fVec.
//	                                       !!
//	                                       !! the value of iflag should not be changed by fcn unless
//	                                       !! the user wants to terminate execution of lmder.
//	                                       !! in this case set iflag to a negative integer.
//	        real(wp), intent(in) :: x(n) !! independent variable vector
//	        real(wp), intent(inout) :: fVec[m] !! value of function at `x`
//	        real(wp), intent(inout) :: fJac[ldfjac, n) !! jacobian matrix at `x`
//
//	        switch (iflag)
//	        case (1)
//	            call ssqfcn(m, n, x, fVec, nProb)
//	            NFEv = NFEv + 1
//	        case (2)
//	            call ssqjac(m, n, x, Fjac, Ldfjac, nProb)
//	            NJEv = NJEv + 1
//	        case default
//	            error stop 'invalid iflag value'
//	        }
//
//	    end subroutine fcn
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
//	       double f
//	        f = real(i, wp)
//	    end function dfloat
//	!*****************************************************************************************
//
//	!*****************************************************************************************
//	!>
//	!  Get Math.expected `x` vectors for each case.
		public static double[] solution(int nProb){
//	    pure function solution(nProb) result(x)
//
//	       // implicit none
//
//	        //integer,intent(in) :: nProb
//	        //real(wp),dimension(:),allocatable :: x
double[] x;

	        switch (nProb) {
	        case(1): x =  new double []  {-0.1000000000000000E+01,-0.1000000000000000E+01,-0.1000000000000000E+01,
	                      -0.1000000000000000E+01,-0.1000000000000000E+01};
return x;
	        case(2): x =  new double []  {-0.9999999999999953E+00,-0.1000000000000005E+01,-0.9999999999999976E+00,
	                      -0.9999999999999956E+00,-0.9999999999999991E+00};
return x;
	        case(3): x =  new double []  {-0.1677968180239693E+03,-0.8339840901198468E+02,0.2211100430795781E+03,
	                      -0.4119920450599233E+02,-0.3275936360479385E+02};
return x;
	        case(4): x =  new double []  {-0.2029999900022674E+02,-0.9649999500113370E+01,-0.1652451975264496E+03,
	                      -0.4324999750056676E+01,0.1105330585100652E+03};
return x;
	        case(5): x =  new double []  {0.1000000000000000E+01,-0.2103615324224772E+03,0.3212042081132130E+02,
	                      0.8113456824980642E+02,0.1000000000000000E+01};
return x;
	        case(6): x =  new double []  {0.1000000000000000E+01,0.3321494858957815E+03,-0.4396851914289522E+03,
	                      0.1636968825825863E+03,0.1000000000000000E+01};
return x;
	        case(7): x =  new double []  {0.1000000000000000E+01,0.1000000000000000E+01};
return x;
	        case(8): x =  new double []  {0.1000000000000000E+01,0.1000000000000000E+01};
return x;
	        case(9): x =  new double []  {0.1000000000000000E+01,0.1000000000000000E+01};
return x;
	        case(10): x =  new double []  {0.1000000000000000E+01,-0.6243301596789443E-17,0.0000000000000000E+00};
return x;
	        case(11): x =  new double []  {0.1000000000000000E+01,0.6563910805155555E-20,0.0000000000000000E+00};
return x;
	        case(12): x =  new double []  {0.1000000000000000E+01,-0.1972152263052530E-29,0.0000000000000000E+00};
return x;
	        case(13): x =  new double []  {0.1652117596168394E-16,-0.1652117596168393E-17,0.2643388153869468E-17,
	                       0.2643388153869468E-17};
return x;
	        case(14): x =  new double []  {0.2016745112510229E-19,-0.2016745112510229E-20,0.3226792180016300E-20,
	                       0.3226792180016300E-20};
return x;
	        case(15): x =  new double []  {0.3226792180016378E-17,-0.3226792180016378E-18,0.5162867488026213E-18,
	                       0.5162867488026213E-18};
return x;
	        case(16): x =  new double []  {0.1141248446549937E+02,-0.8968279137315035E+00};
return x;
	        case(17): x =  new double []  {0.1141300466147456E+02,-0.8967960386859591E+00};
return x;
	        case(18): x =  new double []  {0.1141278178578820E+02,-0.8968051074920677E+00};
return x;
	        case(19): x =  new double []  {0.8241057657583339E-01,0.1133036653471504E+01,0.2343694638941154E+01};
return x;
	        case(20): x =  new double []  {0.8406666738183293E+00,-0.1588480332595655E+09,-0.1643786716535352E+09};
return x;
	        case(21): x =  new double []  {0.8406666738676455E+00,-0.1589461672055184E+09,-0.1644649068577712E+09};
return x;
	        case(22): x =  new double []  {0.1928078104762493E+00,0.1912626533540709E+00,0.1230528010469309E+00,
	                       0.1360532211505167E+00};
return x;
	        case(23): x =  new double []  {0.7286754737686598E+06,-0.1407588031293926E+02,-0.3297779778419661E+08,
	                       -0.2057159419780170E+08};
return x;
	        case(24): x =  new double []  {0.1927984063846549E+00,0.1914736844615448E+00,0.1230924753714115E+00,
	                       0.1361509629062244E+00};
return x;
	        case(25): x =  new double []  {0.5609636471026654E-02,0.6181346346286586E+04,0.3452236346241439E+03};
return x;
	        case(26): x =  new double []  {0.1291118238019656E-10,0.3388670577625951E+05,0.9040289597481549E+03};
return x;
	        case(27): x =  new double []  {-0.1572496150837821E-01,0.1012434882329655E+01,-0.2329917223876743E+00,
	                       0.1260431011028187E+01,-0.1513730313944210E+01,0.9929972729184218E+00};
return x;
	        case(28): x =  new double []  {-0.1572519013866764E-01,0.1012434858601050E+01,-0.2329915458438263E+00,
	                       0.1260429320891618E+01,-0.1513727767065737E+01,0.9929957342632760E+00};
return x;
	        case(29): x =  new double []  {-0.1572470197125875E-01,0.1012434909256583E+01,-0.2329919227616435E+00,
	                       0.1260432929295550E+01,-0.1513733204527069E+01,0.9929990192232205E+00};
return x;
	        case(30): x =  new double []  {-0.1530706441663782E-04,0.9997897039345963E+00,0.1476396349113067E-01,
	                       0.1463423301456550E+00,0.1000821094549679E+01,-0.2617731120708589E+01,
	                       0.4104403139437982E+01,-0.3143612262365286E+01,0.1052626403788335E+01};
return x;
	        case(31): x =  new double []  {-0.1530713348497491E-04,0.9997897039412341E+00,0.1476396297861486E-01,
	                       0.1463423348188786E+00,0.1000821073213774E+01,-0.2617731070847226E+01,
	                       0.4104403076555856E+01,-0.3143612221787109E+01,0.1052626393225985E+01};
return x;
	        case(32): x =  new double []  {-0.1530703652788392E-04,0.9997897039319498E+00,0.1476396369347192E-01,
	                       0.1463423283001695E+00,0.1000821103001098E+01,-0.2617731140510708E+01,
	                       0.4104403164468549E+01,-0.3143612278549893E+01,0.1052626408008490E+01};
return x;
	        case(33): x =  new double []  {-0.6638060467011595E-08,0.1000001644117862E+01,-0.5639322103223269E-03,
	                       0.3478205405036302E+00,-0.1567315040922755E+00,0.1052815177186520E+01,
	                       -0.3247271153392391E+01,0.7288434897798658E+01,-0.1027184824113737E+02,
	                       0.9074113646924774E+01,-0.4541375466623413E+01,0.1012011888539188E+01};
return x;
	        case(34): x =  new double []  {-0.6637102237178952E-08,0.1000001644117874E+01,-0.5639322084840861E-03,
	                       0.3478205404902133E+00,-0.1567315039874129E+00,0.1052815176715104E+01,
	                       -0.3247271152062481E+01,0.7288434895390663E+01,-0.1027184823833726E+02,
	                       0.9074113644905919E+01,-0.4541375465802633E+01,0.1012011888395745E+01};
return x;
	        case(35): x =  new double []  {-0.6638060465675741E-08,0.1000001644117862E+01,-0.5639322102727471E-03,
	                       0.3478205405024820E+00,-0.1567315040810628E+00,0.1052815177127253E+01,
	                       -0.3247271153204238E+01,0.7288434897423372E+01,-0.1027184824066345E+02,
	                       0.9074113646557050E+01,-0.4541375466463498E+01,0.1012011888509363E+01};
return x;
	        case(36): x =  new double []  {0.9999999999999999E+00,0.1000000000000000E+02,0.1000000000000000E+01};
return x;
	        case(37): x =  new double []  {0.2578199266367936E+00,0.2578299767645596E+00};
return x;
	        case(38): x =  new double []  {-0.1159125521567487E+02,0.1320248976116189E+02,-0.4035744915102742E+00,
	                       0.2367363118370096E+00};
return x;
	        case(39): x =  new double []  {-0.1159592742721295E+02,0.1320418669262178E+02,-0.4034173628612446E+00,
	                       0.2367711433960093E+00};
return x;
	        case(40): x =  new double []  {-0.1159025975515195E+02,0.1320206290836587E+02,-0.4036882104139269E+00,
	                       0.2366649366417166E+00};
return x;
	        case(41): x =  new double []  {0.5000000000000000E+00};
return x;
	        case(42): x =  new double []  {0.9817314924683995E+00};
return x;
	        case(43): x =  new double []  {0.9817314852933997E+00};
return x;
	        case(44): x =  new double []  {0.4315366485873882E-01,0.1930916378432737E+00,0.2663285938126948E+00,
	                       0.4999993346289514E+00,0.5000006653710486E+00,0.7336714061873052E+00,
	                       0.8069083621567262E+00,0.9568463351412612E+00};
return x;
	        case(45): x =  new double []  {0.4420534613578274E-01,0.1994906723098809E+00,0.2356191084710600E+00,
	                       0.4160469078925979E+00,0.5000000000000001E+00,0.5839530921074019E+00,
	                       0.7643808915289401E+00,0.8005093276901190E+00,0.9557946538642172E+00};
return x;
	        case(46): x =  new double []  {0.5962026717536065E-01,0.1667087838059440E+00,0.2391710188135143E+00,
	                       0.3988852903458887E+00,0.3988836678710648E+00,0.6011163321289351E+00,
	                       0.6011147096541113E+00,0.7608289811864857E+00,0.8332912161940560E+00,
	                       0.9403797328246394E+00};
return x;
	        case(47): x =  new double []  {0.9794303033498624E+00,0.9794303033498624E+00,0.9794303033498624E+00,
	                       0.9794303033498624E+00,0.9794303033498624E+00,0.9794303033498624E+00,
	                       0.9794303033498624E+00,0.9794303033498624E+00,0.9794303033498624E+00,
	                       0.1205696966501376E+01};
return x;
	        case(48): x =  new double []  {0.9794303033498634E+00,0.9794303033498634E+00,0.9794303033498634E+00,
	                       0.9794303033498634E+00,0.9794303033498634E+00,0.9794303033498634E+00,
	                       0.9794303033498634E+00,0.9794303033498634E+00,0.9794303033498634E+00,
	                       0.1205696966501366E+01};
return x;
	        case(49): x =  new double []  {0.1000000000000002E+01,0.1000000000000002E+01,0.1000000000000002E+01,
	                       0.1000000000000002E+01,0.1000000000000002E+01,0.1000000000000002E+01,
	                       0.1000000000000002E+01,0.1000000000000002E+01,0.1000000000000002E+01,
	                       0.9999999999999840E+00};
return x;
	        case(50): x =  new double []  {0.9977542164428221E+00,0.9977542164428221E+00,0.9977542164428221E+00,
	                       0.9977542164428221E+00,0.9977542164428221E+00,0.9977542164428221E+00,
	                       0.9977542164428221E+00,0.9977542164428221E+00,0.9977542164428221E+00,
	                       0.9977542164428221E+00,0.9977542164428221E+00,0.9977542164428221E+00,
	                       0.9977542164428221E+00,0.9977542164428221E+00,0.9977542164428221E+00,
	                       0.9977542164428221E+00,0.9977542164428221E+00,0.9977542164428221E+00,
	                       0.9977542164428221E+00,0.9977542164428221E+00,0.9977542164428221E+00,
	                       0.9977542164428221E+00,0.9977542164428221E+00,0.9977542164428221E+00,
	                       0.9977542164428221E+00,0.9977542164428221E+00,0.9977542164428221E+00,
	                       0.9977542164428221E+00,0.9977542164428221E+00,0.1067373506715322E+01};
return x;
	        case(51): x =  new double []  {0.1000000000000006E+01,0.1000000000000006E+01,0.1000000000000006E+01,
	                       0.1000000000000006E+01,0.1000000000000006E+01,0.1000000000000006E+01,
	                       0.1000000000000006E+01,0.1000000000000006E+01,0.1000000000000006E+01,
	                       0.1000000000000006E+01,0.1000000000000006E+01,0.1000000000000006E+01,
	                       0.1000000000000006E+01,0.1000000000000006E+01,0.1000000000000006E+01,
	                       0.1000000000000006E+01,0.1000000000000006E+01,0.1000000000000006E+01,
	                       0.1000000000000006E+01,0.1000000000000006E+01,0.1000000000000006E+01,
	                       0.1000000000000006E+01,0.1000000000000006E+01,0.1000000000000006E+01,
	                       0.1000000000000006E+01,0.1000000000000006E+01,0.1000000000000006E+01,
	                       0.1000000000000006E+01,0.1000000000000006E+01,0.1000000000000006E+01,
	                       0.1000000000000006E+01,0.1000000000000006E+01,0.1000000000000006E+01,
	                       0.1000000000000006E+01,0.1000000000000006E+01,0.1000000000000006E+01,
	                       0.1000000000000006E+01,0.1000000000000006E+01,0.1000000000000006E+01,
	                       0.9999999999997748E+00};
return x;
	        case(52): x =  new double []  {0.3754100492440257E+00,0.1935846545431133E+01,-0.1464686767487213E+01,
	                       0.1286753391104403E-01,0.2212270118130737E-01};
return x;
	        case(53): x =  new double []  {0.1309976638100963E+01,0.4315524807599996E+00,0.6336612616028594E+00,
	                       0.5994285609916951E+00,0.7541797682724486E+00,0.9043000823785187E+00,
	                       0.1365799495210074E+01,0.4823731997481072E+01,0.2398684751048711E+01,
	                       0.4568875547914517E+01,0.5675342062730520E+01};
return x;
	        default:
	            x = new double[] {10};
				return x;
	        }
}
//	    end function solution
//	!*****************************************************************************************
//
//	!*****************************************************************************************
//	!
//	!  This subroutine defines the jacobian matrices of eighteen
//	!  nonlinear least squares problems. The problem dimensions are
//	!  as described in the prologue comments of [[ssqfcn]].
		public static void ssqjac(int m, int n, double[] x, double[][] fJac, int ldfJac, int nProb){
//	    subroutine ssqjac(m, n, x, Fjac, Ldfjac, nProb)
//	        implicit none
//
//	        integer, intent(in) :: m !! positive integer input variables. n must not exceed m.
//	        integer, intent(in) :: n !! positive integer input variables. n must not exceed m.
//	        integer, intent(in) :: ldfjac !! a positive integer input variable not less than m
//	                                   !! which specifies the leading dimension of the array fjac.
//	        integer, intent(in) :: nProb !! a positive integer variable which defines the
//	                                  !! number of the problem. nProb must not exceed 18.
//	        real(wp), intent(in) :: x(n) !! an input array of length n.
//	        real(wp), intent(out) :: fJac[ldfjac, n) !! an m by n output array which contains the jacobian
//	                                             !! matrix of the nProb function evaluated at x.

	        double[] v = new double[] {4.0, 2.0, 1.0, 5.0e-1, 
	                                        2.5e-1, 1.67e-1, 1.25e-1, 1.0e-1, 
	                                        8.33e-2, 7.14e-2, 6.25e-2};

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

	        int i, ivar, j, k, mm1, nm1;
	       double div, dx, prod, s2, temp, ti, tmp1, tmp2, tmp3, tmp4, tpi;

//	        fJac[1:m, 1:n) = zero
	       fJac = new double[m][n];
	       for(i = 0;i<m;i++) {
		       for(j = 0;j<n;i++) {
		    	   fJac[i][j] = 0;
		       }
	       }

//	        ! JACOBIAN ROUTINE SELECTOR.

	        switch (nProb){
	        case (2):

	            //! LINEAR FUNCTION - RANK 1.

	            for(j=0;j<n;j++){
	                for(i=0;i<m;i++){
	                    fJac[i][j] = dfloat(i)*dfloat(j);
	                }
	            }
break;

	        case (3):

	            //! LINEAR FUNCTION - RANK 1 WITH ZERO COLUMNS AND ROWS.

	            for(j=0;j<n;j++){
	                for(i=0;i<m;i++){
	                    fJac[i][j] = zero;
	                }
	            }
	            nm1 = n - 1;
	            mm1 = m - 1;
	            if (nm1 >= 2) {
	               for(j = 1;i< mm1;i++) {
	                    for(i = 1;i< mm1;i++) {
	                        fJac[i][j] = dfloat(i - 1)*dfloat(j);
	                    }
	                }
	            }
break;

	        case (4):

	            //! ROSENBROCK FUNCTION.

	            fJac[0][0] = -c20*x[0];
	            fJac[0][1] = ten;
	            fJac[1][0] = -one;
	            fJac[1][1] = zero;
break;

	        case (5):

	            //! HELICAL VALLEY FUNCTION.

	            tpi = eight*Math.atan(one);
	            temp = x[0]*x[0] + x[1]*x[1];
	            tmp1 = tpi*temp;
	            tmp2 = Math.sqrt(temp);
	            fJac[0][0] = c100*x[1]/tmp1;
	            fJac[0][1] = -c100*x[0]/tmp1;
	            fJac[0][2] = ten;
	            fJac[1][0] = ten*x[0]/tmp2;
	            fJac[1][1] = ten*x[1]/tmp2;
	            fJac[1][2] = zero;
	            fJac[2][0] = zero;
	            fJac[2][1] = zero;
	            fJac[2][2] = one;
break;

	        case (6):

	            //! POWELL SINGULAR FUNCTION.

	            for(j = 0;j< 4;j++){
	                for(i=0;i< 4;i++){
	                    fJac[i][j] = zero;
	                }
	            }
	            fJac[0][0] = one;
	            fJac[0][1] = ten;
	            fJac[1][2] = Math.sqrt(five);
	            fJac[1][3] = -fJac[1][2];
	            fJac[2][1] = two*(x[1] - two*x[2]);
	            fJac[2][2] = -two*fJac[2][1];
	            fJac[3][0] = two*Math.sqrt(ten)*(x[0] - x[3]);
	            fJac[3][3] = -fJac[3][0];

break;
	        case (7):

	            //! FREUDENSTEIN AND ROTH FUNCTION.

	            fJac[0][0] = one;
	            fJac[0][1] = x[1]*(ten - three*x[1]) - two;
	            fJac[1][2] = one;
	            fJac[1][3] = x[1]*(two + three*x[1]) - c14;
break;

	        case (8):

	            //! BARD FUNCTION.

	            for(i=0;i< 15;i++){
	                tmp1 = dfloat(i);
	                tmp2 = dfloat(16 - i);
	                tmp3 = tmp1;
	                if (i > 7) {tmp3 = tmp2;}
	                tmp4 = (x[1]*tmp2 + x[2]*tmp3)*(x[1]*tmp2 + x[2]*tmp3);
	                fJac[i][0] = -one;
	                fJac[i][1] = tmp1*tmp2/tmp4;
	                fJac[i][2] = tmp1*tmp3/tmp4;
	            }
break;

	        case (9):

	            //! KOWALIK AND OSBORNE FUNCTION.

	            for(i=0;i< 11;i++){
	                tmp1 = v[i]*(v[i] + x[1]);
	                tmp2 = v[i]*(v[i] + x[2]) + x[3];
	                fJac[i][0] = -tmp1/tmp2;
	                fJac[i][1] = -v[i]*x[0]/tmp2;
	                fJac[i][2] = fJac[i][0]*fJac[i][1];
	                fJac[i][3] = fJac[i][2]/v[i];
	            }
break;

	        case (10):

	            //! MEYER FUNCTION.

	            for(i=0;i< 16;i++){
	                temp = five*dfloat(i) + c45 + x[2];
	                tmp1 = x[1]/temp;
	                tmp2 = Math.exp(tmp1);
	                fJac[i][0] = tmp2;
	                fJac[i][1] = x[0]*tmp2/temp;
	                fJac[i][2] = -tmp1*fJac[i][1];
	            }
break;

	        case (11):

	            //! WATSON FUNCTION.

	            for(i=0;i< 29;i++){
	                div = dfloat(i)/c29;
	                s2 = zero;
	                dx = one;
	                for(j=0;j<n;j++){
	                    s2 = s2 + dx*x[j];
	                    dx = div*dx;
	                }
	                temp = two*div*s2;
	                dx = one/div;
	                for(j=0;j<n;j++){
	                    fJac[i][j] = dx*(dfloat(j - 1) - temp);
	                    dx = div*dx;
	                }
	            }
	            for(j=0;j<n;j++){
	                for(i = 29;i<31;i++){
	                    fJac[i][j] = zero;
	                }
	            }
	            fJac[30][0] = one;
	            fJac[31][1] = -two*x[0];
	            fJac[31][2] = one;
break;

	        case (12):

	            //! BOX 3-DIMENSIONAL FUNCTION.

	            for(i=0;i<m;i++){
	                temp = dfloat(i);
	                tmp1 = temp/ten;
	                fJac[i][0] = -tmp1*Math.exp(-tmp1*x[0]);
	                fJac[i][1] = tmp1*Math.exp(-tmp1*x[1]);
	                fJac[i][2] = Math.exp(-temp) - Math.exp(-tmp1);
	            }
break;

	        case (13):

	            //! JENNRICH AND SAMPSON FUNCTION.

	            for(i=0;i<m;i++){
	                temp = dfloat(i);
	                fJac[i][0] = -temp*Math.exp(temp*x[0]);
	                fJac[i][1] = -temp*Math.exp(temp*x[1]);
	            }
break;

	        case (14):

	            //! BROWN AND DENNIS FUNCTION.

	            for(i=0;i<m;i++){
	                temp = dfloat(i)/five;
	                ti = Math.sin(temp);
	                tmp1 = x[0] + temp*x[1] - Math.exp(temp);
	                tmp2 = x[2] + ti*x[3] - Math.cos(temp);
	                fJac[i][0] = two*tmp1;
	                fJac[i][1] = temp*fJac[i][1];
	                fJac[i][2] = two*tmp2;
	                fJac[i][3] = ti*fJac[i][3];
	            }
break;

	        case (15):

	            //! CHEBYQUAD FUNCTION.

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

	            //! BROWN ALMOST-LINEAR FUNCTION.

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
	                    for(k=0;k<n;k++){
	                        if (k != j) {prod = x[k]*prod;}
	                    }
	                }
	                fJac[n][j] = prod/temp;
	            }
break;

	        case (17):

	            //! OSBORNE 1 FUNCTION.

	            for(i=0;i< 33;i++){
	                temp = ten*dfloat(i - 1);
	                tmp1 = Math.exp(-x[3]*temp);
	                tmp2 = Math.exp(-x[4]*temp);
	                fJac[i][0] = -one;
	                fJac[i][1] = -tmp1;
	                fJac[i][2] = -tmp2;
	                fJac[i][3] = temp*x[1]*tmp1;
	                fJac[i][4] = temp*x[2]*tmp2;
	            }
break;

	        case (18):

	            //! OSBORNE 2 FUNCTION.

	            for(i=0;i<65;i++){
	                temp = dfloat(i - 1)/ten;
	                tmp1 = Math.exp(-x[4]*temp);
	                tmp2 = Math.exp(-x[5]*(temp - x[8])*(temp - x[8]));
	                tmp3 = Math.exp(-x[6]*(temp - x[9])*(temp - x[9]));
	                tmp4 = Math.exp(-x[7]*(temp - x[10])*(temp - x[10]));
	                fJac[i][0] = -tmp1;
	                fJac[i][1] = -tmp2;
	                fJac[i][2] = -tmp3;
	                fJac[i][3] = -tmp4;
	                fJac[i][4] = temp*x[0]*tmp1;
	                fJac[i][5] = x[1]*(temp - x[8])*(temp - x[8])*tmp2;
	                fJac[i][6] = x[2]*(temp - x[9])*(temp - x[9])*tmp3;
	                fJac[i][7] = x[3]*(temp - x[10])*(temp - x[10])*tmp4;
	                fJac[i][8] = -two*x[1]*x[5]*(temp - x[8])*tmp2;
	                fJac[i][9] = -two*x[2]*x[6]*(temp - x[9])*tmp3;
	                fJac[i][10] = -two*x[3]*x[7]*(temp - x[10])*tmp4;
	            }
break;

	        default:

	            //! LINEAR FUNCTION - FULL RANK.

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
//	    end subroutine ssqjac
//	!*****************************************************************************************
//
//	!*****************************************************************************************
//	!>
//	!  This subroutine specifies the standard starting points for the
//	!  functions defined by subroutine [[ssqfcn]]. the subroutine returns
//	!  in x a multiple (factor) of the standard starting point. for
//	!  the 11th function the standard starting point is zero, so in
//	!  this case, if factor is not unity, { the subroutine returns
//	!  the vector  x[j] = factor, j=1,...,n.
		public static void initpt(int n, double[] x, int nProb, double factor){
//	    subroutine initpt(n, x, nProb, factor)
//	        implicit none
//
//	        integer, intent(in) :: n !! a positive integer input variable.
//	        integer, intent(in) :: nProb !! a positive integer input variable which defines the
//	                                  !! number of the problem. nProb must not exceed 18.
//	        real(wp), intent(in) :: factor !! an input variable which specifies the multiple of
//	                                    !! the standard starting point. if factor is unity, no
//	                                    !! multiplication is performed.
//	        real(wp), intent(out) :: x(n) !! an output array of length n which contains the standard
//	                                   !! starting point for problem nProb multiplied by factor.

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
			x = new double[n];
//	        x(1:n) = zero
			for(int i=0;i<n;i++) {
				x[i]=0;
			}

	        //! SELECTION OF INITIAL POINT.

	        switch (nProb){
	        case (4):

	            //! ROSENBROCK FUNCTION.

	            x[0] = -c1;
	            x[1] = one;
break;

	        case (5):

	            //! HELICAL VALLEY FUNCTION.

	            x[0] = -one;
	            x[1] = zero;
	            x[2] = zero;
break;

	        case (6):

	            //! POWELL Math.sinGULAR FUNCTION.

	            x[0] = three;
	            x[1] = -one;
	            x[2] = zero;
	            x[3] = one;
break;

	        case (7):

	            //! FREUDENSTEIN AND ROTH FUNCTION.

	            x[0] = half;
	            x[1] = -two;
break;

	        case (8):

	            //! BARD FUNCTION.

	            x[0] = one;
	            x[1] = one;
	            x[2] = one;
break;

	        case (9):

	            //! KOWALIK AND OSBORNE FUNCTION.

	            x[0] = c2;
	            x[1] = c3;
	            x[2] = c4;
	            x[3] = c3;
break;

	        case (10):

	            //! MEYER FUNCTION.

	            x[0] = c5;
	            x[1] = c6;
	            x[2] = c7;
break;

	        case (11):

	            //! WATSON FUNCTION.

	            for(j=0;j<n;j++){
	                x[j] = zero;
	            }
break;

	        case (12):

	            //! BOX 3-DIMENSIONAL FUNCTION.

	            x[0] = zero;
	            x[1] = ten;
	            x[2] = twenty;
break;

	        case (13):

	            //! JENNRICH AND SAMPSON FUNCTION.

	            x[0] = c8;
	            x[1] = c9;
break;

	        case (14):

	            //! BROWN AND DENNIS FUNCTION.

	            x[0] = twntf;
	            x[1] = five;
	            x[2] = -five;
	            x[3] = -one;
break;

	        case (15):

	            //! CHEBYQUAD FUNCTION.

	            h = one/dfloat(n + 1);
	            for(j=0;j<n;j++){
	                x[j] = dfloat(j)*h;
	            }
break;

	        case (16):

	            //! BROWN ALMOST-LINEAR FUNCTION.

	            for(j=0;j<n;j++){
	                x[j] = half;
	            }
break;

	        case (17):

	            //! OSBORNE 1 FUNCTION.

	            x[0] = half;
	            x[1] = c10;
	            x[2] = -one;
	            x[3] = c11;
	            x[4] = c5;
break;

	        case (18):

	            //! OSBORNE 2 FUNCTION.

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

	            //! LINEAR FUNCTION - FULL RANK OR RANK 1.

	            for(j=0;j<n;j++){
	                x[j] = one;
	            }
break;

	        }

	        //! COMPUTE MULTIPLE OF INITIAL POINT.

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
//	    end subroutine initpt
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

//	        integer,intent(in) :: m !! positive integer input variable. n must not exceed m.
//	        integer,intent(in) :: n !! positive integer input variable. n must not exceed m.
//	        integer,intent(in) :: nProb !! a positive integer input variable which defines the
//	                                    !! number of the problem. nProb must not exceed 18.
//	        real(wp),intent(in) :: x(n) !! an input array of length n.
//	        real(wp),intent(out) :: fVec[m] !! an output array of length m which contains the nProb
//	                                        !! function evaluated at x.

	        double zero  = 0.0;
	        double zp25  = 0.25;
	        double zp5   = 0.5;
	        double one   = 1.0;
	        double two   = 2.0;
	        double five  = 5.0;
	        double eight = 8.0;
	        double ten   = 10.0;
	        double c13   = 13.0;
	        double c14   = 14.0;
	        double c29   = 29.0;
	        double c45   = 45.0;

	        int i, iev, ivar, j, nm1;
	       double div, dx, prod, sum, s1, s2, temp, ti, 
	                    tmp1, tmp2, tmp3, tmp4, tpi;

	        double[] v =  new double[] { 
	             4.0, 2.0, 1.0, 5.0e-1, 2.5e-1, 1.67e-1, 
	             1.25e-1, 1.0e-1, 8.33e-2, 7.14e-2, 6.25e-2};
	        double[] y1 =  new double[] { 
	             1.4e-1, 1.8e-1, 2.2e-1, 2.5e-1, 2.9e-1, 3.2e-1, 
	             3.5e-1, 3.9e-1, 3.7e-1, 5.8e-1, 7.3e-1, 9.6e-1, 
	             1.34, 2.1, 4.39};
	        double[] y2 =  new double[] {1.957e-1, 1.947e-1, 
	             1.735e-1, 1.6e-1, 8.44e-2, 6.27e-2, 4.56e-2, 3.42e-2,  
	             3.23e-2, 2.35e-2, 2.46e-2 };
	        double[] y3 =  new double[] { 
	             3.478e4, 2.861e4, 2.365e4, 1.963e4, 
	             1.637e4, 1.372e4, 1.154e4, 9.744e3, 8.261e3, 7.03e3,   
	             6.005e3, 5.147e3, 4.427e3, 3.82e3, 3.307e3, 2.872e3 };
	        double[] y4 =  new double[] { 8.44e-1,
	             9.08e-1, 9.32e-1, 9.36e-1, 9.25e-1, 9.08e-1, 8.81e-1,  
	             8.5e-1, 8.18e-1, 7.84e-1, 7.51e-1, 7.18e-1, 6.85e-1,   
	             6.58e-1, 6.28e-1, 6.03e-1, 5.8e-1, 5.58e-1, 5.38e-1,   
	             5.22e-1, 5.06e-1, 4.9e-1, 4.78e-1, 4.67e-1, 4.57e-1,   
	             4.48e-1, 4.38e-1, 4.31e-1, 4.24e-1, 4.2e-1, 4.14e-1,   
	             4.11e-1, 4.06e-1 };
	        double[] y5 =  new double[] { 1.366, 
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
	             2.92e-1, 1.62e-1, 9.8e-2, 5.4e-2 };

//	        ! initialize:
//	        fVec(1:m) = zero
	        fVec = new double[m];
	        for(i = 0;i<m;i++) {
	        	fVec[i] = 0;
	        }

//	        ! FUNCTION ROUTINE SELECTOR.

	        switch (nProb){
	        case (2):

	            //! LINEAR FUNCTION - RANK 1.

	            sum = zero;
	            for(j=0;j<n;j++){
	                sum = sum + dfloat(j)*x[j];
	            }
	            for(i=0;i<m;i++){
	                fVec[i] = dfloat(i)*sum - one;
	            };
break;
	        case (3):

	            //! LINEAR FUNCTION - RANK 1 WITH ZERO COLUMNS AND ROWS.

	            sum = zero;
	            nm1 = n - 1;
	            if (nm1 >= 2) {
	                for( j = 1;j< nm1;j++){
	                    sum = sum + dfloat(j)*x[j];
	                }
	            }
	            for(i=0;i<m;i++){
	                fVec[i] = dfloat(i - 1)*sum - one;
	            }
	            fVec[m] = -one;
break;
	        case (4):

	            //! ROSENBROCK FUNCTION.

	            fVec[1] = ten*(x[1] - x[0]*x[0]);
	            fVec[2] = one - x[0];
break;
	        case (5):

	            //! HELICAL VALLEY FUNCTION.

	            tpi = eight*Math.atan(one);
	            tmp1 = Math.signum(x[1]) *Math.abs(zp25); //(zp25, x[1])
	            if (x[0] > zero) {tmp1 = Math.atan(x[1]/x[0])/tpi;}
	            if (x[0] < zero) {tmp1 = Math.atan(x[1]/x[0])/tpi + zp5;}
	            tmp2 = Math.sqrt(x[0]*x[0] + x[1]*x[1]);
	            fVec[1] = ten*(x[2] - ten*tmp1);
	            fVec[2] = ten*(tmp2 - one);
	            fVec[3] = x[2];
break;
	        case (6):

	            //! POWELL Math.sinGULAR FUNCTION.

	            fVec[0] = x[0] + ten*x[1];
	            fVec[1] = Math.sqrt(five)*(x[2] - x[3]);
	            fVec[2] = (x[1] - two*x[2])*(x[1] - two*x[2]);
	            fVec[3] = Math.sqrt(ten)*(x[0] - x[3])*(x[0] - x[3]);
break;
	        case (7):

	            //! FREUDENSTEIN AND ROTH FUNCTION.

	            fVec[0] = -c13 + x[0] + ((five - x[1])*x[1] - two)*x[1];
	            fVec[1] = -c29 + x[0] + ((one + x[1])*x[1] - c14)*x[1];
break;
	        case (8):

	            //! BARD FUNCTION.

	            for(i=0;i< 15;i++){
	                tmp1 = dfloat(i);
	                tmp2 = dfloat(16 - i);
	                tmp3 = tmp1;
	                if (i > 8) tmp3 = tmp2;
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

	            //! MEYER FUNCTION.

	            for(i=0;i< 16;i++){
	                temp = five*dfloat(i) + c45 + x[2];
	                tmp1 = x[1]/temp;
	                tmp2 = Math.exp(tmp1);
	                fVec[i] = x[0]*tmp2 - y3[i];
	            }
break;
	        case (11):

	            //! WATSON FUNCTION.

	            for(i=0;i< 29;i++){
	                div = dfloat(i)/c29;
	                s1 = zero;
	                dx = one;
	                for(j = 1;j< n;j++){
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

	            //! BOX 3-DIMENSIONAL FUNCTION.

	            for(i=0;i<m;i++){
	                temp = dfloat(i);
	                tmp1 = temp/ten;
	                fVec[i] = Math.exp(-tmp1*x[0]) - Math.exp(-tmp1*x[1]) + (Math.exp(-temp) - Math.exp(-tmp1))*x[2];
	            }
break;
	        case (13):

	            // JENNRICH AND SAMPSON FUNCTION.

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

	            // CHEBYQUAD FUNCTION.

	            for(i=0;i<m;i++){
	                fVec[i] = zero;
	            }
break;
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
	                if (iev > 0) fVec[i] = fVec[i] + one/(dfloat(i)*dfloat(i) - one);
	                iev = -iev;
	            }
break;
	        case (16):

	            //! BROWN ALMOST-LINEAR FUNCTION.

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

	            // OSBORNE 1 FUNCTION.

	            for(i=0;i< 33;i++){
	                temp = ten*dfloat(i - 1);
	                tmp1 = Math.exp(-x[3]*temp);
	                tmp2 = Math.exp(-x[4]*temp);
	                fVec[i] = y4[i] - (x[0] + x[1]*tmp1 + x[2]*tmp2);
	            }
break;
	        case (18):

	            //! OSBORNE 2 FUNCTION.

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

	            //! LINEAR FUNCTION - FULL RANK.

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
//
//	    end subroutine ssqfcn
//	!*****************************************************************************************
//
//	!*****************************************************************************************
//	    end program test_lmder
//	!*****************************************************************************************
}
