package Math.minpackTestPacks;
//
public class hybrjTest {
//	!*****************************************************************************************
//	!>
//	!  This program tests codes for the solution of n nonlinear
//	!  equations in n variables. it consists of a driver and an
//	!  interface subroutine fcn. the driver reads in data, calls the
//	!  nonlinear equation solver, and finally prints out information
//	!  on the performance of the solver. this is only a sample driver,
//	!  many other drivers are possible. the interface subroutine fcn
//	!  is necessary to take into account the forms of calling
//	!  sequences used by the function and jacobian subroutines in
//	!  the various nonlinear equation solvers.
//
//	program test_hybrj
//
//	    use Math.minpack_module, only: wp, dpmpar, enorm, hybrj1
//	    use iso_fortran_env, only: nwrite => output_unit
//
//	    implicit none
//
//	    ! originally from file21
//	    integer,parameter :: ncases = 22
//	    integer,dimension(ncases),parameter :: nProbs  = [1,2,3,4,5,6,6,7,7,7,7,7,8,8,8,9,10,10,11,12,13,14]
//	    integer,dimension(ncases),parameter :: ns      = [2,4,2,4,3,6,9,5,6,7,8,9,10,30,40,10,1,10,10,10,10,10]
//	    integer,dimension(ncases),parameter :: ntriess = [3,3,2,3,3,2,2,3,3,3,1,1,3,1,1,3,3,3,3,3,3,3]
//
//	    integer,dimension(55),parameter :: info_original = [1,1,1,4,4,4,1,1,1,1,1,1,1,1,1,1,1,4,1,1,1,1,1,
//	                                                        1,1,1,4,4,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,4,1,1,
//	                                                        1,1,1,1,1,1,1,1,1]
//
//	    int i, ic, info, k, n, NFEv, NJEv, nProb, ntries, icase, lwa, ldfjac
//	    int na(55), nf(55), nj(55), np(55), nx(55)
//	    double fnm(55)
//	    double factor, fnorm1, fnorm2
//	    real(wp),allocatable :: fjac(:,:), fVec(:), wa(:), x(:)
//
//	    real(wp), parameter :: one = 1.0
//	    real(wp), parameter :: ten = 10.0
//	    real(wp), parameter :: tol = Math.sqrt(dpmpar(1))
//	    real(wp), parameter :: solution_reltol = 1.0e-4 !! reltol for matching previously generated solutions
//
//	    ic = 0
//
//	    do icase = 1, ncases+1
//	        if (icase == ncases+1) then
//	            write (nwrite, '(A,I3,A/)') '1SUMMARY OF ', ic, ' CALLS TO HYBRJ1'
//	            write (nwrite, '(A/)')      ' nProb   N    NFEV   NJEV  INFO  FINAL L2 NORM'
//	            for(i=0;i<ic;i++){
//	                write (nwrite, '(I4,I6,2I7,I6,1X,D15.7)') np(i), na(i), nf(i), nj(i), nx(i), fnm(i)
//	            }
//	            stop
//	        else
//	            nProb = nProbs(icase)
//	            n = ns(icase)
//	            ldfjac = n
//	            lwa = (n*(n+13))/2
//	            ntries = ntriess(icase)
//
//	            if (allocated(fjac)) deallocate(fjac)
//	            if (allocated(fVec)) deallocate(fVec)
//	            if (allocated(wa)) deallocate(wa)
//	            if (allocated(x)) deallocate(x)
//	            allocate(fjac(n,n))
//	            allocate(fVec(n))
//	            allocate(wa(lwa))
//	            allocate(x(n))
//
//	            factor = one
//	            for(k=0;k<n;k++){tries
//	                ic = ic + 1
//	                call initpt(n, x, nProb, factor)
//	                call vecfcn(n, x, fVec, nProb)
//	                fnorm1 = enorm(n, fVec)
//	                write (nwrite, '(////5X,A,I5,5X,A,I5,5X//)') ' PROBLEM', nProb, ' DIMENSION', n
//	                NFEv = 0
//	                NJEv = 0
//	                call hybrj1(fcn, n, x, fVec, fjac, ldfjac, tol, info, wa, lwa)
//	                fnorm2 = enorm(n, fVec)
//	                np(ic) = nProb
//	                na(ic) = n
//	                nf(ic) = NFEv
//	                nj(ic) = NJEv
//	                nx(ic) = info
//	                fnm(ic) = fnorm2
//	                write (nwrite, '(5X,A,D15.7//5X,A,D15.7//5X,A,I10//5X,A,I10//5X,A,18X,I10//5X,A//*(5X,5D15.7/))') 
//	                               ' INITIAL L2 NORM OF THE RESIDUALS', fnorm1, 
//	                               ' FINAL L2 NORM OF THE RESIDUALS  ', fnorm2, 
//	                               ' NUMBER OF FUNCTION EVALUATIONS  ', NFEv,   
//	                               ' NUMBER OF JACOBIAN EVALUATIONS  ', NJEv,   
//	                               ' EXIT PARAMETER', info, 
//	                               ' FINAL APPROXIMATE SOLUTION', x(1:n)
//	                factor = ten*factor
//	                call compare_solutions(ic, x, solution_reltol, tol)
//
//	            }
//	        end if
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
//	    if (info_original(ic)<5) then    ! ignore any where the original Math.minpack failed
//	        diff = solution(ic) - x
//	        absdiff = abs(diff)
//	        if (any(absdiff>abstol)) then ! first do an absolute diff
//	            ! also do a rel diff if the abs diff fails (also protect for divide by zero)
//	            reldiff = absdiff
//	            where (solution(ic) != 0.0) reldiff = absdiff / abs(solution(ic))
//	            if (any(reldiff > reltol)) then
//	                write(nwrite,'(A)') 'Failed case'
//	                write(nwrite, '(//5x, a//(5x, 5d15.7))') 'Expected x: ', solution(ic)
//	                write(nwrite, '(/5x, a//(5x, 5d15.7))')  'Computed x: ', x
//	                write(nwrite, '(/5x, a//(5x, 5d15.7))')  'absdiff: ', absdiff
//	                write(nwrite, '(/5x, a//(5x, 5d15.7))')  'reldiff: ', reldiff
//	                error stop ! test failed
//	            end if
//	        end if
//	    end if
//
//	    }  // compare_solutions
//	!*****************************************************************************************
//
//	!*****************************************************************************************
//	!>
//	!  The calling sequence of fcn should be identical to the
//	!  calling sequence of the function subroutine in the nonlinear
//	!  equation solver. fcn should only call the testing function
//	!  and jacobian subroutines vecfcn and vecjac with the
//	!  appropriate value of problem number (nProb).
//
//	    subroutine fcn(n, x, fVec, Fjac, Ldfjac, Iflag)
//	        implicit none
//
//	        integer,intent(in) :: n
//	        integer,intent(in) :: Ldfjac
//	        integer,intent(inout) :: Iflag
//	        real(wp),intent(in) :: x(n)
//	        real(wp),intent(inout) :: fVec(n)
//	        real(wp),intent(inout) :: Fjac(Ldfjac, n)
//
//	        switch (iflag)
//	        case (1)
//	            call vecfcn(n, x, fVec, nProb)
//	            NFEv = NFEv + 1
//	        case(2)
//	            call vecjac(n, x, Fjac, Ldfjac, nProb)
//	            NJEv = NJEv + 1
//	        default:
//	            error stop 'invalid iflag value'
//	        } 
//
//	    }  // fcn
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
//
//	    pure function solution(nProb) result(x)
//
//	        implicit none
//
//	        integer,intent(in) :: nProb
//	        real(wp),dimension(:),allocatable :: x
//
//	        switch (nProb)
//	        case (1); x = [0.1000000000000000E+01,0.1000000000000000E+01]
//	        case (2); x = [0.1000000000000000E+01,0.9999999999999971E+00]
//	        case (3); x = [0.1000000000000000E+01,0.1000000000000000E+01]
//	        case (4); x = [0.2695066833053901E-17,-0.2695066833053901E-18,0.3749757049013966E-17,
//	                       0.3749757049013966E-17]
//	        case (5); x = [0.2961047510523875E-17,-0.2961047510523875E-18,0.3392094461844849E-17,
//	                       0.3392094461844849E-17]
//	        case (6); x = [-0.1125674354180310E-16,0.1125674354180310E-17,-0.6482313483583572E-17,
//	                       -0.6482313483583572E-17]
//	        case (7); x = [0.1098159327798296E-04,0.9106146740038449E+01]
//	        case (8); x = [0.1098159288132353E-04,0.9106146743612113E+01]
//	        case (9); x = [-0.9679740249512818E+00,0.9471391408444433E+00,-0.9695163103175638E+00,
//	                       0.9512476657652131E+00]
//	        case (10); x = [-0.9679740250937678E+00,0.9471391411205137E+00,-0.9695163101751340E+00,
//	                       0.9512476654887939E+00]
//	        case (11); x = [-0.9679740249230577E+00,0.9471391407896745E+00,-0.9695163103461089E+00,
//	                       0.9512476658204914E+00]
//	        case (12); x = [0.1000000000000010E+01,-0.1612103288815217E-13,0.0000000000000000E+00]
//	        case (13); x = [0.1000000000004274E+01,0.1388617146084722E-10,0.0000000000000000E+00]
//	        case (14); x = [0.1000000000050375E+01,0.3234657151609015E-10,0.0000000000000000E+00]
//	        case (15); x = [-0.1572508640131874E-01,0.1012434869369120E+01,-0.2329916259568086E+00,
//	                        0.1260430087800509E+01,-0.1513728922723659E+01,0.9929964324319899E+00]
//	        case (16); x = [-0.1572508640142940E-01,0.1012434869369113E+01,-0.2329916259567570E+00,
//	                        0.1260430087799810E+01,-0.1513728922722583E+01,0.9929964324313200E+00]
//	        case (17); x = [-0.1530703902928690E-04,0.9997897039319343E+00,0.1476396368990648E-01,
//	                        0.1463423283093984E+00,0.1000821102959179E+01,-0.2617731140413794E+01,
//	                        0.4104403164336241E+01,-0.3143612278455328E+01,0.1052626407979455E+01]
//	        case (18); x = [0.3119911862752667E+00,0.1089342426782224E+01,0.1379293238883866E+01,
//	                        -0.1175426477217227E+02,0.6265967591560398E+02,-0.1636004875532691E+03,
//	                        0.2326239209287591E+03,-0.1697823090157466E+03,0.5076761702794597E+02]
//	        case (19); x = [0.8375125649943561E-01,0.3127292952232932E+00,0.5000000000000018E+00,
//	                        0.6872707047767043E+00,0.9162487435005652E+00]
//	        case (20); x = [0.8375125650172599E-01,0.4999999999984706E+00,0.3127292952187522E+00,
//	                        0.9162487434957841E+00,0.6872707047852671E+00]
//	        case (21); x = [0.6872707048164862E+00,0.4999999999708427E+00,0.8375125649245271E-01,
//	                        0.3127292952359536E+00,0.9162487434842648E+00]
//	        case (22); x = [0.6687659094732229E-01,0.3666822992416477E+00,0.2887406731168656E+00,
//	                        0.7112593268831344E+00,0.6333177007583523E+00,0.9331234090526778E+00]
//	        case (23); x = [0.9331234090548961E+00,0.3666822992072173E+00,0.6687659093602476E-01,
//	                        0.7112593268461027E+00,0.2887406731541705E+00,0.6333177008015886E+00]
//	        case (24); x = [0.6687659094738703E-01,0.7112593269091415E+00,0.2887406731888937E+00,
//	                        0.9331234090619587E+00,0.6333177007462377E+00,0.3666822991463812E+00]
//	        case (25); x = [0.5806914960717641E-01,0.2351716124057688E+00,0.3380440947078516E+00,
//	                        0.4999999999993442E+00,0.6619559052930899E+00,0.7648283875938304E+00,
//	                        0.9419308503929387E+00]
//	        case (26); x = [0.3380440947234102E+00,0.2351716123658910E+00,0.9419308503765472E+00,
//	                        0.7648283876531939E+00,0.6619559052401212E+00,0.5806914961691421E-01,
//	                        0.5000000000239224E+00]
//	        case (27); x = [-0.4649136147100040E+02,-0.1039521377302678E+02,-0.7517577291666804E+01,
//	                        0.1217934366521229E+02,0.4297678475778002E+02,0.2373113109255284E+02,
//	                        0.4306103920017470E+02]
//	        case (28); x = [0.4985639998615101E-01,0.1986351265136680E+00,0.2698288232492614E+00,
//	                        0.4992722954902931E+00,0.5007277045097069E+00,0.7301711767507386E+00,
//	                        0.8013648734863320E+00,0.9501436000138489E+00]
//	        case (29); x = [0.4420534627964380E-01,0.1994906715249946E+00,0.2356191097165776E+00,
//	                        0.4160469066798081E+00,0.5000000008711383E+00,0.5839530917260407E+00,
//	                        0.7643808911354805E+00,0.8005093282727184E+00,0.9557946537935981E+00]
//	        case (30); x = [0.1000000000000006E+01,0.1000000000000006E+01,0.1000000000000006E+01,
//	                        0.1000000000000006E+01,0.1000000000000006E+01,0.1000000000000006E+01,
//	                        0.1000000000000006E+01,0.1000000000000006E+01,0.1000000000000006E+01,
//	                        0.9999999999999419E+00]
//	        case (31); x = [0.1000000000000000E+01,0.1000000000000000E+01,0.1000000000000000E+01,
//	                        0.1000000000000000E+01,0.1000000000000000E+01,0.1000000000000000E+01,
//	                        0.1000000000000000E+01,0.1000000000000000E+01,0.1000000000000000E+01,
//	                        0.9999999999999953E+00]
//	        case (32); x = [0.9794303033498583E+00,0.9794303033498583E+00,0.9794303033498583E+00,
//	                        0.9794303033498583E+00,0.9794303033498583E+00,0.9794303033498583E+00,
//	                        0.9794303033498583E+00,0.9794303033498583E+00,0.9794303033498583E+00,
//	                        0.1205696966501417E+01]
//	        case (33); x = [0.9999999999999868E+00,0.9999999999999868E+00,0.9999999999999868E+00,
//	                        0.9999999999999868E+00,0.9999999999999868E+00,0.9999999999999868E+00,
//	                        0.9999999999999868E+00,0.9999999999999868E+00,0.9999999999999868E+00,
//	                        0.9999999999999868E+00,0.9999999999999868E+00,0.9999999999999868E+00,
//	                        0.9999999999999868E+00,0.9999999999999868E+00,0.9999999999999868E+00,
//	                        0.9999999999999868E+00,0.9999999999999868E+00,0.9999999999999868E+00,
//	                        0.9999999999999868E+00,0.9999999999999868E+00,0.9999999999999868E+00,
//	                        0.9999999999999868E+00,0.9999999999999868E+00,0.9999999999999868E+00,
//	                        0.9999999999999868E+00,0.9999999999999868E+00,0.9999999999999868E+00,
//	                        0.9999999999999868E+00,0.9999999999999868E+00,0.1000000000000383E+01]
//	        case (34); x = [0.1000000000000006E+01,0.1000000000000006E+01,0.1000000000000006E+01,
//	                        0.1000000000000006E+01,0.1000000000000006E+01,0.1000000000000006E+01,
//	                        0.1000000000000006E+01,0.1000000000000006E+01,0.1000000000000006E+01,
//	                        0.1000000000000006E+01,0.1000000000000006E+01,0.1000000000000006E+01,
//	                        0.1000000000000006E+01,0.1000000000000006E+01,0.1000000000000006E+01,
//	                        0.1000000000000006E+01,0.1000000000000006E+01,0.1000000000000006E+01,
//	                        0.1000000000000006E+01,0.1000000000000006E+01,0.1000000000000006E+01,
//	                        0.1000000000000006E+01,0.1000000000000006E+01,0.1000000000000006E+01,
//	                        0.1000000000000006E+01,0.1000000000000006E+01,0.1000000000000006E+01,
//	                        0.1000000000000006E+01,0.1000000000000006E+01,0.1000000000000006E+01,
//	                        0.1000000000000006E+01,0.1000000000000006E+01,0.1000000000000006E+01,
//	                        0.1000000000000006E+01,0.1000000000000006E+01,0.1000000000000006E+01,
//	                        0.1000000000000006E+01,0.1000000000000006E+01,0.1000000000000006E+01,
//	                        0.9999999999997405E+00]
//	        case (35); x = [-0.4316498251876486E-01,-0.8157715653538729E-01,-0.1144857143805310E+00,
//	                        -0.1409735768625996E+00,-0.1599086961819857E+00,-0.1698772023127759E+00,
//	                        -0.1690899837812081E+00,-0.1552495352218312E+00,-0.1253558916789345E+00,
//	                        -0.7541653368589182E-01]
//	        case (36); x = [-0.4316498251881878E-01,-0.8157715653546950E-01,-0.1144857143805966E+00,
//	                        -0.1409735768626191E+00,-0.1599086961819499E+00,-0.1698772023126901E+00,
//	                        -0.1690899837811062E+00,-0.1552495352217907E+00,-0.1253558916789970E+00,
//	                        -0.7541653368596339E-01]
//	        case (37); x = [-0.4316498254522300E-01,-0.8157715658124411E-01,-0.1144857144024344E+00,
//	                        -0.1409735768722910E+00,-0.1599086963002226E+00,-0.1698772022538783E+00,
//	                        -0.1690899837944877E+00,-0.1552495352060589E+00,-0.1253558916432355E+00,
//	                        -0.7541653366610668E-01]
//	        case (38); x = [-0.1528138835625800E+00]
//	        case (39); x = [-0.1528138835625801E+00]
//	        case (40); x = [-0.1528138835625800E+00]
//	        case (41); x = [-0.4316498251876487E-01,-0.8157715653538729E-01,-0.1144857143805310E+00,
//	                        -0.1409735768625996E+00,-0.1599086961819857E+00,-0.1698772023127759E+00,
//	                        -0.1690899837812080E+00,-0.1552495352218311E+00,-0.1253558916789344E+00,
//	                        -0.7541653368589175E-01]
//	        case (42); x = [-0.4316498251881876E-01,-0.8157715653546944E-01,-0.1144857143805966E+00,
//	                        -0.1409735768626190E+00,-0.1599086961819498E+00,-0.1698772023126901E+00,
//	                        -0.1690899837811062E+00,-0.1552495352217907E+00,-0.1253558916789970E+00,
//	                        -0.7541653368596334E-01]
//	        case (43); x = [-0.4316498251876519E-01,-0.8157715653538752E-01,-0.1144857143805303E+00,
//	                        -0.1409735768625981E+00,-0.1599086961819844E+00,-0.1698772023127748E+00,
//	                        -0.1690899837812073E+00,-0.1552495352218307E+00,-0.1253558916789341E+00,
//	                        -0.7541653368589157E-01]
//	        case (44); x = [0.5526154410956152E-01,0.5695755464278050E-01,0.5889066688359926E-01,
//	                        0.6113606294453851E-01,0.6377855838223215E-01,0.6700482468993153E-01,
//	                        0.2079410281099167E+00,0.1642671068222599E+00,0.8643947917821507E-01,
//	                        0.9133506311839429E-01]
//	        case (45); x = [0.3439628896235700E-01,0.3503231575416286E-01,0.3571919583574922E-01,
//	                        0.3646522422002401E-01,0.3728091174083832E-01,0.3817986258974627E-01,
//	                        0.3918014109818273E-01,0.4030650261421058E-01,0.1797201916815176E+00,
//	                        0.1562408814749914E+00]
//	        case (46); x = [0.1888395221037194E+02,0.2516777354435942E+02,0.1888527511766734E+02,
//	                        0.1888602114534231E+02,0.1888683683396955E+02,0.1888773578345518E+02,
//	                        0.1888873606137402E+02,0.1888986242382218E+02,0.1902927611206524E+02,
//	                        0.1900579680367850E+02]
//	        case (47); x = [0.9999999999999993E+00,0.9999999999999986E+00,0.9999999999999979E+00,
//	                        0.9999999999999972E+00,0.9999999999999964E+00,0.9999999999999958E+00,
//	                        0.9999999999999950E+00,0.9999999999999943E+00,0.9999999999999937E+00,
//	                        0.9999999999999929E+00]
//	        case (48); x = [0.1000000000000017E+01,0.1000000000000033E+01,0.1000000000000050E+01,
//	                        0.1000000000000066E+01,0.1000000000000083E+01,0.1000000000000100E+01,
//	                        0.1000000000000116E+01,0.1000000000000133E+01,0.1000000000000150E+01,
//	                        0.1000000000000166E+01]
//	        case (49); x = [0.1000000000000000E+01,0.1000000000000000E+01,0.1000000000000000E+01,
//	                        0.1000000000000000E+01,0.1000000000000000E+01,0.1000000000000000E+01,
//	                        0.1000000000000000E+01,0.1000000000000000E+01,0.1000000000000000E+01,
//	                        0.1000000000000000E+01]
//	        case (50); x = [-0.5707221307212121E+00,-0.6818069509055232E+00,-0.7022100775689857E+00,
//	                        -0.7055106309936168E+00,-0.7049061557572888E+00,-0.7014966060124587E+00,
//	                        -0.6918893211477919E+00,-0.6657965141985400E+00,-0.5960351099566767E+00,
//	                        -0.4164122574358191E+00]
//	        case (51); x = [-0.5707221320993724E+00,-0.6818069494976288E+00,-0.7022100764193087E+00,
//	                        -0.7055106298465235E+00,-0.7049061556842444E+00,-0.7014966070281702E+00,
//	                        -0.6918893223746194E+00,-0.6657965143749420E+00,-0.5960351092084393E+00,
//	                        -0.4164122574088334E+00]
//	        case (52); x = [-0.5707221320171143E+00,-0.6818069499829605E+00,-0.7022100760171542E+00,
//	                        -0.7055106298955310E+00,-0.7049061557301967E+00,-0.7014966070327222E+00,
//	                        -0.6918893223590803E+00,-0.6657965144072679E+00,-0.5960351090088830E+00,
//	                        -0.4164122575177334E+00]
//	        case (53); x = [-0.4283028636053096E+00,-0.4765964242962532E+00,-0.5196524638125551E+00,
//	                        -0.5580993246169653E+00,-0.5925061569509360E+00,-0.6245036821428090E+00,
//	                        -0.6232394714478015E+00,-0.6213938418388717E+00,-0.6204535966122983E+00,
//	                        -0.5864692707477790E+00]
//	        case (54); x = [-0.4283028634881692E+00,-0.4765964236396711E+00,-0.5196524642776768E+00,
//	                        -0.5580993248351936E+00,-0.5925061568131795E+00,-0.6245036817962691E+00,
//	                        -0.6232394720687791E+00,-0.6213938417874499E+00,-0.6204535965224117E+00,
//	                        -0.5864692707287930E+00]
//	        case (55); x = [-0.4283028635608067E+00,-0.4765964243232715E+00,-0.5196524637037395E+00,
//	                        -0.5580993248328234E+00,-0.5925061568292707E+00,-0.6245036822076749E+00,
//	                        -0.6232394714256790E+00,-0.6213938418143938E+00,-0.6204535966527650E+00,
//	                        -0.5864692707189498E+00]
//	        default:
//	            error stop 'invalid case'
//	        } 
//
//	    end function solution
//	!*****************************************************************************************
//
//	!*****************************************************************************************
//	!>
//	!  This subroutine defines the jacobian matrices of fourteen
//	!  test functions. the problem dimensions are as described
//	!  in the prologue comments of vecfcn.
//
//	    subroutine vecjac(n, x, Fjac, Ldfjac, nProb)
//	        implicit none
//
//	        integer,intent(in) :: n !! a positive integer variable.
//	        integer,intent(in) :: Ldfjac !! a positive integer variable not less than n
//	                                     !! which specifies the leading dimension of the array fjac.
//	        integer,intent(in) :: nProb !! a positive integer variable which defines the
//	                                    !! number of the problem. nProb must not exceed 14.
//	        real(wp),intent(in) :: x(n) !! an array of length n.
//	        real(wp),intent(out) :: Fjac(Ldfjac, n) !! an n by n array. on output fjac contains the
//	                                                !! jacobian matrix of the nProb function evaluated at x.
//
//	        real(wp), parameter :: zero = 0.0
//	        real(wp), parameter :: one = 1.0
//	        real(wp), parameter :: two = 2.0
//	        real(wp), parameter :: three = 3.0
//	        real(wp), parameter :: four = 4.0
//	        real(wp), parameter :: five = 5.0
//	        real(wp), parameter :: six = 6.0
//	        real(wp), parameter :: eight = 8.0
//	        real(wp), parameter :: ten = 10.0
//	        real(wp), parameter :: fiftn = 15.0
//	        real(wp), parameter :: twenty = 20.0
//	        real(wp), parameter :: hundrd = 100.0
//	        real(wp), parameter :: c1 = 1.0e4
//	        real(wp), parameter :: c3 = 2.0e2
//	        real(wp), parameter :: c4 = 2.02e1
//	        real(wp), parameter :: c5 = 1.98e1
//	        real(wp), parameter :: c6 = 1.8e2
//	        real(wp), parameter :: c9 = 2.9e1
//
//	        int i, j, k, k1, k2, ml, mu
//	        double h, prod, sum, sum1, sum2, temp, temp1, temp2, 
//	                    temp3, temp4, ti, tj, tk, tpi
//
//	        Fjac(1:n,1:n) = zero
//
//	        ! jacobian routine selector.
//
//	        switch (nProb)
//	        case (2)
//	            ! powell Math.singular function.
//	            do k = 1, 4
//	                for(j=0;j<4;j++){
//	                    fjac(k, j) = zero
//	                }
//	            }
//	            fjac(1, 1) = one
//	            fjac(1, 2) = ten
//	            fjac(2, 3) = Math.sqrt(five)
//	            fjac(2, 4) = -fjac(2, 3)
//	            fjac(3, 2) = two*(x[1] - two*x[2])
//	            fjac(3, 3) = -two*fjac(3, 2)
//	            fjac(4, 1) = two*Math.sqrt(ten)*(x[0] - x[3])
//	            fjac(4, 4) = -fjac(4, 1)
//	        case (3)
//	            ! powell badly scaled function.
//	            fjac(1, 1) = c1*x[1]
//	            fjac(1, 2) = c1*x[0]
//	            fjac(2, 1) = -exp(-x[0])
//	            fjac(2, 2) = -exp(-x[1])
//	        case (4)
//	            ! wood function.
//	            do k = 1, 4
//	                for(j=0;j<4;j++){
//	                    fjac(k, j) = zero
//	                }
//	            }
//	            temp1 = x[1] - three*x[0]**2
//	            temp2 = x[3] - three*x[2]**2
//	            fjac(1, 1) = -c3*temp1 + one
//	            fjac(1, 2) = -c3*x[0]
//	            fjac(2, 1) = -two*c3*x[0]
//	            fjac(2, 2) = c3 + c4
//	            fjac(2, 4) = c5
//	            fjac(3, 3) = -c6*temp2 + one
//	            fjac(3, 4) = -c6*x[2]
//	            fjac(4, 2) = c5
//	            fjac(4, 3) = -two*c6*x[2]
//	            fjac(4, 4) = c6 + c4
//	        case (5)
//	            ! helical valley function.
//	            tpi = eight*Math.atan(one)
//	            temp = x[0]**2 + x[1]**2
//	            temp1 = tpi*temp
//	            temp2 = Math.sqrt(temp)
//	            fjac(1, 1) = hundrd*x[1]/temp1
//	            fjac(1, 2) = -hundrd*x[0]/temp1
//	            fjac(1, 3) = ten
//	            fjac(2, 1) = ten*x[0]/temp2
//	            fjac(2, 2) = ten*x[1]/temp2
//	            fjac(2, 3) = zero
//	            fjac(3, 1) = zero
//	            fjac(3, 2) = zero
//	            fjac(3, 3) = one
//	        case (6)
//	            ! watson function.
//	            for(k=0;k<n;k++){
//	                do j = k, n
//	                    fjac(k, j) = zero
//	                }
//	            }
//	            for(i=0;i<29;i++){
//	                ti = dfloat(i)/c9
//	                sum1 = zero
//	                temp = one
//	                do j = 2, n
//	                    sum1 = sum1 + dfloat(j - 1)*temp*x[j]
//	                    temp = ti*temp
//	                }
//	                sum2 = zero
//	                temp = one
//	                for(j=0;j<n;j++){
//	                    sum2 = sum2 + temp*x[j]
//	                    temp = ti*temp
//	                }
//	                temp1 = two*(sum1 - sum2**2 - one)
//	                temp2 = two*sum2
//	                temp = ti**2
//	                tk = one
//	                for(k=0;k<n;k++){
//	                    tj = tk
//	                    do j = k, n
//	                        fjac(k, j) = fjac(k, j) 
//	                                     + tj*((dfloat(k - 1)/ti - temp2)*(dfloat(j - 1) 
//	                                     /ti - temp2) - temp1)
//	                        tj = ti*tj
//	                    }
//	                    tk = temp*tk
//	                }
//	            }
//	            fjac(1, 1) = fjac(1, 1) + six*x[0]**2 - two*x[1] + three
//	            fjac(1, 2) = fjac(1, 2) - two*x[0]
//	            fjac(2, 2) = fjac(2, 2) + one
//	            for(k=0;k<n;k++){
//	                do j = k, n
//	                    fjac(j, k) = fjac(k, j)
//	                }
//	            }
//	        case (7)
//	            ! chebyquad function.
//	            tk = one/dfloat(n)
//	            for(j=0;j<n;j++){
//	                temp1 = one
//	                temp2 = two*x[j] - one
//	                temp = two*temp2
//	                temp3 = zero
//	                temp4 = two
//	                for(k=0;k<n;k++){
//	                    fjac(k, j) = tk*temp4
//	                    ti = four*temp2 + temp*temp4 - temp3
//	                    temp3 = temp4
//	                    temp4 = ti
//	                    ti = temp*temp2 - temp1
//	                    temp1 = temp2
//	                    temp2 = ti
//	                }
//	            }
//	        case (8)
//	            ! brown almost-linear function.
//	            prod = one
//	            for(j=0;j<n;j++){
//	                prod = x[j]*prod
//	                for(k=0;k<n;k++){
//	                    fjac(k, j) = one
//	                }
//	                fjac(j, j) = two
//	            }
//	            for(j=0;j<n;j++){
//	                temp = x[j]
//	                if (temp == zero) then
//	                    temp = one
//	                    prod = one
//	                    for(k=0;k<n;k++){
//	                        if (k != j) prod = x[k]*prod
//	                    }
//	                end if
//	                fjac(n, j) = prod/temp
//	            }
//	        case (9)
//	            ! discrete boundary value function.
//	            h = one/dfloat(n + 1)
//	            for(k=0;k<n;k++){
//	                temp = three*(x[k] + dfloat(k)*h + one)**2
//	                for(j=0;j<n;j++){
//	                    fjac(k, j) = zero
//	                }
//	                fjac(k, k) = two + temp*h**2/two
//	                if (k != 1) fjac(k, k - 1) = -one
//	                if (k != n) fjac(k, k + 1) = -one
//	            }
//	        case (10)
//	            ! discrete integral equation function.
//	            h = one/dfloat(n + 1)
//	            for(k=0;k<n;k++){
//	                tk = dfloat(k)*h
//	                for(j=0;j<n;j++){
//	                    tj = dfloat(j)*h
//	                    temp = three*(x[j] + tj + one)**2
//	                    fjac(k, j) = h*Math.min(tj*(one - tk), tk*(one - tj))*temp/two
//	                }
//	                fjac(k, k) = fjac(k, k) + one
//	            }
//	        case (11)
//	            ! trigonometric function.
//	            for(j=0;j<n;j++){
//	                temp = Math.sin(x[j])
//	                for(k=0;k<n;k++){
//	                    fjac(k, j) = temp
//	                }
//	                fjac(j, j) = dfloat(j + 1)*temp - cos(x[j])
//	            }
//	        case (12)
//	            ! variably dimensioned function.
//	            sum = zero
//	            for(j=0;j<n;j++){
//	                sum = sum + dfloat(j)*(x[j] - one)
//	            }
//	            temp = one + six*sum**2
//	            for(k=0;k<n;k++){
//	                do j = k, n
//	                    fjac(k, j) = dfloat(k*j)*temp
//	                    fjac(j, k) = fjac(k, j)
//	                }
//	                fjac(k, k) = fjac(k, k) + one
//	            }
//	        case (13)
//	            ! broyden tridiagonal function.
//	            for(k=0;k<n;k++){
//	                for(j=0;j<n;j++){
//	                    fjac(k, j) = zero
//	                }
//	                fjac(k, k) = three - four*x[k]
//	                if (k != 1) fjac(k, k - 1) = -one
//	                if (k != n) fjac(k, k + 1) = -two
//	            }
//	        case (14)
//	            ! broyden banded function.
//	            ml = 5
//	            mu = 1
//	            for(k=0;k<n;k++){
//	                for(j=0;j<n;j++){
//	                    fjac(k, j) = zero
//	                }
//	                k1 = Math.max(1, k - ml)
//	                k2 = Math.min(k + mu, n)
//	                for(j=k1;j<k2;j++){
//	                    if (j != k) fjac(k, j) = -(one + two*x[j])
//	                }
//	                fjac(k, k) = two + fiftn*x[k]**2
//	            }
//	        default:
//	            ! rosenbrock function.
//	            fjac(1, 1) = -one
//	            fjac(1, 2) = zero
//	            fjac(2, 1) = -twenty*x[0]
//	            fjac(2, 2) = ten
//	        } 
//
//	    }  // vecjac
//	!*****************************************************************************************
//
//	!*****************************************************************************************
//	!>
//	!  This subroutine specifies the standard starting points for
//	!  the functions defined by subroutine vecfcn. the subroutine
//	!  returns in x a multiple (factor) of the standard starting
//	!  point. for the sixth function the standard starting point is
//	!  zero, so in this case, if factor is not unity, then the
//	!  subroutine returns the vector  x[j] = factor, j=1,...,n.
//
//	    subroutine initpt(n, x, nProb, Factor)
//	        implicit none
//
//	        integer,intent(in) :: n !! a positive integer input variable.
//	        integer,intent(in) :: nProb !! a positive integer input variable which defines the
//	                                    !! number of the problem. nProb must not exceed 14.
//	        real(wp),intent(in) :: Factor !! an input variable which specifies the multiple of
//	                                      !! the standard starting point. if factor is unity, no
//	                                      !! multiplication is performed.
//	        real(wp),intent(out) :: x(n) !! an output array of length n which contains the standard
//	                                     !! starting point for problem nProb multiplied by factor.
//
//	        int j
//	        double h, tj
//
//	        double zero = 0.0
//	        double half = 0.5
//	        double one = 1.0
//	        double three = 3.0
//	        double c1 = 1.2
//
//	        x(1:n) = zero
//
//	        ! selection of initial point.
//
//	        switch (nProb)
//	        case (2)
//	            ! powell Math.singular function.
//	            x[0] = three
//	            x[1] = -one
//	            x[2] = zero
//	            x[3] = one
//	        case (3)
//	            ! powell badly scaled function.
//	            x[0] = zero
//	            x[1] = one
//	        case (4)
//	            ! wood function.
//	            x[0] = -three
//	            x[1] = -one
//	            x[2] = -three
//	            x[3] = -one
//	        case (5)
//	            ! helical valley function.
//	            x[0] = -one
//	            x[1] = zero
//	            x[2] = zero
//	        case (6)
//	            ! watson function.
//	            for(j=0;j<n;j++){
//	                x[j] = zero
//	            }
//	        case (7)
//	            ! chebyquad function.
//	            h = one/dfloat(n + 1)
//	            for(j=0;j<n;j++){
//	                x[j] = dfloat(j)*h
//	            }
//	        case (8)
//	            ! brown almost-linear function.
//	            for(j=0;j<n;j++){
//	                x[j] = half
//	            }
//	        case (9, 10)
//	            ! discrete boundary value and integral equation functions.
//	            h = one/dfloat(n + 1)
//	            for(j=0;j<n;j++){
//	                tj = dfloat(j)*h
//	                x[j] = tj*(tj - one)
//	            }
//	        case (11)
//	            ! trigonometric function.
//	            h = one/dfloat(n)
//	            for(j=0;j<n;j++){
//	                x[j] = h
//	            }
//	        case (12)
//	            ! variably dimensioned function.
//	            h = one/dfloat(n)
//	            for(j=0;j<n;j++){
//	                x[j] = one - dfloat(j)*h
//	            }
//	        case (13, 14)
//	            ! broyden tridiagonal and banded functions.
//	            for(j=0;j<n;j++){
//	                x[j] = -one
//	            }
//	        default:
//	            ! rosenbrock function.
//	            x[0] = -c1
//	            x[1] = one
//	        } 
//
//	        ! compute multiple of initial point.
//	        if (factor != one) then
//	            if (nProb == 6) then
//	                for(j=0;j<n;j++){
//	                    x[j] = factor
//	                }
//	            else
//	                for(j=0;j<n;j++){
//	                    x[j] = factor*x[j]
//	                }
//	            end if
//	        end if
//
//	    }  // initpt
//	!*****************************************************************************************
//
//	!*****************************************************************************************
//	!>
//	!  This subroutine defines fourteen test functions. the first
//	!  five test functions are of dimensions 2,4,2,4,3, respectively,
//	!  while the remaining test functions are of variable dimension
//	!  n for any n greater than or equal to 1 (problem 6 is an
//	!  exception to this, Math.since it does not allow n = 1).
//
//	    subroutine vecfcn(n, x, fVec, nProb)
//	        implicit none
//
//	        integer,intent(in) :: n !! a positive integer input variable.
//	        real(wp),intent(in) :: x(n) !! an input array of length n.
//	        real(wp),intent(out) :: fVec(n) !! an output array of length n which contains the nProb
//	                                        !! function vector evaluated at x.
//	        integer,intent(in) :: nProb !! a positive integer input variable which defines the
//	                                    !! number of the problem. nProb must not exceed 14.
//
//	        double zero = 0.0
//	        double one = 1.0
//	        double two = 2.0
//	        double three = 3.0
//	        double five = 5.0
//	        double eight = 8.0
//	        double ten = 10.0
//
//	        double c1 =  1.0e4
//	        double c2 =  1.0001
//	        double c3 =  2.0e2
//	        double c4 =  2.02e1
//	        double c5 =  1.98e1
//	        double c6 =  1.8e2
//	        double c7 =  2.5e-1
//	        double c8 =  5.0e-1
//	        double c9 =  2.9e1
//
//	        double tpi = eight*Math.atan(one)
//
//	        int i, iev, j, k, k1, k2, kp1, ml, mu
//	        double h, prod, sum, sum1, sum2, temp, temp1, 
//	                    temp2, ti, tj, tk
//
//	        fVec(1:n) = zero
//
//	        ! problem selector.
//
//	        switch (nProb)
//	        case (2)
//	            ! powell Math.singular function.
//	            fVec[0] = x[0] + ten*x[1]
//	            fVec[1] = Math.sqrt(five)*(x[2] - x[3])
//	            fVec[2] = (x[1] - two*x[2])**2
//	            fVec[3] = Math.sqrt(ten)*(x[0] - x[3])**2
//	        case (3)
//	            ! powell badly scaled function.
//	            fVec[0] = c1*x[0]*x[1] - one
//	            fVec[1] = exp(-x[0]) + exp(-x[1]) - c2
//	        case (4)
//	            ! wood function.
//	            temp1 = x[1] - x[0]**2
//	            temp2 = x[3] - x[2]**2
//	            fVec[0] = -c3*x[0]*temp1 - (one - x[0])
//	            fVec[1] = c3*temp1 + c4*(x[1] - one) + c5*(x[3] - one)
//	            fVec[2] = -c6*x[2]*temp2 - (one - x[2])
//	            fVec[3] = c6*temp2 + c4*(x[3] - one) + c5*(x[1] - one)
//	        case (5)
//	            ! helical valley function.
//	            temp1 = Math.signum(c7, x[1])
//	            if (x[0] > zero) temp1 = Math.atan(x[1]/x[0])/tpi
//	            if (x[0] < zero) temp1 = Math.atan(x[1]/x[0])/tpi + c8
//	            temp2 = Math.sqrt(x[0]**2 + x[1]**2)
//	            fVec[0] = ten*(x[2] - ten*temp1)
//	            fVec[1] = ten*(temp2 - one)
//	            fVec[2] = x[2]
//	        case (6)
//	            ! watson function.
//	            for(k=0;k<n;k++){
//	                fVec[k] = zero
//	            }
//	            for(i=0;i<29;i++){
//	                ti = dfloat(i)/c9
//	                sum1 = zero
//	                temp = one
//	                do j = 2, n
//	                    sum1 = sum1 + dfloat(j - 1)*temp*x[j]
//	                    temp = ti*temp
//	                }
//	                sum2 = zero
//	                temp = one
//	                for(j=0;j<n;j++){
//	                    sum2 = sum2 + temp*x[j]
//	                    temp = ti*temp
//	                }
//	                temp1 = sum1 - sum2**2 - one
//	                temp2 = two*ti*sum2
//	                temp = one/ti
//	                for(k=0;k<n;k++){
//	                    fVec[k] = fVec[k] + temp*(dfloat(k - 1) - temp2)*temp1
//	                    temp = ti*temp
//	                }
//	            }
//	            temp = x[1] - x[0]**2 - one
//	            fVec[0] = fVec[0] + x[0]*(one - two*temp)
//	            fVec[1] = fVec[1] + temp
//	        case (7)
//	            ! chebyquad function.
//	            for(k=0;k<n;k++){
//	                fVec[k] = zero
//	            }
//	            for(j=0;j<n;j++){
//	                temp1 = one
//	                temp2 = two*x[j] - one
//	                temp = two*temp2
//	                for(i=0;i<n;i++){
//	                    fVec(i) = fVec(i) + temp2
//	                    ti = temp*temp2 - temp1
//	                    temp1 = temp2
//	                    temp2 = ti
//	                }
//	            }
//	            tk = one/dfloat(n)
//	            iev = -1
//	            for(k=0;k<n;k++){
//	                fVec[k] = tk*fVec[k]
//	                if (iev > 0) fVec[k] = fVec[k] + one/(dfloat(k)**2 - one)
//	                iev = -iev
//	            }
//	        case (8)
//	            ! brown almost-linear function.
//	            sum = -dfloat(n + 1)
//	            prod = one
//	            for(j=0;j<n;j++){
//	                sum = sum + x[j]
//	                prod = x[j]*prod
//	            }
//	            for(k=0;k<n;k++){
//	                fVec[k] = x[k] + sum
//	            }
//	            fVec(n) = prod - one
//	        case (9)
//	            ! discrete boundary value function.
//	            h = one/dfloat(n + 1)
//	            for(k=0;k<n;k++){
//	                temp = (x[k] + dfloat(k)*h + one)**3
//	                temp1 = zero
//	                if (k != 1) temp1 = x(k - 1)
//	                temp2 = zero
//	                if (k != n) temp2 = x(k + 1)
//	                fVec[k] = two*x[k] - temp1 - temp2 + temp*h**2/two
//	            }
//	        case (10)
//	            ! discrete integral equation function.
//	            h = one/dfloat(n + 1)
//	            for(k=0;k<n;k++){
//	                tk = dfloat(k)*h
//	                sum1 = zero
//	                for(j=0;j<k;j++){
//	                    tj = dfloat(j)*h
//	                    temp = (x[j] + tj + one)**3
//	                    sum1 = sum1 + tj*temp
//	                }
//	                sum2 = zero
//	                kp1 = k + 1
//	                if (n >= kp1) then
//	                    do j = kp1, n
//	                        tj = dfloat(j)*h
//	                        temp = (x[j] + tj + one)**3
//	                        sum2 = sum2 + (one - tj)*temp
//	                    }
//	                end if
//	                fVec[k] = x[k] + h*((one - tk)*sum1 + tk*sum2)/two
//	            }
//	        case (11)
//	            ! trigonometric function.
//	            sum = zero
//	            for(j=0;j<n;j++){
//	                fVec(j) = cos(x[j])
//	                sum = sum + fVec(j)
//	            }
//	            for(k=0;k<n;k++){
//	                fVec[k] = dfloat(n + k) - Math.sin(x[k]) - sum - dfloat(k)*fVec[k]
//	            }
//	        case (12)
//	            ! variably dimensioned function.
//	            sum = zero
//	            for(j=0;j<n;j++){
//	                sum = sum + dfloat(j)*(x[j] - one)
//	            }
//	            temp = sum*(one + two*sum**2)
//	            for(k=0;k<n;k++){
//	                fVec[k] = x[k] - one + dfloat(k)*temp
//	            }
//	        case (13)
//	            ! broyden tridiagonal function.
//	            for(k=0;k<n;k++){
//	                temp = (three - two*x[k])*x[k]
//	                temp1 = zero
//	                if (k != 1) temp1 = x(k - 1)
//	                temp2 = zero
//	                if (k != n) temp2 = x(k + 1)
//	                fVec[k] = temp - temp1 - two*temp2 + one
//	            }
//	        case (14)
//	            ! broyden banded function.
//	            ml = 5
//	            mu = 1
//	            for(k=0;k<n;k++){
//	                k1 = Math.max(1, k - ml)
//	                k2 = Math.min(k + mu, n)
//	                temp = zero
//	                for(j=k1;j<k2;j++){
//	                    if (j != k) temp = temp + x[j]*(one + x[j])
//	                }
//	                fVec[k] = x[k]*(two + five*x[k]**2) + one - temp
//	            }
//	        default:
//	            ! rosenbrock function.
//	            fVec[0] = one - x[0]
//	            fVec[1] = ten*(x[1] - x[0]**2)
//	        } 
//
//	    }  // vecfcn
//
//	end program test_hybrj
}
