package Minpack;

import java.util.ArrayList;

public class TestPack {
	
	//Main to run tests
	//Need to finish converting the rest of the tests
	public static void main(String args[]) {
		
	}
	
	public static void chkderTest() {

	}
	
	
	public static void hybrdTest() {
		
		ArrayList<double[]> initialX = new ArrayList<double[]>();
		int ncases = 22;
		double[] nprobs = {1,2,3,4,5,6,6,7,7,7,7,7,8,8,8,9,10,10,11,12,13,14};
		double[] ns =     {2,4,2,4,3,6,9,5,6,7,8,9,10,30,40,10,1,10,10,10,10,10};
		double[] ntriess = {3,3,2,3,3,2,2,3,3,3,1,1,3,1,1,3,3,3,3,3,3,3};
				
		//double[] x0 = {0n};
        double[] x1 = {3,-1,0,1};  // powell singular function.
        initialX.add(x1);
        double[] x2 = {0,1}; // powell badly scaled function.
        initialX.add(x2);
        double[] x3 = {-3,-1,-3,-1}; //wood function
        initialX.add(x3);
        double[] x4 = {-1,0,0}; //helical valley function.
        initialX.add(x4);
        double[] x5 = {0,0,0,0,0,0}; //watson function.
        initialX.add(x5);
        double[] x5a = {0,0,0,0,0,0,0,0,0}; //watson function.
        initialX.add(x5a);
		double[] x6 = {1*1/(5+1),2*1/(5+1),3*1/(5+1),4*1/(5+1),5*1/(5+1)}; //chebyquad function. 5n
		initialX.add(x6);
		 double[] x6a = {1*1/(5+1),2*1/(5+1),3*1/(5+1),4*1/(5+1),5*1/(5+1),6*1/(5+1)}; //chebyquad function. 6n
		initialX.add(x6a);
		 double[] x6b = {1*1/(5+1),2*1/(5+1),3*1/(5+1),4*1/(5+1),5*1/(5+1),6*1/(5+1),7*1/(5+1)}; //chebyquad function.7n
		initialX.add(x6b);
		 double[] x6c = {1*1/(5+1),2*1/(5+1),3*1/(5+1),4*1/(5+1),5*1/(5+1),6*1/(5+1),7*1/(5+1),8*1/(5+1)}; //chebyquad function.8n
		initialX.add(x6c);
		 double[] x6d = {1*1/(5+1),2*1/(5+1),3*1/(5+1),4*1/(5+1),5*1/(5+1),6*1/(5+1),7*1/(5+1),8*1/(5+1),9*1/(5+1)}; //chebyquad function.9n
		initialX.add(x6d);
		 double[] x7 = {.5,.5,.5,.5,.5,.5,.5,.5,.5,.5}; //brown almost-linear function.10n
		initialX.add(x7);
		 double[] x7a = {.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5}; //brown almost-linear function.30n
		initialX.add(x7a);
		 double[] x7b = {.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5}; //brown almost-linear function.40n
		initialX.add(x7b);
		 double[] x8 = {(1/(1/10+1))*((1/(1/10+1))-1),(2/(1/10+1))*((2/(1/10+1))-1),(3/(1/10+1))*((3/(1/10+1))-1),(4/(1/10+1))*((4/(1/10+1))-1),(5/(1/10+1))*((5/(1/10+1))-1),(6/(1/10+1))*((6/(1/10+1))-1),(7/(1/10+1))*((7/(1/10+1))-1),(8/(1/10+1))*((8/(1/10+1))-1),(9/(1/10+1))*((9/(1/10+1))-1),(10/(1/10+1))*((10/(1/10+1))-1)}; //discrete boundary value and integral equation10n
		// functions.
		initialX.add(x8);
		 double[] x9 = {(1/(1/10+1))*((1/(1/10+1))-1)}; //discrete boundary value and integral equation1n
		// functions.
		initialX.add(x9);
		 double[] x9a = {(1/(1/10+1))*((1/(1/10+1))-1),(2/(1/10+1))*((2/(1/10+1))-1),(3/(1/10+1))*((3/(1/10+1))-1),(4/(1/10+1))*((4/(1/10+1))-1),(5/(1/10+1))*((5/(1/10+1))-1),(6/(1/10+1))*((6/(1/10+1))-1),(7/(1/10+1))*((7/(1/10+1))-1),(8/(1/10+1))*((8/(1/10+1))-1),(9/(1/10+1))*((9/(1/10+1))-1),(10/(1/10+1))*((10/(1/10+1))-1)}; //discrete boundary value and integral equation10n
		// functions.
		initialX.add(x9a);
		 double[] x10 = {1/10,1/10,1/10,1/10,1/10,1/10,1/10,1/10,1/10,1/10}; //trigonometric function.10n
		initialX.add(x10);
		 double[] x11 = {1-1/10,1-2/10,1-3/10,1-4/10,1-5/10,1-6/10,1-7/10,1-8/10,1-9/10,1-10/10}; //variably dimensioned function.10n
		initialX.add(x11);
		 double[] x12 = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}; //broyden tri-diagonal and banded functions.10n
		initialX.add(x12);
		 double[] x13 = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}; //broyden tri-diagonal and banded functions.10n
		initialX.add(x13);
        
        //default case
        double[] x14 = {-1.2,1}; //rosenbrock function.
        
//        ! compute multiple of initial point.
//
//        if (factor /= one) then
//            if (nprob == 6) then
//                do j = 1, n
//                    x(j) = factor
//                end do
//            else
//                do j = 1, n
//                    x(j) = factor*x(j)
//                end do
//            end if
//        end if
		ArrayList<SystemOfEquations> problems = new ArrayList<SystemOfEquations>();
		
		
		// Test 1 POWELL SINGULAR FUNCTION.
		SystemOfEquations fX = new SystemOfEquations();
		fX.addFunction(x ->  x[0] + 10.0*x[1]);
		fX.addFunction(x ->  Math.sqrt(5)*(x[2] - x[3]));
		fX.addFunction(x ->  (x[1] - 2*x[2])*(x[1] - 2*x[2]));
		fX.addFunction(x -> Math.sqrt(10.0)*((x[0] - x[3])*(x[0] - x[3])));
		problems.add(fX);
		
//		// Test 2! POWELL BADLY SCALED FUNCTION.
//        fX.addFunction(x -> c1*x(1)*x(2) - one
//        fX.addFunction(x -> exp(-x(1)) + exp(-x(2)) - c2
//        
//	case (4)
//    ! WOOD FUNCTION.
//    temp1 = x(2) - x(1)**2
//    temp2 = x(4) - x(3)**2
//    fX.addFunction(x -> -c3*x(1)*temp1 - (one - x(1))
//    fX.addFunction(x -> c3*temp1 + c4*(x(2) - one) + c5*(x(4) - one)
//    fX.addFunction(x -> -c6*x(3)*temp2 - (one - x(3))
//    fX.addFunction(x -> c6*temp2 + c4*(x(4) - one) + c5*(x(2) - one)
//	case (5)
//    ! HELICAL VALLEY FUNCTION.
//    tpi = eight*atan(one)
//    if (x(1) > zero) then
//        temp1 = atan(x(2)/x(1))/tpi
//    else if (x(1) < zero) then
//        temp1 = atan(x(2)/x(1))/tpi + c8
//    else
//        temp1 = sign(c7, x(2)) ! does this ever happen?
//    end if
//    temp2 = sqrt(x(1)**2 + x(2)**2)
//    fX.addFunction(x -> ten*(x(3) - ten*temp1)
//    fX.addFunction(x -> ten*(temp2 - one)
//    fX.addFunction(x -> x(3)
//case (6)
//    ! WATSON FUNCTION.
//    do k = 1, n
//        Fvec(k) = zero
//    end do
//    do i = 1, 29
//        ti = dfloat(i)/c9
//        sum1 = zero
//        temp = one
//        do j = 2, n
//            sum1 = sum1 + dfloat(j - 1)*temp*x(j)
//            temp = ti*temp
//        end do
//        sum2 = zero
//        temp = one
//        do j = 1, n
//            sum2 = sum2 + temp*x(j)
//            temp = ti*temp
//        end do
//        temp1 = sum1 - sum2**2 - one
//        temp2 = two*ti*sum2
//        temp = one/ti
//        do k = 1, n
//            Fvec(k) = Fvec(k) + temp*(dfloat(k - 1) - temp2)*temp1
//            temp = ti*temp
//        end do
//    end do
//    temp = x(2) - x(1)**2 - one
//    fX.addFunction(x -> = Fvec(1) + x(1)*(one - two*temp)
//    fX.addFunction(x -> = Fvec(2) + temp
//case (7)
//    ! CHEBYQUAD FUNCTION.
//    do k = 1, n
//        Fvec(k) = zero
//    end do
//    do j = 1, n
//        temp1 = one
//        temp2 = two*x(j) - one
//        temp = two*temp2
//        do i = 1, n
//            Fvec(i) = Fvec(i) + temp2
//            ti = temp*temp2 - temp1
//            temp1 = temp2
//            temp2 = ti
//        end do
//    end do
//    tk = one/dfloat(n)
//    iev = -1
//    do k = 1, n
//        Fvec(k) = tk*Fvec(k)
//        if (iev > 0) Fvec(k) = Fvec(k) + one/(dfloat(k)**2 - one)
//        iev = -iev
//    end do
//case (8)
//    ! BROWN ALMOST-LINEAR FUNCTION.
//    sum = -dfloat(n + 1)
//    prod = one
//    do j = 1, n
//        sum = sum + x(j)
//        prod = x(j)*prod
//    end do
//    do k = 1, n
//        Fvec(k) = x(k) + sum
//    end do
//    Fvec(n) = prod - one
//case (9)
//    ! DISCRETE BOUNDARY VALUE FUNCTION.
//    h = one/dfloat(n + 1)
//    do k = 1, n
//        temp = (x(k) + dfloat(k)*h + one)**3
//        temp1 = zero
//        if (k /= 1) temp1 = x(k - 1)
//        temp2 = zero
//        if (k /= n) temp2 = x(k + 1)
//        Fvec(k) = two*x(k) - temp1 - temp2 + temp*h**2/two
//    end do
//case (10)
//    ! DISCRETE INTEGRAL EQUATION FUNCTION.
//    h = one/dfloat(n + 1)
//    do k = 1, n
//        tk = dfloat(k)*h
//        sum1 = zero
//        do j = 1, k
//            tj = dfloat(j)*h
//            temp = (x(j) + tj + one)**3
//            sum1 = sum1 + tj*temp
//        end do
//        sum2 = zero
//        kp1 = k + 1
//        if (n >= kp1) then
//            do j = kp1, n
//                tj = dfloat(j)*h
//                temp = (x(j) + tj + one)**3
//                sum2 = sum2 + (one - tj)*temp
//            end do
//        end if
//        Fvec(k) = x(k) + h*((one - tk)*sum1 + tk*sum2)/two
//    end do
//case (11)
//    ! TRIGONOMETRIC FUNCTION.
//    sum = zero
//    do j = 1, n
//        Fvec(j) = cos(x(j))
//        sum = sum + Fvec(j)
//    end do
//    do k = 1, n
//        Fvec(k) = dfloat(n + k) - sin(x(k)) - sum - dfloat(k)*Fvec(k)
//    end do
//case (12)
//    ! VARIABLY DIMENSIONED FUNCTION.
//    sum = zero
//    do j = 1, n
//        sum = sum + dfloat(j)*(x(j) - one)
//    end do
//    temp = sum*(one + two*sum**2)
//    do k = 1, n
//        Fvec(k) = x(k) - one + dfloat(k)*temp
//    end do
//case (13)
//    ! BROYDEN TRIDIAGONAL FUNCTION.
//    do k = 1, n
//        temp = (three - two*x(k))*x(k)
//        temp1 = zero
//        if (k /= 1) temp1 = x(k - 1)
//        temp2 = zero
//        if (k /= n) temp2 = x(k + 1)
//        Fvec(k) = temp - temp1 - two*temp2 + one
//    end do
//case (14)
//    ! BROYDEN BANDED FUNCTION.
//    ml = 5
//    mu = 1
//    do k = 1, n
//        k1 = max(1, k - ml)
//        k2 = min(k + mu, n)
//        temp = zero
//        do j = k1, k2
//            if (j /= k) temp = temp + x(j)*(one + x(j))
//        end do
//        Fvec(k) = x(k)*(two + five*x(k)**2) + one - temp
//    end do
//case default
//    ! ROSENBROCK FUNCTION.
//    Fvec(1) = one - x(1)
//    Fvec(2) = ten*(x(2) - x(1)**2)
//end select
		
//		Expected Solutions for each test
		ArrayList<double[]> solutions = new ArrayList<double[]>();
		
	        double[]  s0 = {0.1000000000000000E+01,0.1000000000000000E+01};
	        double[]  s1 = {0.1000000000000000E+01,0.1000000000000000E+01};
	        double[]  s2 = {0.1000000000000000E+01,0.1000000000000054E+01};
	        double[]  s3 = {-0.1717232736993598E-17,0.1717232736993598E-18,0.4885791716489701E-17,
	                       0.4885791716489701E-17};
	        double[]  s4 = {0.1137338840805565E-17,-0.1137338840805565E-18,0.1509231876962185E-17,
	                        0.1509231876962185E-17};
	        double[]  s5 = {0.2071578632741476E-17,-0.2071578632741476E-18,0.3365608460018520E-17,
	                        0.3365608460018520E-17};
	        double[]  s6 = {0.1098159327798559E-04,0.9106146740037904E+01};
	        double[]  s7 = {0.1098159288127784E-04,0.9106146743611500E+01};
	        double[]  s8 = {-0.9679740249513952,0.9471391408446626,-0.9695163103174519,
	                        0.9512476657649955};
	        double[]  s9 = {-0.9679740249362032,0.9471391408151544,-0.9695163103329795,
	                        0.9512476657950136};
	        double[]  s10 = {-0.9679740247487649,0.9471391404515958,-0.9695163105216826,
	                        0.9512476661606586};
	        double[]  s11 = {0.1000000000000010E+01,-0.1612103012704913E-13,-0.8125674485239315E-34};
	        double[]  s12 = {0.1000000000004274E+01,0.1388633807362216E-10,0.0000000000000000};
	        double[]  s13 = {0.9999999999999840,0.1238719384388127E-14,0.0000000000000000};
	        double[]  s14 = {-0.1572508640134011E-01,0.1012434869369118E+01,-0.2329916259567960,
	                        0.1260430087800365E+01,-0.1513728922723441E+01,0.9929964324318560};
	        double[]  s15 = {-0.1572508639053690E-01,0.1012434869369907E+01,-0.2329916259623205,
	                        0.1260430087870332E+01,-0.1513728922830754E+01,0.9929964324984680};
	        double[]  s16 = {-0.1530703652147214E-04,0.9997897039319488,0.1476396369354227E-01,
	                        0.1463423282995085,0.1000821103004075E+01,-0.2617731140517609E+01,
	                        0.4104403164477184E+01,-0.3143612278555425E+01,0.1052626408009917E+01};
	        double[]  s17= {-0.1530703713913476E-04,0.9997897039319459,0.1476396369267279E-01,
	                        0.1463423283016577,0.1000821102994390E+01,-0.2617731140495362E+01,
	                        0.4104403164446781E+01,-0.3143612278533631E+01,0.1052626408003180E+01};
	        double[]  s18 = {0.8375125649983552E-01,0.3127292952224503,0.5000000000008663,
	                        0.6872707047760241,0.9162487435008237};
	        double[]  s19 = {0.6872707047770207,0.9162487435004409,0.4999999999996554,
	                        0.8375125649942240E-01,0.3127292952234606};
	        double[]  s20 = {0.5000000001403082,0.6872707047172968,0.8375125651152422E-01,
	                        0.3127292951388482,0.9162487434920225};
	        double[]  s21 = {0.6687659094768218E-01,0.3666822990106460,0.2887406733471217,
	                        0.7112593271119373,0.6333177005294922,0.9331234090531204};
	        double[]  s22 = {0.9331234090539151,0.2887406731191550,0.6687659094608109E-01,
	                        0.7112593268807791,0.3666822992420550,0.6333177007580146};
	        double[]  s23 = {0.3666822992460977,0.6333177007631426,0.7112593268766961,
	                        0.6687659094559376E-01,0.9331234090527033,0.2887406731157666};
	        double[]  s24 = {0.5806914977912385E-01,0.2351716115667845,0.3380440957994889,
	                        0.4999999991484015,0.6619559063505691,0.7648283868041619,
	                        0.9419308505514702};
	        double[]  s25 = {0.3380440947502161,0.4999999997939603,0.2351716123793344,
	                        0.7648283876237569,0.9419308503662596,0.5806914962097134E-01,
	                        0.6619559054655015};
	        double[]  s26 = {-0.3267366079625573E+02,-0.2996209843926172E+02,-0.8587775264169514E+02,
	                        0.2222113097994968E+02,0.5957249137089175E+02,-0.1038025653217158E+01,
	                        0.8600842862942351E+02};
	        double[]  s27 = {0.4985640222318974E-01,0.1986351285003365,0.2698288337443381,
	                        0.4992723176748156,0.5007277255753518,0.7301712224978171,
	                        0.8013649159179719,0.9501436601762751};
	        double[]  s28 = {0.4420534615691015E-01,0.1994906721904692,0.2356191086780681,
	                        0.4160469078466623,0.4999999996232831,0.5839530926716184,
	                        0.7643808911925417,0.8005093278089829,0.9557946538314642};
	        double[]  s29 = {0.1000000000000008E+01,0.1000000000000008E+01,0.1000000000000008E+01,
	                        0.1000000000000008E+01,0.1000000000000008E+01,0.1000000000000008E+01,
	                        0.1000000000000008E+01,0.1000000000000008E+01,0.1000000000000008E+01,
	                        0.9999999999999193};
	        double[]  s30 = {0.1000000000000002E+01,0.1000000000000002E+01,0.1000000000000002E+01,
	                        0.1000000000000002E+01,0.1000000000000002E+01,0.1000000000000002E+01,
	                        0.1000000000000002E+01,0.1000000000000002E+01,0.1000000000000002E+01,
	                        0.9999999999999805};
	        double[]  s31 = {0.1000000000000126E+01,0.1000000000000126E+01,0.1000000000000126E+01,
	                        0.1000000000000126E+01,0.1000000000000126E+01,0.1000000000000126E+01,
	                        0.1000000000000126E+01,0.1000000000000126E+01,0.1000000000000126E+01,
	                        0.9999999999987337};
	        double[]  s32 = {0.1000000000000116E+01,0.1000000000000116E+01,0.1000000000000116E+01,
	                        0.1000000000000116E+01,0.1000000000000116E+01,0.1000000000000116E+01,
	                        0.1000000000000116E+01,0.1000000000000116E+01,0.1000000000000116E+01,
	                        0.1000000000000116E+01,0.1000000000000116E+01,0.1000000000000116E+01,
	                        0.1000000000000116E+01,0.1000000000000116E+01,0.1000000000000116E+01,
	                        0.1000000000000116E+01,0.1000000000000116E+01,0.1000000000000116E+01,
	                        0.1000000000000116E+01,0.1000000000000116E+01,0.1000000000000116E+01,
	                        0.1000000000000116E+01,0.1000000000000116E+01,0.1000000000000116E+01,
	                        0.1000000000000116E+01,0.1000000000000116E+01,0.1000000000000116E+01,
	                        0.1000000000000116E+01,0.1000000000000116E+01,0.9999999999965015};
	        double[]  s33 = {0.9999999999999820,0.9999999999999820,0.9999999999999820,
	                        0.9999999999999820,0.9999999999999820,0.9999999999999820,
	                        0.9999999999999820,0.9999999999999820,0.9999999999999820,
	                        0.9999999999999820,0.9999999999999820,0.9999999999999820,
	                        0.9999999999999820,0.9999999999999820,0.9999999999999820,
	                        0.9999999999999820,0.9999999999999820,0.9999999999999820,
	                        0.9999999999999820,0.9999999999999820,0.9999999999999820,
	                        0.9999999999999820,0.9999999999999820,0.9999999999999820,
	                        0.9999999999999820,0.9999999999999820,0.9999999999999820,
	                        0.9999999999999820,0.9999999999999820,0.9999999999999820,
	                        0.9999999999999820,0.9999999999999820,0.9999999999999820,
	                        0.9999999999999820,0.9999999999999820,0.9999999999999820,
	                        0.9999999999999820,0.9999999999999820,0.9999999999999820,
	                        0.1000000000000726E+01};
	        double[]  s34 = {-0.4316498251876486E-01,-0.8157715653538729E-01,-0.1144857143805310,
	                        -0.1409735768625996,-0.1599086961819857,-0.1698772023127759,
	                        -0.1690899837812081,-0.1552495352218312,-0.1253558916789345,
	                        -0.7541653368589182E-01};
	        double[]  s35 = {-0.4316498251881876E-01,-0.8157715653546944E-01,-0.1144857143805966,
	                        -0.1409735768626190,-0.1599086961819499,-0.1698772023126901,
	                        -0.1690899837811062,-0.1552495352217907,-0.1253558916789970,
	                        -0.7541653368596339E-01};
	        double[]  s36 = {-0.4316498254524553E-01,-0.8157715658128553E-01,-0.1144857144024714,
	                        -0.1409735768723229,-0.1599086963003020,-0.1698772022538506,
	                        -0.1690899837944776,-0.1552495352060321,-0.1253558916432041,
	                        -0.7541653366609002E-01};
	        double[]  s37 = {-0.1528138835625800};
	        double[]  s38 = {-0.1528138835625801};
	        double[]  s39 = {-0.1528138835625801};
	        double[]  s40 = {-0.4316498251876487E-01,-0.8157715653538730E-01,-0.1144857143805310,
	                        -0.1409735768625996,-0.1599086961819857,-0.1698772023127759,
	                        -0.1690899837812081,-0.1552495352218312,-0.1253558916789344,
	                        -0.7541653368589175E-01};
	        double[]  s41 = {-0.4316498251881876E-01,-0.8157715653546946E-01,-0.1144857143805966,
	                        -0.1409735768626190,-0.1599086961819498,-0.1698772023126901,
	                        -0.1690899837811062,-0.1552495352217907,-0.1253558916789970,
	                        -0.7541653368596334E-01};
	        double[]  s42 = {-0.4316498251876519E-01,-0.8157715653538752E-01,-0.1144857143805303,
	                        -0.1409735768625980,-0.1599086961819844,-0.1698772023127748,
	                        -0.1690899837812073,-0.1552495352218307,-0.1253558916789341,
	                        -0.7541653368589159E-01};
	        double[]  s43 = {0.5526115715943968E-01,0.5695713693095779E-01,0.5889020336090119E-01,
	                        0.6113556183519356E-01,0.6377799292665667E-01,0.6700414432471043E-01,
	                        0.2079417258113421,0.1642681193175131,0.8643704095817571E-01,
	                        0.9133212907808361E-01};
	        double[]  s44 = {0.3439628896235289E-01,0.3503231575416022E-01,0.3571919583574593E-01,
	                        0.3646522422001942E-01,0.3728091174083566E-01,0.3817986258974846E-01,
	                        0.3918014109819012E-01,0.4030650261419996E-01,0.1797201916815169,
	                        0.1562408814749922};
	        double[]  s45 = {0.1888395221036672E+02,0.2516777354434988E+02,0.1888527511739392E+02,
	                        0.1888602114554983E+02,0.1888683683452867E+02,0.1888773578345333E+02,
	                        0.1888873606083097E+02,0.1888986242379029E+02,0.1902927611240455E+02,
	                        0.1900579680335179E+02};
	        double[]  s46 = {0.9999999999999992,0.9999999999999984,0.9999999999999977,
	                        0.9999999999999969,0.9999999999999961,0.9999999999999953,
	                        0.9999999999999946,0.9999999999999938,0.9999999999999930,
	                        0.9999999999999922};
	        double[]  s47 = {0.1000000000000000E+01,0.1000000000000000E+01,0.1000000000000000E+01,
	                        0.1000000000000000E+01,0.1000000000000001E+01,0.1000000000000001E+01,
	                        0.1000000000000001E+01,0.1000000000000001E+01,0.1000000000000001E+01,
	                        0.1000000000000001E+01};
	        double[]  s48 = {0.1000000000000000E+01,0.1000000000000000E+01,0.1000000000000000E+01,
	                        0.1000000000000000E+01,0.1000000000000000E+01,0.1000000000000000E+01,
	                        0.1000000000000000E+01,0.1000000000000000E+01,0.1000000000000000E+01,
	                        0.1000000000000000E+01};
	        double[]  s49 = {-0.5707221307212121,-0.6818069509055232,-0.7022100775689857,
	                        -0.7055106309936168,-0.7049061557572888,-0.7014966060124587,
	                        -0.6918893211477919,-0.6657965141985400,-0.5960351099566767,
	                        -0.4164122574358191};
	        double[]  s50 = {-0.5707221320932939,-0.6818069495100820,-0.7022100764111258,
	                        -0.7055106298493696,-0.7049061556844529,-0.7014966070294095,
	                        -0.6918893223739674,-0.6657965143753547,-0.5960351092038981,
	                        -0.4164122574142932};
	        double[]  s51 = {-0.5707221320171143,-0.6818069499829604,-0.7022100760171542,
	                        -0.7055106298955310,-0.7049061557301967,-0.7014966070327222,
	                        -0.6918893223590803,-0.6657965144072677,-0.5960351090088830,
	                        -0.4164122575177334};
	        double[]  s52 = {-0.4283028636053099,-0.4765964242962535,-0.5196524638125549,
	                        -0.5580993246169652,-0.5925061569509362,-0.6245036821428087,
	                        -0.6232394714478015,-0.6213938418388717,-0.6204535966122983,
	                        -0.5864692707477792};
	        double[]  s53 = {-0.4283028634881692,-0.4765964236396713,-0.5196524642776766, -0.5580993248351936,-0.5925061568131795,-0.6245036817962692, -0.6232394720687789,-0.6213938417874499,-0.6204535965224117, -0.5864692707287930};
	        double[]  s54 = {-0.4283028635608067,-0.4765964243232715,-0.5196524637037395, -0.5580993248328234,-0.5925061568292707,-0.6245036822076749, -0.6232394714256790,-0.6213938418143937,-0.6204535966527651 -0.5864692707189498};
	        
	        solutions.add(s1);
	        solutions.add(s2);
	        solutions.add(s3);
	        solutions.add(s4);
	        solutions.add(s5);
	        solutions.add(s6);
	        solutions.add(s7);
	        solutions.add(s8);
	        solutions.add(s9);
	        solutions.add(s10);
	        solutions.add(s11);
	        solutions.add(s12);
	        solutions.add(s13);
	        solutions.add(s14);
	        solutions.add(s15);
	        solutions.add(s16);
	        solutions.add(s17);
	        solutions.add(s18);
	        solutions.add(s19);
	        solutions.add(s20);
	        solutions.add(s21);
	        solutions.add(s22);
	        solutions.add(s23);
	        solutions.add(s24);
	        solutions.add(s25);
	        solutions.add(s26);
	        solutions.add(s27);
	        solutions.add(s28);
	        solutions.add(s29);
	        solutions.add(s30);
	        solutions.add(s31);
	        solutions.add(s32);
	        solutions.add(s33);
	        solutions.add(s34);
	        solutions.add(s35);
	        solutions.add(s36);
	        solutions.add(s37);
	        solutions.add(s38);
	        solutions.add(s39);
	        solutions.add(s40);
	        solutions.add(s41);
	        solutions.add(s42);
	        solutions.add(s43);
	        solutions.add(s44);
	        solutions.add(s45);
	        solutions.add(s46);
	        solutions.add(s47);
	        solutions.add(s48);
	        solutions.add(s49);
	        solutions.add(s50);
	        solutions.add(s51);
	        solutions.add(s52);
	        solutions.add(s53);
	        solutions.add(s54);
	        
	        
	       //run tests
	        double[] xF;
	        xF = Minpack.hybrd1(fX,x1,1.490116119384766E-008);
	        
	       //compare solutions
	        boolean pass = true;
	        double abstol = 0.00000001;
	        double[] reldiff = new double[40];
	        double reltol = 0.0001;
	        double[] absdiff = new double[40];
	        for(int i = 0;i<s4.length;i++) {
	        	absdiff[i] = Math.abs(s4[i] - x1[i]);
	        	if(absdiff[i]>abstol) {
	        		if(s4[i]!=0){
	        			reldiff[i] = absdiff[i]/Math.abs(s4[i]);
	        		}
	        		if(reldiff[i] > reltol) {
	        			System.out.println("Solution has failed for" + "s4");
	        			pass = false;
	        		}
	        		pass = false;
	        	}
	        }
	        if(pass) {
	        	System.out.println("Success at "+"s4");	
	        }
	        pass = true;
	       
	        System.out.println("fin");
		
	}

	public static void hybrjTest() {

	}

	public static void lmderTest() {

	}

	public static void lmdifTest() {

	}

	public static void lmstrTest() {

	}

}
