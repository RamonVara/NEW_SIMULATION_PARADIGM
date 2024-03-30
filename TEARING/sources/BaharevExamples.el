/*-----------------------------------------------------------------------------------------
 LIBRARY: TEARING
 FILE: BaharevExamples
 CREATION DATE: 08/12/2023
-----------------------------------------------------------------------------------------*/
--This EL file contains the Ecosin
FUNCTION REAL f(IN REAL x)
	BODY
		RETURN x
END FUNCTION

COMPONENT BaharevExample1(INTEGER N = 20)
	DATA	
		REAL x_o = 0.1
		REAL x_N1 = 0.1
	DECLS
		BOUND REAL k = 1
		REAL x[N]
	CONTINUOUS	
		f(x_o) + k * f(x[1]) + x[2] =  (k+2)* 0.1
		
		EXPAND (i IN 2, N-1) 
			f(x[i-1]) + k * f(x[i]) + x[i+1] = (k+2) * 0.1
		
		f(x[N-1]) + k * f(x[N]) + x_N1 = (k+2)* 0.1
END COMPONENT

COMPONENT BaharevExample2(INTEGER N = 20)
	DATA	
		REAL x_o = 1
		REAL x_N1 = 1
	DECLS
		BOUND REAL k = 1
		REAL x[N]	
	CONTINUOUS		
		f(x_o) + f(x[1]) + k*x[2] = k+2
		
		EXPAND (i IN 2, N-1) 
				f(x[i-1]) + f(x[i]) + k*x[i+1] = k+2
		
		f(x[N-1]) + f(x[N]) + k*x_N1 = k+2		
END COMPONENT

COMPONENT BaharevExample3(INTEGER N =  300)
	DECLS
		REAL x[N]
	CONTINUOUS
		EXPAND (i IN 1, N-1) 
				f(x[i]) + f(x[N]) = 2
		
		SUM(i IN 1, N-1; f(x[i])) + x[N] = N

END COMPONENT
