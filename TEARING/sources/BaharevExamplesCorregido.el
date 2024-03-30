/*-----------------------------------------------------------------------------------------
 LIBRARY: TEARING
 FILE: BaharevExamples
 CREATION DATE: 08/12/2023
-----------------------------------------------------------------------------------------*/
USE SPARSE
FUNCTION NO_TYPE fsys_Example1(
			IN INTEGER iopt, 
			IN REAL x[], 
			IN REAL rinput[],
			OUT REAL routput[],
			OUT REAL F[], 
			OUT NLEQN_SYSTEM_BASE nleqn_system
			)
	DECLS
		REAL x_o 
		REAL x_N1
		REAL k
		INTEGER N
	BODY
		x_o = rinput[1]
		x_N1 = rinput[2]
		k = rinput[3]
		N = rinput[4]
		IF(iopt == 0) THEN
		   --Definition of the number of unknowns
			nleqn_system.neqn = 20
			--Creation of the unknowns
			FOR(j IN 1, nleqn_system.neqn)
		      nleqn_system.CreateUnk(j, 0., -1e6, 1e6, concatStrings(concatStrings("x[", integerToString(j)),"]"),0)
		   END FOR		
			--Creation of the unknowns
			FOR(j IN 1, nleqn_system.neqn)
		  		nleqn_system.CreateEqn(1,concatStrings("Equation ", integerToString(j)))
			END FOR
		ELSEIF(iopt == 1 OR iopt == 2) THEN 
			--Calculation of the residues	
	   	F[1] =  x_o + k * x[1] + x[2] - (k+2)* 0.1
			nleqn_system.Jacob.stampItem(1, 1, k)
			nleqn_system.Jacob.stampItem(1, 2, 1.)
			FOR(i IN 2, N-1)				
				F[i] = x[i-1] + k * x[i] + x[i+1] - (k+2) * 0.1
				nleqn_system.Jacob.stampItem(i, i-1, 1.)
				nleqn_system.Jacob.stampItem(i, i, k)
				nleqn_system.Jacob.stampItem(i, i+1, 1.)
			END FOR
			F[N] =   x[N-1] + k * x[N] + x_N1 -(k+2) * 0.1
			nleqn_system.Jacob.stampItem(N, N-1, 1.)
			nleqn_system.Jacob.stampItem(N, N, k)
		END IF
		RETURN
END FUNCTION 

COMPONENT BaharevExample1Corr(INTEGER N = 20)
	DATA	
		REAL x_o = 0.1
		REAL x_N1 = 0.1
		INTEGER iter_damp = 0   "Number of initial damped iterations"
		INTEGER iter_max = 40   "Maximum number of iterations"	
		REAL FRAC_damp = 0.1    "Damping factor for the initial damped iterations - fraction of the NR step" 
		REAL tol = 1e-10		   "Tolerance of the NR"
		INTEGER NR_debug = 2    "Debugging level of the Newton-Raphson"
		INTEGER STAMP_debug = 0 "Debug level of the stamping process"
	DECLS
		BOUND REAL k = 1
		REAL x[N] = 0
		REAL rinput[4]
		REAL routput[1]
	OBJECTS
		NLEQN_SYSTEM(MAX_EQN=50, MAX_STAMP = 100) nleqn_system	
	CONTINUOUS	
		rinput[1] = x_o
		rinput[2] = x_N1
		rinput[3] = k
		rinput[4] = N
		SEQUENTIAL
			nleqn_system.NR_debug = NR_debug
			nleqn_system.Jacob.STAMP_debug = STAMP_debug
			nleqn_system.Solve(fsys_Example1, iter_max, iter_damp, FRAC_damp, tol, rinput, routput)  
			nleqn_system.GetUnknown(x)
		END SEQUENTIAL
END COMPONENT

FUNCTION NO_TYPE fsys_Example2(
			IN INTEGER iopt, 
			IN REAL x[], 
			IN REAL rinput[],
			OUT REAL routput[],
			OUT REAL F[], 
			OUT NLEQN_SYSTEM_BASE nleqn_system
			)
	DECLS
		REAL x_o 
		REAL x_N1
		REAL k
		INTEGER N
	BODY
		x_o = rinput[1]
		x_N1 = rinput[2]
		k = rinput[3]
		N = rinput[4]
		IF(iopt == 0) THEN
		   --Definition of the number of unknowns
			nleqn_system.neqn = 20
			--Creation of the unknowns
			FOR(j IN 1, nleqn_system.neqn)
		      nleqn_system.CreateUnk(j, 0., -1e6, 1e6, concatStrings(concatStrings("x[", integerToString(j)),"]"),0)
		   END FOR		
			--Creation of the equations
			FOR(j IN 1, nleqn_system.neqn)
		  		nleqn_system.CreateEqn(1,concatStrings("Equation ", integerToString(j)))
			END FOR
		ELSEIF(iopt == 1 OR iopt == 2) THEN 
			--Calculation of the residues	
	     	F[1] = x_o + x[1] + k*x[2] - (k+2)
			nleqn_system.Jacob.stampItem(1, 1, 1.)
			nleqn_system.Jacob.stampItem(1, 2, k)
			FOR(i IN 2, N-1)				
				F[i] =   x[i-1] + x[i] + k*x[i+1] - (k+2)
				nleqn_system.Jacob.stampItem(i, i-1, 1.)
				nleqn_system.Jacob.stampItem(i, i  , 1.)
				nleqn_system.Jacob.stampItem(i, i+1,  k)
			END FOR
			F[N] =   x[N-1] + x[N] + k*x_N1 -(k+2)
			nleqn_system.Jacob.stampItem(N, N-1, 1.)
			nleqn_system.Jacob.stampItem(N, N,   1.)
		END IF
		RETURN
END FUNCTION

COMPONENT BaharevExample2Corr(INTEGER N = 20)
	DATA	
		REAL x_o  = 1.
		REAL x_N1 = 1.
		INTEGER iter_damp = 0   "Number of initial damped iterations"
		INTEGER iter_max = 40   "Maximum number of iterations"	
		REAL FRAC_damp = 0.1    "Damping factor for the initial damped iterations - fraction of the NR step" 
		REAL tol = 1e-10		   "Tolerance of the NR"
		INTEGER NR_debug = 2    "Debugging level of the Newton-Raphson"
		INTEGER STAMP_debug = 0 "Debug level of the stamping process"
	DECLS
		BOUND REAL k = 15
		REAL x[N] = 0.
		REAL rinput[4]
		REAL routput[1]
	OBJECTS
		NLEQN_SYSTEM(MAX_EQN=50, MAX_STAMP = 100) nleqn_system	
	CONTINUOUS	
		rinput[1] = x_o
		rinput[2] = x_N1
		rinput[3] = k
		rinput[4] = N
		SEQUENTIAL
			nleqn_system.NR_debug = NR_debug
			nleqn_system.Jacob.STAMP_debug = STAMP_debug
			nleqn_system.Solve(fsys_Example2, iter_max, iter_damp, FRAC_damp, tol, rinput, routput)  
			nleqn_system.GetUnknown(x)
		END SEQUENTIAL
END COMPONENT

FUNCTION NO_TYPE fsys_Example3(
			IN INTEGER iopt, 
			IN REAL x[], 
			IN REAL rinput[],
			OUT REAL routput[],
			OUT REAL F[], 
			OUT NLEQN_SYSTEM_BASE nleqn_system
			)
	DECLS
		INTEGER N
	BODY
		N = 300
		IF(iopt == 0) THEN
		   --Definition of the number of unknowns
			nleqn_system.neqn = 300
			nleqn_system.STAMP_debug = 2
			--Creation of the unknowns
			FOR(j IN 1, nleqn_system.neqn)
		      nleqn_system.CreateUnk(j, 0., -1e6, 1e6, concatStrings(concatStrings("x[", integerToString(j)),"]"),0)
		   END FOR		
			--Creation of the equations
			FOR(j IN 1, nleqn_system.neqn)
		  		nleqn_system.CreateEqn(1,concatStrings("Equation ", integerToString(j)))
			END FOR
		ELSEIF(iopt == 1 OR iopt == 2) THEN 
			--Calculation of the residues	
			FOR(i IN 1, N-1)				
				F[i] =    x[i] + x[N] - 2.			
				nleqn_system.Jacob.stampItem(i, i, 1.)
				nleqn_system.Jacob.stampItem(i, N, 1.)
			END FOR
			F[N] =   0.
			FOR(i IN 1, N) 
				F[N] = F[N] + x[i]
				nleqn_system.Jacob.stampItem(N, i, 1.)
			END FOR
			F[N] = F[N] - N
		END IF
		RETURN
END FUNCTION
COMPONENT BaharevExample3Corr(INTEGER N =  300)
	DATA	
		INTEGER iter_damp = 0   "Number of initial damped iterations"
		INTEGER iter_max = 40   "Maximum number of iterations"	
		REAL FRAC_damp = 0.1    "Damping factor for the initial damped iterations - fraction of the NR step" 
		REAL tol = 1e-10		   "Tolerance of the NR"
		INTEGER NR_debug = 2    "Debugging level of the Newton-Raphson"
		INTEGER STAMP_debug = 0 "Debug level of the stamping process"
	DECLS
		BOUND REAL k = 15
		REAL x[N] = 0.
		REAL rinput[1]
		REAL routput[1]
	OBJECTS
		NLEQN_SYSTEM(MAX_EQN=400, MAX_STAMP = 1000) nleqn_system
	INIT
		nleqn_system.NR_debug = NR_debug
		nleqn_system.Jacob.STAMP_debug = STAMP_debug
	CONTINUOUS	
	    rinput[1] = TIME
		SEQUENTIAL
			nleqn_system.Solve(fsys_Example3, iter_max, iter_damp, FRAC_damp, tol, rinput, routput)  
			nleqn_system.GetUnknown(x)
		END SEQUENTIAL
END COMPONENT
