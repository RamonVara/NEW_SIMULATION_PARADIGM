/*-----------------------------------------------------------------------------------------
 LIBRARY: SPARSE_EXAMPLES
 FILE: TestNleqnSystem
 CREATION DATE: 12/11/2023
-----------------------------------------------------------------------------------------*/
USE SPARSE
FUNCTION NO_TYPE fsystem(
			IN INTEGER iopt, 
			IN REAL x[], 
			IN REAL rinput[],
			OUT REAL routput[],
			OUT REAL F[], 
			OUT NLEQN_SYSTEM_BASE nleqn_system
			)
	DECLS
	BODY
		IF(iopt == 0) THEN
		   --Definition of the number of unknowns
			nleqn_system.neqn = 4
			--Creation of the unknowns
			FOR(j IN 1, nleqn_system.neqn)
		      nleqn_system.CreateUnk(j, 0.1, -1e6, 1e6, concatStrings(concatStrings("x[", integerToString(j)),"]"),0)
		   END FOR		
			--Creation of the equations
		   nleqn_system.CreateEqn(1,"-x[1]**2 - x[2]**2 - x[3]**2 + x[4]    = 0")
			nleqn_system.CreateEqn(2," x[1]**2 + x[2]**2 + x[3]**2 + x[4]**2 = 0")
			nleqn_system.CreateEqn(3," x[1] - x[2] = 0")
			nleqn_system.CreateEqn(4," x[2] - x[3] = 0")
		ELSEIF(iopt == 1 OR iopt == 2) THEN 
			--Calculation of the residues	
	   	F[1] =  - x[1]**2 - x[2]**2 - x[3]**2 + x[4]
			F[2] =    x[1]**2 + x[2]**2 + x[3]**2 + x[4]**2 - 1
			F[3] =    x[1] - x[2]
			F[4] =    x[2] - x[3]
			--Calculation of the Jacobian 
			IF(iopt == 2) THEN
				nleqn_system.Jacob.stampItem(1, 1, -2.*x[1])
				nleqn_system.Jacob.stampItem(1, 2, -2.*x[2])
				nleqn_system.Jacob.stampItem(1, 3, -2.*x[3])
				nleqn_system.Jacob.stampItem(1, 4, 1.0)
				nleqn_system.Jacob.stampItem(2, 1, 2.*x[1])
				nleqn_system.Jacob.stampItem(2, 2, 2.*x[2])
				nleqn_system.Jacob.stampItem(2, 3, 2.*x[3])
				nleqn_system.Jacob.stampItem(2, 4, 2.*x[4])
				nleqn_system.Jacob.stampItem(3, 1, 1.0)
				nleqn_system.Jacob.stampItem(3, 2, -1.0)		
				nleqn_system.Jacob.stampItem(4, 2,  1.0)
				nleqn_system.Jacob.stampItem(4, 3, -1.0)
			END IF
		END IF
		RETURN
END FUNCTION 

COMPONENT TestNleqnSystem01
	DATA
		INTEGER iter_damp = 2   "Number of initial damped iterations"
		INTEGER iter_max = 20   "Maximum number of iterations"	
		REAL FRAC_damp = 0.1    "Damping factor for the initial damped iterations - fraction of the NR step" 
		REAL tol = 1e-6		   "Tolerance of the NR"
		INTEGER NR_debug = 2    "Debugging level of the Newton-Raphson"
		INTEGER STAMP_debug = 0 "Debug level of the stamping process"
	DECLS
		DISCR REAL FRAC = 1
		DISCR REAL errorx 
		DISCR REAL errorf
		DISCR	REAL x[20]
		DISCR	REAL rinput[1]
		DISCR REAL routput[1]
	OBJECTS
		NLEQN_SYSTEM(MAX_EQN=100, MAX_STAMP = 200) nleqn_system		
	INIT
		nleqn_system.NR_debug = NR_debug
		nleqn_system.Jacob.STAMP_debug = STAMP_debug
		nleqn_system.Solve(fsystem, iter_max, iter_damp, FRAC_damp, tol, rinput, routput)  
		nleqn_system.GetUnknown(x)
END COMPONENT