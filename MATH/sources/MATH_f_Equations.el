--------------------------------------------------------------------------------
-- EA Internacional 2004   EcosimPro 3.3 Simulation Source Code
--
--    FILE NAME:  MATH_f_Equations.el
--    FILE TYPE:  External objects of the MATH library
--    DESCRIPTION:  Defines a non-linear system of equations solver for the
--                  MATH library
--    NOTES:  Based on EcosimPro MATH library
--    AUTHOR:  Borja Garcia Gutierrez
--    CREATION DATE:  12-Jan-2004
--------------------------------------------------------------------------------


-- External Objects

--------------------------------------------------------------------------------
-- Function EcoNonLinearEqnSolver
--------------------------------------------------------------------------------
-- Purpose: to find a zero of a system of N nonlinear functions in N variables X by a modification of the Powell hybrid method: FVEC (X) = 0.
-- How to use:
/*
INFREV: IN, OUT, INTEGER, reverse communication parameter INFREV is an integer variable hat must be set by the user to a non-positive value initially, in order to indicate the start of the calculation. The value of INFREV on each return indicates the reason for the return. 
If it is 1,3 or 4, then the user's program has to provide new values of FVEC, and then recall the subroutine with INFREV still set to this number. 
If it is 2, then the user's program has to provide new values of FVEC and FJAC, and then recall the subroutine with INFREV set to 2.
         INFREV on return
            = 0: Solver terminated
            = 1: Calculate FVEC at starting point
            = 2: Calculate FVEC and FJAC (only if IOPT=1)
            = 3: Calculate FVEC to compute Jacobian approximation
            = 4: Calculate FVEC to compute next step)
 IOPT  : IN, INTEGER, selection parameter
         IOPT specifies how the Jacobian will be calculated.
         1: the user must supply the Jacobian.
         2: the code will approximate the Jacobian by forward-differencing.
         3: same as IOPT=2, but for the first call (INVREF<=0) the user must supply the Jacobian.
 N     : IN, INTEGER, N.gt.0
         Number of functions and variables.
 X     : IN, OUT, DOUBLE(N)
         On input X must contain an initial estimate of the solution vector.
         On output X contains the values for which FVEC has to be computed in the next iteration step (INFREV=1,2,3 or 4) or the final estimate of the solution vector (INFREV=0).

 FVEC  : IN, DOUBLE(N)
        FVEC has to contain the functions evaluated at X  (INFREV=1,2,3,4).
 FJAC  : IN, DOUBLE(N,N)
         FJAC has to contain the Jakobian (INFREV=2) The Jacobian has to be provided only if IOPT=1.
         If IOPT=2 or 3, FJAC is used as a work array.  It contains the approximation of the Jacobian.
 TOL   : IN, DOUBLE, TOL.ge.0.d0
         Termination occurs when the algorithm estimates that the relative error between X and the solution is at most TOL.
         If TOL=0.D0 is given, a tolerance of sqrt(machine precision) is assumed.
 INFO  : OUT, INTEGER
		 IF INFREV > 0 THEN  Iteration continues
					INFO =  0  Internal step
						 = 10  Internal step completed
		  ELSEIF INFREV = 0 THEN Iteration terminated. Termination status:
					INFO =  0  Improper input parameters (wrong initialization)
					INFO =  1  Algorithm estimates that the relative error between X and the solution is at most TOL  (regular end).
					INFO =  2  Number of calls to FCN has reached or exceeded 100(N+1) for IOPT=1 or 200(N+1) for IOPT=2.
					INFO =  3  TOL is too small.  No further improvement in the approximate solution X is possible.
					INFO =  4  Iteration is not making good progress.
		 ENDIF
 DUM   : OUT, DOUBLE( (N2+13N)/2 )
         work array of length LDUM = 20 + (N2+13N)/2
 LDUM  : IN, INTEGER
         Provided length of the work array DUM. LDUM must not be less than 20 + (N2+13N)/2.
 IDUM  : IN, INTEGER
         work array of length 20.
 EPSFCN: IN, DOUBLE
         Used in determining a suitable step for the forward-difference approximation. This approximation assumes that the relative errors in the functions are of the order of epsfcn. If epsfcn  is less than the machine precision, it is assumed that the relative errors in the functions are of the order of the  machine precision. If iopt=1, then epsfcn can be ignored (treat it as a dummy argument).
*/

"FORTRAN" FUNCTION NO_TYPE EcoNonLinearEqnSolver (     
    OUT INTEGER ILAST,
    IN INTEGER IOPT,
    IN INTEGER N,
    OUT REAL X[],
    OUT REAL FVEC[],
    OUT REAL FJAC[],
    IN REAL tol,
    OUT INTEGER INFO,
    IN REAL DUM[],
    IN INTEGER LDUM,
    IN INTEGER IDUM[],
    IN REAL EPSFCN) IN "MATH.lib"

	
	
--------------------------------------------------------------------------------
"FORTRAN" FUNCTION NO_TYPE EcoMachDep () 
           IN "MATH.lib"