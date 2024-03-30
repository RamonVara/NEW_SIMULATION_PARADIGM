/*-----------------------------------------------------------------------------------------
 LIBRARY: SPARSE_EXAMPLES
 FILE: TestLinearSystem
 CREATION DATE: 22/11/2023
-----------------------------------------------------------------------------------------*/
USE SPARSE
COMPONENT TestLinearSystem01
	DECLS
		DISCR REAL x[3]
		--5*x +   y + 8*z = 46
		--4*x + 2*y       = 12
		--6*x + 7*y + 4*z = 50
		--Solution 
	OBJECTS
		LEQ_SYSTEM(NEQN = 3, NNZEROS = 8) leq_system
	INIT
		leq_system.SetToZero()
		leq_system.Fill_B(1, 46.)
		leq_system.Fill_B(2, 12.)	
		leq_system.Fill_B(3, 50.)
		leq_system.Fill_A(1, 1, 5.)
		leq_system.Fill_A(1, 2, 1.)	
		leq_system.Fill_A(1, 3, 8.)		
		leq_system.Fill_A(2, 1, 4.)
		leq_system.Fill_A(2, 2, 2.)	
		leq_system.Fill_A(3, 1, 6.)
		leq_system.Fill_A(3, 2, 7.)	
		leq_system.Fill_A(3, 3, 4.)	
		leq_system.Solve()
		x[1] = leq_system.GetUnknown(1)
		x[2] = leq_system.GetUnknown(2)
		x[3] = leq_system.GetUnknown(3)
END COMPONENT


