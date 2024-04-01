/*-----------------------------------------------------------------------------------------
 LIBRARY: SPARSE
 FILE: LargeLinearSystem1
 CREATION DATE: 13/01/2022
-----------------------------------------------------------------------------------------*/
USE SPARSE
COMPONENT TestLinearSystem02(INTEGER n = 2000)
--Componente que genera un gran sistema de ecuaciones lineales
--para probar la resolución de ecuaciones lineales
--La dimensión del sistema de ecuaciones generado es n-2
--El vector solución es 
--		x[j] = TIME * 0.1
--
	DATA
		REAL alpha = 10
	DECLS
		REAL x[n], y
	OBJECTS
			LEQ_SYSTEM(NEQN = n-2, NNZEROS = 3*n-8) leq_system
	INIT

	CONTINUOUS
		x[1] = 0.1*TIME
		x[n] = 0.1*TIME
		SEQUENTIAL
			leq_system.SetToZero()
			FOR (i IN 2, n-1)
				IF(i < 3) THEN
					leq_system.Fill_A(i-1,i-1,alpha)
					leq_system.Fill_A(i-1, i, 1.)
					leq_system.Fill_B(i-1, 1.2*TIME-x[1])
				ELSEIF(i < n-1) THEN
					leq_system.Fill_A(i-1, i-2, 1.)
					leq_system.Fill_A(i-1, i-1, alpha)
					leq_system.Fill_A(i-1, i, 1)
					leq_system.Fill_B(i-1, 1.2*TIME)
				ELSE
					leq_system.Fill_A(i-1, i-2, 1.)
					leq_system.Fill_A(i-1, i-1, alpha)
					leq_system.Fill_B(i-1, 1.2*TIME-x[n])
				END IF
			END FOR
			leq_system.Solve()
			y = 1.
		END SEQUENTIAL
		EXPAND (i IN 2, n-1)
				x[i] = y * leq_system.GetUnknown(i-1) 

END COMPONENT
COMPONENT TestLinearSystem02_0010 IS_A TestLinearSystem02 
	DECLS 
		CLOSE n = 10
END COMPONENT

COMPONENT TestLinearSystem02_0100 IS_A TestLinearSystem02 
	DECLS 
		CLOSE n = 100
END COMPONENT

COMPONENT TestLinearSystem02_0500 IS_A TestLinearSystem02 
	DECLS 
		CLOSE n = 500
END COMPONENT

COMPONENT TestLinearSystem02_1000 IS_A TestLinearSystem02 
	DECLS 
		CLOSE n = 1000
END COMPONENT

