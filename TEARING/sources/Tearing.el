/*-----------------------------------------------------------------------------------------
 LIBRARY: TEARING
 FILE: Tearing
 CREATION DATE: 05/06/2023
-----------------------------------------------------------------------------------------*/
FUNCTION REAL Residue(IN REAL x)
	BODY
		RETURN x
END FUNCTION

COMPONENT Test
	DATA
		REAL R1 = 1000
		REAL R2 = 2000.
	DECLS
		BOUND REAL V1 = 1
		BOUND REAL V3 = 0
		REAL V2
		REAL i1 = 0
		REAL i2 = 0
		REAL V2dummy
	CONTINUOUS
		i1 = (V1-V2)/R1
		i2 = (V2-V3)/R2
		Residue(V2) = EQUAL_OUT(i1 - i2)
		Residue(V2) = V2dummy
END COMPONENT