/*-----------------------------------------------------------------------------------------
 LIBRARY: TEARING
 FILE: Algebraic_Loop1
 CREATION DATE: 06/06/2023
-----------------------------------------------------------------------------------------*/
COMPONENT Algebraic_Loop1
   DECLS
       REAL q1
       REAL q2
       BOUND REAL To = 150
       BOUND REAL T2 = 50
       REAL T1
       REAL cond1
       REAL cond2
       TABLE_1D cond_vs_Tave = {{0., 150}, {20., 25}}
		 REAL residue
   CONTINUOUS
       cond1 = linearInterp1D(cond_vs_Tave, 0.5*(To+T1))
       cond2 = linearInterp1D(cond_vs_Tave, 0.5*(T1+T2))
       q1 = cond1 * (To - T1)
       q2 = cond2 * (T1 - T2)
       residue = q1 - q2
END COMPONENT
