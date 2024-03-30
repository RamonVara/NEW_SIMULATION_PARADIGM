USE SPARSE
FUNCTION NO_TYPE fsystem(IN INTEGER iopt, IN REAL x[], IN REAL rdata[], OUT REAL F[], OUT NLEQN_SYSTEM_BASE nleqn_system)
	DECLS
	BODY
		nleqn_system.neqn = 4
		IF(iopt == 0) THEN
				FOR(j IN 1, nleqn_system.neqn)
		   		nleqn_system.CreateUnk(j, 0.1, -1e6, 1e6, concatStrings(concatStrings("x[", integerToString(j)),"]"),0)
				END FOR		
		   	nleqn_system.CreateEqn(1,"-x[1]**2 - x[2]**2 - x[3]**2 + x[4]    = 0")
				nleqn_system.CreateEqn(2," x[1]**2 + x[2]**2 + x[3]**2 + x[4]**2 = 0")
				nleqn_system.CreateEqn(3," x[1] - x[2] = 0")
				nleqn_system.CreateEqn(4," x[2] - x[3] = 0")
		ELSEIF(iopt == 1 OR iopt == 2) THEN 
			nleqn_system.Jacob.setToZero()		
	   	F[1] =  - x[1]**2 - x[2]**2 - x[3]**2 + x[4]
			F[2] =    x[1]**2 + x[2]**2 + x[3]**2 + x[4]**2 - 1
			F[3] =    x[1] - x[2]
			F[4] =    x[2] - x[3]
			IF(iopt == 2) THEN
				nleqn_system.Jacob.setToZero()
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
				nleqn_system.Jacob.stampItem(4, 2, 1.0)
				nleqn_system.Jacob.stampItem(4, 3, -1.0)
				nleqn_system.Jacob.buildCSR()
			END IF
		END IF
		RETURN
END FUNCTION 