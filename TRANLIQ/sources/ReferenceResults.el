COMPONENT ReferenceResults(INTEGER size = 12)
	DATA 
		STRING fileName = "chaudry_c.txt"
		REAL gain[size] = 1 
		REAL bias[size] = 0
	DECLS
		BOOLEAN AlreadyRead = FALSE
		REAL yref[size]
	OBJECTS
		TABLE table_ref_results[size]					
   INIT
		IF(fileName != "") THEN
			--setMySilentMode(TRUE)
			IF(NOT AlreadyRead) THEN
				FOR(i IN 1, size)
					table_ref_results[i].readCols1D(fileName,1,i+1)
				END FOR
				AlreadyRead = TRUE
			END IF
		END IF
	CONTINUOUS
		SEQUENTIAL
		IF(fileName != "") THEN
			FOR(i IN 1, size)
				yref[i] = table_ref_results[i].interpd1D(TIME) * gain[i] + bias[i]
			END FOR
		END IF
		END SEQUENTIAL
END COMPONENT


