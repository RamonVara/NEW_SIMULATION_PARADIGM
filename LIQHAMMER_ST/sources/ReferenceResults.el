/*-----------------------------------------------------------------------------------------
 LIBRARY: LIQHAMMER_ST
 FILE: ReferenceResults
 CREATION DATE: 19/03/2019
-----------------------------------------------------------------------------------------*/
COMPONENT ReferenceResults(INTEGER size = 1)
	DATA 
		FILEPATH fileName = "my_data_file.txt"
		REAL gain[size] = 1 
		REAL bias[size] = 0
	DECLS
		BOOLEAN AlreadyRead = FALSE
		REAL yref[size]
	OBJECTS
		TABLE table_ref_results[size]					
   INIT
		--setMySilentMode(TRUE)
		IF(NOT AlreadyRead) THEN
			FOR(i IN 1, size)
				table_ref_results[i].readCols1D(fileName,1,i+1)
			END FOR
			AlreadyRead = TRUE
		END IF
	CONTINUOUS
		EXPAND(i IN 1, size)
			yref[i] = table_ref_results[i].interpd1D(TIME) * gain[i] + bias[i]
END COMPONENT