/*-----------------------------------------------------------------------------------------
 LIBRARY: SPARSE
 FILE: sparse.el
 CREATION DATE: 30/01/2020
-----------------------------------------------------------------------------------------*/
USE MATH
"C" FUNCTION  NO_TYPE d_ecosim_superlu( 
/* =========================================================================
* FUNCTION: d_ecosim_superlu
* This function is a wrapper for the SuperLU solver in EcosimPro Language.
* ==========================================================================
*
* iopt             = int specifies the operation:
*                      1, performs LU decomposition for the first time
*                      2, performs triangular solve
*                      3, free all the storage in the end
* n                = dimension of the square sparse matrix
* nnz              = number of nonzeros in the sparse matrix
* nrhs             = number of right-hand sides
* val              = double array containing the nonzero entries 
* colind           = int array containing the column indices of the entries 
* rowptr           = int array containing the row start
* b                = double array containing the right-hand side vector (gets overwritten with solution)
* factors          = pointer to LU factors. (If iopt == 1, it is an output and contains the pointer pointing 
*                    to the structure of the factored matrices. Otherwise, it is an input.)
* info             = info flag from SuperLU
*                  = 0: successful exit   
*                  < 0: if info = -i, the i-th argument had an illegal value   
*                  > 0: if info = i, and i is   
*                  <= A->ncol: U(i,i) is exactly zero. The factorization has   
*                       been completed, but the factor U is exactly   
*                       singular, so the solution and error bounds   
*                       could not be computed.   
*                  = A->ncol+1: U is nonsingular, but RCOND is less than machine
*                       precision, meaning that the matrix is singular to
*                       working precision. Nevertheless, the solution and
*                       error bounds are computed because there are a number
*                       of situations where the computed solution can be more
*                       accurate than the value of RCOND would suggest.   
*                  > A->ncol+1: number of bytes allocated when memory allocation
*                       failure occurred, plus A->ncol.
* pivot_threshold  = pivot threshold
* equil            = 1/0 equilibrate matrix/no equilibrate matrix (not implemented yet)
* refine           = 1/0 refine solution/ no refine solution (not implemented yet)
* debug            = 1/0 for debug&performance info/no debug&performance info 
* =========================================================================	*/
	IN INTEGER iopt         "integer that specifies the operation to be performed",
	IN INTEGER n			   "dimension of the square sparse matrix",
	IN INTEGER nnz			   "number of non-zeros entries in the sparse matrix",
	IN INTEGER nrhs 		   "number of righ hand side terms",
	IN REAL val[]			   "values of the non zero items in row compressed order",
	IN INTEGER colwind[]    "column indexes in row compressed format",
	IN INTEGER rowptr[]     "row pointers in row compressed format",
	OUT REAL b[]            "double array containing the right hand side",
	OUT INTEGER factors[4]  "handle pointing to the LU factors",
	OUT INTEGER info        "return flag from superlu, if =0 successful exit",
	IN REAL pivot_threshold "pivot threshold",
	IN INTEGER equil        "1/0 equilibrate matrix/no equilibrate matrix (not implemented yet)",
	IN INTEGER refine 	   "1/0 refine solution/ no refine solution (not implemented yet)",
	IN INTEGER debug        "1/0 for debug&performance info/no debug&performance info" 
	
	) IN  "libsuperlu_5.1.a"



CLASS SPARSE_MATRIX(
		INTEGER MAX_ROW		"Maximum number of rows and columns", 
		INTEGER MAX_STAMP		"Maximum number of stamps to build the Jacobian matrix"
		)
/* =========================================================================
* CLASS: SPARSE_MATRIX
* Purpose: Filling non-zero items by stampation, 
*          Building the CSR representation of the matrix from the stamps 
*          Solving linear equation systems
*			  Performing diverse operations with the matrix: 
*              -Set to zero all the elements of the matrix
*   	         -Get an array with maximum absolute value of each row
*              -Pre-multiply the matrix by a row vector
*              -Post-multiply the matrix by a column vector
*              -Print the matrix
* ==========================================================================
--Example of Usage of the  CLASS SPARSE MATRIX:
USE SPARSE
COMPONENT TestSparseMatrix01
	DATA
		REAL b[5] = { 10.0,  0.25,  4.50, -0.25, 13.5} 
		--independent terms of the linear equation system
		REAL c[5] = {  3.0,  0.50, -0.50,  2.50, -1.0}
		--1D array
	DECLS			
		DISCR REAL x[5] --Solution of the linear equation system
		DISCR REAL m[5] --Maximum absolute values per row
		DISCR REAL d[5] --Result of postmultiplication
		DISCR REAL e[5] --Result of premultiplication
		DISCR REAL error
	OBJECTS
		--Declaration of a sparse matrix with up to 5 rows
		--and up to 20 possible stampings
   	SPARSE_MATRIX(MAX_ROW=10, MAX_STAMP=20)  A 
	INIT
		A.setToZero()
		-- Let the matrix be 5x5 with the values shown below:
		-- A = [
      --        8,   0,   0,   0,   1;
      --        0, 1+2,   0, 4+1,   0;
      --        1,   0,   3,   0,   0;
      --        0,   1,   0,   3,   0;
      --        3,   0,   1,   0, 2-6;
      --     ];
		--Fill the matrix by stamping as follows:
		A.stampItem(1, 1, 8.)
		A.stampItem(2, 2, 1.)
		A.stampItem(2, 2, 2.)
		A.stampItem(3, 3, 3.)
		A.stampItem(4, 4, 3.)
		A.stampItem(5, 5, 2.)
		A.stampItem(5, 5,-6.)
		A.stampItem(1, 5, 1.)
		A.stampItem(2, 4, 4.)
		A.stampItem(2, 4, 1.)
		A.stampItem(3, 1, 1.)
		A.stampItem(4, 2, 1.)
		A.stampItem(5, 1, 3.)
		A.stampItem(5, 3, 1.)
		--
		--After stamping
		--nstamps    = 14
		--stampIrow  = {  1,  2,  2,  3,  4,  5,  5,  1,  2,  2,  3,  4,  5,  5}
		--stampIcol  = {  1,  2,  2,  3,  4,  5,  5,  5,  4,  4,  1,  2,  1,  3}
		--stampIpos  = {  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14}
		--stampValue = { 8., 1., 2., 3., 3., 2.,-6., 1., 4., 1., 1., 1., 3., 1.}
		A.Print()
		--------------------------------------------------------
	   --Construcción de la matriz en formato CSR
	   A.buildCSR()
      --After this method, the stamps are sorted
		--stampIrow    = { 1,   1,   2,   2,   2,   2,   3,   3,   4,   4,   5,   5,   5,   5},
      --stampIcol    = { 1,   5,   2,   2,   4,   4,   1,   3,   2,   4,   1,   3,   5,   5},
      --stampIpos    = { 1,   8,   2,   3,  10,   9,  11,   4,  12,   5,  13,  14,   6,   7}
      --stampValue   = {8.,  1.,  1.,  2.,  1.,  4.,  1.,  3.,  1.,  3.,  3.,  1.,  2., -6.}
		--The CSR format is built
		--nnzeros = 11 
		--csrColInd = {1,    5,   2,   4,   1,   3,   2,   4,   1,   3,   5}
		--csrValue  = {8.,  1.,  3.,  5.,  1.,  3.,  1.,  3.,  3.,  1., -4.}
		--csrRowPtr = {1,    3,   5,   7,   9,  12}
       --The following array indicates the destination position in CSR format of the unsorted stamps
       --stampToCSR = {  1,  3,   3,   6,   8,  11,  11,   2,   4,   4,   5,   7,   9,  10}
		A.Print()
		--------------------------------------------------------
		--Solving the system of linear equations A·x = b
		--where the vector of independent terms
		--------------------------------------------------------
		--b = [
      --      10.0;
      --      0.25;
      --      4.50;
      --     -0.25;
      --      13.5;
      --    ]			
	   A.LinearEqnSolve(b, x, 1., 0)		
		WRITE("\n\nSolution:\n")
		WRITE("    x[5]   = { %lg , %lg , %lg , %lg , %lg }\n", 
		                      x[1], x[2], x[3], x[4], x[5])
		--------------------------------------------------------
		--Find the maximum absolute value of each row
		A.MaxAbsRow(m)
      WRITE("\n\nMax.Abs.Value of each row:\n")
		WRITE("    m[5]   = { %lg , %lg , %lg , %lg , %lg }\n",
		                      m[1], m[2], m[3], m[4], m[5]) 
		--------------------------------------------------------
		--Pre-multiplication by a row vector
		A.preMultiplyByRow(c, e)
		WRITE("\n\nResult of Premultiplicacion:\n")
		WRITE("    e[5]   = { %lg , %lg , %lg , %lg , %lg }\n", 
		                      e[1], e[2], e[3], e[4], e[5])
		--------------------------------------------------------
		---Postmulitplicacion by column vector
		A.postMultiplyByColumn(c, d)
		WRITE("\n\nResult of Post-multiplicacion:\n")
		WRITE("    d[5]   = { %lg , %lg , %lg , %lg , %lg }\n", 
		                      d[1], d[2], d[3], d[4], d[5])

END COMPONENT*/
	DECLS
		BOOLEAN isCSRBuilt = FALSE  "boolean indicating whether the sparse matrix \
		                             has been built in CSR format or not"
		BOOLEAN isSetToZero = TRUE  "boolean indicating if the sparse matrix \
		                             has been reset to zero or not"
		--Stamps representation
		INTEGER nstamp = 0	       "Total number of stamps"
		INTEGER nrow = 0		       "Total number of rows"
		INTEGER ncol = 0		       "Total number of columns"
		INTEGER stampIrow[MAX_STAMP]  = 0	 "Number of row of the stamp"
		INTEGER stampIcol[MAX_STAMP]  = 0	 "Number of column of the stamp"
		INTEGER stampIpos[MAX_STAMP]  = 0	 "Position of the stamp before sorting"
		REAL    stampValue[MAX_STAMP] = 0.   "Real value of the stamp"	
		INTEGER stampToCSR[MAX_STAMP] = 0    "Mapping of stamp contributions in CSR format"
		INTEGER STAMP_debug			          "Debug level of the stamping process (0/1)"
      --CSR representation of the matrix
		INTEGER nnzeros = 0	"Total number of non zero items"
		INTEGER csrColInd[MAX_STAMP]  = 0	 "Columns of the matrix items in CSR format"
		INTEGER csrRowPtr[MAX_ROW+1]  = 0	 "Pointers to the beginning of each row and to the end of last row"
		REAL    csrValue[MAX_STAMP]   = 0.   "Real Values in CSR format"
	METHODS
      ----------------------------------------------------------------------------
		--Method to set to zero the sparse matrix
      ----------------------------------------------------------------------------
		METHOD NO_TYPE setToZero()
			BODY
				isSetToZero = TRUE
				nstamp = 0
				FOR(i IN 1, nnzeros)
					csrValue[i] = 0
				END FOR
				RETURN
		END METHOD
      ----------------------------------------------------------------------------
		--Method to swap two stamps - Auxiliary method for the Quicksort
		----------------------------------------------------------------------------
		METHOD NO_TYPE SwapItems(
			IN INTEGER i "Position of the first stamp", 
			IN INTEGER j "Position of the second stamp")
			DECLS
				INTEGER swap[3]			
				REAL v
			BODY
				swap[1] = stampIrow[i]
				swap[2] = stampIcol[i]
				swap[3] = stampIpos[i]
				v = stampValue[i]
            --   
				stampIrow[i] = stampIrow[j]
				stampIcol[i] = stampIcol[j]
				stampIpos[i] = stampIpos[j]
				stampValue[i] = stampValue[j]				
		      --
				stampIrow[j] = swap[1]
				stampIcol[j] = swap[2]	
				stampIpos[j] = swap[3]
				stampValue[j] = v
				RETURN
		END METHOD
      ----------------------------------------------------------------------------
		--Method two compare two stamp positions of the matrix to find which one 
		--is first in the CSR format		
      ----------------------------------------------------------------------------
		METHOD BOOLEAN Compare(
			IN INTEGER irow1 "row index of first stamp", 
			IN INTEGER icol1 "column index of first stamp", 
			IN INTEGER irow2 "row index of second stamp", 
			IN INTEGER icol2 "column index of second stamp"
			)
			DECLS
				BOOLEAN bval = FALSE
			BODY
				IF (irow1 < irow2  OR (irow1 == irow2 AND icol1 < icol2)) THEN
					bval = TRUE
				END IF
				RETURN bval
		END METHOD
      ----------------------------------------------------------------------------
		--Partition method called by the Quicksort method
		----------------------------------------------------------------------------
		METHOD INTEGER partition(
			IN INTEGER lo "Start index of the partition",
			IN INTEGER hi "Ending index of the partition"
			)
			DECLS
				INTEGER pivot[2]
				INTEGER i
			BODY
				pivot[1] = stampIrow[hi]
				pivot[2] = stampIcol[hi]
				i = lo 
				FOR (j IN lo, hi)
					IF Compare(stampIrow[j], stampIcol[j], pivot[1], pivot[2]) THEN 
						SwapItems(i, j)
						i = i + 1
					END IF
				END FOR
				SwapItems(i, hi)
				RETURN i
		END METHOD		
      ----------------------------------------------------------------------------		
		--QuickSort Algorithm to sort the sparse matrix stamps in CSR format
      ----------------------------------------------------------------------------
		METHOD NO_TYPE QuickSort(
				IN INTEGER lo  "the starting index of the stamp array or the current sub-array", 
				IN INTEGER hi  "the ending index of the stamp array or the current sub-array"
				)
			DECLS
				INTEGER p
			BODY
				IF(lo < hi) THEN
					p = partition(lo, hi)
					QuickSort(lo,  p - 1)
					QuickSort(p + 1, hi)
				END IF
				RETURN
		END METHOD
      ----------------------------------------------------------------------------		
		--Method to build the matrix in CSR format from the unsorted stamps
      ----------------------------------------------------------------------------		
		METHOD NO_TYPE buildCSR()
			DECLS
				INTEGER i
				INTEGER k
				INTEGER i1, i2
			BODY
				isSetToZero = FALSE
				--It checks if the CSR sorting is still valid
				IF(NOT isCSRBuilt) THEN
					ASSERT (nrow == ncol) ERROR  "Non square sparse matrix"
					--Sorting of the stamps for CSR format using QuickSort
					WRITE("=======================\n")
					WRITE("=======================\n")
					WRITE("Ordenación de la Matriz\n")
					WRITE("=======================\n")
					WRITE("=======================\n")
					QuickSort(1, nstamp)
					--Once the stamps have been sorted, the CSR representation is built
					nnzeros = 0
					--
					FOR (i IN 1, nstamp)
						--It detects if the position of the actual stamp is equal to the position 
						--of the previous stamp
						IF(i > 1 AND (stampIrow[i] == stampIrow[i-1]) AND (stampIcol[i] == stampIcol[i-1])) THEN
							csrValue[nnzeros] = csrValue[nnzeros] + stampValue[i] 
						ELSE
							nnzeros = nnzeros + 1
							csrColInd[nnzeros] = stampIcol[i]
							csrValue[nnzeros] = stampValue[i]
							IF(i == 1) THEN	
								i1 = 1
					         i2 = stampIrow[i]
								ASSERT(i1 == i2) ERROR "Empty row in sparse matrix"
								FOR(k IN i1, i2)
								   csrRowPtr[k]= nnzeros
								END FOR	
							ELSEIF(stampIrow[i] != stampIrow[i-1]) THEN
								i1 = stampIrow[i-1]+1
								i2 = stampIrow[i]
								ASSERT(i1 == i2) ERROR "Empty row in sparse matrix"
								FOR(k IN i1, i2)
								   csrRowPtr[k]= nnzeros							
								END FOR								
							END IF
						END IF
						k = stampIpos[i]
						stampToCSR[k] = nnzeros
					END FOR
					csrRowPtr[nrow+1] = nnzeros + 1
					isCSRBuilt = TRUE
				END IF
				RETURN
	   END METHOD
		----------------------------------------------------------------------------
	   --Method to stamp values in any order
		----------------------------------------------------------------------------
	   METHOD NO_TYPE stampItem(
		   IN INTEGER ir	 	"Row index of the stamp", 
		   IN INTEGER ic		"Column index of the stamp", 
		   IN REAL rval		 "Additive value contribution of the stamp"		        
         )
			DECLS
				INTEGER izero
		   BODY
				ASSERT isSetToZero FATAL "Before stamping items the matrix has to be zeroed"
			   nstamp = nstamp + 1
			   ASSERT (nstamp <= MAX_STAMP) FATAL "Maximum number of SPARSE_MATRIX Stamps exceeded. \
			                                Create the SPARSE_MATRIX object with increased MAX_STAMP"

				IF(ir > nrow) THEN
					   ASSERT (ir <= MAX_ROW) FATAL "Maximum number of SPARSE_MATRIX rows exceeded. \
						                       Create the SPARSE_MATRIX object with increased MAX_ROW"
					   nrow = ir
				END IF
				IF(ic > ncol) THEN
					ASSERT (ic <= MAX_ROW) FATAL "Maximum number of SPARSE_MATRIX columns exceeded. \
  	                                   Create the SPARSE_MATRIX object with increased MAX_ROW"
					   ncol = ic
				END IF
				stampIpos[nstamp] = nstamp
				stampIrow[nstamp] = ir
				stampIcol[nstamp] = ic
				stampValue[nstamp] = rval
				IF(STAMP_debug > 0) THEN
					WRITE("Stamp %d into Matrix Position[%d,%d] = %lg \n", nstamp, ir, ic, rval)
				END IF
				--Tty to use previous sorting if the CSR format has been already built
			   IF(isCSRBuilt) THEN
					izero = stampToCSR[nstamp]
					IF(ic == csrColInd[izero] AND izero >= csrRowPtr[ir] AND izero < csrRowPtr[ir+1]) THEN
				   	csrValue[izero] = csrValue[izero] + rval
					ELSE
						isCSRBuilt = FALSE
					END IF
				END IF
			   RETURN
      END METHOD
      ----------------------------------------------------------------------------
		--Method to solve linear equation system
		----------------------------------------------------------------------------
	   METHOD NO_TYPE LinearEqnSolve(
	      IN REAL b[]   "independent terms", 
         OUT REAL x[]  "solution vector",
			IN REAL pivot_threshold "pivot threshold",
			IN INTEGER superlu_debug    "0/1 Debug level for SuperLU"
		   )
		   DECLS
			   INTEGER factors[4]
			   INTEGER info
			   INTEGER superlu_equil =  1
			   INTEGER superlu_refine = 0
		   BODY		
			   ASSERT(isCSRBuilt) FATAL "Method LinearEqnSolve of class SPARSE_MATRIX \
			                          can be called only if CSR format is built" 
			   FOR(i IN 1, nrow)
				   x[i] = b[i]
			   END FOR			
			   --Factorize Matrix into L & U
			   d_ecosim_superlu(1, nrow, nnzeros, 1, csrValue, csrColInd, csrRowPtr, x, factors, info, \
			                 pivot_threshold, superlu_equil, superlu_refine, superlu_debug)

			   IF (info != 0) THEN
         	   WRITE("\n****Error in SUPERLU during Factorization info = %d\n\n", info)
			   END IF			
			   --Solve Linear System
			   d_ecosim_superlu(2, nrow, nnzeros, 1, csrValue, csrColInd, csrRowPtr, x, factors, info, 
			                 pivot_threshold, superlu_equil, superlu_refine, superlu_debug)	
			   IF (info != 0) THEN
         	   WRITE("\n****Error in SUPERLU during Solve info = %d\n\n", info)
			   END IF			
			   --Destroy Matrix
			   d_ecosim_superlu(3, nrow, nnzeros, 1, csrValue, csrColInd, csrRowPtr, x, factors, info, 
			                 pivot_threshold, superlu_equil, superlu_refine, superlu_debug)	
			   IF (info != 0) THEN
         	   WRITE("\n****Error in SUPERLU during Memory deallocation info = %d\n\n", info)
			   END IF	
			   RETURN
      END METHOD
	   ----------------------------------------------------------------------------		
		--Get an array with maximum absolute value of each row
		----------------------------------------------------------------------------
		METHOD NO_TYPE MaxAbsRow(
		   OUT REAL rmax[] "Vector with the absolute maximum of each row"
			)
		   DECLS
		      INTEGER i
				INTEGER j
				INTEGER k, k1, k2
			BODY
			   ASSERT(isCSRBuilt) FATAL "Method MaxAbsRow of class SPARSE_MATRIX \
			                        can be called only if CSR format is built" 
				FOR(i IN 1, nrow)
				   k1 = csrRowPtr[i]
				   k2 = csrRowPtr[i+1]-1
				   rmax[i] = 0.
				   FOR(k IN k1, k2)
					   rmax[i] = max(rmax[i], abs(csrValue[k]))
				   END FOR
				   IF(rmax[i] == 0) THEN
					   WRITE("****WARNING - Singular or Almost Singular Jacobian\n")
				   END IF
			   END FOR
				RETURN
      END METHOD
		----------------------------------------------------------------------------
      -- Method to pre-multiply a row vector by the matrix
      ----------------------------------------------------------------------------
      METHOD NO_TYPE preMultiplyByRow(
         IN REAL X[] "Input Row vector",
         OUT REAL Y[] "Results vector")
         DECLS
            INTEGER i "Loop index for rows"
            INTEGER j "Loop index for columns"
            INTEGER k, k1, k2 "Loop indices for matrix elements"
         BODY
         -- Check if CSR format is built before proceeding
         ASSERT(isCSRBuilt) FATAL "Method preMultiplyByRow of class SPARSE_MATRIX \
                                  can be called only if CSR format is built"
         -- Initialize the result vector to zero
         FOR(j IN 1, ncol)
            Y[j] = 0.
         END FOR
         -- Perform matrix-vector multiplication
         FOR(i IN 1, nrow)
            k1 = csrRowPtr[i]
            k2 = csrRowPtr[i+1]-1
            FOR(k IN k1, k2) 
                j = csrColInd[k]
                -- Accumulate the product of matrix element and vector element
                Y[j] = Y[j] + csrValue[k] * X[i]
            END FOR
         END FOR
         -- Return the result vector
         RETURN
      END METHOD
		----------------------------------------------------------------------------
		--Method to multiply the matrix by a column vector		
		----------------------------------------------------------------------------
		METHOD NO_TYPE postMultiplyByColumn(
				IN REAL X[] "Input Column vector", 
				OUT REAL Y[] "Results vector")
			DECLS
				INTEGER i  "Loop index for rows"
				INTEGER j  "Column index"
				INTEGER k, k1, k2 "Loop indices for matrix elements"
			BODY
			   -- Check if CSR format is built before proceeding
			   ASSERT(isCSRBuilt) FATAL "Method postMultiplyByColumn of class SPARSE_MATRIX \
			                          can be called only if CSR format is built" 
		      -- Initialize the result vector to zero
				FOR(i IN 1, nrow)
					Y[i] = 0.
				END FOR
            -- Perform matrix-vector multiplication
				FOR(i IN 1, nrow)
					k1 = csrRowPtr[i]
					k2 = csrRowPtr[i+1]-1
					FOR(k IN k1, k2) 
						j = csrColInd[k]
                  -- Accumulate the product of matrix element and vector element
						Y[i] = Y[i] + csrValue[k]*X[j]
					END FOR
				END FOR
				-- Return the result vector
				RETURN
		END METHOD		
		----------------------------------------------------------------------------
		--Method to print the sparse matrix
		----------------------------------------------------------------------------		
		METHOD NO_TYPE Print()
			DECLS
			BODY
				IF(isCSRBuilt) THEN
					WRITE("\n\nnnzeros = %10d\n", nnzeros)
					WRITE("nrow     = %10d\n", nrow)
					WRITE("ncol = %10d\n\n", ncol)
					FOR(j IN 1, nnzeros)
						WRITE("csrColInd[%d]= %d \t val[%d] = %14.7g\n", j, csrColInd[j], j, csrValue[j])
					END FOR
					WRITE("\n")
					FOR(j IN 1, nrow+1)
						WRITE("csrRowPtr[%d] = %d \n", j, csrRowPtr[j])
					END FOR
				ELSE
					WRITE("\nnstamp = %10d\n", nstamp)
					WRITE("   Row      Column         Value         Ipos    Ipos CSR\n")
					WRITE("  ======    ======    ==============    ======    ======\n")
					FOR(j IN 1, nstamp)
						WRITE("  %6d    %6d    %14.7lg    %6d    %6d\n", stampIrow[j], stampIcol[j], stampValue[j], stampIpos[j], stampToCSR[j])
					END FOR
				END IF
				RETURN
		END METHOD
END CLASS

CLASS NLEQN_SYSTEM_BASE (
		INTEGER MAX_EQN		"Maximum number of equations", 
		INTEGER MAX_STAMP		"Maximum number of stamps to build the Jacobian matrix"
	) 
	DECLS
		INTEGER neqn = 0							"Number of equations"
		INTEGER iter = 0							"Iteration number of the NR"
		INTEGER NR_debug = 1						"Debug level of the Newton Raphson- range 0/1/2"
		INTEGER STAMP_debug = 0					"Debug level of the stamping process range 0/1" 
		INTEGER superLU_debug = 0				"Debug level of the superLU solver 0/1" 
		REAL Fun[MAX_EQN] = 0.					"Vector of new residues"
		REAL FunOld[MAX_EQN] = 0.				"Vector of old residues"
		REAL Xunk[MAX_EQN] = 0.					"Vector of new values of the unknowns"
		REAL XunkOld[MAX_EQN] = 0.				"Vector of old values of the unknowns"
		REAL XunkMax[MAX_EQN] =  1.e10		"Vector of maximum limits for the unknowns"
		REAL XunkMin[MAX_EQN] = -1.e10		"Vector of minimum limits for the unknowns"
		INTEGER freezeIterUnk[MAX_EQN]=0    "Number of freezed iterations for the unknowns"
		STRING descrUnk[MAX_EQN]		      "String with a description of the unknowns"
		STRING descrEqn[MAX_EQN]		      "String with a description of the equation"
		REAL pivot_threshold = 1.	         "Specifies the threshold used for a diagontal \
		                                     entry to be an acceptable pivot"														 
		REAL tol	= 1.e-7							"Tolerance"										
		BOOLEAN converged = FALSE				"Boolean, if TRUE, the solver has converged"
		CONST REAL MinFracLimit = 0.1
	OBJECTS
		SPARSE_MATRIX (MAX_ROW = MAX_EQN, MAX_STAMP = MAX_STAMP) Jacob		"Jacobian matrix"
	METHODS
		----------------------------------------------------------------------------
		-- Method to set to zero the residues of the functions and the Jacobian
		----------------------------------------------------------------------------
	   METHOD NO_TYPE setToZero()
		   BODY
			   --nstamp = 0
			   FOR (i IN 1, neqn)
				   Fun[i] = 0
			   END FOR
			   Jacob.setToZero()
			   neqn = 0			
			RETURN
	   END METHOD
		----------------------------------------------------------------------------
		--Method to create an unknown
		----------------------------------------------------------------------------		
		METHOD NO_TYPE CreateUnk(
				IN INTEGER i "Unknown number", 
				IN REAL Xo	 "Initial value of the unknown", 
				IN REAL Xmin "Minimum allowable value for the unknown", 
				IN REAL Xmax "Maximum allowable value for the unknown", 
				IN STRING unkn_descrip "Description of the unknown variable", 
				IN INTEGER freezeIter "Number of iterations to freeze the unknown "
				)
			BODY		
				Xunk[i] = Xo
				IF(Xmin != 0. OR Xmax != 0.) THEN
					XunkMax[i] = Xmax
					XunkMin[i] = Xmin
				ELSE
					XunkMax[i] =  1e10
					XunkMin[i] = -1e10
				END IF
				ASSERT(Xunk[i] >= XunkMin[i] AND Xunk[i] <= XunkMax[i]) ERROR "Initial value of unknown is out of allowable range"
				descrUnk[i] = unkn_descrip
				freezeIterUnk[i] = freezeIter
				RETURN
		END METHOD
		----------------------------------------------------------------------------
		--Method to create the comment or description of an equation
		----------------------------------------------------------------------------		
		METHOD NO_TYPE CreateEqn(
				IN INTEGER i "Equation number" , 
				IN STRING eqn_descrip "Description of the equation")
			BODY		
				descrEqn[i] = eqn_descrip
				RETURN
		END METHOD	
		----------------------------------------------------------------------------
	   --Method to update the unknowns after each iteration.
		--It takes into account the maximum and minimum limit on each variable
		----------------------------------------------------------------------------		
		METHOD NO_TYPE UpdateUnknownVec(
			IN REAL FRAC,
			IN REAL XunkOld[], 
			IN REAL dx_nr[], 
			--IN REAL dx_sp[],
			OUT REAL Xunk[],
			OUT REAL errorx[],
			OUT REAL errorx_tot,
			OUT INTEGER iworstx
			)
		   DECLS
		   	REAL fracLimit
			   REAL errorx_max
	         REAL MinFracLimit = 0.5			
		   BODY
				--Calculate the fractional limit of the Newton-Raphson step
			   FOR(j IN 1, neqn)
				   fracLimit = FRAC
               -- Adjust the limit to ensure that neither the maximum nor the minimum
					-- limits are exceeded
				   IF(dx_nr[j] > 0.) THEN
					   fracLimit = 0.99 * (XunkMax[j] - XunkOld[j])/dx_nr[j]
				   ELSEIF(dx_nr[j] < 0.) THEN
					   fracLimit = 0.99 * (XunkMin[j] - XunkOld[j])/dx_nr[j]
				   END IF
				   IF(fracLimit < MinFracLimit) THEN
					   fracLimit = MinFracLimit
				   END IF
				   IF(fracLimit < FRAC) THEN
					   FRAC = fracLimit
				   END IF
			   END FOR

			   FOR(j IN 1, neqn)
				   IF(iter >= freezeIterUnk[j]) THEN
					   Xunk[j] = min(max(XunkOld[j] + FRAC * dx_nr[j], XunkMin[j]), XunkMax[j])
				   ELSE
					   Xunk[j] = Xunk[j]
				   END IF	
			   END FOR
			   errorx_tot = 0.
			   errorx_max = 0.
			   iworstx = 0
			   FOR (j IN 1, neqn)
				   IF(abs(Xunk[j]) > 10*tol) THEN 
					   errorx[j] = (Xunk[j]-XunkOld[j])/abs(XunkOld[j])
				   ELSE
					   errorx[j] =  Xunk[j]-XunkOld[j]
				   END IF
				   IF(abs(errorx[j]) > errorx_max) THEN
					   iworstx = j
					   errorx_max = abs(errorx[j])
				   END IF
					--Accumulate the total variable error
				   errorx_tot = errorx[j]**2 + errorx_tot
			   END FOR
			   errorx_tot = errorx_tot**0.5
			   RETURN
	   END METHOD
		----------------------------------------------------------------------------
		-- Method to evaluate the error of the residues
		-- Errors are scaled by dividing each by the maximum absolute value
		-- in the corresponding row of the Jacobian matrix.
		-- A total error is calculated by composing individual errors
		-- using the Euclidean norm.
		-- Additionally, the worst-scaled residue is identified.
		----------------------------------------------------------------------------
	   METHOD NO_TYPE EvalFunError(
				IN REAL f[],
				IN REAL df[],
				OUT REAL errorf[], 
				OUT REAL errorf_tot, 
				OUT INTEGER iworstf)
		   DECLS
			   REAL errorf_max
		   BODY
				errorf_tot = 0.
				errorf_max = 0.
				iworstf = 0
				FOR(j IN 1, neqn)
				   -- Calculate the scaled error for each equation
					errorf[j] = f[j]/df[j]
					-- Update the worst-scaled residue if needed
					IF(abs(errorf[j]) > errorf_max) THEN
						iworstf = j
						errorf_max = abs(errorf[j])
					END IF
					-- Accumulate the total error
			   	errorf_tot = errorf_tot + errorf[j]**2
				END FOR
				-- Calculate the Euclidean norm of the total error
				errorf_tot = errorf_tot**0.5
				RETURN
      END METHOD
		METHOD NO_TYPE GetUnknown(OUT REAL X[])
		   BODY
			   FOR (i IN 1, neqn)
				   X[i] = Xunk[i]
			   END FOR
			   RETURN
	   END METHOD
END CLASS


TYPEDEF FUNCTION NO_TYPE func_Fsystem(
/*=========================================================================
--Declaration of the function pointer to calculate the residues and
--the Jacobian from the unknowns and the data
===========================================================================*/
		IN INTEGER iopt		"int that specifies the operation",
		IN REAL x[] 		   "vector of unknowns", 
		IN REAL rinput[] 		"Real vector for data input",
		OUT REAL routput[]   "Real vector for output of results" ,
		OUT REAL F[]   		"Residues",
		OUT NLEQN_SYSTEM_BASE nleqn_system  "Base of the non linear equation system"
		) 
/*
* iopt             = Integer specifying the operation:
*                     1, to set the number of equations, create unknowns and equations
*                     2, to calculate the residuals
*                     3, to calculate residuals and functions
* X[]              = Real array containing a vector of unknowns
* rinput[]         = Real array used for communication of inputs
*                    between the calling Solve method and the current function
* routput[]        = Real array used for communication of output
*                    between the calling Solve method and the current function
* F[]              = Real array with the residues of the equations
* nleqn_system     = Object of type NLEQN_SYSTEM_BASE, a parent class 
*                    for NLEQN_SYSTEM. NLEQN_SYSTEM_BASE includes all attributes
*                    and methods of NLEQN_SYSTEM except the SOLVE method.
*                    The SOLVE method calls this function using the "this"
*                    pointer, allowing the function to access the class to define
*                    the number of equations, create unknowns and equations,
*                    and fill the Jacobian by stamping
*/

CLASS NLEQN_SYSTEM IS_A NLEQN_SYSTEM_BASE
/* =========================================================================
* CLASS: NLEQN_SYSTEM
* Purpose: This class is designed for the resolution of 
*          sparse non-linear equation systems within the EcosimPro language, 
*          utilizing the Newton-Raphson Method.
*          It is particularly well-suited for large-scale simulations, 
*          specifically handling sparse systems of equations with a large size, 
*          often reaching thousands of equations.
* ==========================================================================
--Example of Usage of the  CLASS SPARSE MATRIX:
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
         END IF
      END IF
      RETURN
END FUNCTION
COMPONENT TestNleqnSystem01
	DATA
		INTEGER iter_damp = 2  "Number of initial damped iterations"
		INTEGER iter_max = 20  "Maximum number of iterations"	
		REAL FRAC_damp = 0.1   
		REAL tol = 1e-6		  "Tolerance"
		INTEGER NR_debug = 2    "Debugging level of the Newton-Raphson"
		INTEGER STAMP_debug = 0
	DECLS
		DISCR REAL FRAC = 1
		DISCR REAL errorx 
		DISCR REAL errorf
		DISCR	REAL x[20]
		DISCR	REAL rdata[10]
	OBJECTS
		NLEQN_SYSTEM(MAX_EQN=20, MAX_STAMP = 100) nleqn_system		
	INIT
		nleqn_system.NR_debug = NR_debug
		nleqn_system.Jacob.STAMP_debug = STAMP_debug
		nleqn_system.Solve(fsystem, iter_max, iter_damp, FRAC_damp, tol, rdata)  
		nleqn_system.GetUnknown(x)
END COMPONENT*/
	METHODS

----------------------------------------------------------------------------
-- Method to apply the Newton-Raphson (NR) method on a system of non-linear equations
--  Some details of the algorithm are based on two papers:
--  1) Numerical Solution of Constrained Non-Linear Algebraic Equations 
--     by Mordechai Shacham
--  2) A class of Methods for Solving Non Linear Simultaneous Equations by C.G.Broyden
-- 
----------------------------------------------------------------------------
   METHOD NO_TYPE Solve(
      IN FUNC_PTR <func_Fsystem> fsys  "Pointer to the function with the equations",
      IN INTEGER iter_max               "Maximum number of iterations", 
      IN INTEGER iter_damp              "Number of initial damped iterations",
      IN REAL FRAC_damp                 "Multiplier of the NR step during the first iter_damp iterations",
      IN REAL tol_arg                   "Tolerance of the NR",
      IN REAL rinput[]                  "Auxiliary vector for communication of inputs to fsys",
		OUT REAL routput[]                "Auxiliary vector for communication of outputs from fsys"
      )
      DECLS
			REAL FRAC                      "Fractional value for updating the Newton step"
         REAL df_max[MAX_EQN]           "Maximum absolute values of Jacobian rows"
         REAL dx_nr[MAX_EQN]            "Newton step for unknowns"
         REAL dx_sp[MAX_EQN]            "Gradient step for unknowns"
         REAL J_dx_sp[MAX_EQN]          "Jacobian times the gradient step"
         REAL nFun[MAX_EQN]             "Negative values of the system of equations"
         REAL fracLimit                 "Limiting value for fractional update"
         REAL eta                       "Correction factor for reducing the Newton step"
         INTEGER i, j, k, k1, k2        "Loop indices and counters"
         REAL errorx[MAX_EQN]           "Error in unknowns at each iteration"
         REAL errorx_max                "Maximum error in unknowns"
         REAL errorx_tot                "Total error in unknowns"
         REAL errorf[MAX_EQN]           "Error in equations at each iteration"
         REAL errorf_max                "Maximum error in equations"
         REAL errorf_tot                "Total error in equations"
         REAL errorf_tot_old            "Total error in equations at the previous iteration"
         INTEGER iworstx, iworstf       "Indices of the worst errors in unknowns and equations"
		BODY
			iter = 1
			tol = tol_arg
			converged = FALSE
			fsys(0, Xunk, rinput, routput, Fun, this)
			WHILE((iter <= iter_max) AND ((NOT converged)))	
			   IF(iter == 1) THEN
				   IF(NR_debug > 1) THEN
					   WRITE("*********************************************************************************************************************\n")  
					   WRITE("*****************************   List of Unknowns  *******************************************************************\n")  
					   FOR(j IN 1, neqn)
						   WRITE("Unknown  %4d \t %s\n", j, descrUnk[j])
					   END FOR
					   WRITE("*****************************   List of Equations  *******************************************************************\n")  
					   FOR(j IN 1, neqn)
						   WRITE("Equation %4d \t %s\n", j, descrEqn[j])
					   END FOR		
						WRITE("*********************************************************************************************************************\n")  
				   END IF
					Jacob.setToZero()
				   fsys(2, Xunk, rinput, routput, Fun, this) 
					Jacob.buildCSR()
				   Jacob.MaxAbsRow(df_max)
				   EvalFunError(Fun, df_max, errorf, errorf_tot, iworstf)
			   END IF
			   FOR(j IN 1, neqn)
				   XunkOld[j] = Xunk[j]
				   FunOld[j] = Fun[j]
				   nFun[j] =  -Fun[j]
				   errorf_tot_old = errorf_tot
			   END FOR
			   --Calculation of the Newton step
			   Jacob.LinearEqnSolve(nFun, dx_nr, 1., superLU_debug)
			   --Calculation of the gradient step
			   Jacob.preMultiplyByRow(nFun, dx_sp)
		 	   Jacob.postMultiplyByColumn(dx_sp, J_dx_sp)
				FRAC=1
				IF(iter <= iter_damp) THEN
					FRAC = FRAC_damp
			   END IF				
				--Update of the Unknowns taking into account limitations. 
				UpdateUnknownVec(FRAC, XunkOld, dx_nr, Xunk, errorx, errorx_tot, iworstx)
				Jacob.setToZero()
			   fsys(2, Xunk, rinput, routput, Fun, this)
				Jacob.buildCSR()
				EvalFunError(Fun, df_max, errorf, errorf_tot, iworstf)
				--Limitation of the Newton step according to the Original Broyden paper:
				--A class of Methods for Solving Non Linear Simultaneous Equations by C.G.Broyden
				IF(errorf_tot > errorf_tot_old) THEN
					eta = (errorf_tot / errorf_tot_old)**2
					FRAC = FRAC*((1+6.*eta)**0.5 - 1)/(3.*eta)
					--
					UpdateUnknownVec(FRAC, XunkOld, dx_nr, Xunk, errorx, errorx_tot, iworstx)
					Jacob.setToZero()
				   fsys(2, Xunk, rinput, routput, Fun, this)
					Jacob.buildCSR()
					EvalFunError(Fun, df_max, errorf, errorf_tot, iworstf)
				END IF		

				IF(NR_debug > 1) THEN
				   WRITE("*********************************************************************************************************************\n") 
				   WRITE("  No     Xunk. Old       Xunk. New     Dx Increment     Funct. Old      Funct. New    Scal.Unk.Error  Scal.Eqn.Error\n")
				   WRITE("===== =============== =============== =============== =============== =============== =============== ===============\n")
				   FOR(j IN 1, neqn)
					   WRITE(" %4d  %14.6lg  %14.6lg  %14.6lg  %14.6lg  %14.6lg  %14.6lg  %14.6lg\n", j, XunkOld[j], Xunk[j], FRAC*dx_nr[j], FunOld[j], Fun[j], errorx[j], errorf[j])
				   END FOR
				   WRITE("*********************************************************************************************************************\n")  
				END IF
				IF(NR_debug > 0)THEN
					WRITE("*** iter = %3d  Fraction NR= %-6.4lf  Unk.Error= %-10.3lg  Fun.Error= %-10.3lg  \
				      Worst Unk.= %-5d  Worst Eqn= %-5d\n", iter, FRAC, errorx_tot, errorf_tot, iworstx, iworstf)  
				END IF
				ASSERT(neqn == Jacob.nrow AND neqn == Jacob.ncol) FATAL  "Inconsistent dimension of Equation System and Jacobian"
				converged = (iter > 0) AND (errorx_tot < tol) AND (errorf_tot < 10 * tol)
				iter = iter + 1
			END WHILE		
      END METHOD
END CLASS

CLASS LEQ_SYSTEM 	(
		INTEGER NEQN		"Maximum number of equations", 		 
		INTEGER NNZEROS	"Maximum number of non zero matrix coefficients"
	)
/* =========================================================================
* CLASS: LEQ_SYSTEM
* Purpose: This class is designed for the resolution of sparse linear equation 
*          systems within the EcosimPro language, .
* ==========================================================================*/
	DECLS
		INTEGER nnzeros = 0
		INTEGER irow[NNZEROS] = 0
		INTEGER icol[NNZEROS] = 0
		INTEGER ipos[NNZEROS] = 0
		REAL A[NNZEROS] = 0.
		REAL B[NEQN] = 0.
		REAL X[NEQN] = 0.
		INTEGER colind[NNZEROS] = 0
		INTEGER rowptr[NEQN+1] = 0
		REAL superlu_pivot_threshold = 1.
		INTEGER superlu_equil =  0
		INTEGER superlu_refine = 0
		INTEGER superlu_debug =  0
	METHODS	
	METHOD NO_TYPE SetToZero()
		BODY
			nnzeros = 0
			FOR (i IN 1, NEQN)
				B[i] = 0
			END FOR
			FOR (i IN 1, NNZEROS)
				A[i] = 0.
			END FOR
			FOR(i IN 1, NNZEROS)
				irow[i] = 0
				icol[i] = 0
				A[i] = 0.
			END FOR
			RETURN
	END METHOD	
	
	METHOD NO_TYPE Fill_B(IN INTEGER ieqn, IN REAL b)
		BODY 
		B[ieqn] = b
	END METHOD
	
	METHOD NO_TYPE Fill_A(IN INTEGER ieqn, IN INTEGER iunk, IN REAL x)
		BODY
			nnzeros = nnzeros + 1
			ASSERT(nnzeros <= NNZEROS) FATAL "Dimension NNZEROS of LEQ_SYTEM object lower than needed"
			irow[nnzeros] = ieqn
			colind[nnzeros] = iunk
			IF(nnzeros ==1) THEN
				rowptr[irow[nnzeros]] = nnzeros
			END IF
			IF(nnzeros > 1 AND irow[nnzeros] != irow[nnzeros-1]) THEN
				rowptr[irow[nnzeros]]= nnzeros
			END IF

			IF(NNZEROS == nnzeros) THEN
			rowptr[NEQN+1] = nnzeros + 1
			END IF 		
			A[nnzeros] = x		
			RETURN
	END METHOD
	
	
	METHOD NO_TYPE Solve()
		DECLS
			INTEGER ldb
			INTEGER factors[4]
			INTEGER info
			INTEGER colind0[NNZEROS], rowptr0[NEQN+1]
			REAL value[NNZEROS]
			REAL dx[NEQN]
		BODY
		   FOR(i IN 1, NNZEROS)
			   value[i] = A[i]
				colind0[i] = colind[i]
			END FOR
			FOR(i IN 1, NEQN+1)
				rowptr0[i] = rowptr[i]
			END FOR
			FOR(j IN 1, NEQN)
				dx[j] =  B[j]
			END FOR
			
			--Factorize Matrix into L & U
			d_ecosim_superlu(1, NEQN, nnzeros, 1, value, colind0, rowptr0, dx, factors, info, superlu_pivot_threshold, superlu_equil, superlu_refine, superlu_debug)
			IF (info != 0) THEN
         	 WRITE("\n****Error in SUPERLU Factorization info = %d\n\n", info)
			END IF

			--Solve Linear System
			d_ecosim_superlu(2, NEQN, nnzeros, 1, value, colind0, rowptr0, dx, factors, info, superlu_pivot_threshold, superlu_equil, superlu_refine, superlu_debug)	
			IF (info != 0) THEN
         	 WRITE("\n****Error in SUPERLU Solve info = %d\n\n", info)
			END IF			
			--Destroy Matrix
			d_ecosim_superlu(3, NEQN, nnzeros, 1, value, colind0, rowptr0, dx, factors, info, superlu_pivot_threshold, superlu_equil, superlu_refine, superlu_debug)		
			
			FOR(j IN 1, NEQN)
				X[j] =  dx[j] 
			END FOR
	END METHOD
	
	METHOD REAL GetUnknown(IN INTEGER iunk)
		BODY
			RETURN X[iunk]
	END METHOD
END CLASS

