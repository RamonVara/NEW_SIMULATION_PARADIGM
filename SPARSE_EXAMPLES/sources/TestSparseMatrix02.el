/*-----------------------------------------------------------------------------------------
 LIBRARY: SPARSE
 FILE: TestSparseMatrix02
 CREATION DATE: 11/11/2023
-----------------------------------------------------------------------------------------*/
--Example of Usage of the  CLASS SPARSE MATRIX 
USE SPARSE
COMPONENT TestSparseMatrix02
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
		FOR(k IN 0, 10)
			A.setToZero()
			-- Let the matrix be 5x5 with the values shown below:
			-- A = [
      	--       8+k,    0,    0,    0,     1;
      	--         0,  1+2,    0,  4+1,     0;
      	--         1,    0,    3,    0,     0;
      	--         0,    1,    0,    3,     0;
      	--         3,    0,    1,    0, 2-6+k;
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
			A.stampItem(1, 1, k)
			A.stampItem(5, 5, k)
			--The following two stamps change the stamp structure for k=6 y k=7
			IF(k > 5 AND k <=7) THEN
				A.stampItem(k-5, 2, 1.)
			END IF
		   A.Print()
	      A.buildCSR()
		   --------------------------------------------------------
		   --Solving the system of linear equations A·x = b
		   --where the vector of independent terms
		   --------------------------------------------------------	
	      A.LinearEqnSolve(b, x, 1., 1)		
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
		   A.postMultiplyByColumn(x, d)
		   WRITE("\n\nResult of Post-multiplicacion:\n")
		   WRITE("    d[5]   = { %lg , %lg , %lg , %lg , %lg }\n", 
		                      d[1], d[2], d[3], d[4], d[5])
		   error = ((b[1]-d[1])**2 + (b[2]-d[2])**2 + (b[3]-d[3])**2 + (b[4]-d[4])**2 + + (b[4]-d[4])**2)**0.5
		   error = error /(b[1]**2 + b[2]**2 + b[3]**2 + b[4]**2 + b[5]**2)**0.5
		   ASSERT (error < 1e-7) ERROR "Solution of the Linear Equation System has two much error"
		END FOR
END COMPONENT
