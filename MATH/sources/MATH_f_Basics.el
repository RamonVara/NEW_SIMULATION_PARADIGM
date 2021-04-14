--------------------------------------------------------------------------------
-- EA Internacional 2004   EcosimPro 3.3 Simulation Source Code
--
--    FILE NAME:  MATH_f_Basics.el
--    FILE TYPE:  Functions of the MATH library
--    DESCRIPTION:  Defines basic math functions for the MATH library
--    NOTES:  Based on EcosimPro MATH library
--    AUTHOR:  Ramon Perez Vara
--    CREATION DATE:  12-Jan-2004
--------------------------------------------------------------------------------


-- Functions
--------------------------------------------------------------------------------
-- Function max
--------------------------------------------------------------------------------
-- Purpose:
--    To calculate the maximum value of two.
--------------------------------------------------------------------------------
FUNCTION REAL max
   (
   REAL x,
   REAL y
   )

   BODY
      IF (x > y) THEN
         RETURN x
      END IF

      RETURN y

END FUNCTION

--------------------------------------------------------------------------------
-- Function min
--------------------------------------------------------------------------------
-- Purpose:
--    To calculate the minimum value of two.
--------------------------------------------------------------------------------
FUNCTION REAL min
   (
   REAL x,
   REAL y
   )

   BODY
      IF (x < y) THEN
         RETURN x
      END IF

      RETURN y

END FUNCTION

--------------------------------------------------------------------------------
-- Function max_impl
--------------------------------------------------------------------------------
-- Purpose:
--    To calculate the maximum value of two; output also in second argument
--------------------------------------------------------------------------------
FUNCTION REAL max_impl(IN REAL x, OUT REAL z)
BODY
        IF (x > z ) THEN
            z = x
            RETURN x
        END IF
        RETURN z
END FUNCTION


--------------------------------------------------------------------------------
-- Function min_impl
--------------------------------------------------------------------------------
-- Purpose:
--    To calculate the minimum value of two; output also in second argument
--------------------------------------------------------------------------------
FUNCTION REAL min_impl(IN REAL x, OUT REAL z)
BODY
        IF (x < z ) THEN
           	z = x
	     		RETURN x
        END IF
        RETURN z
END FUNCTION

--------------------------------------------------------------------------------
-- Function sign
--------------------------------------------------------------------------------
-- Purpose:
--    To evaluate the sign of a value. The output is as follows:
--        1.    when    x >= 0.
--       -1.    when    x <  0.
--------------------------------------------------------------------------------
FUNCTION REAL sign
   (
   REAL x
   )

   BODY 
      IF (x >= 0.) THEN
         RETURN 1.
      END IF

      RETURN -1.

END FUNCTION

--------------------------------------------------------------------------------
-- Function sign2
--------------------------------------------------------------------------------
-- Purpose:
--    To evaluate the sign of a value with an epsilon. The output is as follows:
--        1.    when    x >= - eps
--       -1.    when    x <  -eps
--------------------------------------------------------------------------------
FUNCTION REAL sign2
   (
   REAL x,
   REAL eps
   )

   BODY 
      IF (x >= eps) THEN
         RETURN 1.
      ELSEIF (x<=-eps) THEN
         RETURN -1.
      ELSE
         RETURN x/eps
      END IF

      RETURN -1.

END FUNCTION

--------------------------------------------------------------------------------
-- Function sign3
--------------------------------------------------------------------------------
-- Purpose:
--    To evaluate the sign of a value with an epsilon. The output is as follows:
--        1.    when    x >= - eps
--       -1.    when    x <  -eps
--------------------------------------------------------------------------------
FUNCTION REAL sign3
   (
   REAL x,
   REAL y,
   REAL eps
   )

   BODY 
      IF (x/y >= eps) THEN
         RETURN 1.
      ELSEIF (x<=-eps) THEN
         RETURN -1.
      ELSE
         RETURN x/eps
      END IF

      RETURN -1.

END FUNCTION
--------------------------------------------------------------------------------
-- Function bound
--------------------------------------------------------------------------------
-- Purpose:
--    To bound or limit values. It should not be used to limit a dynamic
--    variable. The output value is as follows:
--       y = xmin    when    x < xmin
--       y = x       when    xmin <= x <= xmax
--       y = xmax    when    x > xmax
--------------------------------------------------------------------------------
FUNCTION REAL bound
   (
   IN REAL x                      "Input value",
   IN REAL xmin                   "Bottom bound",
   IN REAL xmax                   "Top bound"
   )

   DECLS
      REAL y                      "Output value"

   BODY
      IF (x < xmin) THEN
         y = xmin
      ELSEIF (x > xmax) THEN
         y = xmax
      ELSE
         y = x
      END IF

      RETURN y

END FUNCTION

--------------------------------------------------------------------------------
-- Function sum1D
--------------------------------------------------------------------------------
-- Purpose:
--    Function to calculate the sum of all elements in a 1D array
--------------------------------------------------------------------------------
FUNCTION REAL sum1D
   (
   IN INTEGER Nx                  "Number of array elements (-)",
   IN REAL A[Nx]                  "Array of elements (-)"
   )

   DECLS
      REAL sum = 0.

   BODY
      FOR (i IN 1,Nx)
         sum = sum + A[i]
      END FOR

   RETURN sum

END FUNCTION

--------------------------------------------------------------------------------
-- Function sum2D
--------------------------------------------------------------------------------
-- Purpose:
--    Function to calculate the sum of all elements in a 2D array
--------------------------------------------------------------------------------
FUNCTION REAL sum2D
   (
   IN INTEGER Nx                  "Number of files of matrix (-)",
   IN INTEGER Ny                  "Number of columns of matrix (-)",
   IN REAL A[Nx,Ny]               "Array matrix (-)"
   )

   DECLS
      REAL sum = 0.

   BODY
      FOR (i IN 1,Nx)
         FOR (j IN 1,Ny)
            sum = sum + A[i,j]
         END FOR
      END FOR

   RETURN sum

END FUNCTION

--------------------------------------------------------------------------------
-- Function discr
--------------------------------------------------------------------------------
-- Purpose:
--    Empty function for solve boundary problems with non assigned
--    variables in CONTINUOUS block.
--------------------------------------------------------------------------------
FUNCTION NO_TYPE discr
   (
   OUT REAL y
   )

   BODY
      RETURN

END FUNCTION

--------------------------------------------------------------------------------
-- Function discr_array
--------------------------------------------------------------------------------
-- Purpose:
--    Empty function for solve boundary problems with non assigned
--    array variables in CONTINUOUS block.
--------------------------------------------------------------------------------
FUNCTION NO_TYPE discr_array
   (
   OUT REAL y[]
   )

   BODY
      RETURN

END FUNCTION

--------------------------------------------------------------------------------
-- Function GetRealArrayItem
--------------------------------------------------------------------------------
-- Purpose:
--    Real function to extract an element of a real array.
--------------------------------------------------------------------------------
FUNCTION REAL GetRealArrayItem
   (
   IN REAL input[]     "Real array",
   IN INTEGER i        "Index")

   BODY
        RETURN input[i]
END FUNCTION
