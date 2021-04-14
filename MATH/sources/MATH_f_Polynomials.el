--------------------------------------------------------------------------------
-- EA Internacional 2004   EcosimPro 3.3 Simulation Source Code
--
--    FILE NAME:  MATH_f_Polynomials.el
--    FILE TYPE:  Functions of the MATH library
--    DESCRIPTION:  Defines math function for polynomials for the MATH library
--    NOTES:  Based on EcosimPro MATH library
--    AUTHOR:  Ramon Perez Vara
--    CREATION DATE:  12-Jan-2004
--------------------------------------------------------------------------------


-- Functions
--------------------------------------------------------------------------------
-- Function poly
--------------------------------------------------------------------------------
-- Purpose:
--    To calculate: f(x) = A[1] + A[2] * x + .... + A[n] * x**(n-1).
--------------------------------------------------------------------------------
FUNCTION REAL poly
   (
   IN REAL x                      "Real value to be substituted",
   IN INTEGER n                   "Order of the polynomial plus 1",
   IN REAL A[n]                   "Real array of polynomial coefficients"
   )

   DECLS
      INTEGER j                   "loop counter"

      REAL y                      "Output value"

   BODY
      FOR (j=n, y=0.; j>0; j=j-1)
         y = y * x + A[j]
      END FOR

      RETURN y

END FUNCTION

--------------------------------------------------------------------------------
-- Function poly_int
--------------------------------------------------------------------------------
-- Purpose:
--    To calculate the integral of a polynomial p(x).
--       p(x) = A[1] + A[2] * x + .... + A[n] * x**(n-1)
--       f(x) = A[1] * x + A[2] * x**2 / 2 + .... + A[n] * x**(n) / n
--------------------------------------------------------------------------------
FUNCTION REAL poly_int
   (
   IN REAL x                      "Real value to be substituted",
   IN INTEGER n                   "Order of the polynomial plus 1",
   IN REAL A[n]                   "Real array of polynomial coefficients"
   )

   DECLS
      INTEGER j                   "loop counter"

      REAL y                      "Output value"

   BODY
      FOR (j=n, y=0.; j>0; j=j-1)
         y = y * x + A[j] / j
      END FOR

      y = y * x

      RETURN y

END FUNCTION

--------------------------------------------------------------------------------
-- Function poly_dx_int
--------------------------------------------------------------------------------
-- Purpose:
--    To calculate the integral of a polynomial p(x) divided by x.
--       p(x)     = A[1] + A[2] * x + .... + A[n] * x**(n-1)
--       p(x) / x = A[1] / x + A[2] + .... + A[n] * x**(n-2)
--       f(x)     = A[1] * log(x) + A[2] * x + A[3] * x**2 / 2 + ... + \
--                  A[n] * x**(n-1) / (n-1)
--------------------------------------------------------------------------------
FUNCTION REAL poly_dx_int
   (
   IN REAL x                      "Real value to be substituted",
   IN INTEGER n                   "Order of the polynomial plus 1",
   IN REAL A[n]                   "Real array of polynomial coefficients"
   )

   DECLS
      INTEGER j                   "loop counter"

      REAL y                      "Output value"

   BODY
      FOR (j = n,  y = 0.; j > 1; j = j -   1)
         y = y * x + A[j] / (j - 1.)
      END FOR

      y = y * x + A[1] * log(x)

      RETURN y

END FUNCTION

--------------------------------------------------------------------------------
-- Function fac
--------------------------------------------------------------------------------
-- Purpose:
--    To evaluate the factorial of an integer number.
--------------------------------------------------------------------------------
FUNCTION REAL fac
   (
   IN INTEGER n
   )

   DECLS
      INTEGER k

      REAL y

   BODY
      k = n
      y = 1.

      WHILE (k > 1)
         y = y * k
         k = k - 1
      END WHILE

      RETURN y

END FUNCTION

--------------------------------------------------------------------------------
-- Function PolyMult
--------------------------------------------------------------------------------
-- Purpose:
--    To evaluate the multiplication of two polynomial.
--------------------------------------------------------------------------------
FUNCTION NO_TYPE PolyMult
   (
   INTEGER k,
   REAL a[k],
   INTEGER l,
   REAL b[l],
   OUT INTEGER n,
   OUT REAL c[n]
   )

   BODY
      n = k + l

      FOR (i IN 1,n)
         c[i] = 0.
      END FOR

      FOR (i IN 1,k)
         FOR (j IN 1,l)
            c[i+j-1] = c[i+j-1] + a[i] * b[j]
         END FOR
      END FOR

END FUNCTION
