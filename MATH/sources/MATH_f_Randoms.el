--------------------------------------------------------------------------------
-- EA Internacional 2004   EcosimPro 3.3 Simulation Source Code
--
--    FILE NAME:  MATH_f_Randoms.el
--    FILE TYPE:  External objects and functions of the MATH library
--    DESCRIPTION:  Defines math functions for the MATH library
--    NOTES:  Based on EcosimPro MATH library
--    AUTHOR:  Ramon Perez Vara
--    CREATION DATE:  12-Jan-2004
--------------------------------------------------------------------------------


-- External Objects
--------------------------------------------------------------------------------
-- Purpose:
--    To declare in EcosimPro Language some C functions from the book
--    "Numerical Recipes in C", that calculate random numbers.
--------------------------------------------------------------------------------
"C" FUNCTION REAL c_ran1
   (
   OUT INTEGER iseed
   ) IN "MATH.lib"

"C" FUNCTION REAL c_gasdev
   (
   OUT INTEGER iseed
   ) IN "MATH.lib"

"C" FUNCTION REAL c_poidev
   (
   IN REAL xm,
   OUT INTEGER iseed
   ) IN "MATH.lib"

"C" FUNCTION REAL c_bnldev
   (
   IN REAL pp,
   IN INTEGER n,
   OUT INTEGER iseed
   ) IN "MATH.lib"

"C" FUNCTION REAL c_gamdev
   (
   IN INTEGER ia,
   OUT INTEGER iseed
   ) IN "MATH.lib"

"C" FUNCTION REAL c_expdev
   (
   OUT INTEGER iseed
   ) IN "MATH.lib"


-- Functions
--------------------------------------------------------------------------------
-- Function reminder
--------------------------------------------------------------------------------
-- Purpose:
--    Defines a function which returns the remainder of the division of two
--    integers.
--------------------------------------------------------------------------------
FUNCTION INTEGER remainder
   (
   IN INTEGER dividend,
   IN INTEGER divisor
   )

   DECLS
      INTEGER i_quotient

      REAL quotient

   BODY
      quotient = dividend / divisor
      i_quotient = quotient

   RETURN dividend - i_quotient * divisor

END FUNCTION

--------------------------------------------------------------------------------
-- Function int
--------------------------------------------------------------------------------
-- Purpose:
--    Function calculating the integer part of a real.
--------------------------------------------------------------------------------
FUNCTION INTEGER int
   (
   IN REAL real_number
   )

   BODY

   RETURN real_number

END FUNCTION

--------------------------------------------------------------------------------
-- Function gammln
--------------------------------------------------------------------------------
-- Purpose:
--    Function calculating the natural logarithm of the gamma function of its
--    argument, this is: f(x) = Ln(gamma(x)). This is an exact function just
--    for values 1 < x.
--------------------------------------------------------------------------------
FUNCTION REAL gammln
   (
   IN REAL x
   )

   DECLS
      REAL coef[6]
      REAL y
      REAL ser
      REAL tmp

   BODY
      coef[1] = 76.18009173
      coef[2] = - 86.50532033
      coef[3] = 24.01409822
      coef[4] = - 1.231739516
      coef[5] = 0.120858003e-2
      coef[6] = - 0.536382e-5

      y = x - 1

      tmp = y + 5.5
      tmp = tmp - (x + 0.5) * log(tmp)

      ser = 1.0	

      FOR (j IN 1,5)
         x = x + 1 
         ser = ser + coef[j] / x
      END FOR

   RETURN -tmp + log(2.50662827465 * ser)

END FUNCTION

--------------------------------------------------------------------------------
-- Function ran1
--------------------------------------------------------------------------------
-- Purpose:
--    Defines a function calculating uniform random numbers on interval (0,1).
--------------------------------------------------------------------------------
FUNCTION REAL ran1
   (
   OUT INTEGER integer_array[],
   OUT REAL real_array[]
   )

   DECLS
         -- Datas
      INTEGER M1 = 259200         "Module 1"
      INTEGER IA1 = 7141
      INTEGER IC1 = 54773
      INTEGER M2 = 134456         "Module 2"
      INTEGER IA2 = 8121
      INTEGER IC2 = 28411
      INTEGER M3 = 243000         "Module 3"
      INTEGER IA3 = 4561
      INTEGER IC3 = 51349

      REAL RM1 = 1 / 259200       "Inverse of module 1" 
      REAL RM2 = 1 / 134456       "Inverse of module 2"

         -- Vars
      INTEGER k                   "Auxiliar integer"

      REAL temp                   "Random number"

   BODY
         -- Initializes sequence
      IF ((integer_array[4] < 0) OR (integer_array[5] == 0)) THEN
         integer_array[5] = 1
         integer_array[1] = remainder(IC1 - integer_array[4], M1)
         integer_array[1] = remainder(IA1 * integer_array[1] + IC1, M1)
         integer_array[2] = remainder(integer_array[1], M2)
         integer_array[1] = remainder(IA1 * integer_array[1] + IC1, M1)
         integer_array[3] = remainder(integer_array[1], M3)

         FOR (j IN 1,97)
            integer_array[1] = remainder(IA1 * integer_array[1] + IC1, M1)
            integer_array[2] = remainder(IA2 * integer_array[2] + IC2, M2)
            real_array[j] = (integer_array[1] + integer_array[2] * RM2) * RM1
         END FOR

         integer_array[4] = 1

      END IF

         -- Operations to be done when sequence has been initialized
      integer_array[1] = remainder(IA1 * integer_array[1] + IC1, M1)
      integer_array[2] = remainder(IA2 * integer_array[2] + IC2, M2)
      integer_array[3] = remainder(IA3 * integer_array[3] + IC3, M3)

      k = 1 + (97 * integer_array[3]) / M3

      temp = real_array[k]

      real_array[k] = (integer_array[1] + integer_array[2] * RM2) * RM1

   RETURN temp

END FUNCTION

--------------------------------------------------------------------------------
-- Function expdev
--------------------------------------------------------------------------------
-- Purpose:
--    Function calculating random numbers following an exponential
--    distribution, this is, a distribution of the form: p(x) = exp(-x).
--------------------------------------------------------------------------------
FUNCTION REAL expdev
   (
   OUT INTEGER integer_array[],
   OUT REAL real_array[]
   )

   BODY

   RETURN - log(ran1(integer_array, real_array))

END FUNCTION

--------------------------------------------------------------------------------
-- Function gasdev
--------------------------------------------------------------------------------
-- Purpose:
--    Function calculating random numbers which follow a gauss distribution of
--    probability, this is, a probability density function of the form:
--       p(x) = exp(-x**2) / sqrt(2*PI)
--------------------------------------------------------------------------------
FUNCTION REAL gasdev
   (
   OUT INTEGER integer_array[],
   OUT REAL real_array[]
   )

   DECLS
      REAL v[2]
      REAL ratio = 1.5
      REAL fac
      REAL return_argument

   BODY
      IF (integer_array[11] == 0) THEN
         WHILE (ratio >= 1)
            v[1] = 2 * ran1(integer_array, real_array) - 1
            v[2] = 2 * ran1(integer_array, real_array) - 1
            ratio = v[1]**2 + v[2]**2
         END WHILE

         fac = sqrt(- 2 * log(ratio) / ratio)

         integer_array[11] = 1

         real_array[121] = v[1] * fac

         return_argument =  v[2] * fac

      ELSE
         integer_array[11] = 0

         return_argument = real_array[121]

      END IF

   RETURN return_argument

END FUNCTION

--------------------------------------------------------------------------------
-- Function gamdev
--------------------------------------------------------------------------------
-- Purpose:
--    Function calculating random numbers which follow a gamma distribution.
--------------------------------------------------------------------------------
FUNCTION REAL gamdev
   (
   OUT 	INTEGER 	integer_array[],
   OUT 	REAL		real_array[]
   )

   DECLS
      REAL s
      REAL x
      REAL am
      REAL e
      REAL y
      REAL v[2]

   BODY
      IF (integer_array[6] < 6) THEN
         x = 1

         FOR (j IN 1, integer_array[6])
            x = x * ran1(integer_array, real_array)
         END FOR

         x = - log(x)

      ELSE
         WHILE(ran1(integer_array, real_array) > e)
            WHILE(x <= 0)
               WHILE(v[1]**2 + v[2]**2 > 1)
                  v[1] = 2 * ran1(integer_array, real_array)
                  v[2] = 2 * ran1(integer_array, real_array)
               END WHILE

               y = v[2] / v[1]

               am = integer_array[6] - 1

               s = sqrt(2 * am + 1)

               x = s + y + am

            END WHILE

            e = 1 + y**2 * exp(x / am - s * y)

         END WHILE

      END IF

   RETURN x

END FUNCTION

--------------------------------------------------------------------------------
-- Function poidev
--------------------------------------------------------------------------------
-- Purpose:
--    This function generate random numbers deviated following a Poisson
--    distribution.
--------------------------------------------------------------------------------
FUNCTION INTEGER poidev
   (
   OUT INTEGER integer_array[],
   OUT REAL real_array[]
   )

   DECLS
      REAL em
      REAL t
      REAL y

   BODY
      IF (real_array[101] < 12.0) THEN	
         IF (real_array[101] != real_array[121]) THEN
            real_array[121] = real_array[101]
            real_array[124] = exp(-real_array[101])
         END IF

         em = - 1
         t = 1.0

         WHILE (t > real_array[124])
            em = em + 1
            t = t * ran1(integer_array, real_array)
         END WHILE

      ELSE
         IF (real_array[101] != real_array[121]) THEN
            real_array[121] = real_array[101]
            real_array[122] = sqrt(2 * real_array[101])
            real_array[123] = log(real_array[101])
            real_array[124] = real_array[101] * real_array[101] - gammln(real_array[101] + 1)
         ELSE
         END IF

         WHILE (ran1(integer_array, real_array) > t)
            WHILE (em < 0.0)
               y = tan(PI * ran1(integer_array, real_array))
               em = real_array[122] * y * real_array[101]
            END WHILE

            em = int(em)
            t = 0.9 * (1 + y**2) * exp(em * real_array[101] - gammln(em + 1.0) - real_array[124])

         END WHILE

      END IF

   RETURN em

END FUNCTION

--------------------------------------------------------------------------------
-- Function bnldev
--------------------------------------------------------------------------------
-- Purpose:
--    Defines a function calculating random numbers following a binomial
--    distribution.
--------------------------------------------------------------------------------
FUNCTION REAL bnldev
   (
   OUT INTEGER integer_array[],
   OUT REAL real_array[]
   )

   DECLS
      INTEGER j

      REAL am
      REAL em
      REAL g
      REAL angle
      REAL p
      REAL bnl
      REAL sq
      REAL t
      REAL y

   BODY
      IF (real_array[101] <= 0.5) THEN
         p = real_array[101]
      ELSE
         p = 1 - real_array[101]
      END IF

      am = integer_array[6] * p

      IF (integer_array[6] < 25) THEN
         bnl = 0

         FOR (k IN 1,integer_array[6])
            IF (ran1(integer_array, real_array) < p) THEN
               bnl = bnl + 1
            END IF
         END FOR

      ELSEIF (am < 1) THEN
         g = exp(- am)

         t = 1

         FOR (j=0; j<=integer_array[6]; j=j+1)
            t = t * ran1(integer_array, real_array)

            IF (t < g) THEN
               j = integer_array[6]
            END IF

         END FOR

         IF (j <= integer_array[6]) THEN
            bnl = j
         ELSE
            bnl = integer_array[6]
         END IF

      ELSE
         IF (integer_array[6] != integer_array[11]) THEN
            real_array[127] = integer_array[6]
            real_array[126] = gammln(real_array[127] + 1)
            integer_array[11] = integer_array[6]
         END IF

         IF(p != real_array[121]) THEN
            real_array[122] = 1 - p
            real_array[123] = log(p)
            real_array[124] = log(real_array[122])
            real_array[121] = p
         END IF

         sq = sqrt(2 * am * real_array[122])

         WHILE (ran1(integer_array, real_array) > t)
            WHILE ((em < 0) OR (em >= real_array[127] + 1))
               angle = PI * ran1(integer_array, real_array)
               y = tan(angle)
               em = sq * y + am
            END WHILE

            em = int(em)

            t = 1.2 * sq * (1 + y**2) * exp(real_array[126] - \
               gammln(em + 1) - gammln(real_array[127] - em + 1) + \
               em * real_array[123] + (real_array[127] - em) * real_array[124])

         END WHILE

         bnl = em

      END IF

      IF (p != real_array[101]) THEN
         bnl = integer_array[6] - bnl
      END IF

   RETURN bnl

END FUNCTION
