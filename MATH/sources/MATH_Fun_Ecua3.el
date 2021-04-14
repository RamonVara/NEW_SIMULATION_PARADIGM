--Function ecua3:
-- Purpose: to solve cubic equations
-- Arguments: a[1]   x^3 coefficient
--            a[2]   x^2 coefficient
--            a[3]   x^1 coefficient
--            a[4]   x^0 coefficient
--            *r1,*r2,*r3  pointers to the 3 solutions
--
CONST REAL tol=1.e-12		"Numeric tolerance"

FUNCTION NO_TYPE ecua2
	( 
		IN REAL  a1		"Coefficient of the quadratic term", 
		IN REAL a2		"Coefficient of the linear term", 
		IN REAL a3		"Coefficient of the independent term",
		OUT REAL r1		"First solution of the quadratic equation", 
		OUT REAL r2		"Second solution of the quadratic equation"
	)
"Function to solve a quadratic equation"

   DECLS
	REAL det	"Auxiliary variable"
   BODY
	
      r1 = 1e40
	r2 = 1e40
	IF (abs(a1) > tol) THEN
	
		det = a2*a2 - 4*a1*a3

		IF  (det > 0) THEN
			--Calculation of the root with greatest absolute value
			IF (a2 >=0) THEN
				r1 = (-a2 - sqrt(det))/(2*a1)
			ELSE
				r1 = (-a2 + sqrt(det))/(2*a1)
			END IF	
			--Calculation of the root with smallest absolute value
			r2 = (a3/a1) / r1

            END IF
	ELSE
		--Linear equation
           	r1 = -a3/a2
      END IF
 
	RETURN
END FUNCTION

FUNCTION NO_TYPE ecua3
	(
		IN REAL a[4] 	"Array of the equation coefficients", 
		OUT REAL r1		"First solution of the cubic equation", 
		OUT REAL r2		"Second solution of the cubic equation", 
		OUT REAL r3		"Third solution of the cubic equation"
	)
"Function to solve a cubic equation"
   DECLS
	REAL Q		"Intermediate variable"
	REAL Q3		"Cubic value of Q"
	REAL R		"Intermediate variable"
	REAL theta		"Intermediate variable"
	REAL aux		"Auxiliary variable"

      CONST REAL pi =3.14159265358979	"PI number"
      
      INTEGER nroots 	"Number of roots"
   BODY
	r1 = 1e40
	r2 = 1e40
	r3 = 1e40
	nroots = 0
	-- Cubic Equation
	IF (abs(a[1]) > tol) THEN 

		a[2] = a[2] / a[1]
		a[3] = a[3] / a[1]
		a[4] = a[4] / a[1]

		Q = (a[2]* a[2] - 3 * a[3]) / 9.

		R = (2 * a[2]*a[2]*a[2]- 9.*a[2]*a[3] + 27*a[4]) / 54.

		Q3 = Q*Q*Q

		IF (Q3 > R*R) THEN 
		
			--three real roots
      	 	theta = acos(R/ssqrt(Q3))
			r1 = - 2 * ssqrt(Q) * cos(theta/3.) - a[2]/3
			r2 = - 2 * ssqrt(Q) * cos((theta+2*pi)/3.) - a[2]/3
			r3 = - 2 * ssqrt(Q) * cos((theta+4*pi)/3.) - a[2]/3
		
		ELSE
		
			aux = (ssqrt(R*R-Q3)+abs(R))**(1./3.)
			IF (R >=0) THEN
			    r1 = -(aux + Q/aux)  -a[2]/3.
			ELSE
			    r1 = +(aux + Q/aux)  -a[2]/3.
                  END IF
		END IF
	ELSEIF (abs(a[2]) > tol) THEN
          ecua2(a[2], a[3], a[4], r1, r2)
	ELSEIF (abs(a[3]) >  tol) THEN
           r1 = -a[4] / a[3]
      END IF
	RETURN
END FUNCTION
