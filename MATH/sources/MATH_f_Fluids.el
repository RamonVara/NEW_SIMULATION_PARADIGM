--------------------------------------------------------------------------------
-- EA Internacional 2004   EcosimPro 3.3 Simulation Source Code
--
--    FILE NAME:  MATH_f_Fluids.el
--    FILE TYPE:  Functions of the MATH library
--    DESCRIPTION:  Defines math functions for Fluid Simulation
--    NOTES:  Based on EcosimPro MATH library
--    AUTHOR:  Ramon Perez Vara
--    CREATION DATE:  12-Jan-2004
--------------------------------------------------------------------------------


-----------------------------------------------------------------------
--Function to perform safe division, avoiding division by zero
-- In case of division by zero, it returns the previous value
-----------------------------------------------------------------------
FUNCTION REAL div_safe
	(
		REAL num	"Numerator", 
		REAL den	"Denominator", 
		OUT REAL y	"Result of the division"
	)
"Function to perform safe division, avoiding division by zero. \
In case of division by zero, it returns the previous value."

DECLS
    CONST REAL tol = 1.e-8	"Minimun value allowed for the denominator"
BODY
    IF (abs(den) > tol) THEN
        y = num / den
    END IF

    RETURN y
END FUNCTION
-----------------------------------------------------------------------
--Function to calculate a dead band
--------------------------------------------------------------------
FUNCTION REAL deadband
	(
		REAL x	"Input variable", 
		REAL x_up	"Upper limit of the dead band", 
		REAL x_low	"Lower limit of the dead band"
	)
"Function to calculate a dead band"

    DECLS 
        REAL y	"Result"
    BODY
        IF ( x > x_up OR x < x_low) THEN 
           y = x
        ELSE 
           y =  0
        END IF
        RETURN y
END FUNCTION
-----------------------------------------------------------------------
--Function to calculate the fluid properties of the donor-cell
-----------------------------------------------------------------------
FUNCTION REAL donor_cell 
	(
		REAL flow			"Flow variable", 
		REAL upstream_prop	"Upstream value of the property", 
		REAL downstream_prop	"Downstream value of the property"
	)

"Function to calculate the fluid properties of the donor-cell"

DECLS
    REAL c				"Intermediate variable that indicates the flow direction"
    REAL prop			"Value of the fluid property depending on the flow direction"
    CONST REAL tol = 1.e-10	"Numerical tolerance"
BODY
    c = flow/(abs(flow)+tol)
    prop = 0.5 * ((1.+ c) * upstream_prop + (1.-c) * downstream_prop)
    RETURN prop
END FUNCTION 

----------------------------------------------------------------------------------------
--Function to calculate the coulomb friction 
-----------------------------------------------------------------------
FUNCTION REAL rev_fri 
(   
    REAL v        "Velocity (m/s)",
    REAL fc1      "Direct coulomb friction (N)",
    REAL fc2      "Inverse coulomb friction (N)"
)
DECLS
    REAL ffric    "Friction factor (-)"
BODY
        IF(v <= 0) THEN
          ffric = -fc2*tanh(v*1.e4)
        ELSE
          ffric = -fc1*tanh(v*1.e4)
        END IF

    RETURN ffric
END FUNCTION 

-----------------------------------------------------------------------
--Function to calculate signed power with linear zone near 0
--for exponents less than 1
-----------------------------------------------------------------------
FUNCTION REAL fpow_lt_one
	(
		IN REAL x		"Input number", 
		IN REAL xlam	"Variable to define the linear zone", 
		IN REAL n		"Exponent"
	)
"Function to calculate signed power with linear zone near 0"

	DECLS
		REAL y	"Result of the signed root with linear zone near 0"
	BODY
           ASSERT (n <=1 AND n >= 0) ERROR "Wrong order of the root in function fpow_lt_one"

		IF (abs(x)>xlam OR xlam==0) THEN
			y = MATH.sign(x) * abs(x)**n
		ELSE
			y = (2. - n) * (x / xlam**(1.-n)) +(n -1.) * (x*abs(x) / xlam**(2.-n))
		END IF

		RETURN y
END FUNCTION

-----------------------------------------------------------------------
--Function to calculate signed power with linear zone near 0
--for exponents larger than 1
-----------------------------------------------------------------------
FUNCTION REAL fpow_gt_one
	(
		IN REAL y		"Number that is raised", 
		IN REAL ylam	"Variable to define the linear zone", 
		IN REAL n		"Exponent"
	)
"Function to calculate signed power with linear zone near 0"

	DECLS
		    REAL x		"Result"

               REAL xlam	"Auxiliary variable to define the linear zone"
               REAL a		"Quadratic coefficient of a quadratic equation" 
               REAL b		"Linear coefficient of a quadratic equation"
               REAL c		"Independent coefficient of a quadratic equation"
	BODY
           ASSERT (n>=1) ERROR "Wrong exponent in function fpow_gt_one"

		IF (abs(y) > ylam OR ylam==0) THEN
			x = MATH.sign(y) * abs(y)**n
		ELSE
                xlam = abs(ylam)**n
		      a = (1./n -1.) / xlam**(2.-1./n)
                b = (2. - 1./n) / xlam**(1.-1/n)
                c = -abs(y)
                x = MATH.sign(y)*(-b + sqrt(abs(b**2 - 4 * a * c))) / (2 * a)
		END IF

		RETURN x
END FUNCTION

-----------------------------------------------------------------------
--Function to calculate signed power with linear zone near 0
--for exponents larger than 1 and less than 1
--  y = fpow(x, xlam, n)
--  x = fpow(y, xlam**n, 1/n) 
-----------------------------------------------------------------------
FUNCTION REAL fpow
	(
		IN REAL y		"Number that is raised", 
		IN REAL ylam	"Variable to define the linear zone", 
		IN REAL n		"Exponent"
	)
	DECLS
		   REAL x		"Result"

	BODY
           ASSERT (n>=0) ERROR "Wrong exponent in function fpow"

		IF(n==1) THEN
			x = y
		ELSEIF(n > 1) THEN
			x = fpow_gt_one(y, ylam, n)
		ELSE
			x = fpow_lt_one(y, ylam, n)
		END IF

	RETURN x
END FUNCTION