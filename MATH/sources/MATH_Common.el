--------------------------------------------------------------------------------
-- EA Internacional 2004   EcosimPro 3.3 Simulation Source Code
--
--    FILE NAME:  MATH_Common.el
--    FILE TYPE:  Common elements of the MATH library
--    DESCRIPTION:  Defines common constants for the MATH library
--    NOTES:  Based on EcosimPro MATH library
--    AUTHOR:  Ramon Perez Vara
--    CREATION DATE:  12-Jan-2004
--------------------------------------------------------------------------------

-- Constants
CONST REAL ZERO = 0.                                       "Zero value"

CONST REAL Eps = 1.e-8                                     " "
CONST REAL Inf = 1.e38                                     " "

CONST REAL Small = 1.e-16                                  " "

CONST INTEGER Inf_int = 2147483647                         " "

--------------------------------------------------------------------------------
-- Mathematical and cientific constants
--------------------------------------------------------------------------------
CONST REAL E =     2.7182818284590452353602874713527       "Natural logarithmic base"
CONST REAL GAMMA = 0.5772156649015328606065120900824       "Euler's Constant"
CONST REAL PI =    3.1415926535897932384626433832795       "PI, Archimedes' number"



--------------------------------------------------------------------------------
-- Units conversion constants
--------------------------------------------------------------------------------
CONST REAL DtoR = PI / 180.       "Degree to radian"
CONST REAL RtoD = 180. / PI       "Radian to degrees"

--------------------------------------------------------------------------------
-- Thermal Global Variables
--------------------------------------------------------------------------------
CONST REAL STEFAN = 5.6696e-8 UNITS 				"W/m^2*K^4"	"Stefan constant"
CONST REAL TZERO  = 273.15    UNITS u_K			"Celsius scale shift with respect to Kelvin scale"
CONST REAL TMAX   = 1.e5      UNITS u_K			"Maximum Temperature to report a warning"