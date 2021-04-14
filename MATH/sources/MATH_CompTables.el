--------------------------------------------------------------------------------
-- EA Internacional 2004   EcosimPro 3.3 Simulation Source Code
--
--    FILE NAME:  MATH_CompTables.el
--    FILE TYPE:  Components of the MATH library
--    DESCRIPTION:  Defines a table linear interpolation component for the
--                  MATH library
--    NOTES:  Based on EcosimPro MATH library
--    AUTHOR:  Ramon Perez Vara
--    CREATION DATE:  12-Jan-2004
--------------------------------------------------------------------------------


-- Components
--------------------------------------------------------------------------------
-- Component CompLinearInterp_1D
--------------------------------------------------------------------------------
-- Purpose:
--    To make linear interpolation in tables taking into account
--    discontinuities when points are crossed.
--------------------------------------------------------------------------------
COMPONENT CompLinearInterp1D
   DATA
      TABLE_1D table

      REAL signal_in
      REAL signal_out

   DECLS
      INTEGER i, pi

   DISCRETE
      WHEN cellCrossing1D(table, signal_in, pi, i) THEN
         pi = i
      END WHEN

   CONTINUOUS
      signal_out = cellLinearInterp1D(table, signal_in, pi)

END COMPONENT
