/*-----------------------------------------------------------------------------------------
 LIBRARY: LIQHAMMER_ST 
 FILE: LIQHAMMER_ST
 CREATION DATE: 04/09/2017
 ------------------------------------------------------------------------------------------
-- LIQHAMMER_ST is an EcosimPro Library for solving hydraulic transients by the Method
-- of Characteristics (MOC). It is extremely efficient, accurate and robust due to the used
-- method.
--
-- There are two versions of the Library: 
--    The academic version,  named LiqHammer_ST, and the professional version, named 
--    LiqHammerPro
-- 
-----------------------------------------------------------------------------------------*/
USE MATH
--Gravity constant
CONST REAL Go = 9.80665    UNITS "m/s2"      "Standard Earth Gravity"
--

BOOLEAN sample = FALSE              "Boolean variable to define MOC sampling"
INTEGER istep = 0                   "Time step number"
DISCR REAL TIME_sample  UNITS "s"   "Time of previous step"
DISCR REAL dt_MOC       UNITS "s"   "Time step for Method of Characteristics"

--Enumerative Data Types
ENUM PIPE_HYDRO_METHOD = {Global, MOC, AlgebraicWaterHammer, LaxWendroff}
SET_OF(PIPE_HYDRO_METHOD) HYDRO_METHOD = {MOC, AlgebraicWaterHammer, LaxWendroff}
ENUM FluxLimiter = {Superbee, Vanleer, Vanalbada, Minbee, FirstOrder, HighOrder, UserDef}
ENUM InterpMethod = {WaveSpeedAdjustment, SpaceLine, TimeLine, MinimumPoint, CharacteristicLine, WavePathAdjustment}

--Global variables to transmit solver parameters to the other components--------------------------
HIDDEN ENUM HYDRO_METHOD GlobalMethod                                                           --
HIDDEN ENUM FluxLimiter GlobalFluxLimiter                                                       --
--                                                                                              --
HIDDEN DISCR REAL rho         UNITS "kg/m3"  "Fluid density"                                    -- 
HIDDEN DISCR REAL visc        UNITS "Pa·s"   "Fluid dynamic viscosity"                          -- 
HIDDEN DISCR REAL P_vap       UNITS "MPa"    "Vapor pressure (MPa)"                             --
HIDDEN DISCR REAL EPS_FRIC    UNITS "p.u."   "Linearization constant for the friction term (-)" --
HIDDEN DISCR REAL theta1 = 0.5  UNITS "p.u." "Weighting coefficient in space for implicit pipes"--
HIDDEN DISCR REAL theta2 = 0.6  UNITS "p.u"  "Weighting coefficient in time implicit pipes"     --
--                                                                                              --
HIDDEN DISCR REAL H_vap       UNITS "m"            "Absolute vapor pressure head"               --
HIDDEN DISCR REAL H_atm    UNITS "m"               "Barometric pressure head"                   --
HIDDEN DISCR REAL P_atm    UNITS "MPa (abs)"       "Atmospheric pressure"                       --
HIDDEN DISCR REAL G  = 9.80665    UNITS "m/s2"     "Local Gravity"                              --
--                                                                                              --
HIDDEN BOOLEAN CavitationFlag       "Boolean to enable Discrete vapor cavitation model if TRUE" --
HIDDEN BOOLEAN TransientFriction    "Boolean to enable Brunone´s transient friction if TRUE"    --
--------------------------------------------------------------------------------------------------

CONST INTEGER MaxPipes = 200        "Maximum number of pipes"
CONST INTEGER MaxSectPerPipe = 500  "Maximum number of nodes per pipe"
--
INTEGER npipe = 0    "Total number of pipes in the model"
--
--------------------------------------------------------------------------------
-- Function int
--------------------------------------------------------------------------------
-- Purpose: to convert real value to integer value
--------------------------------------------------------------------------------
FUNCTION INTEGER int(IN REAL x)
   DECLS
      INTEGER i
   BODY
      i = x
      RETURN i
END FUNCTION
--------------------------------------------------------------------------------
-- Function mod
--------------------------------------------------------------------------------
-- Purpose: to return the reminder after integer division
--------------------------------------------------------------------------------
FUNCTION INTEGER mod(IN INTEGER a, IN INTEGER d)
   DECLS
      INTEGER k
   BODY
      k = a /d
      RETURN a - k * d
END FUNCTION
--------------------------------------------------------------------------------
-- Function dabs
--------------------------------------------------------------------------------
-- Purpose:
--    To obtain the absolute value of a real number
--------------------------------------------------------------------------------
FUNCTION REAL dabs(IN REAL x)
DECLS
   REAL y
BODY
   y = x
   IF (x < 0) THEN
      y = -x
   END IF
   RETURN y
END FUNCTION
--------------------------------------------------------------------------------
-- Function GELB --Translated from Fortran to EcosimPro Language
--------------------------------------------------------------------------------
-- Purpose:
--    To solve a system of simultaneous linear equation with a coefficient matrix
--    of band structure
--------------------------------------------------------------------------------
--Usage of FUNCTION GELB
--   GELB(R,A,M,N,MUD,MLD,EPS,IER)
--Description of parameters
--   R      - M by N right hand side matrix (destroyed);
--            on return R contains the solution of the equations.
--   A      - M by M coefficient matrix with band structure
--            (destroyed).
--   M      - the number of equations in the system.
--   N      - the number of right hand side vectors.
--   MUD    - the number of upper codiagonals (that means
--            codiagonals above main diagonal).
--   MLD    - the number of lower codiagonals (that means
--            codiagonals below main diagonal).
--   EPS    - an input constant which is used as relative
--            tolerance for test on loss of significance.
--   IER    - resulting error parameter coded as follows
--            IER=0  - no error
--            IER=-1 - no result because of wrong input parameters
--                     M,MUD,MLD or because of pivot element
--                     at any elimination step equal to 0.
--            IER=K  - warning due to possible loss of significance
--                     indicated at elimination step K+1,
--                     where pivot element was less than or equal
--                     to the tolerance EPS times absolutely 
--                     greatest element of matrix A.
--Remarks
--   Band matrix A is assumed to be stored rowwise in the first
--   ME successive storage locations of totally needed MA
--   storage locations, where
--     MA=M*MC-ML*(ML+1)/2   and   ME=MA-MU*(MU+1)/2    with
--     MC=MIN(M,1+MUD+MLD),  ML=MC-1-MLD,  MU=MC-1-MUD.
--   Right hand side matrix R is assumed to be stored columnwise
--   in N*M successive storage locations. On return solution
--   matrix R is stored columnwise too.
--   Input parameters M, MUD, MLD, should satisfy the following
--   restrictions     MUD not less than zero
--                    MLD not less than zero
--                    MUD+MLD not greater than 2*M-2.
--   No action besides error message IER=-1 takes place if these
--   restrictions are not satisfied.
--   The procedure gives results if the restrictions on input
--   parameters are satisfied  and if pivot elements at all
--   elimination steps are different from 0. However warning
--   IER=X - if given - indicates possible loss of significance.
--   In case of a well scaled matrix A and appropriate tolerance
--   EPS, IER=X may be interpreted that matrix A has the rank K.
--   No warning is given if matrix A has no lower codiagonal.
--
--   Solution is done by means of GAUSS elimination with
--   column pivoting only, in order to preserve band structure
--   in remaining coefficient matrices.
--
FUNCTION NO_TYPE GELB(
      OUT REAL R[],
      OUT REAL A[],
      IN INTEGER M,
      IN INTEGER N,
      IN INTEGER MUD,
      IN INTEGER MLD,
      IN REAL EPS,
      OUT INTEGER IER)
   DECLS
      INTEGER I, IC,ID, IDST, II, ILR
      INTEGER J, JJ
      INTEGER KST
      INTEGER MC, MU, ML, MR, MZ, MA, NM
      REAL PIV, TB, TOLX
   BODY
   --Test on wrong input parameters
      IF(MLD < 0)THEN
         IER = -1
         RETURN
      END IF
      IF(MUD < 0) THEN
         IER = -1
         RETURN
      END IF
      MC=1+MLD+MUD
      IF(MC+1-M-M > 0 ) THEN
         IER = -1
         RETURN
      END IF
      --
      --     Prepare integer parameters
      --       MC=number of columns in matrix A
      --       MU=number of zeros to be inserted in the first row of matrix A
      --       ML=number of missing elements in last row of matrix A
      --       MR=index of last row in matrix A with MC elements
      --       MZ=total number of zeros to be inserted in matrix A
      --       MA=total number of storage locations necessary for matrix A
      --        NM=number of elements in matrix R
      IF (MC-M > 0) THEN
         MC=M
      END IF
      MU=MC-MUD-1
      ML=MC-MLD-1
      MR=M-ML
      MZ=(MU*(MU+1))/2
      MA=M*MC-(ML*(ML+1))/2
      NM=N*M
      --
      --     Move elements backward and search for absolutely greatest element
      --    (not necessary in case of a matrix without lower codiagonals)
      IER=0
      PIV=0.
      IF(MLD > 0) THEN
         JJ=MA
         J=MA-MZ
         KST=J
         FOR (K IN 1,KST)
            TB=A[J]
            A[JJ]=TB
            TB=abs(TB)
            IF(TB-PIV > 0) THEN 
               PIV=TB
            END IF
            J=J-1
            JJ=JJ-1
         END FOR
         --
         -- Insert zeros in first MU rows (not necessary in case MZ=0)
         IF (MZ > 0) THEN
            JJ=1
            J=1+MZ
            IC=1+MUD  
            FOR (I IN 1,MU)
               FOR (K IN 1,MC)
                  A[JJ]=0.
                  IF(K-IC <= 0) THEN
                     A[JJ]=A[J]
                     J = J+1
                  END IF
                  JJ=JJ+1
               END FOR
               IC=IC+1
            END FOR
         END IF
      END IF
         --
         --    Generate test value for singularity
         --
         TOLX =EPS*PIV
         --
         --     Start decomposition loop
         KST=1
         IDST=MC
         IC=MC-1
         FOR (K IN 1,M)
            IF(K-MR-1 > 0) THEN
               IDST=IDST-1
            END IF
            ID=IDST
            ILR=K+MLD
            IF(ILR-M > 0) THEN
               ILR=M
            END IF
            II=KST
            --
            --   Pivot search in first column (row indexes from I=K up to I=ILR)
            PIV=0.
            FOR(I IN K,ILR)
               TB=abs(A[II])
               IF (TB-PIV > 0) THEN
                  PIV=TB
                  J=I
                  JJ=II
               END IF
               IF(I-MR > 0) THEN
                  ID=ID-1
               END IF
               II=II+ID
            END FOR
            --
            --    Test on singularity
            IF(PIV <0) THEN
               IER = -1
               RETURN
            END IF
            IF(IER == 0) THEN
               IF(PIV-TOLX <=0 ) THEN
                  IER=K-1
               END IF
            END IF
            PIV=1./A[JJ]
            --
            --    Pivot row reduction and row interchange in right hand side R
            ID=J-K
            FOR(I=K; I <=NM; I = I+M)
               II=I+ID
               TB=PIV*R[II]
               R[II]=R[I]
               R[I]=TB
            END FOR
            --
            --  Pivot row reduction and row interchange in coefficient matrix A
            II=KST
            J=JJ+IC
            FOR(I IN JJ,J)
               TB=PIV*A[I]
               A[I]=A[II]
               A[II]=TB 
               II=II+1
            END FOR
            --
            --     Element reduction
            IF(K-ILR < 0) THEN
               ID=KST
               II=K+1
               MU=KST+1
               MZ=KST+IC
               FOR(I IN II,ILR)
                  --     In matrix A
                  ID=ID+MC
                  JJ=I-MR-1
                  IF(JJ > 0) THEN
                     ID=ID-JJ
                  END IF
                  PIV=-A[ID]
                  J=ID+1
                  FOR(JJ IN MU,MZ)
                     A[J-1]=A[J]+PIV*A[JJ]
                     J=J+1
                  END FOR
                  A[J-1]=0.
                  --
                  -- In matrix R
                  J=K
                  FOR(JJ = I; JJ <= NM; JJ = JJ+M)
                     R[JJ]=R[JJ]+PIV*R[J]
                     J=J+M
                  END  FOR
               END FOR
            END IF
            KST=KST+MC  
            IF(ILR-MR >= 0) THEN
               IC=IC-1
            END IF
            ID=K-MR
            IF (ID > 0) THEN
               KST=KST-ID
            END IF
         END FOR
         --     End of decomposition loop
         --
         --     Back substitution
         IF (MC-1 > 0) THEN
            IC=2
            KST=MA+ML-MC+2
            II=M
            FOR (I IN 2,M) 
               KST=KST-MC
               II=II-1
               J=II-MR
               IF (J >0) THEN
                  KST=KST+J
               END IF
               FOR(J=II; J <=NM; J = J+M)
                  TB=R[J]
                  MZ=KST+IC-2
                  ID=J
                  FOR(JJ IN KST,MZ)
                     ID=ID+1
                     TB=TB-A[JJ]*R[ID]
                  END FOR
                  R[J]=TB
               END FOR
               IF (IC-MC < 0) THEN
                  IC=IC+1
               END IF
            END FOR
         END IF      
      RETURN
END FUNCTION

--------------------------------------------------------------------------------
-- Class PipeVars
--------------------------------------------------------------------------------
-- Purpose:
--    To represent the pipe variables in a global common area
--------------------------------------------------------------------------------
CLASS PipeVars
   DECLS
      REAL l         UNITS "m"      "Pipe length (m)"
      REAL dia       UNITS "m"      "Pipe inside diameter (m)"
      REAL a         UNITS "m3/s"   "Wave speed (m/s)"
      REAL fric      UNITS ""       "Friction factor"
      --
      ENUM InterpMethod interp_method
      REAL A         UNITS "m2"     "Pipe flow area"
      REAL Cr        UNITS "p.u."   "Courant number (-)"
      REAL a_mod     UNITS "m/s"    "Modified wave speed"
      REAL tt        UNITS "s"      "Wave travel time in pipe"
      REAL rnr       UNITS ""       "Exact number of reaches in pipe"
      INTEGER nr     UNITS ""       "Actual number of reaches in pipe"
      REAL zeta      UNITS ""       "Degree of space line interpolation"   
      REAL w[4]      UNITS ""       "Interpolation weights for points (i-1,t) (i-1,t-dt) (i,t) and (i,t-dt)"
      REAL error     UNITS "%"      "Discretization error"
   METHODS
   
      METHOD NO_TYPE FillVars(
            IN REAL l_in      UNITS "m"   "Pipe length", 
            IN REAL dia_in    UNITS "m"   "Pipe inside diameter", 
            IN REAL a_in      UNITS "m/s" "Wave speed",  
            IN REAL fric_in   UNITS ""    "friction factor")
         BODY
            l = l_in
            dia = dia_in
            a = a_in
            fric = fric_in
            --
            tt = l/a
            A = 0.25 * PI *dia*+2
      END METHOD
      
      METHOD NO_TYPE GetPars(
            OUT INTEGER nst                  "Number of sections", 
            OUT REAL a_act    UNITS "m/s"    "Actual wave speed after adjustment", 
            OUT REAL A        UNITS "m2"     "Flow area", 
            OUT REAL B        UNITS "s/m2"   "Impedance parameter", 
            OUT REAL R        UNITS "s2/m5"  "Friction parameter(s2/m5)", 
            OUT REAL Courant  UNITS ""       "Courant number", 
            OUT REAL zp       UNITS ""       "Degree of space line interpolation",
            OUT REAL wp[4]    UNITS ""       "Interpolation weights for points (i-1,t) (i-1,t-dt) (i,t) and (i,t-dt)",
            OUT ENUM InterpMethod m_interp   "Interpolation method")
         DECLS
            REAL dx UNITS "m"    "distance between nodes"
         BODY
            nst = nr + 1 --number of sections equal to number of reaches + 1
            a_act = a_mod
            dx = l/nr
            A = 0.25 * PI * dia**2
            B = a_mod /(G * A)
            R = fric*dx/(2*G*dia*A**2)
            FOR(j IN 1, 4)
               wp[j] = w[j]
            END FOR
            zp = zeta
            Courant = Cr
            m_interp = interp_method
      END METHOD
      METHOD NO_TYPE Write()
         BODY
            WRITE("nr = %d Courant = %g  wp[1] = %g  wp[2] = %g  wp[3] = %g  wp[4] = %g  zp = %g\n",  nr, Cr, w[1], w[2], w[3], w[4], zeta)
      END METHOD
END CLASS

PipeVars pipe[MaxPipes]
----------------------------------------------------------------------------------
-- Function FlexDiscretization
------------------------------------------------------------------------------------
-- Purpose:
--    Function that implements the flexible discretization algorithm described in
--    the following refernce
--
--    Karney, B. W., and Ghidaoui, M. S. (1997), ‘"Flexible discretization algorithm
--    for fixed-grid MOC in pipelines." J. Hydraul. Eng., 123(11), pp. 1004–1011.
------------------------------------------------------------------------------------
FUNCTION BOOLEAN FlexDiscretization(
            IN ENUM HYDRO_METHOD global_method        "Global method",
            IN ENUM InterpMethod interp_method        "Interpolation Method", 
            IN REAL PSIMAX             UNITS "%"      "User specified limit on maximum allowable wave speed adjustment",
            IN REAL Cr_TL              UNITS "p.u."   "Parameter used to set a time line threshold for the courant number",
            IN REAL dt                 UNITS "s"      "Integration step",
            OUT INTEGER ntot_sect                     "Total number of sections",
            OUT REAL error_max         UNITS "%"      "Maximum discretization error",
            OUT REAL error_avg         UNITS "%"      "Average discretization error")
   DECLS
      INTEGER nrtot  "Total number of reaches"
      REAL PSI       "Allowable shift in Courant number (-)"
      REAL xi        "Degree of time line interpolation"
   BODY
      ntot_sect = 0     --Total number of sections or nodes
      nrtot = 0      --Total number of reaches
      error_max = 0
      error_avg = 0
      FOR (i IN 1, npipe)
         pipe[i].rnr = pipe[i].tt / dt          --exact (real) number of reaches
         pipe[i].nr = pipe[i].rnr + 0.5         --rounded (integer) number of reaches 
         pipe[i].Cr = pipe[i].nr / pipe[i].rnr  --Courant number

         IF(pipe[i].Cr > 1.000001 /(1 - 0.01*PSIMAX)) THEN 
            IF(pipe[i].nr >= 2) THEN
               pipe[i].nr = pipe[i].nr - 1
               pipe[i].Cr = pipe[i].nr / pipe[i].rnr
            ELSE
               --Discretization does not fulfill Cr <= 1 
               RETURN FALSE
            END IF
         END IF
         --
         PSI = 0.01 * PSIMAX * pipe[i].Cr --allowable shift in Courant number
         IF(abs(1. - pipe[i].Cr)< PSI) THEN 
            pipe[i].Cr = 1
         ELSEIF(interp_method == WaveSpeedAdjustment) THEN
            --Discretization does not satisfy required wave speed adjustment
            RETURN FALSE
         ELSEIF(pipe[i].Cr <= Cr_TL + PSI) THEN
            IF((pipe[i].Cr - PSI) > 0.5) THEN
               pipe[i].Cr = pipe[i].Cr - PSI
            ELSEIF((pipe[i].Cr + PSI) > 0.5) THEN 
               pipe[i].Cr = 0.5
            END IF
         ELSE
            pipe[i].Cr = pipe[i].Cr + PSI --increase Courant number as much as possible
         END IF
         --Obtain modified wave speed
         pipe[i].a_mod = pipe[i].Cr * pipe[i].l / (pipe[i].nr * dt)
         nrtot = nrtot + pipe[i].nr
         ntot_sect = ntot_sect + pipe[i].nr + 1
         pipe[i].zeta = 0.0
         xi = 0
         pipe[i].interp_method = interp_method
         IF (interp_method == WaveSpeedAdjustment OR interp_method == WavePathAdjustment) THEN
            pipe[i].zeta = 0
            pipe[i].Cr = 1
         ELSEIF(interp_method == SpaceLine AND pipe[i].Cr > Cr_TL) THEN
            pipe[i].zeta = 1 - pipe[i].Cr
         ELSEIF(interp_method == MinimumPoint AND pipe[i].Cr > Cr_TL) THEN
            pipe[i].zeta = (1. - pipe[i].Cr) / (1. + pipe[i].Cr**2)
            xi = pipe[i].Cr * pipe[i].zeta
         ELSEIF(interp_method == CharacteristicLine AND pipe[i].Cr > Cr_TL) THEN
            pipe[i].zeta = 0.5 * (1. - pipe[i].Cr)
            xi = (1. - pipe[i].Cr) / (2. * pipe[i].Cr)
         ELSE --TimeLine Interpolation or Cr[i]< Cr_TL
            pipe[i].interp_method = TimeLine
            xi = (1.0 - pipe[i].Cr)/pipe[i].Cr
         END IF
         pipe[i].w[1] = (1. - pipe[i].zeta) * (1 - xi)
         pipe[i].w[2] = (1. - pipe[i].zeta) * xi
         pipe[i].w[3] = pipe[i].zeta * (1 - xi)
         pipe[i].w[4] = pipe[i].zeta * xi
         pipe[i].error = abs(1 - (pipe[i].a_mod/pipe[i].a)) 
         IF (pipe[i].Cr > Cr_TL) THEN
            pipe[i].error = pipe[i].error + abs(1 - pipe[i].Cr)
         ELSE
            pipe[i].error = pipe[i].error + abs(0.5-pipe[i].Cr)
         END IF
         pipe[i].error = pipe[i].error * 100.
         error_max = max(error_max, pipe[i].error) 
         error_avg = error_avg + pipe[i].error * pipe[i].nr
      END FOR
      error_max = error_max 
      error_avg = (error_avg / nrtot) 
      RETURN TRUE
END FUNCTION

COMPONENT SolverControl "MOC Solver Control" 
   DATA
      --Solver data
      ENUM HYDRO_METHOD global_method = MOC "Global integration method"
      ENUM InterpMethod interp_method = WaveSpeedAdjustment "Sectioning Method for MOC application"
      INTEGER nrsp_min = 1    RANGE 1, MaxSectPerPipe "Minimum number of reaches in shortest pipe"
      INTEGER nrsp_max = 10   RANGE 1, MaxSectPerPipe "Maximum number of reaches in shortest pipe"
      REAL PSIMAX = 10        UNITS "%"      "User specified limit on maximum allowable wave speed adjustment"
      REAL perc_incr = 0.1    UNITS "%"      "Percentage increment for searching the optimum time step"
      REAL Cr_TL = 0.55       UNITS "p.u."   "Parameter used to set a time line threshold for the courant number"
      REAL eps_fric = 0.81    UNITS "p.u."   "Linearization constant for the friction term"
      BOOLEAN cavitation_flag = TRUE         "Boolean flag to activate discrete cavitation model"
      BOOLEAN transient_friction = FALSE     "Boolean to activate transient friction calculation if true"
      REAL theta1 = 0.5  UNITS "p.u."        "Weighting coefficient in space for implicit pipes"--
      REAL theta2 = 0.6  UNITS "p.u"         "Weighting coefficient in time implicit pipes"     --    
      --Fluid data
      REAL rho  = 998.        UNITS "kg/m3"  "Fluid density"
      REAL visc = 0.01        UNITS "Pa·s"   "Fluid dynamic viscosity"
      REAL P_vap = 2338.8e-6  UNITS "MPa"    "Vapor pressure "
      -----------------------------------------------------------------------
      --System data
      REAL P_atm = 0.101325    UNITS "MPa"   RANGE 0., 10000.     "Atmospheric pressure"
      REAL G = 9.80665         UNITS "m/s2"  RANGE 0.01, 10000.   "Local gravity"
   DECLS

      INTEGER nrtot                       "Total number of reaches"
      INTEGER ntot_sect                   "Total number of sections"
      --
      DISCR REAL dt_optim     UNITS "s"   "Optimum time step"
      INTEGER nrsp_optim                  "Number of reaches in pipe with shortest travel time for optimum discretization"
      INTEGER ntot_sect_optim             "Number of reaches in pipe with shortest travel time for optimum discretization"
      DISCR REAL error_max_optim
      DISCR REAL error_avg_optim
      --
      DISCR REAL cpu = 0      UNITS "s"    "CPU TIME"

      HIDDEN BOOLEAN valid_discretization 
      HIDDEN BOOLEAN found 
      
      HIDDEN INTEGER ip_max               "Number of the pipe with maximum travel time"
      HIDDEN INTEGER ip_min               "Number of the pipe with minimum travel time"
      HIDDEN INTEGER nrsp                 "Number of reaches in pipe with shortes travel time"
      HIDDEN INTEGER n_max             "Maximum index for nrsp search"
      HIDDEN INTEGER n_min             "Minimum index for nrsp search"
      --
      HIDDEN DISCR REAL PSI               "Allowable shift in Courant number (-)"
      HIDDEN DISCR REAL dt_local         UNITS "s" "Local time step"
      HIDDEN DISCR REAL dt_max            UNITS "s"   "Maximum time step of the search"
      HIDDEN DISCR REAL dt_min            UNITS "s"   "Minimum time step of the search"
      HIDDEN DISCR REAL dt_min_prev       UNITS "s"   "Previous minimum time step of the search"
      HIDDEN DISCR REAL error_max         UNITS "%"   "Maximum discretization error"
      HIDDEN DISCR REAL error_avg         UNITS "%"   "Average discretization error"
      HIDDEN DISCR REAL tt_max            UNITS "s"   "Maximum wave travel time"
      HIDDEN DISCR REAL tt_min            UNITS "s"   "Minimum wave travel time"
      HIDDEN DISCR REAL tt_tot            UNITS "s"   "Total travel time"
   INIT PRIORITY -100
      --Low priority INIT block to postpone its execution after 
      --the INIT block of the pipe componets that initialize the pipe class
      
      --Reset cpu time to start the cpu count
      resetTime()
      
      --Assign data needed by the components to global variables for their transmission
      LIQHAMMER_ST.rho = rho
      LIQHAMMER_ST.visc = visc
      LIQHAMMER_ST.H_vap = P_vap * 1e6 / (rho * G)
      LIQHAMMER_ST.H_atm = P_atm * 1e6 / (rho * G)
      LIQHAMMER_ST.P_atm = P_atm
      LIQHAMMER_ST.EPS_FRIC = eps_fric
      LIQHAMMER_ST.CavitationFlag = cavitation_flag
      LIQHAMMER_ST.TransientFriction = transient_friction
      LIQHAMMER_ST.GlobalMethod = global_method
      LIQHAMMER_ST.theta1 = theta1
      LIQHAMMER_ST.theta2 = theta2
      --
      
      found = FALSE
      --Determination of minimum wave travel time
      tt_min = pipe[1].tt
      tt_max = pipe[1].tt
      tt_tot = pipe[1].tt
      ip_max = 1
      ip_min = 1
      FOR (i IN 2, npipe)
         tt_tot = tt_tot + pipe[i].tt
         IF (pipe[i].tt < tt_min) THEN
            tt_min = pipe[i].tt
            ip_min = i
         END IF
         IF (pipe[i].tt > tt_max) THEN
            tt_max = pipe[i].tt
            ip_max = i
         END IF
      END FOR     

      ASSERT (nrsp_min <= nrsp_max) FATAL "Wrong data in component of SolverControl type, nrsp_max < nrsp_min\n"
      --
      nrsp = nrsp_min
      WHILE(NOT valid_discretization AND nrsp <= nrsp_max)
         IF(nrsp > nrsp_min) THEN 
            dt_min_prev = dt_min
         END IF
         dt_min = max((tt_min / nrsp) / (1 + PSIMAX/100.), tt_min/(nrsp + 0.5))
         dt_max = min((tt_min / nrsp) / (1 - PSIMAX/100.), tt_min/(nrsp - 0.5))
         IF (nrsp > nrsp_min) THEN
            dt_max = min(dt_max, dt_min_prev)
         END IF
         
         n_max = ((tt_min/nrsp)/dt_min-1)*100/perc_incr 
         n_min = ((tt_min/nrsp)/dt_max -1)*100/perc_incr
         
         error_max_optim  =  1e6
         --Loop to search the optimum time increment
         FOR(j IN  n_min, n_max)
            dt_local = (tt_min/nrsp) / (1 +  j * 0.01 * perc_incr)         
            valid_discretization = FlexDiscretization(global_method, interp_method, PSIMAX, Cr_TL, dt_local, ntot_sect, error_max, error_avg)         
            IF(valid_discretization == TRUE) THEN
               found = TRUE
               IF(error_max < error_max_optim) THEN
                  error_max_optim = error_max
                  dt_optim = dt_local
               END IF
            END IF
         END FOR
         nrsp = nrsp + 1
      END WHILE
      
      IF(found == TRUE) THEN
         --MOC time increment is assigned the optimum time step
         dt_MOC = dt_optim
         FlexDiscretization(global_method, interp_method, PSIMAX, Cr_TL, dt_optim, ntot_sect, error_max_optim, error_avg_optim)
      ELSE
         STOP "*===== FATAL SECTIONING ERROR: NO VALID DISCRETIZATION HAS BEEN FOUND. TRY TO INCREASE nrsp_max or to increase PSIMAX====*" 
      END IF
      sample = TRUE
      istep = 0
      IF(CINT <  dt_MOC) THEN
         WRITE("*===== Comunication interval changed to %lg because it was lower than MOC time step\n", dt_MOC)
         CINT = dt_MOC 
      END IF
   DISCRETE
         WHEN (sample) THEN
            sample = FALSE AFTER 0
            TIME_sample = TIME
            sample = TRUE AFTER dt_MOC
            istep = istep + 1 AFTER 0
            cpu = getTime()
         END WHEN
END COMPONENT
--------------------------------------------------------------------------
--Fluid Port Type
--------------------------------------------------------------------------
PORT fluid
    SUM REAL Q          UNITS "m3/s"   "Volume flow through the port"
    EQUAL REAL H        UNITS "m"      "Piezometric head at the port"
    DISCR EQUAL REAL z  UNITS "m"      "Elevation at the port"    
END PORT  
--------------------------------------------------------------------------
FUNCTION NO_TYPE FluxSourceCalc(IN REAL c, IN REAL g, IN REAL dia, IN REAL A, IN REAL f, IN REAL u[2], OUT REAL flux[2], OUT REAL source[2])
      BODY
         flux[1] = (c**2/(g*A)) * u[2]
         flux[2] = A * g * u[1]
         source[1] = 0
         source[2] = -f * u[2] * dabs(u[2]) / (2 * dia  *  A)
         RETURN
END FUNCTION

DISCR REAL cpu_moc_end_nodes
DISCR REAL cpu_awh_end_nodes
DISCR REAL cpu_moc_inner_nodes

INTEGER ncall_moc_end_nodes
INTEGER ncall_awh_end_nodes
INTEGER ncall_moc_inner_nodes

CLASS MOC_Nodes
   DECLS
      REAL HN[MaxSectPerPipe]       UNITS "m"      "Head in section i at time t + dt"
      REAL QUN[MaxSectPerPipe]      UNITS "m3/s"   "Volume flow in section i at time t + dt"
      REAL QDN[MaxSectPerPipe]      UNITS "m3/s"   "Volume flow in section i at time t + dt"
      REAL H[MaxSectPerPipe]        UNITS "m"      "Head in section i at time t"
      REAL QU[MaxSectPerPipe]       UNITS "m3/s"   "Volume flow in section i at time t"
      REAL QD[MaxSectPerPipe]       UNITS "m3/s"   "Volume flow in section i at time t"
      REAL HO[MaxSectPerPipe]       UNITS "m"      "Head in section i at time t - dt"
      REAL QUO[MaxSectPerPipe]      UNITS "m3/s"   "Volume flow in section i at time t - dt"    
      REAL VCAV[MaxSectPerPipe]     UNITS "m3"     "Volume of vapor cavity in section i at time t + dt "
      REAL QDO[MaxSectPerPipe]      UNITS "m3/s"   "Volume flow in section i at time t - dt"
      REAL P[MaxSectPerPipe]        UNITS "MPa"    "Static pressure"  
      REAL z[MaxSectPerPipe]        UNITS "m"      "Elevation of the pipe nodes"
      BOOLEAN icav[MaxSectPerPipe]  "Cavitatiton flag"
      --Variables for applying algebraic water hammer
      REAL H_in[MaxSectPerPipe+1]   UNITS "m"      "Circular Buffer to store piezometric head at inlet"
      REAL Q_in[MaxSectPerPipe+1]   UNITS "m3/s"   "Circular Buffer to store flow at inlet"
      REAL H_out[MaxSectPerPipe+1]  UNITS "m"      "Circular Buffer to store piezometric head at outlet"
      REAL Q_out[MaxSectPerPipe+1]  UNITS "m3/s"   "Circular Buffer to store flow at outlet"
      --Interpolated points for applying characteristic equation at fist and last node
      REAL HS1                      UNITS "m"      "Interpolated Head at point S for 1st node"
      REAL QS1                      UNITS "m3/s"   "Interpolated Flow at point S for 1st node"
      REAL HRn                      UNITS "m"      "Interpolated Head at point R for last node"
      REAL QRn                      UNITS "m3/s"   "Interpolated Flow at point R for last node"
      INTEGER i
   METHODS
   --Method to initialize the nodes
   METHOD NO_TYPE InitNodes(
         IN INTEGER nst                "Number of sections", 
         IN REAL z_in   UNITS "m"      "Inlet elevation", 
         IN REAL z_out  UNITS "m"      "Outlet elevation",
         IN REAL R      UNITS "s2/m5"  "Friction parameter(s2/m5)", 
         IN REAL Ho     UNITS "m"      "Initial head at pipe inlet", 
         IN REAL Qo     UNITS "m3/s"   "Initial volume flow")
      BODY
         HO[1] = Ho
         H[1]  = Ho
         HN[1] = Ho
         QUO[1] = Qo
         QU[1]  = Qo
         QUN[1] = Qo
         QDO[1] = Qo
         QD[1]  = Qo
         QDN[1] = Qo
         z[1] = z_in
         icav[1] = FALSE
         FOR (i IN 2, nst)       
            HO[i] = HO[i-1] - R * Qo * dabs(Qo)
            H[i]  = H[i-1] - R * Qo * dabs(Qo)
            HN[i] = HN[i-1] - R * Qo * dabs(Qo)
            --
            QUO[i]  = Qo
            QU[i]  = Qo
            QUN[i] = Qo
            QDO[i]  = Qo
            QD[i]  = Qo
            QDN[i] = Qo
            icav[i] = FALSE
            z[i] = z_in + (z_out - z_in) * (i-1) / (nst-1)
         END FOR
      END METHOD
      --Method to apply the MOC (method of characteristics) at the inner nodes
      METHOD NO_TYPE MOC_InnerNodes(
            IN INTEGER nst    "Number of sections", 
            IN REAL f_in_H    UNITS "m"      "Piezometric head at inlet",
            IN REAL f_in_Q    UNITS "m3/s"   "Volume flow at inlet", 
            IN REAL f_out_H   UNITS "m"      "Piezometric head at outlet", 
            IN REAL f_out_Q   UNITS "m3/s"   "Volume flow at outlet",   
            IN REAL B         UNITS "s/m2"   "Impedance parameter", 
            IN REAL R         UNITS "s2/m5"  "Friction parameter(s2/m5)", 
            IN REAL dia       UNITS "m"      "Inside diameter", 
            IN REAL A         UNITS "m2"     "Flow area", 
            IN REAL zeta      UNITS ""       "Degree of space line interpolation",  
            IN REAL w[4]      UNITS ""       "Interpolation weights for points (i-1,t) (i-1,t-dt) (i,t) and (i,t-dt)",
            OUT REAL P_max    UNITS "MPa"    "Pipe maximum pressure", 
            OUT REAL P_min    UNITS "MPa"    "Pipe minimum pressure"
            )
         DECLS
            REAL cpu UNITS "s"      "CPU time"
            REAL HR  UNITS "m"      "Piezometric head at point R"
            REAL QR  UNITS "m3/s"   "Volume flow at point R"
            REAL HS  UNITS "m"      "Piezometric head at point S"
            REAL QS  UNITS "m3/s"   "Volume flow at point S"
            REAL Bp  UNITS "s/m2"   "Linear coeff. of C+ compatibiity equation"
            REAL Cp  UNITS "m"      "Const. coeff. of C+ compatibiity equation"
            REAL Bm  UNITS "s/m2"   "Linear coeff. of C- compatibiity equation"
            REAL Cm  UNITS "m"      "Const. coeff. of C- compatibiity equation"
            REAL C   UNITS ""       "Vardy's shear decay coefficient"
            REAL k   UNITS ""       "Decay coefficient of transient friction model"
            REAL Re  UNITS ""       "Reynolds number"
         BODY
            cpu = getTime()
            ncall_moc_inner_nodes = ncall_moc_inner_nodes + 1
            
            --Nodal update
            IF (istep > 0) THEN

            FOR(i IN 1, nst)
               HO[i] = H[i]
               QUO[i] = QU[i]
               QDO[i] = QD[i]
            END FOR
            H[1] = f_in_H
            H[nst] = f_out_H
            QD[1] = f_in_Q
            QU[nst] = f_out_Q
            FOR(i IN 2, nst-1)
               H[i] = HN[i]
               QU[i] = QUN[i]
               QD[i] = QDN[i]
            END FOR
            END IF
      FOR (i IN 2, nst-1)
         IF(TransientFriction == TRUE) THEN
            --Hamid Shamloo, Reyhaneh Nooroz and Maryam Mousavifard 
            --"A review of one-dimensional unsteady friction models for transient pipe flow"
            Re = rho * 0.5 * dabs(QUO[i] + QDO[i]) * dia / (A * visc)
            IF(Re < 1616.393) THEN
               C = 4.76e-3
            ELSE
               C = 7.41 / Re ** (log10(14.3/Re**0.05))
            END IF
            k = 0.5 * C**0.5
         ELSE 
            C = 0
            k = 0
         END IF
         
         HR = H[i-1] * w[1] + HO[i-1] * w[2] + H[i] * w[3] + HO[i] * w[4] 
         QR = QD[i-1] * w[1] + QDO[i-1] * w[2] + QU[i] * w[3] + QUO[i] * w[4] 
         HS = H[i+1] * w[1] + HO[i+1] * w[2] + H[i] * w[3] + HO[i] * w[4]  
         QS = QU[i+1] * w[1] + QUO[i+1] * w[2] + QD[i] * w[3] + QDO[i] * w[4] 
         
         Cp = HR + QR * (B - R * (1.-zeta) * dabs(QR) * (1. - EPS_FRIC))  \
                  +  0.5 * k * B * (QUO[i]  - sign(QR)* dabs(QS - QR))

         Bp = B * (1 + 0.5 * k) + R * (1 - zeta) * dabs(QR) * EPS_FRIC

         Cm = HS - QS * (B -  R * (1.-zeta) * dabs(QS) * (1. - EPS_FRIC)) \
                  - 0.5 * k * B * (QDO[i]  - sign(QS) * dabs(QS - QR))

         Bm = B * (1 + 0.5 * k) + R * (1 - zeta) * dabs(QS) * EPS_FRIC
         --Calculation of the new heads and flows at innner point 
         --in Time+DT
         IF(icav[i] == FALSE) THEN 
            HN[i] = (Cp * Bm  + Cm * Bp) / (Bp + Bm) 
            QUN[i] = (Cp - Cm) / (Bp + Bm)
            QDN[i] = QUN[i]
            IF (LIQHAMMER_ST.CavitationFlag == TRUE AND HN[i] < z[i] - H_atm + H_vap) THEN
               icav[i] = TRUE
               HN[i] = z[i] - H_atm + H_vap
               QUN[i] = (Cp - HN[i]) / Bp
               QDN[i]  = (HN[i] - Cm) / Bm
               VCAV[i] =  0.5 * dt_MOC *(QDN[i] - QUN[i])
            END IF
         ELSE
            HN[i] = z[i] - H_atm + H_vap
            QUN[i] = (Cp - HN[i]) / Bp
            QDN[i]  = (HN[i] - Cm) / Bm
            VCAV[i] = VCAV[i] + 0.5 * dt_MOC *(QDN[i] + QD[i] - QUN[i] - QU[i])
            IF(VCAV[i] < 0.) THEN
               VCAV[i] = 0.
               icav[i] = FALSE
               HN[i] = (Cp * Bm  + Cm * Bp)/(Bp + Bm)
               QUN[i] = (Cp - HN[i]) / Bp
               QDN[i] = QUN[i]
            END IF
         END IF
      END FOR
      FOR (i IN 1, nst)
         P[i] = rho * G * H[i] / 1e6
         P_max = max(P[i], P_max)
         P_min = min(P[i], P_min)
      END FOR

         cpu = getTime() - cpu
         cpu_moc_inner_nodes = cpu_moc_inner_nodes + cpu
      END METHOD
         
      METHOD NO_TYPE AWH_EndNodes(
            IN INTEGER nst                "Number of sections", 
            IN REAL chi    UNITS ""       "Degree of time interpolation",
            IN REAL B      UNITS "s/m2"   "Impedance parameter", 
            IN REAL R      UNITS "s2/m5"  "Friction parameter(s2/m5)", 
            OUT REAL Bm1   UNITS "s/m2"   "Linear coeff. of C- compatibiity equation at 1st section", 
            OUT REAL Cm1   UNITS "m"      "Constant coeff. of C- compatibiity equation at 1st section",
            OUT REAL Bpn   UNITS "s/m2"   "Linear coeff. of C+ compatibiity equation at last section", 
            OUT REAL Cpn   UNITS "m"      "Constant coeff. of C+ compatibiity equation at last section"
            )
      DECLS
         INTEGER i,  i1
         INTEGER nr
         REAL chia
         REAL cpu
         --Recommended Linearization term for algebraic water hammer is zero
         CONST REAL eps_fric = 0  
      BODY
         cpu = getTime()
         ncall_awh_end_nodes = ncall_awh_end_nodes + 1
         nr = nst - 1

         i = mod(istep, nst+1) + 1
         IF(i < nst+1) THEN
            i = i+1
         ELSE
            i = 1
         END IF
         IF (chi >= 0) THEN
            IF(i != nst+1) THEN
               i1 = i+1 
            ELSE
               i1 = 1
            END IF
            HS1 = (1-chi)*H_out[i] + chi*H_out[i1]
            QS1 = (1-chi)*Q_out[i] + chi*Q_out[i1]
            --Head and volume flow interpolated in R
            HRn = chi*H_in[i1] + (1-chi)*H_in[i]
            QRn = chi*Q_in[i1] + (1-chi)*Q_in[i]
         ELSE 
            IF (i != 1) THEN
               i1 = i-1
            ELSE
               i1 = nst+1
            END IF
            chia = dabs(chi)
            HS1 = (1-chia)*H_out[i] + chia*H_out[i1]
            QS1 = (1-chia)*Q_out[i] + chia*Q_out[i1]
            --Head and volume flow interpolated in R
            HRn = chia*H_in[i1] + (1-chia)*H_in[i]
            QRn = chia*Q_in[i1] + (1-chia)*Q_in[i]
         END IF
         --Coefficients of the characteristic equation at the first node
         Cm1 =  HS1 - QS1 * (B - nr * R * dabs(QS1) * (1. - eps_fric))
         Bm1 = B + chi * R * dabs(QS1) * eps_fric
         --Coefficients of the characteristic equation at the last node
         Cpn =  HRn + QRn * (B - nr * R * dabs(QRn) * (1. - eps_fric))
         Bpn = B + nr * R * dabs(QRn) * eps_fric
         cpu = getTime() - cpu
         cpu_awh_end_nodes = cpu_awh_end_nodes + cpu
      END METHOD
      
      METHOD NO_TYPE MOC_EndNodes(
            IN INTEGER nst    "Number of sections in pipe", 
            IN REAL chi    UNITS ""       "Degree of space line interpolation", 
            IN REAL B      UNITS "s/m2"   "Impedance parameter", 
            IN REAL R      UNITS "s2/m5"  "Friction parameter(s2/m5)", 
            IN REAL D      UNITS "m"      "Diameter", 
            IN REAL A      UNITS "m2"     "Flow area",  
            OUT REAL Bm1   UNITS "s/m2"   "Linear coeff. of C- compatibiity equation at 1st section", 
            OUT REAL Cm1   UNITS "m"      "Constant coeff. of C- compatibiity equation at 1st section",
            OUT REAL Bpn   UNITS "s/m2"   "Linear coeff. of C+ compatibiity equation at last section", 
            OUT REAL Cpn   UNITS "m"      "Constant coeff. of C+ compatibiity equation at last section"
            )
         DECLS
            REAL C
            REAL k
            REAL Re "Reynolds number at first and last node"
            REAL cpu
         BODY
            ncall_moc_end_nodes = ncall_moc_end_nodes + 1
            cpu = getTime()
            --Interpolation of Point S at 1st node
            HS1 = (1-chi)*H[1] + chi*H[2]
            QS1 = (1-chi)*QD[1] + chi*QU[2]
            --Interpolation of Point R at last node
            HRn = chi*H[nst-1] + (1-chi)*H[nst]
            QRn = chi*QD[nst-1] + (1-chi)*QU[nst]

            IF(TransientFriction == TRUE) THEN
               --Coefficients of the characteristic equation at the first node
               Re = dabs(QD[1]) * D / (A * visc)
               IF(Re < 1616.393) THEN
                  C = 4.76e-3
               ELSE
                  C = 7.41 / Re ** (log10(14.3/Re**0.05))
               END IF
               k = 0.5 * C**0.5
               Cm1 =  HS1 - QS1 * (B - chi * R * dabs(QS1) * (1. - EPS_FRIC)) \
                     -  k * B * (QD[1]  - sign(QS1) * dabs(QS1 - QD[1]))
               Bm1 = B * (1 + k) + chi * R * dabs(QS1) * EPS_FRIC    
               --Coefficients of the characteristic equation at the last node
               Re = dabs(QU[nst]) * D / (A * visc)
               IF(Re < 1616.393) THEN
                  C = 4.76e-3
               ELSE
                  C = 7.41 / Re ** (log10(14.3/Re**0.05))
               END IF
               k = 0.5 * C**0.5
               Cpn =  HRn + QRn * (B - chi * R * dabs(QRn) * (1. - EPS_FRIC)) \
                     +  k * B * (QU[nst]  - sign(QRn)* dabs(QU[nst] - QRn))
               Bpn =  B * (1+ k ) + chi * R * dabs(QRn) * EPS_FRIC      
            ELSE 
               --Coefficients of the characteristic equation at the first node
               Cm1 =  HS1 - QS1 * (B - chi * R * dabs(QS1) * (1. - EPS_FRIC)) 
               Bm1 = B  + chi * R * dabs(QS1) * EPS_FRIC
               --Coefficients of the characteristic equation at the last node
               Cpn =  HRn + QRn * (B - chi * R * dabs(QRn) * (1. - EPS_FRIC))
               Bpn = B + chi * R * dabs(QRn) * EPS_FRIC
            END IF
            cpu = getTime() - cpu
            cpu_moc_end_nodes = cpu_moc_end_nodes + cpu
            RETURN
      END METHOD
END CLASS

ABSTRACT COMPONENT AbstractPipe  "Absract Pipe Component"
   PORTS
      IN  LIQHAMMER_ST.fluid f_in      "fluid inlet"
      OUT LIQHAMMER_ST.fluid f_out  "fluid outlet"
   DATA

      REAL l=600.       UNITS "m"      "Pipe length"
      REAL dia=0.5      UNITS "m"      "Internal pipe diameter"
      REAL f = 0.018    UNITS "-"      "Friction pressure loss coefficient"
      REAL a = 1200.    UNITS "m/s"    "Fluid speed of sound"
      REAL Ho = 150     UNITS "m"      "Initial head at pipe inlet"
      REAL Qo = 0.477   UNITS "m3/s"   "Initial flow"
END COMPONENT

COMPONENT Pipe IS_A AbstractPipe "Pipe using explicit MOC or Algebraic water hammer"
   DATA
      ENUM PIPE_HYDRO_METHOD method = Global "Solution method"
   DECLS
      ENUM PIPE_HYDRO_METHOD actual_method = Global
      ENUM InterpMethod m_interp
      INTEGER ipipe        "Pipe internal number" 
      DISCR REAL a_mod     UNITS "m/s"    "Adjusted wave speed"
      DISCR REAL A         UNITS "m2"     "Pipe cross-sectional area"
      DISCR REAL B         UNITS "s/m2"   "Impedance parameter"
      DISCR REAL Cr        UNITS ""       "Courant mumber (-)"

      INTEGER nst          "Number of sections"
      INTEGER nr           "Number of reaches"


      DISCR REAL P_max =  -1.e10 UNITS "MPa" "Maximum Pressure in pipe"
      DISCR REAL P_min =   1.e10 UNITS "MPa" "Minimum pressure in pipe"
      
      --HIDDEN Variables
      HIDDEN INTEGER i                          "Position to be updated in the queue for water hammer algebraic calculation"
      HIDDEN DISCR REAL R        UNITS "s2/m5"  "Friction parameter(s2/m5)"
      HIDDEN DISCR REAL w[4]     UNITS ""       "Interpolation weights for points (i-1,t) (i-1,t-dt) (i,t) and (i,t-dt)"
      HIDDEN DISCR REAL zp       UNITS ""       "Degree of space line interpolation"   
      HIDDEN REAL Cm1            UNITS "m"      "Const. coeff. of C- compatibiity equation in first node"
      HIDDEN REAL Bm1            UNITS "s/m2"   "Linear coeff. of C- compatibiity equation in first node"
      HIDDEN REAL Cpn            UNITS "m"      "Const. coeff. of C+ compatibiity equation in last node"
      HIDDEN REAL Bpn            UNITS "s/m2"   "Linear coeff. of C+ compatibiity equation in last node"    

   OBJECTS
      MOC_Nodes nodes
   INIT PRIORITY -50
      IF (ipipe == 0) THEN
         npipe = npipe + 1 
         ipipe = npipe  
      END IF
      pipe[ipipe].FillVars(l, dia, a, f)
   DISCRETE
      WHEN(sample) THEN
         --Initialization at time step 0
         IF(istep == 0) THEN
            IF(method == Global) THEN
               actual_method = GlobalMethod
            ELSE
               actual_method = method
            END IF
            pipe[ipipe].GetPars(nst, a_mod, A, B, R, Cr, zp, w, m_interp)
            nr = nst-1
            nodes.InitNodes(nst, f_in.z, f_out.z , R, Ho, Qo)
            FOR(j IN 1, nst+1)
               nodes.H_in[j] = Ho
               nodes.Q_in[j] = Qo
               nodes.H_out[j] = Ho - nr * R * Qo * dabs(Qo)
               nodes.Q_out[j] = Qo
            END FOR
         END IF
         --End of initialization------
         IF (actual_method == MOC) THEN
               nodes.MOC_InnerNodes(nst, f_in.H, f_in.Q, f_out.H, f_out.Q, B, R, dia, A, zp, w, P_max, P_min)
         ELSEIF(actual_method == AlgebraicWaterHammer AND istep > 0) THEN
               --PRINT("sample - algebraic waterhammer")
               i = mod(istep, nst+1) + 1
               nodes.H_in[i] = f_in.H 
               nodes.Q_in[i] = f_in.Q  
               nodes.H_out[i] = f_out.H  
               nodes.Q_out[i] = f_out.Q
         END IF
      END WHEN
   CONTINUOUS  
      --Head and volume flow interpolate
      SEQUENTIAL
         IF(actual_method == MOC OR actual_method == LaxWendroff) THEN
            nodes.MOC_EndNodes(nst, Cr * (TIME - TIME_sample)/dt_MOC, B,  R, dia, A, Bm1, Cm1, Bpn, Cpn)
         ELSEIF(actual_method == AlgebraicWaterHammer) THEN 
            nodes.AWH_EndNodes(nst,  Cr * (TIME - TIME_sample)/dt_MOC - ((1./Cr)-1)*nr, B, R, Bm1, Cm1, Bpn, Cpn) 
         END IF
      END SEQUENTIAL
      -- Interpolated characteristics equation for pipe ports
      f_in.H = Cm1 + Bm1 * f_in.Q
      
      f_out.H = Cpn - Bpn * f_out.Q
END COMPONENT

CLASS NODES
   DECLS
      REAL H   UNITS "m"      "Piezometric head"
      REAL Q   UNITS "m3/s"   "Volume flow"
END CLASS

COMPONENT PipeImpl IS_A AbstractPipe (INTEGER nr = 2)
   DECLS
      REAL dx
      REAL HN[nr+1],QN[nr+1] 
      HIDDEN REAL c1, c2, c3, c4, d1, d2, d3, d4
      REAL BM[10*nr+6]
      REAL R[6*nr]
      REAL A            UNITS "m2"  "Flow Area"
      REAL dt_local
      CONST INTEGER i = 1
      HIDDEN INTEGER ia = 0
      HIDDEN INTEGER irow = 0
      HIDDEN INTEGER ma 
      HIDDEN INTEGER mc
      HIDDEN CONST INTEGER mup = 2
      HIDDEN CONST INTEGER mld = 2
      HIDDEN INTEGER ier = 0
   OBJECTS
      NODES node[nr+1]
   INIT
      A = PI * dia**2 /4.
      f_in.H = Ho
      f_out.H = Ho - f * (l/dia) * Qo * abs(Qo) /(2. * A**2 * G)
      f_in.Q = Qo
      f_out.Q = Qo
      FOR(i IN 1, nr+1)
         node[i].H  = Ho - (i-1) * f * (l/dia) * Qo * abs(Qo) /(2. * A**2 * G)/nr
         node[i].Q = Qo
      END FOR
   DISCRETE
      WHEN (sample) THEN 
         IF(istep > 0) THEN
         FOR(i IN 1, nr+1)
            node[i].H = HN[i]
            node[i].Q = QN[i]
         END FOR
         END IF
      END WHEN
   CONTINUOUS
      A = PI * dia**2 /4.
      dx = l/nr
      dt_local = TIME - TIME_sample + 1e-9
      SEQUENTIAL
         ia = 0
         irow = 0
         
         FOR(i IN 1, nr)
            d1    = 2.*(1-theta1) - theta2*((node[i].Q + node[i+1].Q) / A) * (dt_local/dx) + f * dt_local * EPS_FRIC * abs(node[i].Q + node[i+1].Q)/(4.*dia*A)
            d2    = 2.*  theta1   + theta2*((node[i].Q + node[i+1].Q) / A) * (dt_local/dx) + f * dt_local * EPS_FRIC * abs(node[i].Q + node[i+1].Q)/(4.*dia*A)
            d3    = 2. * theta2 * G * A * dt_local /dx
            d4    = 2.*(1-theta2)*G*A*(dt_local/dx)*(node[i+1].H - node[i].H) - 2.* (theta1*node[i+1].Q + (1.-theta1)*node[i].Q) \
                     + ((node[i].Q + node[i+1].Q)/A)*(1-theta2)* (node[i+1].Q - node[i].Q) * (dt_local/dx) \
                     + f *(1-EPS_FRIC) * dt_local * (node[i].Q + node[i+1].Q)*abs(node[i].Q + node[i+1].Q)/(4.*dia*A)
              
            c1   = (a**2)* theta2 *(dt_local/dx)
            c2   = (1-theta1) * G * A  - 0.5 * G * theta2* (node[i].Q + node[i+1].Q)*(dt_local/dx)
            c3   =     theta1 * G * A +  0.5 * G * theta2* (node[i].Q + node[i+1].Q)*(dt_local/dx)
            c4  = -G * A *(theta1*node[i+1].H + (1-theta1)*node[i].H) + \
                  0.5 * G * (node[i].Q + node[i+1].Q) * (1-theta2) * (node[i+1].H - node[i].H)* (dt_local/dx)\
                  +(a**2)*(1-theta2)*(node[i+1].Q - node[i].Q) * (dt_local/dx)
            --
            mc = min(2*nr, 1 + mup + mld)
            ma = 2 * nr * mc  - (mc - mld -1) * (mc - mld)/2
               IF(nr == 1) THEN
                  irow = 1
                  ia = 0
                  BM[ia+1] = d1
                  BM[ia+2] = d2
                  BM[ia+3] = 0  --d3
                  R[1] = -d4
                  R[1+2*nr] =  d3
                  R[1+4*nr] = -d3
                  irow = 2
                  ia = ia + 3
                  BM[ia+1]= -c1
                  BM[ia+2]=  c1
                  BM[ia+3]=  0 --c3
                  BM[ia+4]=  0
                  R[2] = -c4
                  R[2+2*nr] = -c2
                  R[2+4*nr] = -c3
               ELSE
                  IF(i == 1) THEN
                     irow = 1
                     ia = 0
                     BM[ia+1] = d1
                     BM[ia+2] = d2
                     BM[ia+3] = d3
                     R[1] = -d4
                     R[1+2*nr] =  d3
                     R[1+4*nr] =  0
                     irow = 2
                     ia = ia + 3
                     BM[ia+1]= -c1
                     BM[ia+2]=  c1
                     BM[ia+3]=  c3
                     BM[ia+4]=  0
                     R[2] = -c4
                     R[2+2*nr] = -c2
                     R[2+4*nr] =  0
                     ia = ia + 4
                  ELSEIF(i < nr) THEN
                     irow = irow + 1
                     BM[ia+1] = 0
                     BM[ia+2] =  d1
                     BM[ia+3] = -d3
                     BM[ia+4] =  d2
                     BM[ia+5] =  d3
                     R[irow] = -d4
                     R[irow+2*nr] = 0
                     R[irow+4*nr] = 0
                     ia = ia + 5
                     irow = irow + 1
                     BM[ia+1] =  -c1
                     BM[ia+2] =  c2
                     BM[ia+3] =  c1
                     BM[ia+4] =  c3
                     BM[ia+5] =  0
                     R[irow] = -c4
                     R[irow+2*nr] = 0
                     R[irow+4*nr] = 0
                     ia = ia + 5       
                  ELSE
                     irow = irow + 1
                     BM[ia+1] = 0
                     BM[ia+2] =  d1
                     BM[ia+3] = -d3
                     BM[ia+4] =  d2
                     R[irow] = -d4
                     R[irow+2*nr] =  0
                     R[irow+4*nr] = -d3
                     ia = ia +4
                     irow = irow + 1
                     BM[ia+1] =  -c1
                     BM[ia+2] =   c2
                     BM[ia+3] =   c1
                     R[irow] =  -c4
                     R[irow+2*nr] =  0
                     R[irow+4*nr] =  -c3
                     ia = ia+3
                  END IF
               END IF
            END FOR
            IF (nr >=1) THEN
               GELB(R, BM, 2*nr, 3, mup, mld, 1e-6, ier)
               IF(ier!=0) THEN
                  WRITE ("*** TROUBLE IN SOLUTION OF LINEAR EQUATION SYSTEM IN IMPLICIT PIPE  IER= %d\n", ier)
               END IF
            END IF
      END SEQUENTIAL
      
      EXPAND_BLOCK(nr == 0) 
            f_in.H = f_out.H
            f_in.Q = f_out.Q
      END EXPAND_BLOCK
      
      EXPAND_BLOCK(nr == 1) 
         d1*f_in.Q - d3*f_in.H + d2*f_out.Q + d3*f_out.H = -d4
        -c1*f_in.Q + c2*f_in.H + c1*f_out.Q + c3*f_out.H = -c4
      END EXPAND_BLOCK
      EXPAND_BLOCK(nr >1)
         f_in.Q = R[1] + R[1+2*nr]*f_in.H + R[1+4*nr]*f_out.H
         f_out.Q = R[2*nr] + R[2*nr+2*nr] * f_in.H + R[2*nr+4*nr] * f_out.H
         EXPAND (j IN 2, nr)
            QN[j] = R[2*(j-1)] + R[2*(j-1)+2*nr] * HN[1] + R[2*(j-1)+4*nr] * HN[nr+1]
         EXPAND (j IN 2, nr)
            HN[j] = R[2*j-1] + R[2*j-1+2*nr] * HN[1] + R[2*j-1+4*nr] * HN[nr+1]
      END EXPAND_BLOCK
      f_in.H = HN[1]
      f_in.Q = QN[1]
      f_out.H = HN[nr+1]
      f_out.Q = QN[nr+1]
END COMPONENT

ENUM TankType = {Atmospheric, Pressurized}

COMPONENT Tank "Tank component"
   PORTS
      OUT LIQHAMMER_ST.fluid f_out
   DATA
      ENUM TankType type = Atmospheric          "Type of tank"
      REAL H_tank = 100.      UNITS "m"         "Elevation of the level in the tank" 
      REAL z = 0.             UNITS "m"         "Elevation of the outlet connection"
      REAL P_gas = 0.101325   UNITS "MPa (abs)" "Gas pressure"
   INIT
      f_out.z = z
   CONTINUOUS
      f_out.H = IF (type == Atmospheric) 
                     H_tank
                ELSE  --(type == Pressurized)
                     H_tank + (P_gas - LIQHAMMER_ST.P_atm) * 1e6 / (rho * G) + z 
                  
END COMPONENT

COMPONENT ExitValve
   PORTS
      IN LIQHAMMER_ST.fluid f_in
   DATA
      REAL Cd_A = 0.009    UNITS "m2"  "Effective exit area of the valve"
      REAL z  = 0.         UNITS "m"   "Valve elevation" 
   DECLS
      BOUND REAL pos       UNITS "p.u."   "Valve position"
      REAL Q               UNITS "m3/s"   "Flow leaving the valve"
      REAL V_bubble        UNITS "m3"     "Bubble volume"
      --Hidden variables
      HIDDEN DISCR REAL Q_old             UNITS "m3/s"   "Growth rate of bubble at previous time step"    
      HIDDEN DISCR REAL V_bubble_old = 0  UNITS "m3"     "Bubble volume at previous time step"
      HIDDEN DISCR REAL TIME_old = 0      UNITS "s"      "Time of previous step"
      HIDDEN INTEGER icav = 0                            "Cavitation flag"          
   INIT
      f_in.z = z
      V_bubble = 0
      V_bubble_old = 0
      Q_old = 0
      TIME_old = TIME
      icav = 0
   DISCRETE
      --Event to detect bubble formation  
      WHEN (LIQHAMMER_ST.CavitationFlag == TRUE AND f_in.H - z < - H_atm + H_vap) THEN
         icav = 1
         Q_old = Q - f_in.Q
         TIME_old = TIME
      END WHEN
      
      WHEN (sample) THEN
         IF (icav == 1) THEN
            V_bubble_old = V_bubble
            Q_old = Q - f_in.Q
            TIME_old = TIME
         END IF
      END WHEN
      --Event to detect bubble colapse    
      WHEN (V_bubble < 0) THEN
         V_bubble = 0
         V_bubble_old = 0
         icav = 0
      END WHEN
   CONTINUOUS
      Q = Cd_A * pos * sqrt(2.* G ) * fsqrt((f_in.H - z), 0.001)
      
   IMPL(f_in.Q) 0 = (1-icav) * (f_in.Q - Q) \
                  + icav * (f_in.H -z + H_atm - H_vap)
                  
       V_bubble  = ZONE (icav == 0) 0
                    OTHERS V_bubble_old + 0.5 * ((Q  - f_in.Q) + Q_old) * (TIME - TIME_old)
END COMPONENT



COMPONENT Flow_in "Assigned inlet flow condition"
   PORTS
      OUT fluid f_out   "Outlet fluid port"
   DATA
      REAL  z  = 0.   UNITS "m"  "Valve elevation"
   DECLS
      BOUND REAL Q = 0     UNITS "m3/s"   "Assigned inlet flow"
      REAL V_bubble = 0    UNITS "m3"     "Bubble volume"
      --Hidden Variables
      HIDDEN DISCR REAL Q_old             UNITS "m3/s"   "Growth rate of bubble at previous time step"
      HIDDEN DISCR REAL V_bubble_old      UNITS "m3"     "Bubble volume at previous time step"  
      HIDDEN DISCR REAL TIME_old          UNITS "s"      "Time of previous step"
      HIDDEN INTEGER    icav                             "Cavitation flag"
   INIT
      f_out.z = z
      V_bubble = 0
      V_bubble_old = 0
      Q_old = 0
      TIME_old = TIME
      icav = 0
   DISCRETE
      --Event to detect bubble formation 
      WHEN (LIQHAMMER_ST.CavitationFlag == TRUE AND f_out.H - z < - H_atm + H_vap) THEN
         icav = 1
         Q_old = f_out.Q - Q
         TIME_old = TIME
      END WHEN
      
      WHEN (sample) THEN
         IF (icav == 1) THEN
            V_bubble_old = V_bubble
            Q_old = f_out.Q - Q
            TIME_old = TIME
         END IF
      END WHEN
      
      --Event to detect bubble colapse
      WHEN (V_bubble < 0) THEN
         V_bubble = 0
         V_bubble_old = 0
         Q_old = 0
         icav = 0
      END WHEN
            
   CONTINUOUS
        0 = (1-icav) * (Q - f_out.Q) + icav * (f_out.H -z + H_atm - H_vap)
        V_bubble  = ZONE (icav == 0) 0
                    OTHERS V_bubble_old + 0.5 * ((f_out.Q - Q) + Q_old) * (TIME - TIME_old)
END COMPONENT

COMPONENT Flow_out   "Assigned outlet flow condition"
   PORTS
      IN fluid f_in  "Inlet fluid port"
   DATA
      REAL  z  = 0.  UNITS "m" "Valve elevation"
   DECLS
      BOUND REAL Q = 0     UNITS "m3/s"   "Assigned outlet flow"
      REAL V_bubble = 0    UNITS "m3"     "Bubble volume"
      --Hidden Variables
      HIDDEN DISCR REAL Q_old             UNITS "m3/s"   "Growth rate of bubble at previous time step"    
      HIDDEN DISCR REAL V_bubble_old      UNITS "m3"     "Bubble volume at previous time step"  
      HIDDEN DISCR REAL TIME_old          UNITS "s"      "Time of previous step"
      HIDDEN INTEGER    icav                             "Cavitation flag"
   INIT
      f_in.z = z
      V_bubble = 0.
      V_bubble_old = 0.
      Q_old = 0.
      TIME_old = TIME
      icav = 0.
   DISCRETE
      --Event to detect bubble formation  
      WHEN (LIQHAMMER_ST.CavitationFlag == TRUE AND f_in.H - z < - H_atm + H_vap) THEN
         icav = 1
         Q_old = Q - f_in.Q
         TIME_old = TIME
      END WHEN
      
      WHEN (sample) THEN
         IF (icav == 1) THEN
            V_bubble_old = V_bubble
            Q_old = Q - f_in.Q
            TIME_old = TIME
         END IF
      END WHEN
      
      --Event to detect bubble colapse    
      WHEN (V_bubble < 0) THEN
         V_bubble = 0
         V_bubble_old = 0
         icav = 0
      END WHEN
            
    CONTINUOUS
      --Be careful with ibubble, sometimes the initial value is not consistant
        0 = (1-icav) * (Q - f_in.Q) + icav * (f_in.H -z + H_atm - H_vap)
        V_bubble  = ZONE (icav == 0) 0
                    OTHERS V_bubble_old +  0.5 * ((Q - f_in.Q) + Q_old) * (TIME - TIME_old)
END COMPONENT

ABSTRACT COMPONENT Junction_In_Out     
   PORTS
        IN LIQHAMMER_ST.fluid f_in        "Inlet fluid port"
        OUT LIQHAMMER_ST.fluid f_out      "Outlet fluid port"
   DATA
      REAL z_in  = 0.                     UNITS "m"      "Inlet elevation"
		BOOLEAN outlet_elev_same_as_inlet = TRUE  		"Outlet elevation, same as inlet elevation"
		REAL z_out = 0.							UNITS "m"		"Outlet elevation"
   DECLS
      REAL Q                      		   UNITS "m3/s"   "Flow through junction"
      --
      REAL V_bubble_in                    UNITS "m3"     "Volume of imlet bubble"
      HIDDEN DISCR REAL Q_old_in          UNITS "m3/s"   "Growth rate of bubble at inlet" 
      HIDDEN DISCR REAL V_bubble_old_in   UNITS "m3"     "Volume of inlet bubble at previous time step"
      HIDDEN DISCR REAL TIME_old_in       UNITS "s"      "Time of previous step at inlet"
      HIDDEN INTEGER icav_in                       "Cavitation flag at inlet"
      --
      REAL V_bubble_out                   UNITS "m3"     "Volume of imlet bubble"
      HIDDEN DISCR REAL Q_old_out         UNITS "m3/s"   "Growth rate of bubble at inlet" 
      HIDDEN DISCR REAL V_bubble_old_out  UNITS "m3"     "Volume of inlet bubble at previous time step"
      HIDDEN DISCR REAL TIME_old_out      UNITS "s"      "Time of previous step at inlet"
      HIDDEN INTEGER    icav_out                         "Cavitation flag at inlet" 
   INIT
		IF(outlet_elev_same_as_inlet) THEN
      		f_in.z = z_in
			f_out.z = z_in
		ELSE
			f_in.z = z_in
			f_out.z = z_out
		END IF
      --Initialization of variables for bubble in inlet side
      V_bubble_in = 0.
      V_bubble_old_in = 0.
      Q_old_in = 0.
      icav_in = 0.
      --Initialization of variables for bubble in outlet side
      V_bubble_out = 0.
      V_bubble_old_out = 0.
      Q_old_out = 0.
      icav_out = 0.
      --
      TIME_old_in = TIME
   DISCRETE 
      --Events to detect bubble formation    
      WHEN (LIQHAMMER_ST.CavitationFlag == TRUE AND f_in.H - f_in.z < - H_atm + H_vap) THEN
         icav_in = 1
         Q_old_in = Q - f_in.Q
         TIME_old_in = TIME
      END WHEN
      
      WHEN (LIQHAMMER_ST.CavitationFlag == TRUE AND f_out.H - f_out.z < - H_atm + H_vap) THEN
         icav_out = 1
         Q_old_in = f_out.Q - Q
         TIME_old_out = TIME
      END WHEN
      
      WHEN (sample) THEN
         IF (icav_in == 1) THEN
            V_bubble_old_in = V_bubble_in
            Q_old_in = Q - f_in.Q
         END IF
         IF (icav_out == 1) THEN
            V_bubble_old_out = V_bubble_out
            Q_old_out = f_out.Q - Q
         END IF
         TIME_old_in = TIME
         TIME_old_out = TIME
      END WHEN    
      --Events to detect bubble colapse      
      WHEN (V_bubble_in < 0) THEN
         V_bubble_in = 0
         V_bubble_old_in = 0
         icav_in = 0
      END WHEN

      WHEN (V_bubble_out < 0) THEN
         V_bubble_out = 0
         V_bubble_old_out = 0
         icav_out = 0
      END WHEN
   CONTINUOUS

IMPL(f_in.Q)   0 = (1-icav_in) * (Q - f_in.Q) + icav_in * (f_in.H - f_in.z + H_atm - H_vap)
IMPL(f_out.Q)  0 = (1-icav_out) * (f_out.Q - Q) + icav_out * (f_out.H - f_out.z + H_atm - H_vap)

        V_bubble_in  = ZONE (icav_in == 0) 0
                    OTHERS V_bubble_old_in +  0.5 * ((Q - f_in.Q) + Q_old_in) * (TIME - TIME_old_in)
        V_bubble_out  = ZONE (icav_out == 0) 0
                    OTHERS V_bubble_old_out +  0.5 * ((f_out.Q -Q) + Q_old_out) * (TIME - TIME_old_out)
                    
END COMPONENT

COMPONENT Valve IS_A Junction_In_Out "Two ways valve"

    DATA
      REAL A =  0.009      UNITS "m2"  "Valve reference area"
      REAL k_loss = 1      UNITS ""    "Loss coefficient"
    DECLS
      BOUND REAL pos     UNITS "p.u."  "Valve position"
   CONTINUOUS
         Q = A * pos * sqrt(2. * G / k_loss ) * fsqrt((f_in.H - f_out.H), 0.001)               
END COMPONENT

COMPONENT NonReturnValve IS_A Junction_In_Out "Non return valve with ideal behavior"

   DATA
      REAL A =  0.009      UNITS "m2"  "Valve reference area"
      REAL k_loss = 1      UNITS ""    "Loss coefficient for forward flow"
   CONTINUOUS  
      Q =  IF (f_in.H >  f_out.H) A * sqrt(2. * G / k_loss ) * fsqrt((f_in.H - f_out.H), 0.001)
                        ELSE 0 
END COMPONENT

ENUM PumpTransient = {Trip, Start, ConstantSpeed}

COMPONENT Pump IS_A Junction_In_Out "Four quadrants pump"
   DATA
      REAL tdh_r = 10.        UNITS "m"      "Total dynamic head at rated conditions"
      REAL Q_r = 20.          UNITS "m3/s"   "Volume flow at at rated conditions"
      REAL n_r = 1490.        UNITS "rpm"    "Pump speed at rated conditions"
      REAL eff_r = 0.84       UNITS "p.u."   "Efficiency at rated conditions"            
        --
      REAL Ns  = 25     		 UNITS "SI"  "Specific speed of the pump in metric units, n_r in rpm, Q_r in m3/s and tdh_r in m"		
      INTEGER N_stages = 1          "Number of stages, only used to calculate the informative specific speed"
      INTEGER N_suctions = 1			"Number of suctions, = 1 if single suction, = 2 if double suction, only used to calculate the informative specific speed"
      BOOLEAN User_Curves = FALSE   "TRUE to use the curves defined by the user"     
       
        --User Specified 1-D Tables
      TABLE_1D wh_vs_theta  = \
         {{ 0.0000000, 0.0713998, 0.1427997, 0.2141995, 0.2855993, 0.3569992, 0.4283990, 
            0.4997988, 0.5711987, 0.6425985, 0.7139983, 0.7853982, 0.8567980, 0.9281978, 
            0.9995977, 1.0709975, 1.1423973, 1.2137972, 1.2851970, 1.3565968, 1.4279967, 
            1.4993965, 1.5707963, 1.6421962, 1.7135960, 1.7849958, 1.8563957, 1.9277955, 
            1.9991953, 2.0705952, 2.1419950, 2.2133948, 2.2847947, 2.3561945, 2.4275943, 
            2.4989942, 2.5703940, 2.6417938, 2.7131937, 2.7845935, 2.8559933, 2.9273932, 
            2.9987930, 3.0701928, 3.1415927, 3.2129925, 3.2843923, 3.3557922, 3.4271920, 
            3.4985918, 3.5699917, 3.6413915, 3.7127913, 3.7841912, 3.8555910, 3.9269908, 
            3.9983907, 4.0697905, 4.1411903, 4.2125901, 4.2839900, 4.3553898, 4.4267896, 
            4.4981895, 4.5695893, 4.6409891, 4.7123890, 4.7837888, 4.8551886, 4.9265885, 
            4.9979883, 5.0693881, 5.1407880, 5.2121878, 5.2835876, 5.3549875, 5.4263873, 
            5.4977871, 5.5691870, 5.6405868, 5.7119866, 5.7833865, 5.8547863, 5.9261861, 
            5.9975860, 6.0689858, 6.1403856, 6.2117855, 6.2831853},
            { 0.634,     0.643,     0.646,     0.640,     0.629,     0.613,     0.595,           
              0.575,     0.552,     0.533,     0.516,     0.505,     0.504,     0.510,
              0.512,     0.522,     0.539,     0.559,     0.580,     0.601,     0.630,
              0.662,     0.692,     0.722,     0.753,     0.782,     0.808,     0.832,       
              0.857,     0.879,     0.904,     0.930,     0.959,     0.996,     1.027,
              1.060,     1.090,     1.124,     1.165,     1.204,     1.238,     1.258,
              1.271,     1.282,     1.288,     1.281,     1.260,     1.225,     1.172,     
              1.107,     1.031,     0.942,     0.842,     0.733,     0.617,     0.500,     
              0.368,     0.240,     0.125,     0.011,    -0.102,    -0.168,    -0.255,
             -0.342,    -0.423,    -0.494,    -0.556,    -0.620,    -0.655,    -0.670,
             -0.670,    -0.660,    -0.655,    -0.640,    -0.600,    -0.570,    -0.520,
             -0.470,    -0.430,    -0.360,    -0.275,    -0.160,    -0.040,     0.130,     
              0.295,     0.430,     0.550,     0.620,     0.634 }}    "wh vs theta "
      TABLE_1D wbeta_vs_theta =   \
         {{ 0.0000000, 0.0713998, 0.1427997, 0.2141995, 0.2855993, 0.3569992, 0.4283990, 
            0.4997988, 0.5711987, 0.6425985, 0.7139983, 0.7853982, 0.8567980, 0.9281978, 
            0.9995977, 1.0709975, 1.1423973, 1.2137972, 1.2851970, 1.3565968, 1.4279967, 
            1.4993965, 1.5707963, 1.6421962, 1.7135960, 1.7849958, 1.8563957, 1.9277955, 
            1.9991953, 2.0705952, 2.1419950, 2.2133948, 2.2847947, 2.3561945, 2.4275943, 
            2.4989942, 2.5703940, 2.6417938, 2.7131937, 2.7845935, 2.8559933, 2.9273932, 
            2.9987930, 3.0701928, 3.1415927, 3.2129925, 3.2843923, 3.3557922, 3.4271920, 
            3.4985918, 3.5699917, 3.6413915, 3.7127913, 3.7841912, 3.8555910, 3.9269908, 
            3.9983907, 4.0697905, 4.1411903, 4.2125901, 4.2839900, 4.3553898, 4.4267896, 
            4.4981895, 4.5695893, 4.6409891, 4.7123890, 4.7837888, 4.8551886, 4.9265885, 
            4.9979883, 5.0693881, 5.1407880, 5.2121878, 5.2835876, 5.3549875, 5.4263873, 
            5.4977871, 5.5691870, 5.6405868, 5.7119866, 5.7833865, 5.8547863, 5.9261861, 
            5.9975860, 6.0689858, 6.1403856, 6.2117855, 6.2831853},
          {  -0.684,    -0.547,    -0.414,    -0.292,    -0.187,    -0.105,    -0.053,
             -0.012,     0.042,     0.097,     0.156,     0.227,     0.300,     0.371,     
              0.444,     0.522,     0.596,     0.672,     0.738,     0.763,     0.797,   
              0.837,     0.865,     0.883,     0.886,     0.877,     0.859,     0.838,     
              0.804,     0.758,     0.703,     0.645,     0.583,     0.520,     0.454,     
              0.408,     0.370,     0.343,     0.331,     0.329,     0.338,     0.354, 
              0.372,     0.405,     0.450,     0.486,     0.520,     0.552,     0.579,     
              0.603,     0.616,     0.617,     0.606,     0.582,     0.546,     0.500, 
              0.432,     0.360,     0.288,     0.214,     0.123,     0.037,    -0.053,
             -0.161,    -0.248,    -0.314,    -0.372,    -0.580,    -0.740,    -0.880,
             -1.000,    -1.120,    -1.250,    -1.370,    -1.490,    -1.590,    -1.660,
             -1.690,    -1.770,    -1.650,    -1.590,    -1.520,    -1.420,    -1.320,
             -1.230,    -1.100,    -0.980,    -0.820,    -0.684}} "wbeta vs theta"

      ENUM PumpTransient pump_transient = Trip  "Type of pump transient"
      REAL I = 1                 UNITS "kg·m2"  "Moment of Inertia"
      REAL TIME_trip = 1         UNITS "s"      "Time of pump trip or shutdown"
      REAL n_o_trip = 1490.      UNITS "rpm"    "Initial pump speed previous to the trip time"
      REAL TIME_start = 1        UNITS "s"      "Time of pump start"
      REAL DTIME_start = 5       UNITS "s"      "Duration of the pump start ramp"
      REAL n_start_end = 1490.   UNITS "rpm"    "Speed at the end of the start ramp"
      REAL n_constant = 1490.    UNITS "rpm"    "Constant speed"
   DECLS
      REAL alpha              "Adimensional speed"
      REAL beta               "Adimensional torque"
      REAL h                  "Adimensional head"
      REAL n                  UNITS "rpm" "Pump speed"
      REAL tdh                UNITS "m"   "Total dynamic head"
      REAL torque_r           UNITS "N·m" "Torque at rated conditions"
      REAL torque             UNITS "N·m" "Brake torque" 
      REAL theta              "Variable defined as PI + atan2(v, alpha)"
      REAL v                  "Adimensional volume flow"
      REAL wh                 "Dimensionless turbomachine characteristics defined as h / (v**2 + alpha**2)"
      REAL wbeta              "Dimensionless turbomachine characteristics defined as beta / (v**2 + alpha**2)"
      REAL torque_motor       UNITS "N·m"    "Motor torque"
      DISCR REAL wh_array[89]
      DISCR REAL wbeta_array[89]
      --Auxiliary variables to manage the pump curves
      INTEGER m
      HIDDEN REAL a0, a1
      HIDDEN REAL b0, b1
      DISCR REAL Ns_info            UNITS "SI"  "Informative specific speed of the pump in metric units, n_r in rpm, Q_r in m3/s and tdh_r in m"
      DISCR REAL n_old              UNITS "rpm" "Rotational speed at previous time step"
      DISCR REAL torque_old         UNITS "N.m" "Brake torque at previous time step"
      DISCR REAL torque_motor_old   UNITS "N·m" "Motor torque at previous time step"
      DISCR REAL TIME_old
      CONST REAL dx = 2. * PI / 88
      ---------------------------------------------------------------------------------------
      -- Two D Table of dimensionless head (H/Hr) / ((Q/Qr)**2 + (N/Nr)**2) vs theta
      ---------------------------------------------------------------------------------------
      HIDDEN CONST TABLE_2D wh_vs_theta_Ns = \
       {  { 25, 147, 261 },
          { 0.0000000, 0.0713998, 0.1427997, 0.2141995, 0.2855993, 0.3569992, 0.4283990, 
            0.4997988, 0.5711987, 0.6425985, 0.7139983, 0.7853982, 0.8567980, 0.9281978, 
            0.9995977, 1.0709975, 1.1423973, 1.2137972, 1.2851970, 1.3565968, 1.4279967, 
            1.4993965, 1.5707963, 1.6421962, 1.7135960, 1.7849958, 1.8563957, 1.9277955, 
            1.9991953, 2.0705952, 2.1419950, 2.2133948, 2.2847947, 2.3561945, 2.4275943, 
            2.4989942, 2.5703940, 2.6417938, 2.7131937, 2.7845935, 2.8559933, 2.9273932, 
            2.9987930, 3.0701928, 3.1415927, 3.2129925, 3.2843923, 3.3557922, 3.4271920, 
            3.4985918, 3.5699917, 3.6413915, 3.7127913, 3.7841912, 3.8555910, 3.9269908, 
            3.9983907, 4.0697905, 4.1411903, 4.2125901, 4.2839900, 4.3553898, 4.4267896, 
            4.4981895, 4.5695893, 4.6409891, 4.7123890, 4.7837888, 4.8551886, 4.9265885, 
            4.9979883, 5.0693881, 5.1407880, 5.2121878, 5.2835876, 5.3549875, 5.4263873, 
            5.4977871, 5.5691870, 5.6405868, 5.7119866, 5.7833865, 5.8547863, 5.9261861, 
            5.9975860, 6.0689858, 6.1403856, 6.2117855, 6.2831853  
          },
          { { 0.634,     0.643,     0.646,     0.640,     0.629,     0.613,     0.595,           
              0.575,     0.552,     0.533,     0.516,     0.505,     0.504,     0.510,
              0.512,     0.522,     0.539,     0.559,     0.580,     0.601,     0.630,
              0.662,     0.692,     0.722,     0.753,     0.782,     0.808,     0.832,       
              0.857,     0.879,     0.904,     0.930,     0.959,     0.996,     1.027,
              1.060,     1.090,     1.124,     1.165,     1.204,     1.238,     1.258,
              1.271,     1.282,     1.288,     1.281,     1.260,     1.225,     1.172,     
              1.107,     1.031,     0.942,     0.842,     0.733,     0.617,     0.500,     
              0.368,     0.240,     0.125,     0.011,    -0.102,    -0.168,    -0.255,
             -0.342,    -0.423,    -0.494,    -0.556,    -0.620,    -0.655,    -0.670,
             -0.670,    -0.660,    -0.655,    -0.640,    -0.600,    -0.570,    -0.520,
             -0.470,    -0.430,    -0.360,    -0.275,    -0.160,    -0.040,     0.130,     
              0.295,     0.430,     0.550,     0.620,     0.634     
            },
            {-0.690,    -0.599,    -0.512,    -0.418,    -0.304,    -0.181,    -0.078,     
             -0.011,     0.032,     0.074,     0.130,     0.190,     0.265,     0.363,     
              0.461,     0.553,     0.674,     0.848,     1.075,     1.337,     1.629, 
              1.929,     2.180,     2.334,     2.518,     2.726,     2.863,     2.948,     
              3.026,     3.015,     2.927,     2.873,     2.771,     2.640,     2.497,     
              2.441,     2.378,     2.336,     2.288,     2.209,     2.162,     2.140,     
              2.109,     2.054,     1.970,     1.860,     1.735,     1.571,     1.357,     
              1.157,     1.016,     0.927,     0.846,     0.744,     0.640,     0.500,              
              0.374,     0.191,     0.001,    -0.190,    -0.384,    -0.585,    -0.786,     
             -0.972,    -1.185,    -1.372,    -1.500,    -1.940,    -2.160,    -2.290,
             -2.350,    -2.350,    -2.230,    -2.200,    -2.130,    -2.050,    -1.970,
             -1.895,    -1.810,    -1.730,    -1.600,    -1.420,    -1.130,    -0.950,     
             -0.930,    -0.950,    -1.000,    -0.920,    -0.690 
            },
            {-2.230,    -2.000,    -1.662,    -1.314,    -1.089,    -0.914,    -0.760,
             -0.601,    -0.440,    -0.284,    -0.130,     0.055,     0.222,     0.357,     
              0.493,     0.616,     0.675,     0.680,     0.691,     0.752,     0.825,     
              0.930,     1.080,     1.236,     1.389,     1.548,     1.727,     1.919,
              2.066,     2.252,     2.490,     2.727,     3.002,     3.225,     3.355,              
              3.475,     3.562,     3.604,     3.582,     3.540,     3.477,     3.321,     
              3.148,     2.962,     2.750,     2.542,     2.354,     2.149,     1.909, 
              1.702,     1.506,     1.310,     1.131,     0.947,     0.737,     0.500,   
              0.279,     0.082,    -0.112,    -0.300,    -0.505,    -0.672,    -0.797, 
             -0.872,    -0.920,    -0.949,    -0.960,    -1.080,    -1.300,    -1.500,
             -1.700,    -1.890,    -2.080,    -2.270,    -2.470,    -2.650,    -2.810, 
             -2.950,    -3.040,    -3.100,    -3.150,    -3.170,     -3.170,   -3.130, 
             -3.070,    -2.960,    -2.820,    -2.590,    -2.230 
            }
         }    
      }
      ---------------------------------------------------------------------------------------
      -- Two D Table of dimensionless head (T/Tr) / ((Q/Qr)**2 + (N/Nr)**2) vs theta
      ---------------------------------------------------------------------------------------  
      HIDDEN CONST TABLE_2D wbeta_vs_theta_Ns = \
       {  { 25, 147, 261 },
          {    0.0000000, 0.0713998, 0.1427997, 0.2141995, 0.2855993, 0.3569992, 0.4283990, 
            0.4997988, 0.5711987, 0.6425985, 0.7139983, 0.7853982, 0.8567980, 0.9281978, 
            0.9995977, 1.0709975, 1.1423973, 1.2137972, 1.2851970, 1.3565968, 1.4279967, 
          1.4993965, 1.5707963, 1.6421962, 1.7135960, 1.7849958, 1.8563957, 1.9277955, 
          1.9991953, 2.0705952, 2.1419950, 2.2133948, 2.2847947, 2.3561945, 2.4275943, 
          2.4989942, 2.5703940, 2.6417938, 2.7131937, 2.7845935, 2.8559933, 2.9273932, 
          2.9987930, 3.0701928, 3.1415927, 3.2129925, 3.2843923, 3.3557922, 3.4271920, 
          3.4985918, 3.5699917, 3.6413915, 3.7127913, 3.7841912, 3.8555910, 3.9269908, 
          3.9983907, 4.0697905, 4.1411903, 4.2125901, 4.2839900, 4.3553898, 4.4267896, 
          4.4981895, 4.5695893, 4.6409891, 4.7123890, 4.7837888, 4.8551886, 4.9265885, 
          4.9979883, 5.0693881, 5.1407880, 5.2121878, 5.2835876, 5.3549875, 5.4263873, 
          5.4977871, 5.5691870, 5.6405868, 5.7119866, 5.7833865, 5.8547863, 5.9261861, 
          5.9975860, 6.0689858, 6.1403856, 6.2117855, 6.2831853  
          },
          {{  -0.684,    -0.547,    -0.414,    -0.292,    -0.187,    -0.105,    -0.053,
             -0.012,     0.042,     0.097,     0.156,     0.227,     0.300,     0.371,     
              0.444,     0.522,     0.596,     0.672,     0.738,     0.763,     0.797,   
              0.837,     0.865,     0.883,     0.886,     0.877,     0.859,     0.838,     
              0.804,     0.758,     0.703,     0.645,     0.583,     0.520,     0.454,     
              0.408,     0.370,     0.343,     0.331,     0.329,     0.338,     0.354, 
              0.372,     0.405,     0.450,     0.486,     0.520,     0.552,     0.579,     
              0.603,     0.616,     0.617,     0.606,     0.582,     0.546,     0.500, 
              0.432,     0.360,     0.288,     0.214,     0.123,     0.037,    -0.053,
             -0.161,    -0.248,    -0.314,    -0.372,    -0.580,    -0.740,    -0.880,
             -1.000,    -1.120,    -1.250,    -1.370,    -1.490,    -1.590,    -1.660,
             -1.690,    -1.770,    -1.650,    -1.590,    -1.520,    -1.420,    -1.320,
             -1.230,    -1.100,    -0.980,    -0.820,    -0.684
          },
          {  -1.420,    -1.328,    -1.211,    -1.056,    -0.870,    -0.677,    -0.573,
             -0.518,    -0.380,    -0.232,    -0.160,     0.000,     0.118,     0.308,
              0.442,     0.574,     0.739,     0.929,     1.147,     1.370,     1.599,     
              1.839,     2.080,     2.300,     2.480,     2.630,     2.724,     2.687,     
              2.715,     2.688,     2.555,     2.434,     2.288,     2.110,     1.948,              
              1.825,     1.732,     1.644,     1.576,     1.533,     1.522,     1.519,     
              1.523,     1.523,     1.490,     1.386,     1.223,     1.048,     0.909,     
              0.814,     0.766,     0.734,     0.678,     0.624,     0.570,     0.500,     
              0.407,     0.278,     0.146,     0.023,    -0.175,    -0.379,    -0.585,
             -0.778,    -1.008,    -1.277,    -1.560,    -2.070,    -2.480,    -2.700,
             -2.770,    -2.800,    -2.800,    -2.760,    -2.710,    -2.640,    -2.540,
             -2.440,    -2.340,    -2.240,    -2.120,    -2.000,    -1.940,    -1.900,
             -1.900,    -1.850,    -1.750,    -1.630,    -1.420 
          },
          {
             -2.260,    -2.061,    -1.772,    -1.465,    -1.253,    -1.088,    -0.921,
             -0.789,    -0.632,    -0.457,    -0.300,    -0.075,     0.052,     0.234,     
              0.425,     0.558,     0.630,     0.621,     0.546,     0.525,     0.488,     
              0.512,     0.660,     0.850,     1.014,     1.162,     1.334,     1.512,     
              1.683,     1.886,     2.105,     2.325,     2.580,     2.770,     2.886,              
              2.959,     2.979,     2.962,     2.877,     2.713,     2.556,     2.403,     
              2.237,     2.080,     1.950,     1.826,     1.681,     1.503,     1.301,     
              1.115,     0.960,     0.840,     0.750,     0.677,     0.604,     0.500,     
              0.352,     0.161,    -0.040,    -0.225,    -0.403,    -0.545,    -0.610,
             -0.662,    -0.699,    -0.719,    -0.730,    -0.810,    -1.070,    -1.360,
             -1.640,    -1.880,    -2.080,    -2.270,    -2.470,    -2.650,    -2.810,
             -2.950,    -3.040,    -3.100,    -3.150,    -3.170,    -3.200,    -3.160,
             -3.090,    -2.990,    -2.860,    -2.660,    -2.260  
         }
         }    
      } 
      INTEGER ip, ipa, index
   INIT
		IF(User_Curves == FALSE) THEN
      	WRITE ("Default Homologous Curves of Pump Automatically Generated for an Specific Speed = %g\n", Ns)
      	Ns_info = n_r * (Q_r/N_suctions)**0.5 / (tdh_r/N_stages)**0.75
			WRITE("Specified specific speed Ns = %lg  Calculated specific speed = %lg\n", Ns, Ns_info)
			IF (abs(Ns_info-Ns)/Ns_info > 0.1) THEN
				WRITE("*************************************************************************************************\n")
				WRITE("******WARNING specified specific speed differ by more than 10%% from calculated specific speed ***\n")
				WRITE("*************************************************************************************************\n")
			END IF
		END IF
      FOR (j IN 1, 89)
         theta = 2*PI*(j-1)/88.
         IF (User_Curves == FALSE) THEN
            wh_array[j] =     splineInterp2D(wh_vs_theta_Ns, Ns, theta)
            wbeta_array[j] =  splineInterp2D(wbeta_vs_theta_Ns, Ns, theta)
         ELSE
            wh_array[j] =     splineInterp1D(wh_vs_theta, theta)
            wbeta_array[j] =  splineInterp1D(wbeta_vs_theta, theta)
         END IF
         WRITE("%d \t %g \t %g \t %g \n", j, theta, wh_array[j], wbeta_array[j])
      END FOR
      TIME_old =  TIME
      IF (pump_transient == Trip) THEN
         n = n_o_trip
         n_old = n_o_trip
      ELSEIF(pump_transient == Start) THEN
         n = 0
         n_old = 0
      ELSEIF(pump_transient == ConstantSpeed) THEN
         n = n_constant
         n_old = n_constant
      END IF
         alpha = n/n_r
         v = Q/Q_r
         IF(Q != 0 OR n != 0.) THEN
            theta = atan2(v,alpha)+ PI
         ELSE
            theta = 0
         END IF
   DISCRETE
      WHEN (sample) THEN       
         n_old = n AFTER 0
         torque_motor_old = torque_motor
         torque_old = torque 
         TIME_old = TIME  
      END WHEN
   CONTINUOUS
      --following statement cannot go in the INIT block, because the density is calculated
      --in the INIT block of the Solver_Control que tiene baja prioridad.
      
      torque_r =  rho * G * tdh_r * Q_r / eff_r / (n_r * PI/30.)
      Q = v  * Q_r
IMPL(alpha)    n = alpha * n_r

      theta = PI + atan2(v , alpha )
      SEQUENTIAL 
         m = theta/dx + 1
         m = min(m, 88)
         a1 =(wh_array[m+1]-wh_array[m])/dx
         a0 = wh_array[m+1]-a1*m*dx      
         b1 =(wbeta_array[m+1]-wbeta_array[m])/dx
         b0 = wbeta_array[m+1]-b1*m*dx    
         wh  = a1*theta + a0  
         wbeta = b1*theta + b0
         h = wh * (alpha**2 + v**2)
         beta = wbeta * (alpha**2 + v**2)
      END SEQUENTIAL
      ---------------------------------------------------------------
      -- Calculation of total dynamic head (tdh) and instantaneous torque (torque)
      ---------------------------------------------------------------
      tdh = tdh_r * h 
      torque = torque_r * beta 
IMPL(v)     LIQHAMMER_ST.G * f_out.H = LIQHAMMER_ST.G * f_in.H + tdh * LIQHAMMER_ST.Go
      torque_motor = 0
      SEQUENTIAL
         IF (pump_transient == Trip) THEN
            IF(TIME <  TIME_trip) THEN
               n = n_o_trip
            ELSE  
               n = n_old +  0.5*(((torque_motor_old -  torque_old)+ (torque_motor  - torque))/I) * (30./PI) * (TIME - TIME_old) 
            END IF
         ELSEIF (pump_transient == Start) THEN
            IF (TIME < TIME_start) THEN
               n = 0
            ELSEIF (TIME < TIME_start + DTIME_start) THEN
               n = (TIME - TIME_start)/DTIME_start * n_start_end
            ELSE
               n = n_start_end
            END IF
         ELSEIF (pump_transient == ConstantSpeed) THEN
            n = n_constant
         END IF
      END SEQUENTIAL

END COMPONENT

COMPONENT Collector
   PORTS
      IN LIQHAMMER_ST.fluid f_in          "Fluid inlet port"
   DATA
      REAL z = 0  UNITS "m"   "Collector elevation"
   DECLS
      REAL V_bubble = 0
      HIDDEN DISCR REAL Q_in_old
      HIDDEN DISCR REAL V_bubble_old = 0
      HIDDEN DISCR REAL TIME_old = 0
      HIDDEN INTEGER icav = 0
   INIT
      f_in.z = z
      V_bubble = 0
      icav = 0
   DISCRETE
      --Event to detect bubble formation  
      WHEN (LIQHAMMER_ST.CavitationFlag == TRUE AND f_in.H - z < - H_atm + H_vap) THEN
         icav = 1
         Q_in_old = f_in.Q
         TIME_old = TIME
      END WHEN
      
      WHEN (sample) THEN
         IF (icav == 1) THEN
            V_bubble_old = V_bubble
            Q_in_old = f_in.Q
            TIME_old = TIME
         END IF
      END WHEN

      --Event to detect bubble colapse    
      WHEN (V_bubble < 0) THEN
         V_bubble = 0
         V_bubble_old = 0
         icav = 0
      END WHEN
            
    CONTINUOUS
      --Be careful with ibubble, sometimes the initial value is not consistant
IMPL()  0 = (1-icav) * f_in.Q + icav * (f_in.H -z + H_atm - H_vap)
        V_bubble  = ZONE (icav == 0) 0
                    OTHERS V_bubble_old -  0.5 * (f_in.Q + Q_in_old) * (TIME - TIME_old)
END COMPONENT





