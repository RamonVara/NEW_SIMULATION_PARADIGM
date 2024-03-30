/*-----------------------------------------------------------------------------------------
 LIBRARY: THERMAL_NEW
 FILE: Thermal_NSP  
 CREATION DATE: 27/03/2024
-----------------------------------------------------------------------------------------*/

FUNCTION NO_TYPE tridiag(
--Function implementing the Thomas algorithm to solve tridiagonal linear equation systems
      IN REAL a[],
      IN REAL b[],
      IN REAL c[],
      IN REAL r[],
      OUT REAL u[],
      IN INTEGER n)
   DECLS
      REAL bet
      REAL gam[n]
      INTEGER j
   BODY
      bet=b[1]
      ASSERT (abs(bet) > 0) FATAL  "Tridiag function fails"
      u[1]=r[1]/bet
      FOR (j IN 2,n) --Decomposition and forward substitution.
         gam[j]=c[j-1]/bet
         bet=b[j]-a[j]*gam[j]
         ASSERT (abs(bet) > 0) FATAL  "Tridiag function fails"
         u[j]=(r[j]-a[j]*u[j-1])/bet
      END FOR
      FOR(k IN 1, n-1)
         j = n - k
         u[j]=u[j]-gam[j+1]*u[j+1]
      END FOR
      RETURN 
END FUNCTION

CLASS WallClass (INTEGER MAX_NODES)
   DECLS
      INTEGER icount = 0
      INTEGER nodes
      REAL k
      REAL rho
      REAL cp
      REAL e
      REAL A
      REAL dx
      REAL Time0
      REAL T0[MAX_NODES]
      REAL T[MAX_NODES]
      REAL alpha
      CONST REAL theta = 0.5
      REAL a[MAX_NODES]
      REAL b[MAX_NODES]
      REAL c[MAX_NODES]
      REAL d[MAX_NODES]
   METHODS
      METHOD NO_TYPE Initialize(IN INTEGER nodes_cmp, IN REAL k_cmp, 
                   IN REAL rho_cmp, IN REAL cp_cmp, IN REAL e_cmp, 
                   IN REAL A_cmp, IN REAL To, IN REAL Time)
         BODY
            icount = 0
            nodes = nodes_cmp
            k = k_cmp
            rho = rho_cmp
            cp = cp_cmp
            e = e_cmp
            A = A_cmp
            dx = e/(nodes-1)
            Time0 = Time
            FOR (i IN 1, nodes)
               T0[i] = To
               T[i] = To
            END FOR
            RETURN
      END METHOD
      METHOD BOOLEAN Update(IN REAL Time)
         BODY
            icount = 0
            FOR (i IN 1, nodes)
               T0[i] = T[i]
            END FOR
            Time0 = Time
            RETURN TRUE
      END METHOD
      METHOD NO_TYPE InnerPoints(IN REAL Time, IN REAL Tk_in, IN REAL Tk_out)
         DECLS
            REAL r         
         BODY
            alpha = k /(rho*cp)
            r = theta * (Time-Time0) * alpha / dx**2
            a[1] = 0
            b[1] = 1
            c[1] = 0
            d[1] = Tk_in
            FOR(j IN 2, nodes-1)
               a[j] = -r
               b[j] = (1. + 2.*r)
               c[j] = -r
               d[j] = T0[j] + (1-theta)*((Time-Time0) * alpha  / dx**2)\
                       *(T0[j-1] - 2*T0[j] + T0[j+1])
            END FOR
            a[nodes] = 0
            b[nodes] = 1
            c[nodes] = 0
            d[nodes] = Tk_out
            tridiag(a,b,c,d,T,nodes)
         RETURN
      END METHOD
      METHOD REAL q_in(IN REAL Time, IN REAL Tk_in, IN REAL Tk_out)
         DECLS
            REAL q_in
         BODY
            icount = icount + 1
            IF(icount - icount/2*2 == 1) THEN
               InnerPoints(Time, Tk_in, Tk_out)
            END IF
               q_in = k*A*(T[1]- T[2])/dx
               RETURN q_in
      END METHOD
      METHOD REAL q_out(IN REAL Time, IN REAL Tk_in, IN REAL Tk_out)
         DECLS
            REAL q_out
         BODY
            icount = icount + 1
            IF(icount - icount/2*2 == 1) THEN
               InnerPoints(Time, Tk_in, Tk_out)
            END IF
            q_out = k*A*(T[nodes-1] - T[nodes])/dx
            RETURN q_out
      END METHOD
END CLASS

COMPONENT Wall_NSP "Thermal 1D transient wall using the new approach"
   PORTS
      IN thermal tp_in
      OUT thermal tp_out
   DATA
      INTEGER nodes = 20  "Number of nodes"
      REAL A = 1          UNITS "m2" "Wall Area"
      REAL e = 0.1        "Wall Thickness"
      REAL rho = 1000.    "Densiy of wall material"
      REAL cp = 500       "Specific heat of wall material"
      REAL k = 50.        "Thermal conductivity of wall material"      
      REAL To = 290.      "Initial Wall Temperature"
   OBJECTS
      WallClass(MAX_NODES=1000) wall
   INIT
      wall.Initialize(nodes, k, rho, cp, e, A, To, TIME) 
   DISCRETE
      ASSERT(wall.Update(TIME)) WARNING "This message should never appear"
   CONTINUOUS
      tp_in.q = wall.q_in(TIME, tp_in.Tk, tp_out.Tk)
      tp_out.q =  wall.q_out(TIME, tp_in.Tk, tp_out.Tk)
END COMPONENT
