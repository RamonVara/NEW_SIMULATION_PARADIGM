/*-----------------------------------------------------------------------------------------
 LIBRARY: THERMAL_NEW
 FILE: Thermal_Classic
 CREATION DATE: 27/03/2024
-----------------------------------------------------------------------------------------*/
PORT thermal
	EQUAL REAL Tk = 290.
	SUM REAL q
END PORT
COMPONENT Wall (INTEGER nodes = 20 ) "Thermal 1D transient wall using the classical approach"
   PORTS
      IN thermal tp_in
      OUT thermal tp_out
   DATA
      REAL A = 1          "Wall Area"
      REAL e = 0.1        "Wall Thickness"
      REAL rho = 1000.    "Densiy of wall material"
      REAL cp = 500       "Specific heat of wall material"
      REAL k = 50.        "Thermal conductivity of wall material"      
      REAL To = 290.      "Initial Wall Temperature"
	DECLS
		DISCR REAL dx
		REAL T[nodes]
   INIT
		dx = e/(nodes-1)
      FOR(j IN 1, nodes)
			T[j] = To
		END FOR
   CONTINUOUS
		tp_in.Tk = T[1]
		tp_out.Tk = T[nodes]
      tp_in.q = k*A*(T[1]-T[2])/dx
      tp_out.q =  k*A*(T[nodes-1]-T[nodes])/dx
		EXPAND(j IN 2, nodes-1)
			T[j]'= (k/(rho*cp))*(T[j-1]-2*T[j]+T[j+1])/dx**2
END COMPONENT
