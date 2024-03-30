/*-----------------------------------------------------------------------------------------
 LIBRARY: HYDRAULIC_ST_CLASSIC
 FILE: hydro_classic.el
 CREATION DATE: 03/11/2023
-----------------------------------------------------------------------------------------*/
PORT Hydraulic
   SUM REAL Q
	EQUAL REAL H
END PORT

COMPONENT Pipe
   PORTS
      IN Hydraulic h_in
      OUT Hydraulic h_out
   DATA
      REAL L = 600      UNITS "m"    "Pipe length"
      REAL D = 0.5      UNITS "m"    "Pipe inside diametert"
      REAL CHW = 120.   UNITS "-"    "Hazen Williams Coefficient"
      REAL Qo = 0.47744 UNITS "m3/s" "Initial flow"
   DECLS
      REAL H1     UNITS "m"    "Piezometric head at inlet"
      REAL H2     UNITS "m"    "Piezometric head at outlet"
      REAL Q      UNITS "m3/s" "Volume flow"
      REAL Qmin
		 CONST REAL Vmin = 1.e-4
   INIT
		 Q = Qo
   CONTINUOUS
	  	 Qmin = 0.785398163*D**2 * Vmin
		 0 = h_in.H - h_out.H - 10.67 * (L /(CHW**1.852 * D**4.8704))* Q * (abs(Q)+Qmin)**0.852
		 H1 = h_in.H
		 H2 = h_out.H
		 Q = h_in.Q
		 Q = h_out.Q
END COMPONENT

COMPONENT Tank
   PORTS
      OUT Hydraulic h_out
   DATA
		REAL A = 5
      REAL Ho = 100.  UNITS "m" "Tank piezometric head"
      REAL Qo = 0.1  UNITS "m3/s" "Initial guess of the tank outflow"
   DECLS
	   REAL H 
      REAL Q 
	INIT
		H = Ho
	CONTINUOUS
		 H' = - Q/A
		 h_out.H = H
		 Q = h_out.Q
END COMPONENT

COMPONENT Outflow
   PORTS
      IN Hydraulic h_in
   DATA
      REAL Q = 0.    UNITS "m3/s" "Outflow"
   DECLS
      REAL H   UNITS "Piezometric head at the connected node"
	CONTINUOUS
		h_in.Q = Q 
		H = h_in.H
END COMPONENT 


COMPONENT Solver
   DATA 
      REAL Ho = 100. UNITS "m" "Initial guess for the node piezometric head"
      INTEGER iter_damp = 5    "Number of initial damped iteration in NR algorithm"    
      REAL FRAC_damp = 0.1     "Damping factor for the first iter_damp iterations"
      INTEGER iter_max = 20    "Maximum number of iterations in NR algorithm"
      INTEGER NR_debug = 1     "Debug level for Newton-Raphson"
      INTEGER STAMP_debug = 0  "Debug level for Stamping process"
   DECLS
   INIT PRIORITY  0
END COMPONENT 
