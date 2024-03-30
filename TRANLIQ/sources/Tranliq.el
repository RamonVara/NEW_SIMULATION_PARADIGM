
LIBRARY TRANLIQ
USE MATH
CONST REAL g = 9.806
PORT Hydraulic
   SUM REAL Q     UNITS "m3/s"
   EQUAL REAL H   UNITS "m"      
END PORT

ENUM PipeType = {RR, CR, RC, CC}

COMPONENT Pipe(INTEGER n = 100, ENUM PipeType type =RR)
   PORTS
      IN Hydraulic h_in
      OUT Hydraulic h_out
   DATA
      REAL L = 600      UNITS "m"    "Pipe length"
      REAL D = 0.5      UNITS "m"    "Pipe inside diametert"
      REAL a = 1200.    UNITS "m/s"  "Wave speed"
      REAL f = 0.018    UNITS "-"    "Friction factor"
      REAL kdamp = 0.15 UNITS "-"    "Damping factor for artificial viscosity"
      REAL Ho = 150     UNITS "m"    "Initial elevation at inlet section"
      REAL Qo = 0.47744 UNITS "m3/s" "Initial flow"
   DECLS
      DISCR REAL A      UNITS "m"    "Flow area"
      DISCR REAL dx     UNITS "m"    "Nodal length"
      REAL H[n]         UNITS "m"    "Array of piezometric heads"
      REAL Q[n+1]       UNITS "m3/s" "Array of flows"
      REAL W[n]         UNITS "m"    "Array of artificial viscosity head"
      DISCR REAL km     UNITS "-"    "Multiplier of momentum terms"
   INIT
      A = 0.25 * PI * D**2
      dx = L/n
      IF(type == RR) THEN
         km = 1.
      ELSEIF(type == CR OR type == RC) THEN
         km = n/(n-0.5)
      ELSEIF(type == CC) THEN
         km = n/(n-1)
      END IF
      FOR(i IN 1, n)
         IF(type == RR) THEN
            H[i] = Ho - ((i-1+0.5)/n) * f * (L/D) * Qo*abs(Qo)/(2*g*A**2)
         ELSEIF(type == CC) THEN
            H[i] = Ho - ((i-1)/(n-1)) * f * (L/D) * Qo*abs(Qo)/(2*g*A**2)
         ELSEIF(type == RC) THEN
            H[i] = Ho - ((i-1+0.5)/(n-0.5))* f * (L/D) * Qo*abs(Qo)/(2*g*A**2)
         ELSEIF(type == CR) THEN
            H[i] = Ho - ((i-1)/(n-0.5)) * f * (L/D) * Qo*abs(Qo)/(2*g*A**2)
         END IF
      END FOR
      FOR(i IN 1, n+1)
         Q[i] = Qo
      END FOR
      h_in.Q = Qo
      h_out.Q = Qo
   CONTINUOUS
      h_in.Q = Q[1]
      h_out.Q = Q[n+1]     
      EXPAND_BLOCK(type == RR) 
         km*Q[1]' =   - A * g * (H[1] + W[1] - h_in.H)/ (dx/2) \
                       - km*f*Q[1]*abs(Q[1])/(2*D*A)
         km*Q[n+1]' = - A * g * (h_out.H - H[n] - W[n])/ (dx/2) \
                       - km*f*Q[n+1]*abs(Q[n+1])/(2*D*A)
      END EXPAND_BLOCK     
      EXPAND_BLOCK(type == CR) 
         h_in.H = H[1]
         km*Q[n+1]' =  - A * g * (h_out.H - H[n] - W[n])/ (dx/2) \
                        - km*f*Q[n+1]*abs(Q[n+1])/(2*D*A)
      END EXPAND_BLOCK     
      EXPAND_BLOCK(type == RC) 
         h_out.H = H[n]
         km*Q[1]' = - A * g * (H[1] + W[1] - h_in.H)/ (dx/2) \
                     - km*f*Q[1]*abs(Q[1])/(2*D*A)
      END EXPAND_BLOCK  
      EXPAND_BLOCK(type==CC)
         h_in.H = H[1]
         h_out.H = H[n]
      END EXPAND_BLOCK     
      EXPAND (i IN 1, n)
         H[i]' = -a**2/(g*A) * (Q[i+1] - Q[i])/dx
      EXPAND (i IN 2, n)
         km*Q[i]' = - A * g * (H[i]+W[i]-H[i-1]-W[i-1])/dx \
                     - km * f*Q[i]*abs(Q[i])/(2*D*A) 
      EXPAND (i IN 1, n)
         W[i] = - kdamp*a*(Q[i+1]-Q[i])/(g*A)
END COMPONENT

COMPONENT Tank
   PORTS
      OUT Hydraulic h_out
   DATA
      REAL H = 150. UNITS "m" "Tank Level"
   CONTINUOUS
      h_out.H = H
END COMPONENT

COMPONENT Collector
   PORTS
      IN Hydraulic h_in
   CONTINUOUS
      h_in.Q = 0
END COMPONENT

COMPONENT Flow
   PORTS
      IN Hydraulic h_in
   DATA
   DECLS
      BOUND REAL Q UNITS "m3/s" "Specified flow"
   CONTINUOUS
      h_in.Q = Q
END COMPONENT

COMPONENT ValveExit
   PORTS
      IN Hydraulic h_in
   DATA
      REAL Cd = 1    UNITS "-"  "Discharge coefficient"
      REAL A = 0.009 UNITS "m2" "Nozzle area"
   DECLS
      CONST REAL dHlam = 0.001 UNITS "m" "Pressure diff. for laminar flow"
      BOUND REAL pos UNITS "-"  "Valve Position between 0 & 1"
   CONTINUOUS
      h_in.Q = Cd*A*pos*fsqrt(2*g*h_in.H, dHlam)
END COMPONENT
