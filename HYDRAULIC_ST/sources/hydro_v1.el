/*-----------------------------------------------------------------------------------------
 LIBRARY: HYDRO_STEADY
 FILE: hydro_v1.el
 CREATION DATE: 19/10/2023
-----------------------------------------------------------------------------------------*/
USE SPARSE
--Flag to control the part of the INIT block executed
BOOLEAN AssignSolution = FALSE

--Maximum Dimensions
CONST INTEGER MAX_PIPE = 1000
CONST INTEGER MAX_TANK = 100
CONST INTEGER MAX_FLOW = 100

CONST INTEGER MAX_EQN = 2 * MAX_PIPE + 1 + MAX_TANK
--Each pipe makes 5 stamps in the Jacobian and
--each tank makes 2 stamps
CONST INTEGER MAX_STAMP = 5 * MAX_PIPE + 2 * MAX_TANK

--Counters
INTEGER npipe = 0	--Counter for the number of pipes
INTEGER ntank = 0 --Counter for the number of tanks
INTEGER nflow = 0 --Counter for the number of outlet flow conditions
INTEGER node  = 0 --Counter for the number of nodes

PORT Hydraulic
	EQUAL INTEGER inode  
END PORT

CLASS PipeClass
	DECLS
		INTEGER nd1              "Node at inlet"
		INTEGER nd2              "Node at outlet"
      REAL L      UNITS "m"    "Pipe length"
      REAL D		UNITS "m"    "Pipe inside diametert"
      REAL CHW    UNITS "-"    "Hazen Williams Coefficient"
      REAL Qo     UNITS "m3/s" "Initial guess flow"
		REAL H1     UNITS "m"    "Piezometric head at inlet"
		REAL H2     UNITS "m"    "Piezometric head at outlet"
		REAL Q      UNITS "m3/s" "Flow through pipe"
END CLASS

PipeClass P[MAX_PIPE]

CLASS TankClass
	DECLS
		INTEGER nod           "Node at tank outlet"
      REAL H   UNITS "m"    "Tank Piezometriz Head"
      REAL Q   UNITS "m3/s" "Tank outflow"
		REAL Qo  UNITS "m3/s" "Initial guess for pipe flow"
END CLASS

TankClass T[MAX_TANK]

CLASS OutflowClass
	DECLS
		INTEGER nod           "Node at tank outlet"
      REAL H   UNITS "m"    "Tank Piezometriz Head"
		REAL Q   UNITS "m3/s" "Initial guess for pipe flow"
END CLASS

OutflowClass OF[MAX_FLOW]
	
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
		INTEGER ipipe = 0              "Pipe identifier"
		DISCR REAL H1 		UNITS "m"    "Piezometric head at inlet node"
		DISCR REAL H2 		UNITS "m"    "Piezometric head at outlet node"
		DISCR REAL T1 		UNITS "m"    "Temperature at inlet node"
		DISCR REAL T2 		UNITS "m"    "Temperature at outlet node"
		DISCR REAL Q 		UNITS "m3/s" "Volume flow"
   INIT PRIORITY 100
		--Obtain topology and fill global data
		IF(ipipe == 0) THEN
			npipe = npipe+1
         ASSERT (npipe <= MAX_PIPE) FATAL "Maximum number of pipe exceeded. \
			          Recompile library with increased Maximum Dimensions"
			ipipe = npipe
			IF(h_in.inode == 0) THEN
            node = node + 1
            h_in.inode = node
         END IF
			IF(h_out.inode == 0) THEN
				node = node + 1
				h_out.inode = node
			END IF
			P[ipipe].nd1 = h_in.inode
			P[ipipe].nd2 = h_out.inode
			P[ipipe].L =  L
			P[ipipe].D  = D
			P[ipipe].CHW = CHW
			P[ipipe].Qo = Qo
		END IF
		--Assign Solution
		IF(AssignSolution==TRUE) THEN
			H1 = P[ipipe].H1 
			H2 = P[ipipe].H2 
			Q =  P[ipipe].Q  
		END IF
END COMPONENT

COMPONENT Tank
	PORTS
		OUT Hydraulic h_out
	DATA
		REAL H = 100.	UNITS "m" "Tank piezometric head"
		REAL Qo = 0.1	UNITS "m3/s" "Initial guess of the tank outflow"
	DECLS
		INTEGER itank = 0           "Tank identifier"
		DISCR REAL Q   UNITS "m3/s" "Tank outflow"
	INIT PRIORITY 90
		--Obtain topology and fill global data
      IF(itank == 0) THEN
         ntank = ntank + 1 
			ASSERT (ntank <= MAX_TANK) FATAL "Maximum number of tanks exceeded. Recompile library with increased Maximum Dimensions"
         itank = ntank
         IF(h_out.inode == 0) THEN
            node = node + 1
            h_out.inode = node
         END IF
			T[itank].nod = h_out.inode
			T[itank].H = H
			T[itank].Qo = Qo
      END IF
		--Assign Solution
		IF(AssignSolution) THEN
		   Q = T[itank].Q
		END IF
END COMPONENT

COMPONENT Outflow
	PORTS
		IN Hydraulic h_in
	DATA
		REAL Q = 0. 	  UNITS "m3/s" "Outflow"
	DECLS
		INTEGER iflow = 0           "Outflow identifier"
		DISCR REAL H      UNITS "m" "Piezometric head at the connected node"
	INIT PRIORITY 90
      IF(iflow == 0) THEN
         nflow = nflow + 1 
			ASSERT (nflow <= MAX_FLOW) FATAL "Maximum number of flows exceeded. Recompile library with increased Maximum Dimensions"
         iflow = nflow
         IF(h_in.inode == 0) THEN
            node = node + 1
            h_in.inode = node
         END IF
			OF[iflow].nod = h_in.inode
			OF[iflow].Q = Q
      END IF
		IF(AssignSolution) THEN
		   H = OF[iflow].H
		END IF
END COMPONENT

FUNCTION NO_TYPE Fhydro_system
	(
		IN INTEGER iopt        "option: 0 Initialize Equations, \
		                                1 Calculate residuals, \
												  2 Calculate Jacobian",       
		IN REAL Xunk[]         "vector of unknowns", 
		IN REAL rinput[]       "real parameter array for input (not used)",
		OUT REAL routput[]     "real parameter array for output (not used)",
	   OUT REAL F[]   		  "vector of residues",
      OUT NLEQN_SYSTEM_BASE  nleqn_system "Non Linear Equation System"
	)
	DECLS
		INTEGER nd1
		INTEGER nd2
		INTEGER iunkq
		REAL H1, H2, Q
		REAL Qmin
		CONST REAL Vmin = 1e-4  "Min. velocity to avoid singular Jacobian at zero flow"
	BODY
		IF(iopt == 0) THEN
			--Initialize the equation number, the unknowns and the equation descriptions
			nleqn_system.neqn = node + npipe + ntank 
			--Initialization of the unknowns
			FOR(j IN 1, node)
				nleqn_system.CreateUnk(j, rinput[1], -1e6, 1e6, concatStrings(concatStrings("Hnode[", integerToString(j)),"]"),0)
			END FOR
			FOR(j IN 1, npipe)
				nleqn_system.CreateUnk(node+j, P[j].Qo, -1e6, 1e6, concatStrings(concatStrings("Qpipe[", integerToString(j)),"]"),0)
			END FOR
			FOR(j IN 1, ntank)
				nleqn_system.CreateUnk(node+npipe+j, T[j].Qo, -1e6, 1e6, concatStrings(concatStrings("Qtank[", integerToString(j)),"]"),0)
			END FOR
		END IF
		IF(iopt == 1 OR iopt == 2) THEN
			FOR(j IN 1, nleqn_system.neqn)
				F[j] = 0
			END FOR
			nleqn_system.Jacob.setToZero()
			--Pipe Equations
			FOR(j IN 1, npipe)
				nd1 = P[j].nd1
				nd2 = P[j].nd2
				iunkq = node+j
				P[j].H1 = Xunk[nd1]
				P[j].H2 = Xunk[nd2]
				P[j].Q  = Xunk[iunkq]
            --Conservation of mass in node nd1
				F[nd1] = F[nd1] - P[j].Q
				nleqn_system.Jacob.stampItem(nd1, iunkq, -1)	
            --Conservation of mass in node nd2
				F[nd2] = F[nd2] + P[j].Q
				nleqn_system.Jacob.stampItem(nd2, iunkq,  1)			
				--Momentum equation in pipe j
				Qmin = 0.785398163*P[j].D**2 * Vmin
		   	F[node+j] = P[j].H1 - P[j].H2 - 10.67 * (P[j].L /(P[j].CHW**1.852*\
                     P[j].D**4.8704))* P[j].Q * (abs(P[j].Q)+Qmin)**0.852
							
         	nleqn_system.Jacob.stampItem(node+j, nd1,  1.)
         	nleqn_system.Jacob.stampItem(node+j, nd2, -1.)
         	nleqn_system.Jacob.stampItem(node+j, iunkq, -  10.67 * (P[j].L /(P[j].CHW**1.852*\
                     P[j].D**4.8704))*((abs(P[j].Q)+Qmin)**0.852 \    
                     +0.852*abs(P[j].Q)*(abs(P[j].Q)+Qmin)**(0.852-1)))
			END FOR
			--Tank Equations
			FOR(j IN 1, ntank) 
				nd2 = T[j].nod
				T[j].Q = Xunk[node+npipe+j]
				--Conservation of mass in node nd2
				F[nd2] = F[nd2] + T[j].Q 
				nleqn_system.Jacob.stampItem(nd2, node+npipe+j, 1)
				--Piezometric head of nd2 in node nd2
				F[node+npipe+j] = Xunk[nd2] - T[j].H
				nleqn_system.Jacob.stampItem(node+npipe+j, nd2, 1.)
			END FOR
			--Flow Equations
			FOR(j IN 1, nflow)
				nd1 = OF[j].nod
				OF[j].H = Xunk[nd1]
				--Conservation of mass in node nd1
				F[nd1] = F[nd1] - OF[j].Q
			END FOR
			nleqn_system.Jacob.buildCSR()
		END IF
		RETURN
END FUNCTION

COMPONENT Solver
	DATA
	   REAL Ho = 100. UNITS "m" "Initial guess for the node piezometric head"
		REAL tol = 1.e-7         "Tolerance of the Newton-Rapshon"
		INTEGER iter_damp = 5    "Number of initial damped iteration in NR algorithm"		
		REAL FRAC_damp = 0.1		 "Damping factor for the first iter_damp iterations"
		INTEGER iter_max = 20    "Maximum number of iterations in NR algorithm"
		INTEGER NR_debug = 1     "Debug level for Newton-Raphson"
		INTEGER STAMP_debug = 0	 "Debug level for Stamping process"
	DECLS
		HIDDEN DISCR REAL rinput[1]
		HIDDEN DISCR REAL routput[1]
	OBJECTS
		NLEQN_SYSTEM(MAX_EQN=MAX_EQN, MAX_STAMP = MAX_STAMP) nleqn_system
	INIT PRIORITY	0
		IF(NOT AssignSolution) THEN
			nleqn_system.NR_debug = NR_debug
			nleqn_system.Jacob.STAMP_debug = STAMP_debug
			--The initial guess for the node piezometric head is passed through  
			--the auxiliary rinput vector
			rinput[1] = Ho
			nleqn_system.Solve(Fhydro_system, iter_max, iter_damp, FRAC_damp, tol, rinput, routput)  
			AssignSolution = TRUE
			EXEC_INIT()
		END IF
END COMPONENT

