/*-----------------------------------------------------------------------------------------
 LIBRARY: HYDRAULIC_ST
 FILE: hydro_v4.el
 CREATION DATE: 19/10/2023
-----------------------------------------------------------------------------------------*/
USE SPARSE
BOOLEAN Transient = FALSE
BOOLEAN AssignSolution = FALSE
BOOLEAN UpdateState = FALSE
ENUM INIT_COMMANDS = {GetTopology_FillData, Assign_Solution}
--Maximum Dimensions
CONST INTEGER MAX_PIPE = 1000
CONST INTEGER MAX_TANK = 100
CONST INTEGER MAX_FLOW = 100

CONST INTEGER MAX_EQN = 2 * MAX_PIPE + 1 + MAX_TANK
--Each pipe makes 5 stamps in the Jacobian and
--each tank makes 2 stamps
CONST INTEGER MAX_STAMP = 5 * MAX_PIPE + 2 * MAX_TANK

--Counnters
INTEGER npipe = 0	--Counter for the number of pipes
INTEGER ntank = 0 --Counter for the number of tanks
INTEGER nflow = 0 --Counter for the number of outlet flow conditions
INTEGER node  = 0 --Counter for the number of nodes

DISCR REAL DT
DISCR REAL Time0
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
		INTEGER nod           	 "Node at tank outlet"
		REAL A      UNITS "m2"   "Flow area"
		REAL Ho		UNITS "m"	 "Initial piezometric head"
      REAL H   	UNITS "m"    "Tank piezometriz head"
		REAL H_old	UNITS "m"	 "Tank piezometriz head at previous time step"
      REAL Q   	UNITS "m3/s" "Tank outflow"
      REAL Q_old 	UNITS "m3/s" "Tank outflow at previous time step"
		REAL Qo  	UNITS "m3/s" "Initial guess for pipe flow"
END CLASS

TankClass T[MAX_TANK]

CLASS OutflowClass
	DECLS
		INTEGER nod           "Node at tank outlet"
      REAL H   UNITS "m"    "Tank Piezometriz Head"
		REAL Q   UNITS "m3/s" "Initial guess for pipe flow"
END CLASS

OutflowClass OF[MAX_FLOW]

FUNCTION NO_TYPE InitPipe(
		IN REAL L,
		IN REAL D,
		IN REAL CHW,
		IN REAL Qo,
		OUT INTEGER ipipe,
		OUT INTEGER inode1,
		OUT INTEGER inode2,
		OUT REAL H1, 
		OUT REAL H2,
		OUT REAL Q)
	BODY
		--Obtain topology and fill global data
		IF(ipipe == 0) THEN
			npipe = npipe+1
         ASSERT (npipe <= MAX_PIPE) FATAL "Maximum number of pipe exceeded. \
			              Recompile library with increased Maximum Dimensions"
			ipipe = npipe
			IF(inode1 == 0) THEN
            node = node + 1
            inode1 = node
         END IF
			IF(inode2 == 0) THEN
				node = node + 1
				inode2 = node
			END IF
			P[ipipe].nd1 = inode1
			P[ipipe].nd2 = inode2
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
END FUNCTION

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
		INTEGER ipipe = 0
		DISCR REAL H1 		UNITS "m"    "Piezometric head at inlet"
		DISCR REAL H2 		UNITS "m"    "Piezometric head at outlet"
		DISCR REAL Q 		UNITS "m3/s" "Volume flow"
   INIT PRIORITY 100
		InitPipe(L, D, CHW, Qo, ipipe, h_in.inode, h_out.inode, H1, H2, Q)
END COMPONENT

FUNCTION NO_TYPE InitTank(
		IN REAL A,
		IN REAL Ho, 
		IN REAL Qo, 
		OUT INTEGER itank,
		OUT INTEGER inode,
		OUT REAL H,
		OUT REAL Q
		)
	BODY
		--Obtain topology and fill global data
      IF(itank == 0) THEN
         ntank = ntank + 1 
			ASSERT (ntank <= MAX_TANK) FATAL "Maximum number of tanks exceeded. \
			               Recompile library with increased Maximum Dimensions"
         itank = ntank
         IF(inode == 0) THEN
            node = node + 1
            inode = node
         END IF
			T[itank].nod = inode
			T[itank].A = A
			T[itank].Ho = Ho
			T[itank].Qo = Qo
      END IF
		--Assign Solution
		IF(AssignSolution) THEN
			H = T[itank].H
		   Q = T[itank].Q
		END IF
		--IF (UpdateState) THEN 
		--	T[itank].H_old  = T[itank].H
		--	T[itank].Q_old  = T[itank].Q
		--END IF
		RETURN
END FUNCTION

COMPONENT Tank
	PORTS
		OUT Hydraulic h_out
	DATA
		REAL A  = 1000. UNITS "m2"   "Surface area"
		REAL Ho = 100.	 UNITS "m"    "Tank piezometric head"
		REAL Qo = 0.1	 UNITS "m3/s" "Initial guess of the tank outflow"
	DECLS
		INTEGER itank = 0
		DISCR REAL Q   UNITS "m3/s" "Tank outflow"
		DISCR REAL H   UNITS "m3/s" "Tank piezometric head"
	INIT PRIORITY 90
		InitTank(A, Ho, Qo, itank, h_out.inode, H, Q)
END COMPONENT

FUNCTION NO_TYPE InitOutflow(
		IN REAL Q, 
		OUT INTEGER iflow,
		OUT INTEGER inode,
		OUT REAL H
	)
	BODY
		--Obtain topology and fill global data
      IF(iflow == 0) THEN
         nflow = nflow + 1 
			ASSERT (nflow <= MAX_FLOW) FATAL "Maximum number of OutFlows exceeded. \
			               Recompile library with increased Maximum Dimensions"
         iflow = nflow
         IF(inode == 0) THEN
            node = node + 1
            inode = node
         END IF
			OF[iflow].nod = inode
			OF[iflow].Q = Q
      END IF
		--Assign solution
		IF(AssignSolution) THEN
		   H = OF[iflow].H
		END IF
		RETURN 
END FUNCTION

COMPONENT Outflow
	PORTS
		IN Hydraulic h_in
	DATA
		REAL Q = 0. 	UNITS "m3/s" "Outflow"
	DECLS
		INTEGER iflow = 0
		DISCR REAL H   UNITS "Piezometric head at the connected node"
	INIT PRIORITY 90
		InitOutflow(Q, iflow, h_in.inode, H)
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
		REAL Time
	BODY
		IF(iopt == 0 AND NOT Transient) THEN
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
			Time = rinput[2]
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
				P[j].Q  = Xunk[j+node]
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
				T[j].H = Xunk[nd2]
				--Conservation of mass in node nd2
				F[nd2] = F[nd2] + T[j].Q 
				nleqn_system.Jacob.stampItem(nd2, node+npipe+j, 1)
				--Piezometric head of nd2 in node nd2
				IF(NOT Transient) THEN
					F[node+npipe+j] = T[j].H - T[j].Ho
					nleqn_system.Jacob.stampItem(node+npipe+j, nd2, 1.)
					T[j].H_old = T[j].Ho
				ELSE
					F[node+npipe+j] = T[j].H + 0.5*T[j].Q * (Time-Time0) / T[j].A - T[j].H_old + 0.5*T[j].Q_old * (Time-Time0) / T[j].A 
					nleqn_system.Jacob.stampItem(node+npipe+j, nd2, 1.)
					nleqn_system.Jacob.stampItem(node+npipe+j, node+npipe+j, .5*T[j].Q * (Time-Time0) / T[j].A)
				END IF
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
		REAL dt = 0.1            UNITS "s" "Time Step" 
	DECLS
		HIDDEN REAL rinput[2]
		HIDDEN REAL routput[1]
		BOOLEAN IsTimeStep = TRUE 
	OBJECTS
		NLEQN_SYSTEM(MAX_EQN=MAX_EQN, MAX_STAMP = MAX_STAMP) nleqn_system
	INIT PRIORITY	0
		IF(NOT AssignSolution) THEN
		  -- DT = dt
			nleqn_system.NR_debug = NR_debug
			nleqn_system.Jacob.STAMP_debug = STAMP_debug
			--The initial guess for the node piezometric head is passed through  
			--the auxiliary rdata vector
			rinput[1] = Ho
			nleqn_system.Solve(Fhydro_system, iter_max, iter_damp, FRAC_damp, tol, rinput, routput)  
			AssignSolution = TRUE
			UpdateState = TRUE
			EXEC_INIT()
			UpdateState = FALSE
			Time0 = TIME
			IsTimeStep = TRUE AFTER dt
			Transient = TRUE
		END IF

	DISCRETE
		WHEN(IsTimeStep == TRUE) THEN
			rinput[2] = TIME
			--nleqn_system.Solve(Fhydro_system, iter_max, iter_damp, FRAC_damp, tol, rinput, routput)
			AssignSolution = TRUE
			UpdateState = TRUE
			FOR(j IN 1, ntank)
				T[j].H_old  = T[j].H
				T[j].Q_old  = T[j].Q
			END FOR
			EXEC_INIT()
			IsTimeStep = FALSE 
			WRITE("dt=%lg\n", dt)
			IsTimeStep = TRUE AFTER dt
			UpdateState = FALSE
			Time0 = TIME
		END WHEN
	CONTINUOUS
		SEQUENTIAL
			rinput[2] = TIME
			nleqn_system.Solve(Fhydro_system, iter_max, iter_damp, FRAC_damp, tol, rinput, routput)	
			AssignSolution = TRUE
			EXEC_INIT()
		END SEQUENTIAL
END COMPONENT

