/*-----------------------------------------------------------------------------------------
 LIBRARY: NEW_TRANSINENT APPROACH
 FILE: Inline
 CREATION DATE: 02/01/2024
-----------------------------------------------------------------------------------------*/
COMPONENT Analytic
	DATA
		REAL a = -1.  "Data, equation to be solved x' = a*x"
		REAL h = 0.1  "Integration step"
		REAL xo = 10. "Initial value"
	DECLS
		REAL x
	CONTINUOUS
		x = xo * exp(a*TIME)
END COMPONENT
/*-----------------------------------------------------------------------------------------
  COMPONENT Inline1 implements an inline integration in the discrete block
  This implementacion raises discontinuities at the update times
-----------------------------------------------------------------------------------------*/  
COMPONENT Inline1 
"Component to solve x' = a*x by inline integration - Discrete Version"
	DATA
		REAL a = -1.  "Data, equation to be solved x' = a*x"
		REAL h = 0.25 "Integration step"
		REAL xo = 10. "Initial value"
	DECLS
		DISCR REAL x  "Unknown Variable"
		DISCR REAL x_old "Previous stored value"
		BOOLEAN IsTimeStep = FALSE "Boolean flag to indicate end of inline time step" 
	INIT
		x_old = xo
		x = xo
		IsTimeStep = TRUE AFTER h
	DISCRETE
		WHEN(IsTimeStep) THEN
			IsTimeStep = FALSE
			IsTimeStep = TRUE AFTER h
			x_old = x
			x = x_old * (1 + 0.5*a*h) / (1. - 0.5*a*h)
		END WHEN
END COMPONENT 
/*-----------------------------------------------------------------------------------------
  COMPONENT Inline2 implements an inline integration in the continuous block
  This implementacion results in a continuous integration
-----------------------------------------------------------------------------------------*/  
COMPONENT Inline2
"Component to solve x' = a*x by inline integration - Continuous Version"
	DATA
		REAL a = -1.  "Data, equation to be solved x' = a*x"
		REAL h = 0.25  "Integration step"
		REAL xo = 10. "Initial value"
	DECLS
		REAL dt
		REAL x
		DISCR REAL x_old
		DISCR REAL TIME_old 
		BOOLEAN IsTimeStep
	INIT
		x_old = xo
		TIME_old = TIME
		x = xo
		IsTimeStep = TRUE
	DISCRETE
		WHEN(IsTimeStep) THEN
			IsTimeStep = FALSE
			IsTimeStep = TRUE AFTER h
			x_old = x
			TIME_old = TIME
		END WHEN
	CONTINUOUS
		dt = (TIME-TIME_old)
		x = x_old * (1 + 0.5*a*dt) / (1. - 0.5*a*dt)
END COMPONENT 
/*-----------------------------------------------------------------------------------------
   Class for hiding the equation
-----------------------------------------------------------------------------------------*/  
CLASS SimpleEquation
	DECLS
		REAL x_old   "Stored state"
	METHODS
		METHOD NO_TYPE updateState(IN REAL x)
			BODY 
				x_old = x
		END METHOD
		METHOD REAL increment(IN REAL a, IN REAL dt)
			BODY
				RETURN  x_old * (1 + 0.5*a*dt) / (1. - 0.5*a*dt)
		END METHOD
END CLASS
/*-----------------------------------------------------------------------------------------
   Component Inline3. Similar to Inline2 but hiding the equation details in a class
-----------------------------------------------------------------------------------------*/ 
COMPONENT Inline3
"Component to solve x' = a*x by inline integration - Continuous Version & State Hiding"
	DATA
		REAL a = -1.  "Data, equation to be solved x' = a*x"
		REAL h = 0.25  "Integration step"
		REAL xo = 10. "Initial value"
	DECLS
		REAL dt
		REAL x
		DISCR REAL TIME_old 
		BOOLEAN IsTimeStep
	OBJECTS
		SimpleEquation eqn
	INIT
		eqn.updateState(xo)
		TIME_old = TIME
		IsTimeStep = TRUE
	DISCRETE
		WHEN(IsTimeStep) THEN
			IsTimeStep = FALSE
			IsTimeStep = TRUE AFTER h
			eqn.updateState(x)
			TIME_old = TIME
		END WHEN
	CONTINUOUS
		dt = (TIME-TIME_old)
		x =eqn.increment(a, dt)
END COMPONENT 

ENUM Method = {ForwardEuler, BackwardEuler, Trapezoidal}
/*-----------------------------------------------------------------------------------------
  COMPONENT Inline01. Very similar to Inline1, but offering four integration methods
  This implementacion raises discontinuities at the update times
-----------------------------------------------------------------------------------------*/ 
COMPONENT Inline01
"Component to solve x' = a*x by inline integration - Improved Discrete Version"
--it offers 4 integration methods: Forward Euler, Backward Euler, Trapezoidal 
--and second order Gear
DATA
		REAL a = -1.  "Data, equation to be solved x' = a*x"
		REAL h = 0.25  "Integration step"
		ENUM Method method = Trapezoidal
		REAL xo = 10. "Initial value"
	DECLS
		DISCR REAL x
		DISCR REAL x_old_1, x_old_2
		BOOLEAN IsTimeStep = FALSE
	INIT
		x_old_1 = xo
		x_old_2 = xo
		x = xo
		IsTimeStep = TRUE AFTER h
	DISCRETE
		WHEN(IsTimeStep) THEN
			IsTimeStep = FALSE
			IsTimeStep = TRUE AFTER h
			x_old_2 = x_old_1
			x_old_1 = x
			IF(method==ForwardEuler) THEN
				x = x_old_1 * (1+a*h)
			ELSEIF(method == BackwardEuler) THEN
				x = x_old_1 /(1-a*h)		
			ELSEIF(method == Trapezoidal) THEN
				x = x_old_1 * (1 + 0.5*a*h) / (1. - 0.5*a*h)
			END IF
		END WHEN
END COMPONENT 
/*-----------------------------------------------------------------------------------------
  COMPONENT Inline02. Very similar to Inline2, but offering four integration methods
  This implementacion results in a continuous integration
 -----------------------------------------------------------------------------------------*/ 
COMPONENT Inline02
"Component to solve x' = a*x by inline integration - Improved Discrete Version"
--it offers 4 integration methods: Forward Euler, Backward Euler, Trapezoidal 
--and second order Gear
	DATA
		REAL a = -1.  "Data, equation to be solved x' = a*x"
		REAL h = 0.25  "Integration step"
		ENUM Method method = Trapezoidal
		REAL xo = 10. "Initial value"
	DECLS
		REAL x
		DISCR REAL x_old_1, x_old_2
		DISCR REAL TIME_old
		BOOLEAN IsTimeStep
		REAL dt
	INIT
		x_old_1 = xo
		x_old_2 = xo
		TIME_old = TIME
		x = xo
		IsTimeStep = TRUE
	DISCRETE
		WHEN(IsTimeStep) THEN
			IsTimeStep = FALSE
			IsTimeStep = TRUE AFTER h
			x_old_2 = x_old_1
			x_old_1 = x
			TIME_old = TIME
		END WHEN
	CONTINUOUS
		SEQUENTIAL
			dt = (TIME-TIME_old)
			IF(method==ForwardEuler) THEN
				x = x_old_1 * (1+a*dt)
			ELSEIF(method == BackwardEuler) THEN
				x = x_old_1 /(1-a*dt)		
			ELSEIF(method == Trapezoidal) THEN
				x = x_old_1 * (1 + 0.5*a*dt) / (1. - 0.5*a*dt)
			END IF
		END SEQUENTIAL
END COMPONENT 

/*-----------------------------------------------------------------------------------------
   Class for hiding the equation , but offering four integration methods
-----------------------------------------------------------------------------------------*/  
CLASS SimpleEquation0
	DECLS
		REAL x_old_1   "Stored state at previous time step"
		REAL x_old_2   "Stored state at two previous time steps"
	METHODS
		METHOD NO_TYPE updateState(IN REAL x)
			BODY 
			   x_old_2 = x_old_1
				x_old_1 = x
		END METHOD
		METHOD REAL increment(IN REAL a, IN REAL dt, IN ENUM Method method)
			DECLS 
				REAL x
			BODY
				IF(method==ForwardEuler) THEN
					x = x_old_1 * (1+a*dt)
				ELSEIF(method == BackwardEuler) THEN
					x = x_old_1 /(1-a*dt)		
				ELSEIF(method == Trapezoidal) THEN
					x = x_old_1 * (1 + 0.5*a*dt) / (1. - 0.5*a*dt)	
				END IF
				RETURN  x
		END METHOD
END CLASS

COMPONENT Inline03
"Component to solve x' = a*x by inline integration - Continuous Version & State Hiding"
	DATA
		REAL a = -1.  "Data, equation to be solved x' = a*x"
		REAL h = 0.25  "Integration step"
		REAL xo = 10. "Initial value"
 		ENUM Method method = Trapezoidal
	DECLS
		REAL dt
		REAL x
		DISCR REAL TIME_old 
		BOOLEAN IsTimeStep
	OBJECTS
		SimpleEquation0 eqn
	INIT
		eqn.updateState(xo)
		TIME_old = TIME
		IsTimeStep = TRUE
	DISCRETE
		WHEN(IsTimeStep) THEN
			IsTimeStep = FALSE
			IsTimeStep = TRUE AFTER h
			eqn.updateState(x)
			TIME_old = TIME
		END WHEN
	CONTINUOUS
		dt = (TIME-TIME_old)
		x =eqn.increment(a, dt, method)
END COMPONENT 

COMPONENT AllInlineComps
	DECLS
		REAL error
	TOPOLOGY
		Analytic analytic
		Inline1 inline1
		Inline2 inline2
		Inline3 inline3
		Inline01 inline01(method = Trapezoidal)
		Inline02 inline02(method = Trapezoidal)
		Inline03 inline03(method = Trapezoidal)
END COMPONENT		