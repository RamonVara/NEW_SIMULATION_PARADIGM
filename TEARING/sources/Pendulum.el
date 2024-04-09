COMPONENT PendulumTransient
	DATA
		REAL m = 1     
		REAL L = 1      
		REAL g = 9.8   
	DECLS
		BOUND REAL Fx
		REAL x, y
		REAL T
	CONTINUOUS
		x**2 + y**2 = L**2
		m*x'' = Fx - T * x/L
		m*y'' = -m*g + T * y/L 
END COMPONENT

COMPONENT PendulumSteady
	DATA
		REAL m = 1      
		REAL L = 1      
		REAL g = 9.8    
	DECLS
	   BOUND REAL Fx
		REAL x, y
		REAL T
	CONTINUOUS
		x**2 + y**2 = L**2
		0 = Fx  - T * x/L
		0 = -m*g - T * y/L 
END COMPONENT

