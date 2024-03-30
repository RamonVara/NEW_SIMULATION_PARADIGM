
COMPONENT Pendulum
	DATA
		REAL m = 1      UNITS "kg"   "Pendulum mass"
		REAL L = 1      UNITS "m"    "Pendulum length"
		REAL g = 9.8    UNITS "m2/s" "Gravity"
		REAL Fx_o = 9.8 UNITS "N"    "Initial force"
	DECLS
		REAL x, y
		REAL T
	CONTINUOUS
		x**2 + y**2 = L**2
		0 = Fx_o - T * x/L
		0 = -m*g - T * y/L 
END COMPONENT

COMPONENT Pendulum1
	DATA
		REAL m = 1      UNITS "kg"   "Pendulum mass"
		REAL L = 1      UNITS "m"    "Pendulum length"
		REAL g = 9.8    UNITS "m2/s" "Gravity"
	DECLS
		REAL Fx
		--REAL Fy
		REAL x, y
		REAL T
	CONTINUOUS
		x**2 + y**2 = L**2
		m*x'' = Fx - T * x/L
		m*y'' = -m*g + T * y/L 
END COMPONENT