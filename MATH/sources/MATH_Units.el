/*-----------------------------------------------------------------------------------------
 LIBRARY: MATH
 FILE: MATH_Units
 AUTHOR: EA
 COMPANY: EA
 DESCRIPTION: Variables to define Units for variables of libraries 
 CREATION DATE: jueves, marzo 22 2012
-----------------------------------------------------------------------------------------*/
PRIVATE CONST STRING no_units						= "-"			"Adimensional"

PRIVATE CONST STRING u_m								=	"m" 		"metre" 
PRIVATE CONST STRING u_kg								=	"kg" 		"kilograme" 
PRIVATE CONST STRING u_s								=	"s" 		"second" 
PRIVATE CONST STRING u_A								=	"A" 		"Ampere" 
PRIVATE CONST STRING u_K								=	"K" 		"Kelvin" 
PRIVATE CONST STRING u_mol							=	"mol" 	"mole" 
PRIVATE CONST STRING u_cd								=	"cd" 		"candela" 
PRIVATE CONST STRING u_rad							=	"rad" 	"radian" 
PRIVATE CONST STRING u_sr								=	"sr" 		"steradian" 

PRIVATE CONST STRING u_Hz								=	"Hz" 		"Hertz" 
PRIVATE CONST STRING u_1_s							=	"1/s" 	"inverse second"
PRIVATE CONST STRING u_1_h							=	"1/h" 	"inverse hour"
PRIVATE CONST STRING u_1_s2 							= 	"1/s^2" 	"squared inverse second"
PRIVATE CONST STRING u_N								=	"N" 		"Newton" 
PRIVATE CONST STRING u_Pa								=	"Pa" 		"Pascal" 
PRIVATE CONST STRING u_1_Pa 							= 	"1/Pa" 	"inverse Pascal"
PRIVATE CONST STRING u_J								=	"J" 		"Joule" 
PRIVATE CONST STRING u_W								=	"W" 		"Watt" 
PRIVATE CONST STRING u_Coulomb						=	"C" 		"Coulomb" 
PRIVATE CONST STRING u_V								=	"V" 		"Volt" 
PRIVATE CONST STRING u_Farad							=	"F" 		"Farad"
PRIVATE CONST STRING u_Ohm							=	"Ohm" 	"Ohm" 
PRIVATE CONST STRING u_S								=	"S" 		"Siemens" 
PRIVATE CONST STRING u_Wb								=	"Wb" 		"Weber" 
PRIVATE CONST STRING u_T								=	"T" 		"Tesla" 
PRIVATE CONST STRING u_H								=	"H" 		"Henry" 
PRIVATE CONST STRING u_lm								=	"lm" 		"lumen" 
PRIVATE CONST STRING u_lx								=	"lx" 		"lux" 
PRIVATE CONST STRING u_Bq								=	"Bq" 		"Becquerel" 
PRIVATE CONST STRING u_Gy								=	"Gy" 		"Gray" 

--[T]
PRIVATE CONST STRING u_min							=	"min" 	"minute"
PRIVATE CONST STRING u_h								=	"h" 		"hour"

--[T·Theta]
PRIVATE CONST STRING u_Ks								=	"K*s"		"kelvin - second"

--[T2·Theta·L2/(M·L)]
PRIVATE CONST STRING u_m2C_W			= "m2·degC/W"					"Square metre - Celsius degree per watt"

--[M]
PRIVATE CONST STRING u_g								=	"g" 		"gram"
PRIVATE CONST STRING u_mg								=	"mg" 		"miligram"
PRIVATE CONST STRING u_lbm							=	"lbm" 	"pound-mass"

PRIVATE CONST STRING u_kmol							=	"kmol" 	"kilomole"

PRIVATE CONST STRING u_pct							=	"%" 		"percentage"
PRIVATE CONST STRING u_ppm							=	"ppm" 	"parts per million"


--[Theta]
PRIVATE CONST STRING u_C								=	"degC" 		"Celsius degree"
PRIVATE CONST STRING u_F								=	"degF" 		"Fahrenheit degree"
PRIVATE CONST STRING u_R								=	"degR" 		"Reamour degree"

PRIVATE CONST STRING u_1_K 							= 	"1/K" 	
PRIVATE CONST STRING u_1_C								=	"1/degC" 	"Inverse Celsius"

--[Theta/L]
PRIVATE CONST STRING u_K_m							= 	"K/m" 	"Kelvin per meter"

--[Theta/T]
PRIVATE CONST STRING u_C_s							=  "degC/s"		"Celsius per second"
PRIVATE CONST STRING u_K_s							=  "K/s"		"Kelvin per second"

-- [1/(Theta·T)]
PRIVATE CONST STRING u_1_hC						=	"1/(hºC)" 			"inverse hours-Celsius degree "

-- [1/(Theta2·T)]
PRIVATE CONST STRING u_1_hC2					=	"1/(hºC2)" 			"inverse hours-suare Celsius degree "

--[L]
PRIVATE CONST STRING u_km								=	"km" 		"kilometre"
PRIVATE CONST STRING u_cm								=	"cm" 		"centimetre"
PRIVATE CONST STRING u_mm								=	"mm" 		"milimetre"
PRIVATE CONST STRING u_ft								=	"ft" 		"foot"
PRIVATE CONST STRING u_in								=	"in" 		"inch"
--[L2]
PRIVATE CONST STRING u_m2								=	"m^2" 		"square metre"
PRIVATE CONST STRING u_km2							=	"km^2" 		"square kilometre"
PRIVATE CONST STRING u_cm2							=	"cm^2" 		"square centimetre"
PRIVATE CONST STRING u_mm2							=	"mm^2" 		"square milimetre"
PRIVATE CONST STRING u_ft2							=	"ft^2" 		"square foot"
PRIVATE CONST STRING u_in2							=	"in^2" 		"square inch"
--[L3]
PRIVATE CONST STRING u_m3								=	"m^3" 		"cubic metre"
PRIVATE CONST STRING u_dm3							=	"dm^3" 		"cubic decimetre"
PRIVATE CONST STRING u_cm3							=	"cm^3" 		"cubic centimetre"
PRIVATE CONST STRING u_ft3							=	"ft^3" 		"cubic foot"
PRIVATE CONST STRING u_in3							=	"in^3" 		"cubic inch"

PRIVATE CONST STRING u_1_m							=	"1/m"			"reciprocal meter"

--[L/T]
PRIVATE CONST STRING u_m_s							=	"m/s" 		"metres per second"
PRIVATE CONST STRING u_km_h							=	"km/h" 		"kilometres per hour"
PRIVATE CONST STRING u_cm_s							=	"cm/s"		"centimetres per second"
PRIVATE CONST STRING u_mm_s							=	"mm/s" 		"milimetres per second"
PRIVATE CONST STRING u_ft_s							=	"ft/s" 		"feet per second"
PRIVATE CONST STRING u_ft_min						=	"ft/min"		"feet per minute"
PRIVATE CONST STRING u_ft_h							=	"ft/h" 		"feet per hour"
PRIVATE CONST STRING u_in_s							=	"in/s" 		"inches per second"
PRIVATE CONST STRING u_in_min						=	"in/min"		"inches per minute"

PRIVATE CONST STRING u_m3_sm2						=	"m3/(s·m^2)" 		"cubic metres per second and square metre"
PRIVATE CONST STRING u_ft3_sft2						=	"ft3/(s·ft^2)" 	"cubic feet per second and square feet"
PRIVATE CONST STRING u_ft3_minft					=	"ft3/(min·ft^2)" 	"cubic feet per minute and square feet"

--[L/T2]
PRIVATE CONST STRING u_m_s2							=	"m/s^2" 			"metres per square second"
PRIVATE CONST STRING u_cm_s2							=	"cm/s^2" 		"centimetres per square second"
PRIVATE CONST STRING u_ft_s2							=	"ft/s^2" 		"feet per square second"

--[L2/T]
PRIVATE CONST STRING u_m2_s							=	"m^2/s" 				"square metres per second"
PRIVATE CONST STRING u_mm2_s							=	"mm^2/s" 			"square milimetres per second"
PRIVATE CONST STRING u_m3_sm							=	"m3/(s·m)" 			"metres per second and metre"
PRIVATE CONST STRING u_ft2_s							=	"ft^2/s" 			"square feet per second"
PRIVATE CONST STRING u_in2_s							=	"in^2/s" 			"square inches per second"
PRIVATE CONST STRING u_UKgal_hin					=	"UKgal/(h·in)" 	"Imperial gallon per hour and inch"
PRIVATE CONST STRING u_UKgal_hft					=	"UKgal/(h·ft)" 	"Imperial gallon per hour and foot"

--[L2/(T2·Theta)]
PRIVATE CONST STRING u_m2_s2K 						= 	"m^2/(s^2·K)" 		"squared meter per squared second and Kelvin"
PRIVATE CONST STRING u_J_Wh 						= 	"J/(W·h)" 			"Joules per Watt and hour"

--[L3/T]
PRIVATE CONST STRING u_m3_s							=	"m^3/s" 				"cubic metres per second"
PRIVATE CONST STRING u_m3_h							=	"m^3/h" 				"cubic metres per hour"
PRIVATE CONST STRING u_dm3_s							=	"dm^3/s" 			"cubic decimetres per second"
PRIVATE CONST STRING u_l_s							=	"l/s" 				"litres per second"
PRIVATE CONST STRING u_ft3_s							=	"ft^3/s" 			"cubic feet per second"
PRIVATE CONST STRING u_ft3_min						=	"ft^3/min" 			"cubic feet per minute"
PRIVATE CONST STRING u_ft3_h							=	"ft^3/h" 			"cubic feet per hour"
PRIVATE CONST STRING u_cm3_s						= 	"cm^3/s" 			"cubic centimeters per second"


--[L3/T]
PRIVATE CONST STRING u_m3_m3							=	"m^3/m^3" 			"cubic metres per cubic metre"
PRIVATE CONST STRING u_dm3_dm3						=	"dm^3/dm^3" 		"cubic decimetres per cubic decimetre"
PRIVATE CONST STRING u_cm3_cm3						=	"cm^3/cm^3" 		"cubic centimetres per cubic centimetre"
PRIVATE CONST STRING u_ft3_ft3						=	"ft^3/ft^3" 		"cubic feet per cubic foot"

--[L3/M]
PRIVATE CONST STRING u_m3_kg							=	"m^3/kg" 			"cubic metres per kilogram"
PRIVATE CONST STRING u_ft3_lbm						=	"ft^3/lbm" 			"cubic feet per pound-mass"
PRIVATE CONST STRING u_m3_g							=	"m^3/g" 				"cubic metres per gram"
PRIVATE CONST STRING u_dm3_kg						=	"dm^3/kg" 			"cubic decimetres per kilogram"

--[L3/N]
PRIVATE CONST STRING u_m3_kmol						=	"m^3/kmol" 			"cubic metres per kilomole"
PRIVATE CONST STRING u_dm3_kmol						=	"dm^3/kmol"			"cubic decimetres per kilomole"
PRIVATE CONST STRING u_ft3_lbmol					=	"ft^3/lbmol" 		"cubic feet per pound-mole"


--[M/L3]
PRIVATE CONST STRING u_kg_m3							=	"kg/m^3" 			"kilograms per cubic metre"
PRIVATE CONST STRING u_g_cm3							=	"g/cm^3" 			"grams per cubic centimetre"
PRIVATE CONST STRING u_g_m3							=	"g/m^3" 				"grams per cubic metre"
PRIVATE CONST STRING u_mg_m3							=	"mg/m^3" 			"miligrams per cubic metre"
PRIVATE CONST STRING u_lbm_ft3						=	"lbm/ft^3" 			"pounds-mass per cubic foot"

--[M/(L3·Theta)]
PRIVATE CONST STRING u_kg_m3C 					=	 "kg/(m^3·degC)" 			"kilogram per cubic meter and Celsius degree"
PRIVATE CONST STRING u_kg_m3K 						= "kg/(m^3·K)" 				"kilogram per cubic meter and Kelvin"

--[M/M]
PRIVATE CONST STRING u_kg_kg							=	"kg/kg" 				"kilograms per kilogram"
PRIVATE CONST STRING u_g_kg							=	"g/kg" 				"grams per kilogram"
PRIVATE CONST STRING u_mg_kg							=	"mg/kg" 				"miligrams per kilogram"

--[M/(M·T)]
PRIVATE CONST STRING u_kg_skg						= "kg/(kg_total·s)"		"kinetic equation of each reaction (kg/kg_total s)"


--[M/N]
PRIVATE CONST STRING u_kg_kmol						= "kg/kmol" 				"kilograms per kilomol"
PRIVATE CONST STRING u_kg_mol 						= "kg/mol" 					"kilogram per mol"
PRIVATE CONST STRING u_g_mol 						= 	"g/mol"				 "grams per mol"





--[M/T]
PRIVATE CONST STRING u_kg_s							=	"kg/s" 			"kilograms per second"
PRIVATE CONST STRING u_kg_h							=	"kg/h" 			"kilograms per hour"
PRIVATE CONST STRING u_lbm_s							=	"lbm/s" 			"pounds-mass per second"
PRIVATE CONST STRING u_lbm_h							=	"lbm/h" 			"pounds-mass per hour"
PRIVATE CONST STRING u_g_h								=	"g/h" 				"grams per hour"		

--[M/T2]
PRIVATE CONST STRING u_kg_s2 						= 	"kg/(s^2)" 				"kilogram squared second"



--[M/TL]
PRIVATE CONST STRING u_kg_sm							=	"kg/(s·m)" 			"kilograms per metre per second"
PRIVATE CONST STRING u_lbm_sft						=	"lbm/(s·ft)" 		"pounds-mass per foot per second"
PRIVATE CONST STRING u_lbm_hft						=	"lbm/(h·ft)" 		"pounds-mass per foot per hour"

--[M/TL2]
PRIVATE CONST STRING u_kg_sm2						=	"kg/(s·m^2)" 			"kilograms per square metre per second"
PRIVATE CONST STRING u_lbm_sft2						=	"lbm/(s·ft^2)" 		"pounds-mass per square foot per second"
PRIVATE CONST STRING u_lbm_hft2						=	"lbm/(h·ft^2)" 		"pounds-mass per square foot per hour"

--[M/(T·L3)]

PRIVATE CONST STRING u_kg_sm3						= "kg/(s·m^3)" 			"kilograms per cubic metre per second"

PRIVATE CONST STRING u_kg_Ns						= "kg/(N·s)" 			"kilograms per Newton and second"


--[M/L]
PRIVATE CONST STRING u_kg_m							=	"kg/m" 				"kilograms per metre"
PRIVATE CONST STRING u_lbm_ft						=	"lbm/ft" 			"pounds-mass per foot"

--[M/L2]
PRIVATE CONST STRING u_kg_m2							=	"kg/m^2" 			"kilograms per square metre"
PRIVATE CONST STRING u_lbm_ft2						=	"lbm/ft^2" 			"pounds-mass per square foot"

--[M/(L2·Theta)]
PRIVATE CONST STRING u_kg_m2K 						= 	"kg/(m^2·K)" 			"kilogram per squared meter and Kelvin"


--[ML/T]
PRIVATE CONST STRING u_kgm_s							=	"kg·m/s" 			"kilograms - metre per second"
PRIVATE CONST STRING u_lbmft_s						=	"lbm·ft/s" 			"pounds-mass - foot per second"

--[ML2]
PRIVATE CONST STRING u_kgm2							=	"kg·m^2" 			"kilograms - square metre"
PRIVATE CONST STRING u_lbmft2						=	"lbm·ft^2" 			"pounds-mass - square foot"


--[N/T]
PRIVATE CONST STRING u_kmol_s						=	"kmol/s" 			"kilomole per second"
PRIVATE CONST STRING u_kmol_h						=	"kmol/h" 			"kilomole per hour"
PRIVATE CONST STRING u_lbmol_s						=	"lbmol/s"	 		"pounds-mole per second"
PRIVATE CONST STRING u_lbmol_h						=	"lbmol/h" 			"pounds-mole per hour"

--[N/L3]
PRIVATE CONST STRING u_mol_m3						=	"mol/m^3" 			"mole per cubic metre"
PRIVATE CONST STRING u_kmol_m3						=	"kmol/m^3" 			"kilomole per cubic metre"
PRIVATE CONST STRING u_lbmol_ft3					=	"lbmol/ft^3" 		"pounds-mole per cubic foot"

--[N/N]
PRIVATE CONST STRING u_kmol_kmol					=	"kmol/kmol" 			"Molar fraction expressed in kilomole"
PRIVATE CONST STRING u_mol_mol 						=	"mol/mol"  				"Molar fraction expressed in mole"
--[M·L/T]
PRIVATE CONST STRING u_kgm2_s						= 	"kg·m^2/s" 			"kilograms-squared meters per second"

--[M·L/T2]
PRIVATE CONST STRING u_kgf							=	"kgf" 				"kilogram-force"
PRIVATE CONST STRING u_dyn							=	"dyn" 				"Dyne"
PRIVATE CONST STRING u_lbf						=	 "lbf" 			"pounds-force"

--[M·L/(T2·Theta)]
PRIVATE CONST STRING u_kcal_C				=	"kcal/degC" 			"kilocalories per celsius degree"				

--[T2/L2]
PRIVATE CONST STRING u_g_Wh						=	"g/(W·h)"			"grams per Watt and hour"			

--[T2·M/(M·L·T)]
PRIVATE CONST STRING u_kg_kgfh						= "kg/(kgf·h)" 			"kilograms per kilogram-force and hour"

--[T2·L2·M/(M·L·T)]
PRIVATE CONST STRING u_kg_Pas 						= 	"kg/(Pa·s)" 		"kilogram per Pascal and second"


--[M·L2/T2]
PRIVATE CONST STRING u_Nm								=	"N·m" 			"Newton - metres"
PRIVATE CONST STRING u_kgfm							=	"kgf·m" 			"kilograms-force - metres"

--[M·L2/(T2·L)]
PRIVATE CONST STRING u_N_m							=	"N/m" 			"Newton per metre"
PRIVATE CONST STRING u_mN_m							=	"mN/m" 			"miliNewton per metre"
PRIVATE CONST STRING u_dyn_cm						=	"dyn/cm" 		"Dyne per centimetre"

--[M·L/(T2·L2)]
PRIVATE CONST STRING u_MPa							=	"MPa" 			"MegaPascal"
PRIVATE CONST STRING u_kPa							=	"kPa" 			"kiloPascal"
PRIVATE CONST STRING u_atm							=	"atm" 			"Atmosphere"
PRIVATE CONST STRING u_bar							=	"bar" 			"Bar"
PRIVATE CONST STRING u_mmHg							=	"mmHg" 			"milimetres of mercury, torr"
PRIVATE CONST STRING u_N_m2 							= "N/m^2" 				"Newton per squared meter"

--[M·L2/(T2·M·Theta)]
PRIVATE CONST STRING u_J_kgC							=	"J/(kg·degC)" 			"Joules per kilogram - Celsius"

--[M·L2/(T2·N·Theta)]
PRIVATE CONST STRING u_J_kmolK						=	"J/(kmol·K)" 			"Joules per kilomole - Kelvin"

--[M·L2/(T3·Theta)]
PRIVATE CONST STRING u_W_C 							=	"W/degC" 					"Watts per Celsius"

--[M·L2/(T3·L·Theta)]
PRIVATE CONST STRING u_W_mC							=	"W/(m·degC)" 				"Watts per metre - Celsius"

--[M·L2/(T3·L2·Theta)]
PRIVATE CONST STRING u_W_m2C 						= 	"W/(m^2·degC)" 			"Watts per square metre - Celsius"
--[M·L/(T3·L2)]
PRIVATE CONST STRING u_Pa_s 							= 	"Pa/s" 			"Pascal per second"
PRIVATE CONST STRING u_bar_s							= "bar/s"					"bar per second"


--[M·L·T/(T2·L2·L3)]
PRIVATE CONST STRING 	u_Pas_cm3 						= 	"Pa·s/cm^3" 			"Pascal - second per cubic centimeter" 
--[(T2·L2·L3/(M·L)]
PRIVATE CONST STRING 	u_cm3_Pa 						= "cm^3/Pa" 				"Cubic centimeter per Pascal"


--[M/(T·L)]
PRIVATE CONST STRING u_Pas							=	"Pa·s" 			"Pascal - second"
PRIVATE CONST STRING u_dyns_cm2						=	"dyn·s/cm^2" 	"dyne - second per square centimetre"
PRIVATE CONST STRING u_cP								=	"cP" 				"centiPoise"

--[M·L/(T2·L3)]
PRIVATE CONST STRING u_kPa_m							=	"kPa/m" 			"kiloPascals per metre"
PRIVATE CONST STRING u_psi_ft						=	"psi/ft" 		"Pounds per square inch per foot"

-- [M·L/(Theta·T3)]
PRIVATE CONST STRING u_kcal_hC				=	"kcal/(h·degC)" 		"kilocalories per hour and celsius degree"		

-- [M·L/(Theta2·T3)]
CONST STRING u_kcal_hC2				=	"kcal/(h·degC2)" 		"kilocalories per hour and square celsius degree"

--[M·L·M3/(T2·L2·N)]
PRIVATE CONST STRING u_Pam3_mol							= "Pa·m^3/mol"		"activation energy of the Arrhenius equation"


--Energía:
PRIVATE CONST STRING u_kJ								=	"kJ" 			"kiloJoule"
PRIVATE CONST STRING u_kWh							=	"kW·h" 		"kiloWatt - hour"
PRIVATE CONST STRING u_kcal							=	"kcal" 		"kilocalories"
PRIVATE CONST STRING u_Btu							=	"Btu" 		"British Thermal Unit"
PRIVATE CONST STRING u_therm							=	"therm" 		"Therm"

--[M·L2/(T2·L3)]
PRIVATE CONST STRING u_J_m3							=	"J/m^3" 			"Joule per cubic metre"
PRIVATE CONST STRING u_kJ_m3							=	"kJ/m^3" 		"kiloJoule per cubic metre"
PRIVATE CONST STRING u_kcal_m3						=	"kcal/m^3" 		"kilocalories per cubic metre"
PRIVATE CONST STRING u_Btu_ft3						=	"Btu/ft^3" 		"Btu per cubic foot"

--[M·L2/(T2·N)]
PRIVATE CONST STRING u_J_mol							=	"J/mol" 			"Joules per mole"
PRIVATE CONST STRING u_J_kmol						=	"J/kmol"					"Joules per kmole"
PRIVATE CONST STRING u_kJ_kmol						=	"kJ/kmol" 		"kiloJoules per kilomole"
PRIVATE CONST STRING u_kcal_gmol					=	"kcal/gmol" 	"kilocalories per gram-mole"
PRIVATE CONST STRING u_Btu_lbmol					=	"Btu/lbmol" 	"Btu per pound-mole"

--[M·L2/(T2·Theta)]
PRIVATE CONST STRING u_J_K							=	"J/K" 		"Joule per Kelvin"
PRIVATE CONST STRING u_J_C							=	"J/degC"		"Joules per Celsius"

--[M·L2/(T2·N·Theta)]
PRIVATE CONST STRING u_J_molK						=	"J/(mol·K)" 		"Joules per mole - Kelvin"
PRIVATE CONST STRING u_kJ_kmolK						=	"kJ/(kmol·K)" 		"kiloJoules per kilomole - Kelvin"
PRIVATE CONST STRING u_kJ_kmolC 					= "kJ/(kmol·degC)" 		"kiloJoule per kilomol and Celsius degree"
PRIVATE CONST STRING u_kJ_kgC 						= "kJ/(kg·degC)" 		"kiloJoule per kilogram and Celsius degree"

PRIVATE CONST STRING u_cal_gmolC					=	"cal/(gmol·degC)" 	"calories per gram-mole Celsius"
PRIVATE CONST STRING u_Btu_lbmolF					=	"Btu/(lbmol·degF)" 	"Btu per pound-mole Fahrenheit"


--[M·L2/(T2·M)]
PRIVATE CONST STRING u_J_kg							=	"J/kg" 			"Joules per kilogram"
PRIVATE CONST STRING u_kJ_kg							=	"kJ/kg" 			"kiloJoules per kilogram"	
PRIVATE CONST STRING u_cal_g							=	"cal/g" 			"calories per gram"

--[M·T2/(M·L2)]
PRIVATE CONST STRING u_g_kcal				=	"g/kcal" 			"grams per kilocalorie"		

--[M·L2/(T2·M·Theta)]
PRIVATE CONST STRING u_J_kgK							=	"J/(kg·K)" 			"Joules per kilogram - Kelvin"
PRIVATE CONST STRING u_kJ_kgK						=	"kJ/(kg·K)" 		"kiloJoules per kilogram - Kelvin"
PRIVATE CONST STRING u_cal_gK						=	"cal/(g·K)" 		"calories per gram - Kelvin"
PRIVATE CONST STRING u_kcal_kgC						=	"kcal/(kg·degC)" 	"kilocalories per kilogram - Celsius"


--[M·L2·T/(T2·M·Theta)]
PRIVATE CONST STRING u_kJs_kgC						= "kJ·s/(kg·degC)" 		"PRIVATE CONSTant for temperature calculation"

--[M·L2/(T2·L2)]
PRIVATE CONST STRING u_J_cm2							=	"J/cm^2" 			"Joules per square centimetre"
PRIVATE CONST STRING u_lbfft_in2					=	"lbf·ft/in^2" 		"foot-pound force per square inch"


--[M·L2/T3]
PRIVATE CONST STRING u_MW								=	"MW" 				"Magawatts"
PRIVATE CONST STRING u_kW								=	"kW" 				"kilowatts"
PRIVATE CONST STRING u_kcal_h						=	"kcal/h"			"kilocalories per hour"
PRIVATE CONST STRING u_hhp							=	"hhp" 			"hydraulic horsepower"
PRIVATE CONST STRING u_CV								=	"CV" 				"Cheval vapeur"
PRIVATE CONST STRING u_Btu_s							=	"Btu/s" 			"Btu per second"
PRIVATE CONST STRING u_Btu_min						=	"Btu/min" 		"Btu per minute"
PRIVATE CONST STRING u_Btu_h							=	"Btu/h" 			"Btu per hour"
PRIVATE CONST STRING u_HP							= "HP" 					"horse power"
PRIVATE CONST STRING u_kJ_s						=" kJ/s"			"heat transfer"
PRIVATE CONST STRING u_kJ_h						= "kJ/h"		"supplied heat"


--[M·L/(T3·L)]
PRIVATE CONST STRING u_N_ms 							= "N/(m·s)" 		"Newton per metre and second"

--[M·L·T/(T2·L)]
PRIVATE CONST STRING u_Ns_m 							= 	"N·s/m" 			"Newton second per metre"


--[M·L2/T3·L2]
PRIVATE CONST STRING u_W_m2							=	"W/m^2" 				"Watts per square metre"
PRIVATE CONST STRING u_kW_m2							=	"kW/m^2" 			"kilowatts per square metre"
PRIVATE CONST STRING u_cal_scm2						=	"cal/(s·cm^2)" 	"calories per square centimetre per second"
PRIVATE CONST STRING u_cal_hcm2						=	"cal/(h·cm^2)" 	"calories per square centimetre per hour"
PRIVATE CONST STRING u_Btu_sft2						=	"Btu/(s·ft^2)" 	"Btu per square foot per second"
PRIVATE CONST STRING u_Btu_hft2						=	"Btu/(h·ft^2)" 	"Btu per square foot per hour"

--[M·L3/T3·L2·Theta]
PRIVATE CONST STRING u_W_mK						=	"W/(m·K)" 					"Watts per metre Kelvin"
PRIVATE CONST STRING u_calcm_scm2C					=	"cal·cm/(s·cm^2·degC)" 	"calories - centimetre per square centimetre per second - Celsius"
PRIVATE CONST STRING u_Btuft_hft2F					=	"Btu·ft/(h·ft^2·degF)" 	"Btu - foot per square foot per hour - Fahrenheit"
PRIVATE CONST STRING u_kW_m2C						= 	"kW/(m^2·degC)"  		 "overall heat transfer coefficient"

--[M·L2/T3·L3]
PRIVATE CONST STRING u_kW_m3							=	"kW/m^3" 				"kilowatts per cubic metre"
PRIVATE CONST STRING u_hp_ft3						=	"hp/ft^3" 				"horsepower per cubic foot"
PRIVATE CONST STRING u_Btu_sft3						=	"Btu/(s·ft^3)" 		"Btu per cubic metre per second"

--[M·L2/T3·L2·Theta]
PRIVATE CONST STRING u_kW_m2K						=	"kW/(m^2·K)" 				"kilowatts per square metre -kelvin"
PRIVATE CONST STRING u_kcal_hm2C					=	"kcal/(h·m^2·degC)" 		"kilocalorias per square metre per hour Celsius"
PRIVATE CONST STRING u_Btu_hft2F					=	"Btu/(h·ft^2·degF)"			"Btu per square foot per hour Fahrenheit"
PRIVATE CONST STRING u_W_m2K 						= 	"W/(m^2·K)" 			"Watts per square metre Kelvin"

--[M·L2/T3·L3·Theta]
PRIVATE CONST STRING u_kW_m3K						=	"kW/(m^3·K)" 				"kilowatts per cubic metre -kelvin"
PRIVATE CONST STRING u_Btu_sft3F					=	"Btu/(s·ft^3·degF)"			"Btu per cubic foot per second - Fahrenheit"

PRIVATE CONST STRING u_W_m2sr						=	"W/(m^2·sr)" 				"watts per square metre - steradian"
PRIVATE CONST STRING u_W_sr							=	"W/sr" 						"watts per steradian"

--[M·L2/(T3·Theta)]
PRIVATE CONST STRING u_W_K 							= "W/K" 							"Watts per Kelvin"
PRIVATE CONST STRING u_kW_C  						=	"kW/degC"							"kilowatts per Celsius"

--[M·L2/(T3·M·Theta)]
PRIVATE CONST STRING u_W_kgK						= "W/(kg·K)"				"watts per kilogram - kelvin"

--[M·L2/(T3·M)]
PRIVATE CONST STRING u_kJ_skg						= "kJ/(s·kg)"				"specific energy flow"


--[T/L]
PRIVATE CONST STRING u_s_m 							= "s/m" 							"second per metre"
PRIVATE CONST STRING u_g_kNs						= 	"g/(kN·s)" 			"grams per kiloNewton and second"
--[T/M]
PRIVATE CONST STRING u_s_kg						= 	"s/kg" 				"seconds per kilogram"

--
PRIVATE CONST STRING u_Cm2h_kcal					=	"degC·m^2·h/kcal" 			"degree Celsius square metre hour per kilocalorie"
PRIVATE CONST STRING u_Km2_kW						=	"K·m^2/kW"					"Kelvin square metre per kilowatt"

--
PRIVATE CONST STRING u_deg							=	"deg" 						"degree"
PRIVATE CONST STRING u_deg_s							=	"deg/s" 				"degrees per second"
PRIVATE CONST STRING u_rad_s							=	"rad/s" 					"radians per second"
PRIVATE CONST STRING u_rad_s2						=	"rad/s^2" 				"radians per second squared"
PRIVATE CONST STRING u_rpm							=	"rpm" 					"revolutions per minute"
PRIVATE CONST STRING u_rps							=	"rps" 					"revolutions per second"
PRIVATE CONST STRING u_rpm_s							=	"rpm/s" 					"revolutions per minute per second"
PRIVATE CONST STRING u_Nm_rad 						= "N·m/rad" 				"Newton metre per radian"
PRIVATE CONST STRING u_Nms_rad 						= "N·m·s/rad" 				"Newton metre second per radian"
PRIVATE CONST STRING u_rad_m 						= "rad/m" 					"radian per metre"

--Eléctricas
PRIVATE CONST STRING u_A_m							=	"A/m" 				"Ampere per metre"
PRIVATE CONST STRING u_A_m2							=	"A/m^2"				"Ampere per square metre"
PRIVATE CONST STRING u_Am2							=	"A·m^2" 				"Ampere - square metre"
PRIVATE CONST STRING u_A_skmol						=	"A/(kmol·s)"			"Ampere per kilomole - second"
PRIVATE CONST STRING u_A_smol						=	"A/(s·mol)"				"Ampere per mole - second"
PRIVATE CONST STRING u_C_m2							=	"C/m^2" 				"Coulombs per square metre"
PRIVATE CONST STRING u_C_m3							=	"C/m^3"				"Coulombs per cubic metre"
PRIVATE CONST STRING u_C_mm3							=	"C/mm^3" 			"Coulombs per cubic milimetre"

PRIVATE CONST STRING u_mV								=	"mV" 					"milivolts"
PRIVATE CONST STRING u_1_V							=	"1/V"					"Inverse volts"
PRIVATE CONST STRING u_V_m							=	"V/m"					"Volts per metre"
PRIVATE CONST STRING u_V_K 							= 	"V/K"					"Volts per Kelvin"

PRIVATE CONST STRING u_Farad_m							=	"F/m" 				"Farad per metre"

PRIVATE CONST STRING u_H_m							=	"H/m" 				"Henry per metre"

PRIVATE CONST STRING u_cd_m2							=	"cd/m2" 				"candela per square metre"

PRIVATE CONST STRING u_S_m							=	"S/m" 				"Siemens per metre"

PRIVATE CONST STRING u_mWb							=	"mWb" 				"miliweber"

PRIVATE CONST STRING u_Ohm_m							= "Ohm/m"				"Ohm per meter"
PRIVATE CONST STRING u_Ohmm 							= "Ohm·m"				"Ohm meter"

PRIVATE CONST STRING u_bits 							= "bits"					"bits"

PRIVATE CONST STRING u_AU								= "AU" 				"Astronomical Unit"

PRIVATE CONST STRING u_AU_day2							= "AU/day^2" 			"Astronomical Units per squared day"

PRIVATE CONST STRING u_AU3_day2							= "AU^3/day^2" 			"Cubic Astronomical Units per squared day"

PRIVATE CONST STRING u_knots							= "knots" 				"knots"


