-- Generated automatically by - EcosimPro - 6.2.0 


COMPONENT Wall_NA_3  

PORTS



TOPOLOGY
   Wall_NSP layer1(
	      nodes = 5,
         rho = 1000,
         cp = 500,
         k = 0.001,  -- Non default value.
         A = 1,
         e = 0.001,  -- Non default value.
         To = 290)
   Wall_NSP layer2(
			nodes = 5,
         rho = 1000,
         cp = 500,
         k = 0.02,  -- Non default value.
         A = 1,
         e = 0.02,  -- Non default value.
         To = 290)
   Wall_NSP layer3(
         nodes = 5,
         rho = 1000,
         cp = 500,
         k = 0.001,  -- Non default value.
         A = 1,
         e = 0.001,  -- Non default value.
         To = 290)

   CONNECT layer1.tp_out TO layer2.tp_in
   CONNECT layer2.tp_out TO layer3.tp_in


END COMPONENT