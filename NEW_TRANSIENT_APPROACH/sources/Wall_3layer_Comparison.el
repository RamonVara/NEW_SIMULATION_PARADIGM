-- Generated automatically by - EcosimPro - 6.2.0 


COMPONENT Wall_3layer_Comparison 

PORTS
TOPOLOGY
   Wall_NSP layer1_nsp(
	      nodes = 5,
         rho = 1000,
         cp = 500,
         k = 0.001,  -- Non default value.
         A = 1,
         e = 0.001,  -- Non default value.
         To = 290)
   Wall_NSP layer2_nsp(
			nodes = 5,
         rho = 1000,
         cp = 500,
         k = 0.02,  -- Non default value.
         A = 1,
         e = 0.02,  -- Non default value.
         To = 290)
   Wall_NSP layer3_nsp(
         nodes = 5,
         rho = 1000,
         cp = 500,
         k = 0.001,  -- Non default value.
         A = 1,
         e = 0.001,  -- Non default value.
         To = 290)
   CONNECT layer1_nsp.tp_out TO layer2_nsp.tp_in
   CONNECT layer2_nsp.tp_out TO layer3_nsp.tp_in
	 
   Wall_Classic (nodes=5) layer1_classic(
         rho = 1000,
         cp = 500,
         k = 0.001,  -- Non default value.
         A = 1,
         e = 0.001,  -- Non default value.
         To = 290)
   Wall_Classic (nodes=5) layer2_classic(
         rho = 1000,
         cp = 500,
         k = 0.02,  -- Non default value.
         A = 1,
         e = 0.02,  -- Non default value.
         To = 290)
   Wall_Classic (nodes=5) layer3_classic(
         rho = 1000,
         cp = 500,
         k = 0.001,  -- Non default value.
         A = 1,
         e = 0.001,  -- Non default value.
         To = 290)
   CONNECT layer1_classic.tp_out TO layer2_classic.tp_in
   CONNECT layer2_classic.tp_out TO layer3_classic.tp_in
END COMPONENT