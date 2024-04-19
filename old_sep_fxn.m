
% function cost = cost_separation_system(P_flowrates, F_steam, R_ethane)
% 	global MOLMASS_METHANE MOLMASS_HYDROGEN MOLMASS_ETHANE MOLMASS_ETHYLENE ...
% 		 MOLMASS_PROPANE MOLMASS_BUTANE YR_PER_SEC
% 	global T_SEPARATION R PRESS_RXTR R ...
% 	 MAX_OPEX MAX_TFCI MAX_CAPEX G_PER_KT MOLMASS_WATER 

% 	% Product flow rate indicies 
% 	HYDROGEN = 1;
% 	METHANE = 2;
% 	ETHYLENE = 3;
% 	PROPANE = 4;
% 	BUTANE = 5;

% 	% Feed flow rate index
% 	ETHANE = 6;
 
% 	% SEPARATION_EFFICIENCY_FACTOR = 30;
% 	T = T_SEPARATION; % [ K ]

% 	%Using compositions from ASPEN
% 	%Component mole flow rate out of rxtr over total mole flow rate out of reactor
% 	% Mol fractions out of the reactoor

% 	% (mol / s) = (kt / yr) * (g / kt) * (mol / g) * (yr / s)
% 	P_flowrates(METHANE) = P_flowrates(METHANE) * G_PER_KT * (1/MOLMASS_METHANE) * YR_PER_SEC;
% 	P_flowrates(HYDROGEN) = P_flowrates(HYDROGEN) * G_PER_KT * (1/MOLMASS_HYDROGEN) * YR_PER_SEC;
% 	R_ethane = R_ethane * G_PER_KT * (1/MOLMASS_ETHANE) * YR_PER_SEC;
% 	P_flowrates(ETHYLENE) = P_flowrates(ETHYLENE) * G_PER_KT * (1/MOLMASS_ETHYLENE) * YR_PER_SEC;
% 	P_flowrates(PROPANE) = P_flowrates(PROPANE) * G_PER_KT * (1/MOLMASS_PROPANE) * YR_PER_SEC;
% 	P_flowrates(BUTANE) = P_flowrates(BUTANE) * G_PER_KT * (1/MOLMASS_BUTANE) * YR_PER_SEC; % Add this line for butane
% 	F_steam = F_steam * G_PER_KT * (1/MOLMASS_WATER) * YR_PER_SEC;

% 	%CONVERT TO_MOLES

% 	P_tot = sum(P_flowrates(HYDROGEN:BUTANE)) + F_steam + R_ethane;

% 	z_methane = P_flowrates(METHANE) / P_tot;
% 	z_hydrogen = P_flowrates(HYDROGEN) / P_tot;
% 	z_ethane = R_ethane / P_tot;
% 	z_ethylene = P_flowrates(ETHYLENE) / P_tot;
% 	z_propane = P_flowrates(PROPANE) / P_tot;
% 	z_butane = P_flowrates(BUTANE) / P_tot;
% 	z_water = F_steam / P_tot;
	
% 	%Mol fractions leaving each separation system (refer to Isa's drawing in GN)
% 	% leaving sep 1
% 	x_water = 1;

% 	% leaving sep 4
% 	x_ethane = 1;
% 	x_ethylene = 1;

% 	% leaving sep 2
% 	x_butane = 0.9995;
% 	x_propane = 1 - x_butane;

% 	% leaving sep 5 (PSA)
% 	x_methane = 3.01809372499680e-004;
% 	x_hydrogen = 1;
% 		% ?? How should I implement the PSA toggle switch on this
	
% 	%Pressures of PSA system [bar]
% 	P_in = PRESS_RXTR;
% 	P_H2 = 10;				% [ bar ]
% 	P_ME = 1;				% [ bar ]
% 		% These outlet pressures are constant for PSA system. DONT change 
	
% 	%Using flow rates from ASPEN [NOTE: FOR MATLAB USE THE VALUES FROM THE
% 	%SOLN_TABLE. WE USED THESE AS EXPECTED COSTS)

% 	% Flowrates of each exiting stream from the sep system	
% 	F_water = F_steam; 									% mol/s
% 	F_LPG = P_flowrates(BUTANE) + P_flowrates(PROPANE); % (mol / s)
% 	F_ethylene = P_flowrates(ETHYLENE);					% (mol / s)
% 	F_ethane = R_ethane;						% (mol / s)
% 	F_H2 = P_flowrates(HYDROGEN);						% (mol / s)
% 	F_ME = P_flowrates(METHANE);						% (mol / s); 

% 	%(J/s) =    (mol/s) * (J/mol K) * (T) 
% 	W_min_Sep_System = F_water*R*T*log(x_water/z_water) + ...
% 					F_LPG*R*T*log(x_propane/z_propane + ...
% 								  x_butane/z_butane) + ...
% 					F_ethylene*R*T*log(x_ethylene/z_ethylene) + ...
% 					F_ethane*R*T*log(x_ethane/z_ethane) + ...
% 					R*T*( ... 
% 						F_H2*log(P_H2/P_in)+ ...
% 						F_H2*log(x_hydrogen/z_hydrogen) +...
% 						F_ME*log(x_methane/z_methane) +...
% 						F_ME*log(P_ME/P_in)...
% 						);

% 	lamdba_min = 20;
% 	lambda_max = 50;	
% 	cost_energy = 3;		% ( $ / GJ )

% 	if MAX_OPEX
% 	%($/yr)             =   (J/s)     * (GJ/J) * (Work Efficiency) *($/GJ)* (s/yr)
% 		opex =  W_min_Sep_System*1e-9 * lambda_max * cost_energy * 30.24e6;
% 	else
% 		opex = W_min_Sep_System*1e-9 * lamdba_min * cost_energy * 30.24e6;
% 	end

% 	if MAX_CAPEX
% 	%($) 				 = ($/W)    (Efficiency) * (J/s) 
% 		capex = 1 * lambda_max * W_min_Sep_System;
% 	else
% 		capex = 0.5 * lamdba_min * W_min_Sep_System;
% 	end

% 	cost = 2.5 * capex ;
	
% end
