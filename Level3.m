clc; clear; close all;
global S1_MIN S1_MAX S1_POINTS;
global S2_MIN S2_MAX S2_POINTS;
global INVALID_FLOWRATE;
global Fethyl_S1S2_plotOpt;
global F_ETHYLENE;
global M_C2H6;
global MT_PER_KT G_PER_KT GJ_PER_KJ;
global VALUE_ETHANE VALUE_ETHYLENE VALUE_H2_CHEM;
global COST_RATES_STEAM;
global VALUE_H2_FUEL VALUE_CH4_FUEL VALUE_C3H6_FUEL VALUE_C4H8_FUEL;
global VALUE_NATGAS_FUEL VALUE_NUM2OIL_FUEL;
global COST_CO2 COST_WASTESTREAM;
global ENTHALPY_PROPANE ENTHALPY_BUTANE;
global MOLMASS_PROPANE MOLMASS_BUTANE;
global PROFIT_S1S2_OPT;
global HEAT_CAPACITY_ETHANE;
global HEAT_FORMATION_ETHANE;
global STEAM_30C STEAM_50C STEAM_100C STEAM_200C STEAM_500C STEAM_750C;

% DESIGN PARAMETERS_____________________________________________________________
STEAM_TO_FEED_RATIO = 0.6; %0.6 to 1.0

% Note: The primary units of this script are ... 
% Mass			kta
% Energy		GJ
% Pressure 		Bar
% Temperature 	Celcius

% CONSTANTS________________________________________________________________

% Plotting 
S1_MIN = 0.05;
S1_MAX = 0.95;
S1_POINTS = 40;
S2_MIN = 0.05;
S2_MAX = 0.95;
S2_POINTS = 40;
INVALID_FLOWRATE = 0;
Fethyl_S1S2_plotOpt = { ...
	'S_1 Selectivity', ...
	'S_2 Selectivity', ...
	'Ethylene Flowrate [kta]',...
	'P_ethylene_VS_S1_S2.jpg'};
PROFIT_S1S2_OPT = { ...
	'S_1 Selectivity', ...
	'S_2 Selectivity', ...
	'Annual Profit [$ MM USD]',...
	'P_ethylene_VS_S1_S2.jpg'}; 

% Design params
F_ETHYLENE = 200;			% [	kta ]

% Molar mass
M_C2H6 = 28.05;			% [ g / mol ]

% Unit conversions 
MT_PER_KT = 1000;		% [ kt / MT ]
G_PER_KT = 10^9;		% [ g / kt ]
GJ_PER_KJ = 10^-6;		% [ GJ / kJ ]
KG_PER_KT = 10^3;		% [ kg / MT ]
KJ_PER_GJ = 10^-6;		% [ kJ / GJ ]
MT_PER_G = 10^6;		% [ MT / g ]

% Economic | Chemicals
VALUE_ETHANE = 200;		% [ $ / MT ]
VALUE_ETHYLENE = 900;	% [ $ / MT ]
VALUE_H2_CHEM = 1400;	% [ $ / MT ]

% Economic | Steam 
% psia Temp[C] $/kg kJ/kg
COST_RATES_STEAM = [
    30  121		2.38  2213;
    50  138		3.17  2159;
    100 165		4.25  2067;
    200 194		5.32  1960;
    500 242		6.74  1755;
    750 266		7.37  1634
];

% Economic | Fuel
VALUE_H2_FUEL = 3;			% [ $ / GJ ]
VALUE_CH4_FUEL = 3;			% [ $ / GJ ]
VALUE_C3H6_FUEL = 3;		% [ $ / GJ ]
VALUE_C4H8_FUEL = 3;		% [ $ / GJ ]
VALUE_NATGAS_FUEL = 3;		% [ $ / GJ ]
VALUE_NUM2OIL_FUEL = 4.5;	% [ $ / US Gallon ]


% Thermodynamics | Heats of Formation (at 25C)
HEAT_FORMATION_ETHANE = -83.8;			% [ kJ / mol K ] reference Temp = std
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C74840&Units=SI&Mask=1EFF
HEAT_FORMATION_METHANE = -74.87;		% [ kJ / mol K ] reference Temp = std
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C74828&Mask=1
HEAT_FORMATION_ETHYLENE = 52.47;		% [ kJ / mol K ] reference Temp = std 
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C74851&Mask=1
HEAT_FORMATION_HYDROGEN = 0; 			% [ kJ / mol K ] reference Temp = std 
HEAT_FORMATION_PROPANE = -104.7;		% [ kJ / mol K ] reference Temp = std 
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C74986&Mask=1
HEAT_FORMATION_BUTANE = -125.6;			% [ kJ / mol K ] reference Temp = std 
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C106978&Mask=1

% Thermodynamics | Enthalpy of Reactions
ENTHALPY_RXN_1 = HEAT_FORMATION_HYDROGEN + HEAT_FORMATION_ETHYLENE ...
										- HEAT_FORMATION_ETHANE;
ENTHALPY_RXN_2 = HEAT_FORMATION_METHANE + HEAT_FORMATION_PROPANE ...
										- 2 * HEAT_FORMATION_ETHANE; 
ENTHALPY_RXN_3 = HEAT_FORMATION_ETHANE - HEAT_FORMATION_ETHANE ...
										- HEAT_FORMATION_ETHYLENE;

% Thermodynamics | Enthalpy of combustion of gas at standard conditions
ENTHALPY_METHANE = 890;					% [ kJ / mol ]	
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C74828&Mask=1
ENTHALPY_PROPANE = 2219.2;				% [ kJ / mol ]
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C74986&Mask=1
ENTHALPY_BUTANE = 2877.5;				% [ kJ / mol ]
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C106978&Mask=1
HEAT_CAPACITY_ETHANE = 52.71 * 10^-3;	% [ kJ / mol K ] Reference Temp = 300K 
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C74840&Units=SI&Mask=1EFF

% Chemical | Combustion Stochiometery 
CO2_TO_METHANE_COMBUSTION_STOICH = 1;
CO2_TO_PROPANE_COMBUSTION_STOICH = 3;
CO2_TO_BUTANE_COMBUSTION_STOICH = 4;

% Economics | Enviormental
TAX_CO2_PER_MT = 125;				% [ $ / MT ]
TAX_CO2_PER_GJ_METHANE = KJ_PER_GJ * (1 / ENTHALPY_METHANE) * CO2_TO_METHANE_COMBUSTION_STOICH * MT_PER_G * TAX_CO2_PER_MT;
TAX_CO2_PER_GJ_PROPANE = KJ_PER_GJ * (1 / ENTHALPY_PROPANE) * CO2_TO_PROPANE_COMBUSTION_STOICH * MT_PER_G * TAX_CO2_PER_MT;
TAX_CO2_PER_GJ_BUTANE = KJ_PER_GJ * (1 / ENTHALPY_BUTANE) * CO2_TO_BUTANE_COMBUSTION_STOICH * MT_PER_G * TAX_CO2_PER_MT;

% Economics | Post-Tax Value of different fuel sources
EFFECTIVE_VALUE_METHANE_FUEL = VALUE_H2_FUEL - TAX_CO2_PER_GJ_METHANE;
EFFECTIVE_VALUE_PROPANE_FUEL = VALUE_C3H6_FUEL - TAX_CO2_PER_GJ_PROPANE;
EFFECTIVE_VALUE_BUTANE_FUEL = VALUE_C4H8_FUEL - TAX_CO2_PER_GJ_BUTANE;

% Chemical | Molar Mass
MOLMASS_PROPANE = 44.0956;				% [ g / mol ]
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C74986&Mask=1
MOLMASS_BUTANE = 58.1222;				% [ g / mol ]
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C106978&Mask=1
MOLMASS_ETHANE = 30.0690;				% [ g / mol ]
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C74840

% Chemical | Steam Choice indicies
STEAM_COST_ROW = 3;
STEAM_30C = 1;
STEAM_50C = 2;
STEAM_100C = 3;
STEAM_200C = 4;
STEAM_500C = 5;
STEAM_750C = 6;

% Flow rate Indicies | For the flowrates(i) array
HYDROGEN = 1; 
METHANE = 2;
ETHYLENE = 3;
PROPANE = 4;
BUTANE = 5;

% SYSTEM OF EQUARTIONS (EXTENT OF RXN)_____________________________________

A = @(s1, s2)...
	[s1-1	,s1		,s1+1;
     s2		,s2-1	,s2;
     1		,2		,1	];
b = [0;		0;		F_ETHYLENE];

% FUNCTIONS | FLOWRATE_____________________________________________________

P_HYDROGEN = @(xi_1)			xi_1;
P_METHANE = @(xi_2)				xi_2;
P_ETHYLENE = @(xi_1, xi_3)		xi_1 - xi_3;
P_PROPANE = @(xi_2)				xi_2;
P_BUTANE = @(xi_3)				xi_3;

F_ETHANE = @(xi_1, xi_2, xi_3)	xi_1 + 2 * xi_2 + xi_3;

% FUNCTIONS | VALIDATION___________________________________________________

flowrates_valid = @( flowrates ) all(flowrates >= 0);

% FUNCTIONS | ECONOMICS____________________________________________________

% Inputs are in [ kta ] outputs are in [ $ ]
value_ethane = @(P_ethane) P_ethane * MT_PER_KT * VALUE_ETHANE;
value_ethylene = @(P_ethylene) P_ethylene * MT_PER_KT * VALUE_ETHYLENE;
value_h2_chem = @(P_h2_chem) P_h2_chem * MT_PER_KT * VALUE_H2_CHEM;

cost_steam = @(F_steam, steam_rate) F_steam * KG_PER_KT * steam_rate;
% cost_feed = @(F_ethane) F_ethane * VALUE_ETHANE

% FUNCTIONS | THEROMODYNAMICS______________________________________________
% Input: [ kta ] Output: [ GJ ]
% kta * (g / kta) * (mol / g) * (kJ / mol K) * (kJ / mol) * (GJ / kJ)
heat_ethane = @(F_ethane, T0, Tf) F_ethane * G_PER_KT * (1 / MOLMASS_ETHANE) * HEAT_CAPACITY_ETHANE * GJ_PER_KJ * (Tf - T0);
% Input: [ kta ] 	Output: [ GJ ]
heat_rxn1 = @(xi_1) xi_1 * ENTHALPY_RXN_1;
heat_rxn2 = @(xi_2) xi_2 * ENTHALPY_RXN_2;
heat_rxn3 = @(xi_3) xi_3 * ENTHALPY_RXN_3; 
heat_rxn = @(xi) heat_rxn1(xi(1)) + heat_rxn2(xi(2)) + heat_rxn3(xi(3)); 

% SCRIPT___________________________________________________________________

% Iterates through each value of selectivities S1 and S2 to find the economic
% potential for different reaction conditions 
s1_domain = linspace(S1_MIN, S1_MAX, S1_POINTS);
s2_domain = linspace(S2_MIN, S2_MAX, S2_POINTS);
[s1_mesh, s2_mesh] = meshgrid(s1_domain, s2_domain);
ethylene_flowrates = (s1_mesh + s2_mesh) .* 0;
profit = (s1_mesh + s2_mesh) .* 0;	
T_reactor = 800;				% [ C ]
P_reactor = 3;					% [ Bar ]
T_ethane_feed = 25;				% [ C ]

i = 1;
for s1 = s1_domain
	for s2 = s2_domain
		
		% Solve for extents of reaction
		xi = A(s1, s2) \ b;

		% Calculate the flow rates of each species
		P_hydrogen = P_HYDROGEN(xi(1));
		P_methane = P_METHANE(xi(2));
		P_ethylene = P_ETHYLENE(xi(1), xi(3));
		P_propane = P_PROPANE(xi(2));
		P_butane = P_BUTANE(xi(3));	
		F_ethane = F_ETHANE(xi(1), xi(2), xi(3));
		flowrates = [ P_hydrogen, P_methane, P_ethylene, P_propane, P_butane ];
	
		if (flowrates_valid(flowrates))


			% Store all etylene polymer output in a DS so it can be plotted
			ethylene_flowrates(i) = P_ethylene;

% 			% Calculate the heat flux needed to keep reactor isothermal 
			heat_flux = 0;
			F_steam = STEAM_TO_FEED_RATIO * F_ethane;
			heat_flux = heat_flux + heat_ethane(P_ethylene, T_ethane_feed, T_reactor);
% 			heat_flux = heat_flux + heat_steam(F_steam, STEAM_50C, P_reactor, T_reactor); 
			heat_flux = heat_flux + heat_rxn(xi)
% 
% 			% Use the heat flux to calculate the fuel cost	
			combusted_fuel_flow_rates = fuel_combustion(heat_flux, flowrates)
			% combusted_fuel_flow_rates = flowrates * 0;
% 
			% Determine how much of the product streams were combusted to keep the reactor isothermal	
			% Assume: no hydrogen is combusted
			combusted_methane = combusted_fuel_flow_rates(METHANE);
			combusted_propane = combusted_fuel_flow_rates(PROPANE);
			combusted_butane = combusted_fuel_flow_rates(BUTANE);
% 
% 			% VALUE CREATED | Primary Products
			profit(i) = profit(i) + value_ethylene(P_ethylene);
			profit(i) = profit(i) + value_h2_chem(P_hydrogen); % Assume no H2 combusted

			% VALUE CREATED | Non-combusted fuels 
% 			profit(i) = profit(i) + value_methane(P_methane - combusted_methane);
% 			profit(i) = profit(i) + value_propane(P_propane - combusted_propane);
% 			profit(i) = profit(i) + value_butane(P_butane - combusted_butane);	
% % 			
			% COSTS INCURRED
% 			profit(i) = profit(i) - tax_C02(combusted_fuel_flowrates);
			
			profit(i) = profit(i) - cost_steam(F_steam, COST_RATES_STEAM(STEAM_COST_ROW,STEAM_50C));
			profit(i) = profit(i) - value_ethane(F_ethane);
% 			profit(i) = profit(i) - cost_natural_gas_fuel(heatflux, combusted_fuel_flow_rates);
% 			% Assume no #2 Fuel Oil is used
% 			F_waste = 0; % ??????????????????????????
% 			profit(i) = profit(i) - cost_waste_stream(F_steam, F_waste)

		else
			profit(i) = INVALID_FLOWRATE;
			ethylene_flowrates(i) = INVALID_FLOWRATE;
		end 
		i = i + 1;
	end 
end 

profit = profit ./ 10^6; % Convert to Millions of dollars 
plot_contour(s1_mesh, s2_mesh, ethylene_flowrates, Fethyl_S1S2_plotOpt);
plot_contour(s1_mesh, s2_mesh, profit, PROFIT_S1S2_OPT);
disp("Function completed running")

% HELPER FUNCTIONS | PLOTTING______________________________________________

function z = plot_contour(x, y, z, options)
	% Unpack options 
	x_label = options{1};
	y_label = options{2};
	plt_title = options{3};
	plt_saveName = options{4};

	hold on 
	figure
    [C, h] = contourf(x, y, z); % Create filled contours
    clabel(C, h, 'FontSize', 10, 'Color', 'k', 'LabelSpacing', 200); % Customize label properties
	xlabel(x_label);
	ylabel(y_label);
	title(plt_title);
	saveas(gcf, plt_saveName);
	hold off
end

% function value = value_LPG(P_propane, P_butane)
% 	% ASSUMPTION : VALUE OF LPG IS JUST SUM VALUE OF PROPANE + BUTANE BEING
% 	% COMBUSTED
% 	% kt * (g / kt) * (mol / g) * (kJ / mol) * (GJ / KJ) * ($ / GJ)
% 	global G_PER_KT MOLMASS_PROPANE ENTHALPY_PROPANE GJ_PER_KJ ...
% 		VALUE_C3H6_FUEL MOLMASS_BUTANE ENTHALPY_BUTANE VALUE_C4H8_FUEL;

% 	prop_val = P_propane * G_PER_KT * (1 / MOLMASS_PROPANE) * ...
% 							ENTHALPY_PROPANE * GJ_PER_KJ * VALUE_C3H6_FUEL;
% 	but_val = P_butane * G_PER_KT * (1 / MOLMASS_BUTANE) * ...
% 						ENTHALPY_BUTANE * GJ_PER_KJ * VALUE_C4H8_FUEL;
% 	value = prop_val + but_val;
% end

% function cost = kkam(flowrates)
% 	cost = 0;
% end

function combusted_fuel_flowrates = fuel_combustion(heat_flux, flowrates)
	% Longest Chain Hydrocarbons are cheapest to combust
	remaining_heat_flux = 0;


	combusted_fuel_flowrates = flowrates * 0;

end

function heat = heat_steam(F_steam, STEAM_50C, P_reactor, T_reactor)

	heat = 0;

end















