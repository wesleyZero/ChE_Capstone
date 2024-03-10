
% Clear the console
clc; 
% Close all the windows
close all;
% Clear Workspace Variables
clear;

global S1_MIN S1_MAX S1_POINTS;
global S2_MIN S2_MAX S2_POINTS;
global INVALID_FLOWRATE;
global Fethyl_S1S2_plotOpt;
global MT_PER_KT G_PER_KT GJ_PER_KJ;
global VALUE_ETHANE VALUE_ETHYLENE VALUE_HYDROGEN_CHEM;
global COST_RATES_STEAM;
global VALUE_HYDROGEN_FUEL VALUE_METHANE_FUEL VALUE_PROPANE_FUEL VALUE_BUTANE_FUEL;
global VALUE_NATGAS_FUEL VALUE_NUM2OIL_FUEL;
global ENTHALPY_PROPANE ENTHALPY_BUTANE;
global MOLMASS_PROPANE MOLMASS_BUTANE;
global PROFIT_S1S2_OPT;
global HEAT_CAPACITY_ETHANE;
global HEAT_FORMATION_ETHANE;
global STEAM_30PSIA STEAM_50PSIA STEAM_100PSIA STEAM_200PSIA STEAM_500PSIA STEAM_750PSIA;
global HYDROGEN METHANE ETHYLENE PROPANE BUTANE;
global ENTHALPY_METHANE ENTHALPY_PROPANE ENTHALPY_BUTANE HEAT_CAPACITY_ETHANE;
global KT_PER_G KG_PER_KT KJ_PER_GJ MT_PER_G ENTHALPY_NAT_GAS MOLMASS_ETHANE...
	MOLMASS_ETHYLENE MOLMASS_NATGAS;
global MT_CO2_PER_KT_METHANE MT_CO2_PER_KT_PROPANE MT_CO2_PER_KT_BUTANE ...
	MT_CO2_PER_KT_NATURALGAS;
global TAX_CO2_PER_MT;
global STEAM_PRESSURE_COL STEAM_TEMP_COL;
global MOLMASS_METHANE MOLMASS_WATER BAR_PER_PSIA;
global C_TO_K HEAT_CAPACITY_WATER;
global R k1_f k1_r k2 k3 R_2 C_TO_K YR_PER_SEC SEC_PER_YR MOLMASS_HYDROGEN
global PSA_TOGGLE ENTHALPY_HYDROGEN

% USER NOTES____________________________________________________________________

% Note: The primary (high level) units of this script are ... 
% Mass			kta
% Energy		GJ
% Pressure 		Bar
% Temperature 	Celcius
% Moles			Moles
% Value			Dollars

% [ __ ] THIS MEANS DIMENSIONLESS UNITS

% USER INPUTS | DESIGN PARAMETERS_______________________________________________

% Product
P_ETHYLENE_DES = 200;			% [	kta ]
	% Note! This design parameter's units are changed prior to the matrix def 

YEARS_IN_OPERATION = 10 ;

% USER INPUTS | GLOBAL CONSTANTS___________________________________________

% USER INPUTS | 3D PLOT, CONTOUR, LVL 2 & 3 CALCS_________________________

% Reactor Conditions | 3D PLOT & CONTOUR PLOT (S1 S2) && THE LVL3 CALCS
STEAM_TO_FEED_RATIO_MOLS = 0.6;		% [ __ ] 0.6 to 1.0
TEMP_RXTR = 825;				% [ C ] 
PRESS_RXTR = 2;					% [ Bar ]  
TEMP_ETHANE_FEED = 25;			% [ C ]
CONVERSION =  0.17053;			% [ __ ] % Level 2 & 3 Calculations 
USERINPUT_S1 = 0.96971 ;		% [ __ ] % Level 2 & 3 Calculations 
USERINPUT_S2 = 0.00011843; 		% [ __ ] % Level 2 & 3 Calculations 
STEAM_CHOICE = 1;
% 	STEAM_30PSIA = 1;
% 	STEAM_50PSIA = 2;
% 	STEAM_100PSIA = 3;
% 	STEAM_200PSIA = 4;
% 	STEAM_500PSIA = 5;
% 	STEAM_750PSIA = 6;
		% % Steam 
		% % [ psia Temp[C] $/MT  kJ/kg ]
		% COST_RATES_STEAM = [
		%     30  121		2.38  2213;
		%     50  138		3.17  2159;
		%     100 165		4.25  2067;
		%     200 194		5.32  1960;
		%     500 242		6.74  1755;
		%     750 266		7.37  1634
		% ];

% Plotting | 3D PLOT & CONTOUR PLOT (S1 S2)
NUM_POINTS = 10^4; 

% USER INPUTS | RXTR TABLE PARAMETERS_______________________________________

% Reactor Script Parameters | RXTR TABLE OUTPUT
V_MIN = 0.1;					% [ L ]
V_MAX = 10^4;					% [ L ]
NUM_V_POINTS = 20;				% [ __ ]

P_MIN = 2;						% [ Bar ]
P_MAX = 5;						% [ Bar ]
NUM_P_POINTS = 2;				% [ __ ]

T_MIN = 775;					% [ Celcius ]
T_MAX = 825;					% [ Celcius ]
NUM_T_POINTS = 2; 				% [ __ ]

STEAM_MIN = 0.6;				% [ __ ]
STEAM_MAX = 1.0;				% [ __ ]
NUM_STEAM_POINTS = 2;			% [ __ ]

% Table Overrides | RXTR TABLE OUTPUT
T_P_OVERRIDE = true;		
	T_OVERRIDE = 825;			%[C]
	P_OVERRIDE = 2;				%[Bar]
	STEAM_MR_OVERRIDE = 0.6;%	[__]

% Output fuel costs 
CONSOLE_OUTPUT_EFFECTIVE_VALUE_FUELS = true;

% Output the level 2 and 3 calculations 
OUTPUT_LVL3_FLOWRATES_TO_CONSOLE = true;
	SANITY_CHECK_CALCULATIONS = true;

% Plot the 3D and Contour plot's
CALCULATE_ALL_SELECTIVITIES = true;
	PLOT_ECON_3D = true;
	PLOT_ECON_COUNTOUR = true;

% Output the Reactor Design tables
CALCULATE_REACTOR_FLOWS = true;

% PSA Toggle switch
PSA_TOGGLE = true;

% Do you want to add the work of the compressor to the heat flux of heating 
% the steam from the temp it's avilable at, to the temp of the reactor?
ADD_COMPRESSOR_WORK_TO_STEAM_HEATFLUX = true;

% Zeolite and waste stream
% zeo 1.2 - 2.2 wt% absobtion = max of zeolite (g/g)

% NOTE SEARCH FOR "??" TO SEE MY ASSUMPTIONS AND OTHER NOTES IN THE CODE
 
% WORK OF THE COMPRESSOR HAS NOT BEEN IMPLEMENTED
% THE STEAM TO FEED RATIO LIKELY HAS UNIT ISSUES OF (g/g) vs (mol/mol)
% 		I think I implemented both

%_______________________________________________________________________________
% DON'T TOUCH ANYTHING BELOW THIS LINE
%_______________________________________________________________________________


% CONSTANTS | PLOTTING_____________________________________________________

CONSOLE_SECTION_DIVIDER = ...
	"_____________________________________________________________________";
S1_MIN = 0.01;
S1_MAX = 1.00;
S1_POINTS = NUM_POINTS ^ (1/2) ;
S2_MIN = 0.01;
S2_MAX = 1.00;
S2_POINTS = NUM_POINTS ^ (1/2);
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

% CONSTANTS | UNITS________________________________________________________

% Mass
MT_PER_KT = 10^3;		% [ MT / kt ]

G_PER_KT = 10^9;		% [ g / kt ]
KT_PER_G = 10^-9;		% [ kt / g ] 

KG_PER_KT = 10^6;		% [ kg / MT ]

MT_PER_G = 10^-6;		% [ MT / g ]

% Energy
GJ_PER_KJ = 10^-6;		% [ GJ / kJ ]
KJ_PER_GJ = 10^6;		% [ kJ / GJ ]

% Temperature		
C_TO_K = 273.15;		% [ C -> K ]
% Value
MMDOLLA_PER_DOLLA = 10^-6;	% [ $ MM / $]
DOLLA_PER_MMDOLLA = 10^6;	% [ $ / $ MM ]

% Pressure
BAR_PER_PSIA = 0.0689476;	% [ Bar / Psia ]

% Time
YR_PER_SEC = 1 / (3.154 * 10^7);	% [ yr / s ]
SEC_PER_YR = 3.154 * 10^7;			% [ s / yr ]

% Volumes 
M3_PER_L = 0.001;

% CONSTANTS | CHEMICAL_____________________________________________________

% Chemical | Molar Mass
MOLMASS_HYDROGEN = 2.01588;				% [ g / mol ]
	% Source :  https://webbook.nist.gov/cgi/cbook.cgi?ID=1333-74-0
MOLMASS_METHANE = 16.0425;				% [ g / mol ]
	% Source :  https://webbook.nist.gov/cgi/cbook.cgi?ID=74-82-8
MOLMASS_WATER = 18.015;						% [ g / mol ]
	% source : https://pubchem.ncbi.nlm.nih.gov/compound/Water
MOLMASS_CO2 = 44.01;					% [ g / mol ]
	% Source : https://pubchem.ncbi.nlm.nih.gov/compound/Carbon-dioxide-water
MOLMASS_PROPANE = 44.0956;				% [ g / mol ]
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C74986&Mask=1
MOLMASS_BUTANE = 58.1222;				% [ g / mol ]
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C106978&Mask=1
MOLMASS_ETHANE = 30.0690;				% [ g / mol ]
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C74840
MOLMASS_ETHYLENE = 28.0532;				% [ g / mol ]
	% Source = https://webbook.nist.gov/cgi/cbook.cgi?ID=74-85-1&Type=IR-SPEC&Index=QUANT-IR,20
MOLMASS_NATGAS = 16.04;					% [ g / mol ]
	% ASSUMING NATURAL GAS IS ALL METHANE

% Chemical | Combustion Stochiometery 
CO2_TO_METHANE_COMBUSTION_STOICH = 1;
CO2_TO_PROPANE_COMBUSTION_STOICH = 3;
CO2_TO_BUTANE_COMBUSTION_STOICH = 4;
C02_TO_NATGAS_COMBUSTION_STOICH = CO2_TO_METHANE_COMBUSTION_STOICH;
	% Natural gas is asuumed to be entirely methane

% CONSTANTS | THERMODYNAMICS_______________________________________________

% Gas Constant 
R = 8.314;								% [ J / mol K ]
R_2 = 0.0831446261815324;				% [ L bar / K mol ]

% Heat capacities 
HEAT_CAPACITY_WATER = 33.79 * 10^-3;	% [ kJ / mol K ] Ref Temp = 298K
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C14940637&Mask=1&Type=JANAFG&Table=on
HEAT_CAPACITY_ETHANE = 52.71 * 10^-3;	% [ kJ / mol K ] Reference Temp = 300K 
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C74840&Units=SI&Mask=1EFF

% Heats of Formation (at 25C)
HEAT_FORMATION_ETHANE = -83.8;			% [ kJ / mol ] reference Temp = std
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C74840&Units=SI&Mask=1EFF
HEAT_FORMATION_METHANE = -74.87;		% [ kJ / mol ] reference Temp = std
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C74828&Mask=1
HEAT_FORMATION_ETHYLENE = 52.47;		% [ kJ / mol ] reference Temp = std 
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C74851&Mask=1
HEAT_FORMATION_HYDROGEN = 0; 			% [ kJ / mol ] reference Temp = std 
HEAT_FORMATION_PROPANE = -104.7;		% [ kJ / mol ] reference Temp = std 
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C74986&Mask=1
HEAT_FORMATION_BUTANE = -125.6;			% [ kJ / mol ] reference Temp = std 
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C106978&Mask=1

% Enthalpy of combustion (std conditions)
ENTHALPY_HYDROGEN = 286;
	% Source : https://chem.libretexts.org/Courses/University_of_Kentucky/UK%3A_General_Chemistry/05%3A_Thermochemistry/5.3%3A_Enthalpy
ENTHALPY_METHANE = 890;					% [ kJ / mol ]	
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C74828&Mask=1
ENTHALPY_PROPANE = 2219.2;				% [ kJ / mol ]
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C74986&Mask=1
ENTHALPY_BUTANE = 2877.5;				% [ kJ / mol ]
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C106978&Mask=1
ENTHALPY_NAT_GAS = ENTHALPY_METHANE; 
	% Source : https://afdc.energy.gov/fuels/natural_gas_basics.html#:~:text=Natural%20gas%20is%20an%20odorless,used%20in%20the%20United%20States.
	% Natural gas is mostly methane, so assumed to be 100% methane in the calcs

% Enthalpy of Reactions [ kJ / extent rxn]
ENTHALPY_RXN_1 = HEAT_FORMATION_HYDROGEN + HEAT_FORMATION_ETHYLENE ...
										- HEAT_FORMATION_ETHANE;
ENTHALPY_RXN_2 = HEAT_FORMATION_METHANE + HEAT_FORMATION_PROPANE ...
										- 2 * HEAT_FORMATION_ETHANE; 
ENTHALPY_RXN_3 = HEAT_FORMATION_ETHANE - HEAT_FORMATION_ETHANE ...
										- HEAT_FORMATION_ETHYLENE;
% CONSTANTS | ECONOMICS____________________________________________________

%  Chemicals
VALUE_ETHANE = 200;				% [ $ / MT ]
VALUE_ETHYLENE = 900;			% [ $ / MT ]
VALUE_HYDROGEN_CHEM = 1400;		% [ $ / MT ]

% Steam 
% [ psia Temp[C] $/MT  kJ/kg ]
COST_RATES_STEAM = [
    30  121		2.38  2213;
    50  138		3.17  2159;
    100 165		4.25  2067;
    200 194		5.32  1960;
    500 242		6.74  1755;
    750 266		7.37  1634
];

% Accessing the Steam P,T Data 
	STEAM_PRESSURE_COL = 2;
	STEAM_TEMP_COL = 1;
	STEAM_COST_COL = 3;
	STEAM_30PSIA = 1;
	STEAM_50PSIA = 2;
	STEAM_100PSIA = 3;
	STEAM_200PSIA = 4;
	STEAM_500PSIA = 5;
	STEAM_750PSIA = 6;

% Economic | Fuel
VALUE_HYDROGEN_FUEL = 3;			% [ $ / GJ ]
VALUE_METHANE_FUEL = 3;				% [ $ / GJ ]
VALUE_PROPANE_FUEL = 3;				% [ $ / GJ ]
VALUE_BUTANE_FUEL = 3;				% [ $ / GJ ]
VALUE_NATGAS_FUEL = 3;				% [ $ / GJ ]
VALUE_NUM2OIL_FUEL = 4.5;			% [ $ / US Gallon ]

% Economics | Enviormental
TAX_CO2_PER_MT = 125;				% [ $ / MT ]

% [$ / GJ]	= 1GJ(basis) * (KJ / GJ)   * (mol gas / KJ) *           (mol CO2 / mol gas)          *  (g / mol C02)*(MT / g) * ($ / MT)
TAX_CO2_PER_GJ_METHANE = KJ_PER_GJ * (1 / ENTHALPY_METHANE) * CO2_TO_METHANE_COMBUSTION_STOICH * MOLMASS_CO2 * MT_PER_G * TAX_CO2_PER_MT;
TAX_CO2_PER_GJ_PROPANE = KJ_PER_GJ * (1 / ENTHALPY_PROPANE) * CO2_TO_PROPANE_COMBUSTION_STOICH * MOLMASS_CO2 * MT_PER_G * TAX_CO2_PER_MT;
TAX_CO2_PER_GJ_BUTANE = KJ_PER_GJ * (1 / ENTHALPY_BUTANE) * CO2_TO_BUTANE_COMBUSTION_STOICH * MOLMASS_CO2 * MT_PER_G * TAX_CO2_PER_MT;
TAX_CO2_PER_GJ_NATGAS = TAX_CO2_PER_GJ_METHANE; %

% Chemistry | MT of C02 per KT of Fuel used 
% (MT CO2) = 1KT(basis) * (g / KT) * (mol gas/ g gas) * 
MT_CO2_PER_KT_METHANE = G_PER_KT * (1/MOLMASS_METHANE) *...
	...	% (mol CO2 / mol gas) * 		(g CO2 / mol CO2) * (MT / g) 
	CO2_TO_METHANE_COMBUSTION_STOICH * MOLMASS_CO2 * MT_PER_G;
MT_CO2_PER_KT_PROPANE = G_PER_KT * (1/MOLMASS_PROPANE) *...
	CO2_TO_PROPANE_COMBUSTION_STOICH * MOLMASS_CO2 * MT_PER_G;
MT_CO2_PER_KT_BUTANE = G_PER_KT * (1/MOLMASS_BUTANE) *...
	CO2_TO_BUTANE_COMBUSTION_STOICH * MOLMASS_CO2 * MT_PER_G;
MT_CO2_PER_KT_NATURALGAS = MT_CO2_PER_KT_METHANE;

% FUNCTIONS | FLOWRATE_____________________________________________________

P_ETHYLENE = P_ETHYLENE_DES;
P_ETHYLENE_DES = P_ETHYLENE_DES * (1 / MOLMASS_ETHYLENE);	
P_PROPANE = @(s1, s2) 		(s2 / s1 *P_ETHYLENE_DES) * ...
										MOLMASS_PROPANE;
P_BUTANE = @(s1, s2)		(P_ETHYLENE_DES*(1/(2*s1) - s2/s1 - 1/2)) * ...
										MOLMASS_BUTANE;
F_ETHANE = @(s1, s2)		(P_ETHYLENE_DES / s1) * ...
										MOLMASS_ETHANE;
P_METHANE = @(s1, s2) (s2 / s1 * P_ETHYLENE_DES) * ...
										MOLMASS_METHANE;
P_HYDROGEN = @(s1, s2) (P_ETHYLENE_DES * ((1/(2*s1) - s2/s1 + 1/2))) * ...
										MOLMASS_HYDROGEN;

% FUNCTIONS | EXTENT OF REACTION___________________________________________

% Returns molar flowrates [ mol / yr ]
get_xi = @(flowrates) [ flowrates(HYDROGEN) * G_PER_KT / MOLMASS_HYDROGEN, ...
					flowrates(PROPANE) * G_PER_KT / MOLMASS_PROPANE, ...
					flowrates(BUTANE) * G_PER_KT / MOLMASS_BUTANE ];

% FUNCTIONS | VALIDATION___________________________________________________

flowrates_valid = @( flowrates ) all(flowrates >= 0);

% FUNCTIONS | ECONOMICS____________________________________________________

% ($ / yr) =                  (kta) *   (MT / KT) * ($ / MT)
value_ethane = @(P_ethane) P_ethane * MT_PER_KT * VALUE_ETHANE;
value_ethylene = @(P_ethylene) P_ethylene * MT_PER_KT * VALUE_ETHYLENE;
value_h2_chem = @(P_h2_chem) P_h2_chem * MT_PER_KT * VALUE_HYDROGEN_CHEM;
value_methane = @(P_methane) P_methane * MT_PER_KT * VALUE_METHANE_FUEL;
value_propane = @(P_propane) P_propane * MT_PER_KT * VALUE_PROPANE_FUEL;
value_butane = @(P_butane) P_butane * MT_PER_KT * VALUE_BUTANE_FUEL;

% ($ / yr) = 							(kta) *  (MT / kt) * ($ / MT)
cost_steam = @(F_steam, steam_rate) F_steam * MT_PER_KT * steam_rate;

% FUNCTIONS | THEROMODYNAMICS______________________________________________
% (GJ / yr) =                          (kta) *  (g / KT)   * (mol gas/ g gas)    * (kJ / mol K)          * (GJ / KJ) * (K)
heat_ethane = @(F_ethane, T0, Tf) F_ethane * G_PER_KT * (1 / MOLMASS_ETHANE) * HEAT_CAPACITY_ETHANE * GJ_PER_KJ * (Tf - T0);

% (GJ / yr)  =      (mol / yr) * (kJ / mol) *   (GJ / kJ)
heat_rxn1 = @(xi_1) xi_1 * ENTHALPY_RXN_1 * GJ_PER_KJ;
heat_rxn2 = @(xi_2) xi_2 * ENTHALPY_RXN_2 * GJ_PER_KJ;
heat_rxn3 = @(xi_3) xi_3 * ENTHALPY_RXN_3 * GJ_PER_KJ;
heat_rxn = @(xi) heat_rxn1(xi(1)) + heat_rxn2(xi(2)) + heat_rxn3(xi(3)); 

% FUNCTIONS | RATE CONTANTS________________________________________________

% T is [ Kelvin ]  R is [ J / mol K ]
k1_f = @(T) (4.652 * 10^13) * exp( (-273000 / (R * (T ))));
k1_r = @(T) (9.91 * 10^8) * exp( (-137800 / (R * (T ))));
k2 = @(T) (4.652 * 10^11) * exp( (-273000 / (R * (T ))));
k3 = @(T) (7.083 * 10^13) * exp( (-252600 / (R * (T ))));


% DESIGN PARAMS____________________________________________________________
STEAM_TO_FEED_RATIO_MASS = (MOLMASS_WATER / MOLMASS_ETHANE) * STEAM_TO_FEED_RATIO_MOLS;


% SCRIPT___________________________________________________________________

% Economics | Post-Tax Value of different fuel sources 
if (CONSOLE_OUTPUT_EFFECTIVE_VALUE_FUELS)
	disp(" [ $ / GJ ] ")
	EFFECTIVE_VALUE_HYDROGEN_FUEL = VALUE_HYDROGEN_FUEL
	EFFECTIVE_VALUE_METHANE_FUEL = VALUE_METHANE_FUEL + TAX_CO2_PER_GJ_METHANE
	EFFECTIVE_VALUE_PROPANE_FUEL = VALUE_PROPANE_FUEL + TAX_CO2_PER_GJ_PROPANE
	EFFECTIVE_VALUE_BUTANE_FUEL = VALUE_BUTANE_FUEL + TAX_CO2_PER_GJ_BUTANE
	EFFECTIVE_VALUE_NAT_GAS_FUEL = VALUE_NATGAS_FUEL + TAX_CO2_PER_GJ_NATGAS
% 	EFFECTIVE_VALUE_NUM2_FUEL = VALUE_NATGAS_FUEL + TAX_CO2_PER_GJ_NUM2;

end

if (OUTPUT_LVL3_FLOWRATES_TO_CONSOLE)

	% Calculate the flow rates of each species (kta)
	P_hydrogen = P_HYDROGEN(USERINPUT_S1, USERINPUT_S2);
	P_methane = P_METHANE(USERINPUT_S1, USERINPUT_S2);
	P_ethylene = P_ETHYLENE;
	P_propane = P_PROPANE(USERINPUT_S1, USERINPUT_S2);
	P_butane = P_BUTANE(USERINPUT_S1, USERINPUT_S2);
	F_ethane = F_ETHANE(USERINPUT_S1, USERINPUT_S2);
	P_flowrates = [ P_hydrogen, P_methane, P_ethylene, P_propane, P_butane ];	

	disp(CONSOLE_SECTION_DIVIDER)
	if (flowrates_valid(P_flowrates))
		
		fprintf("Flowrates for the reactor given that s1 = %f, s2 = %f conv = %f\n\n", ...
			USERINPUT_S1, USERINPUT_S2, CONVERSION)

		disp(CONSOLE_SECTION_DIVIDER)
		disp("Level 2 Flowrates  in / out of the entire plant [ kt / yr ]")
		P_hydrogen 
		P_methane 
		P_ethylene
		P_propane 
		P_butane 

		disp("Fresh Feed Flowrate")
		F_ethane
		
		disp(CONSOLE_SECTION_DIVIDER)
		disp("Level 3 Flowrates [ kt / yr ] ")

		disp("Recycle Stream Flowrate")

		R_ethane = F_ethane * ((1-CONVERSION) / (CONVERSION))
		% R_ethane = (P_ethylene/USERINPUT_S1) * ((1-CONVERSION)/CONVERSION)

		disp("Reactor Flowrates")

		F_ethane_into_reactor = R_ethane + F_ethane
		
		if SANITY_CHECK_CALCULATIONS
			disp(CONSOLE_SECTION_DIVIDER)
			disp("Sanity Checking the Calculations")
			Conservation_of_mass = F_ethane - sum(P_flowrates)
			if Conservation_of_mass
				fprintf("WARNING : YOU ARE NOT CONSERVING MASS\n\n")
			end
		end
	else
		disp("ERROR : Selectivities S1 S2 chosen are not physically possible")
	end
end

% SCRIPT | PLOTTING_____________________________________________________________

if (CALCULATE_ALL_SELECTIVITIES)
	disp(CONSOLE_SECTION_DIVIDER)
	disp("Calculating all selectivities... ")
	% Iterates through each value of selectivities S1 and S2 to find the economic
	% potential for different reaction conditions 
	s1_domain = linspace(S1_MIN, S1_MAX, S1_POINTS);
	s2_domain = linspace(S2_MIN, S2_MAX, S2_POINTS);
	[s1_mesh, s2_mesh] = meshgrid(s1_domain, s2_domain);
	% All flowrates are initialized as matricies of zeros
	ethylene_flowrates = (s1_mesh + s2_mesh) .* 0;
	hydrogen_flowrates = (s1_mesh + s2_mesh) .* 0;
	methane_flowrates = (s1_mesh + s2_mesh) .* 0;
	ethylene_flowrates = (s1_mesh + s2_mesh) .* 0;
	propane_flowrates = (s1_mesh + s2_mesh) .* 0;
	butane_flowrates = (s1_mesh + s2_mesh) .* 0;
	ethane_flowrates = (s1_mesh + s2_mesh) .* 0;
	
	profit = (s1_mesh + s2_mesh) .* 0;

	% Flow rate Indicies | For the flowrates(i) array
	HYDROGEN = 1; 
	METHANE = 2;
	ETHYLENE = 3;
	PROPANE = 4;
	BUTANE = 5;

	i = 1;
	for s1 = s1_domain
		for s2 = s2_domain
			
			P_hydrogen = P_HYDROGEN(s1, s2);
			P_methane = P_METHANE(s1, s2);
			P_ethylene = P_ETHYLENE;
			P_propane = P_PROPANE(s1, s2);
			P_butane = P_BUTANE(s1, s2);
			F_ethane = F_ETHANE(s1, s2);

			P_flowrates = [ P_hydrogen, P_methane, P_ethylene, P_propane, P_butane ];
		
			if (flowrates_valid(P_flowrates))

				% Store for plotting (kta)
				hydrogen_flowrates(i) = P_HYDROGEN(s1, s2);
				methane_flowrates(i) = P_METHANE(s1, s2);
				ethylene_flowrates(i) = P_ETHYLENE;
				propane_flowrates(i) = P_PROPANE(s1, s2);
				butane_flowrates(i) = P_BUTANE(s1, s2);
				ethane_flowrates(i) = F_ETHANE(s1, s2);
				
				xi = [];
				% Calculate the heat flux needed to keep reactor isothermal 
				heat_flux = 0;
				xi = get_xi(P_flowrates);
				F_steam = STEAM_TO_FEED_RATIO_MASS * F_ethane;
				heat_flux = heat_flux + heat_ethane(F_ethane, TEMP_ETHANE_FEED, TEMP_RXTR);
				heat_flux = heat_flux + heat_steam(F_steam, STEAM_CHOICE, PRESS_RXTR, TEMP_RXTR) ;
				heat_flux = heat_flux + heat_rxn(xi);

	% 			% Use the heat flux to calculate the fuel cost	
				[combusted_fuel_flow_rates, heat_flux_remaining] = fuel_combustion(heat_flux, P_flowrates);

				% Calculate how much natural gas you needed to combust
				F_natural_gas = natgas_combustion(heat_flux_remaining);

				% Determine how much of the product streams were combusted to keep the reactor isothermal	

				combusted_hydrogen = combusted_fuel_flow_rates(HYDROGEN);
				combusted_methane = combusted_fuel_flow_rates(METHANE);
				combusted_propane = combusted_fuel_flow_rates(PROPANE);
				combusted_butane = combusted_fuel_flow_rates(BUTANE);

	% 			% VALUE CREATED | Primary Products
				profit(i) = profit(i) + value_ethylene(P_ethylene);
				profit(i) = profit(i) + value_h2_chem(P_hydrogen - combusted_hydrogen);

				% VALUE CREATED | Non-combusted fuels 
				% profit(i) = profit(i) + value_methane(P_methane - combusted_methane);
					% ?? I don't think you can sell methane. IH - need to
                    % determine energy requirements for compressors +
                    % separation + cooling (will likely need to purchase
                    % Nat Gas)
				profit(i) = profit(i) + value_propane(P_propane - combusted_propane);
				profit(i) = profit(i) + value_butane(P_butane - combusted_butane);	

				% COSTS INCURRED
				profit(i) = profit(i) - tax_C02(combusted_fuel_flow_rates, F_natural_gas);
 				profit(i) = profit(i) - cost_steam(F_steam, COST_RATES_STEAM(STEAM_CHOICE, STEAM_COST_COL));
				profit(i) = profit(i) - value_ethane(F_ethane);
				profit(i) = profit(i) - cost_natural_gas_fuel(F_natural_gas);
% 				profit(i) = profit(i) - cost_reactor(P_flowrates, P, T, S1, S2);	
	% 			profit(i) = profit(i) - cost_waste_stream(F_steam, F_waste)

			else
				profit(i) = INVALID_FLOWRATE;
				ethylene_flowrates(i) = INVALID_FLOWRATE;
			end 
			i = i + 1;
		end 
	end 

	profit = profit ./ 10^6; % Convert to Millions of dollars 
	profit(profit < 0) = 0; % remove irrelvant data
	
	if (PLOT_ECON_COUNTOUR)
		disp("Plotting EP Contour Map")
		plot_contour(s1_mesh, s2_mesh, profit, PROFIT_S1S2_OPT);
	end
	if (PLOT_ECON_3D)
		disp("Plotting 3D EP Surface Function")
		plot_3D(s1_mesh, s2_mesh, profit, PROFIT_S1S2_OPT);
	end

	% Prepare the array of flow rate matrices
	% flowRatesArray = {hydrogen_flowrates, methane_flowrates, ethylene_flowrates, propane_flowrates, butane_flowrates, ethane_flowrates};

	% Call the function with the desired row
	% plotFlowRatesForRow(4, flowRatesArray); % To plot the first row across all matrices
end


% SCRIPT | REACTOR _____________________________________________________________

T_RANGE = linspace(T_MIN, T_MAX, NUM_T_POINTS);
P_RANGE = linspace(P_MIN, P_MAX, NUM_P_POINTS);
STEAM_RANGE = linspace(STEAM_MIN, STEAM_MAX, NUM_STEAM_POINTS);
V_RANGE = [V_MIN, V_MAX]; % WARNING THESE ARE IN LITERS
% H2 Methane Ethane Propane Butane Ethylene 
F_INTIAL_COND = [ 0; 0; 0; 0; 0; 10]; % These are in kta

	% Product flow rate indicies 
	HYDROGEN = 1;
	METHANE = 2;
	ETHYLENE = 3;
	PROPANE = 4;
	BUTANE = 5;

	% Feed flow rate index
	ETHANE = 6;

if (CALCULATE_REACTOR_FLOWS)
	disp("Reactor Script ")
	for T_i = T_RANGE
		for P_i = P_RANGE
			for MR_S_i = STEAM_RANGE
	
				% override the T_i and P_i with user input 
				if T_P_OVERRIDE
					disp("WARNING: OVERRIDE HAS BEEN ACTIVATED")
					T_i = T_OVERRIDE;
					P_i = P_OVERRIDE;
					MR_S_i = STEAM_MR_OVERRIDE;
				end

				fprintf("\n\nT = %f [C], P = %f [bar] MR = %f [__]\n", T_i, P_i, MR_S_i)

				% Setup the PFR Design Equations 
				
				% BASIS CALCULATIONS_______________________________________________________________

				% CONVERT TO MOLES_________________________________________________________________
				% Convert all of the initial conditions to mol / s
				% (mol / s) =			(kt / yr) *               (g / kt) *   ( mol / g )        * ( yr / s)
				F_INTIAL_COND(METHANE) = F_INTIAL_COND(METHANE) * G_PER_KT * (1/MOLMASS_METHANE) * YR_PER_SEC;
				F_INTIAL_COND(HYDROGEN) = F_INTIAL_COND(HYDROGEN) * G_PER_KT * (1/MOLMASS_HYDROGEN) * YR_PER_SEC;
				F_INTIAL_COND(ETHANE) = F_INTIAL_COND(ETHANE) * G_PER_KT * (1/MOLMASS_ETHANE) * YR_PER_SEC;
				F_INTIAL_COND(ETHYLENE) = F_INTIAL_COND(ETHYLENE) * G_PER_KT * (1/MOLMASS_ETHYLENE) * YR_PER_SEC;
				F_INTIAL_COND(PROPANE) = F_INTIAL_COND(PROPANE) * G_PER_KT * (1/MOLMASS_PROPANE) * YR_PER_SEC;
				
				% Calculate the molar flow rate of the steam
				% mol/s = __    * mol / s
				F_steam = MR_S_i * F_INTIAL_COND(ETHANE);
				
				% Solve the system ODE's 
				%	(L, mol / s)           (L, mol/s, Celcius, Bar, mol/s)
				odes = @(V, F) reactionODEs(V, F, T_i, P_i, F_steam);
				[V_soln_ODE, F_soln_ODE] = ode45(odes, V_RANGE, F_INTIAL_COND); 
	
				% Calculate the conversion
				conversion = (F_INTIAL_COND(ETHANE) - F_soln_ODE(:, ETHANE)) / F_INTIAL_COND(ETHANE);

				% put handles length of the solution and the initial ethane flow
				len = length(F_soln_ODE(:, 1));
				F_ethane_initial = ones(len, 1) * F_INTIAL_COND(ETHANE);

				% Calculate the Selectivities, for each row (aka V_rxtr) 
				select_1 = (F_soln_ODE(:, ETHYLENE) ) ./ (F_ethane_initial - F_soln_ODE(:, ETHANE));
				select_2 = (F_soln_ODE(:, PROPANE) ) ./ (F_ethane_initial - F_soln_ODE(:, ETHANE));
				
				% Calculate the inlet volumetric flow rate 
				% (L / s) ???????????????
				P_sum = F_soln_ODE(:, HYDROGEN:BUTANE);
				% Turn these constants into vectors to operation is valid
				F_steam = ones(length(P_sum(:,1)), 1) .* F_steam;
				% put handles on terms, to make the code readable
				sum_flowrates_into_reactor = F_INTIAL_COND(ETHANE) + F_steam;
				% Calculate the flow rate into the reactor
				q0 = (R_2  * (T_i + C_TO_K) / P_i) .* sum_flowrates_into_reactor;
					% This is F.30 in the 'Design PFR Algorithm Appendix'
			
				% PLANT CALCULATIONS_______________________________________________________________

				% Calculate the the flowrates of the plant sized reactor given S1, S2 from ODE's 
				F_ethane = [];
				P_ethylene = [];
				for row = 1:length(select_1)
					% mol / s          = (kt / yr)   * (g / kt)  * (mol / g)            * (yr / s)
					P_ethylene(row, 1) = P_ETHYLENE .* G_PER_KT .* (1/MOLMASS_ETHYLENE) * YR_PER_SEC;
				end
				
				% Calculate the scaling factor of the plant, from the basis
				% mol / mol = ...
				scaling_factor = P_ethylene(:, 1) ./ F_soln_ODE(:, ETHYLENE);
				
				% Calculate the volume of the plant sized reactor
				% L / s = ( L / s )         * (   (mol / s )  ) / (    (mol / s)         )   
				%         BASIS             *     PLANT_FLOW    /      BASIS_FLOW
				V_plant = V_soln_ODE(:, 1) .* scaling_factor;

				% cost of the reactor
				cost_rxt_vec = zeros(size(V_plant));
				for row = 1:length(V_plant)	
					% ( $ )  
					cost_rxt_vec(row) = cost_reactor(V_plant(row,1) * M3_PER_L);
					cost_rxt_vec(row) = cost_rxt_vec(row) / YEARS_IN_OPERATION;
				end
				
				% inlet flow of the plant scaled reactor
				q0_plant = q0(:, 1) .* scaling_factor;
					% Eqn F.35 in 'Design PFR Algorithm Appendix' 

				% Scaling all of the molar flowrates to the size of the plant
				F_soln_ODE(:, METHANE) = F_soln_ODE(:, METHANE) .* scaling_factor;
				F_soln_ODE(:, HYDROGEN) = F_soln_ODE(:, HYDROGEN) .* scaling_factor;
				F_soln_ODE(:, ETHANE) = F_soln_ODE(:, ETHANE) .* scaling_factor;
				F_soln_ODE(:, ETHYLENE) = F_soln_ODE(:, ETHYLENE) .* scaling_factor;
				F_soln_ODE(:, BUTANE) = F_soln_ODE(:, BUTANE) .* scaling_factor;
				F_soln_ODE(:, PROPANE) = F_soln_ODE(:, PROPANE) .* scaling_factor;

				% CONVERT BACK TO MASS__________________________________________________________
				
				% convert back to kta
				% kt / yr =  mol / s    * g / mol         * kt / g   * s / yr
				F_soln_ODE(:, METHANE) = F_soln_ODE(: ,METHANE) * MOLMASS_METHANE * KT_PER_G * SEC_PER_YR;
				F_soln_ODE(:, ETHANE) = F_soln_ODE(:, ETHANE) * MOLMASS_ETHANE * KT_PER_G * SEC_PER_YR;	
				F_soln_ODE(:, HYDROGEN) = F_soln_ODE(:, HYDROGEN) * MOLMASS_HYDROGEN * KT_PER_G * SEC_PER_YR;
				F_soln_ODE(:, ETHYLENE) = F_soln_ODE(:, ETHYLENE) * MOLMASS_ETHYLENE * KT_PER_G * SEC_PER_YR;
				F_soln_ODE(:, BUTANE) = F_soln_ODE(:, BUTANE) * MOLMASS_BUTANE * KT_PER_G * SEC_PER_YR;
				F_soln_ODE(:, PROPANE) = F_soln_ODE(:, PROPANE) * MOLMASS_PROPANE * KT_PER_G * SEC_PER_YR;

								% Check if you're conserving mass
				conserv_mass = zeros(length(F_soln_ODE(:,1)), 1);

				% ECONOMIC CALCULATIONS____________________________________________________________
				profit = zeros(length(F_soln_ODE(:,1)), 1);
				for i = 1:length(F_soln_ODE(:, 1))
					
					% DEBUGGING
					if i > 500
						disp("")
					end 
					% P_flowrates = [ P_hydrogen, P_methane, P_ethylene, P_propane, P_butane ];
					P_flowrates = F_soln_ODE(i , HYDROGEN:BUTANE);
					
					P_hydrogen = P_flowrates(HYDROGEN);
					P_methane = P_flowrates(METHANE);
					P_ethylene = P_flowrates(ETHYLENE);
					P_propane = P_flowrates(PROPANE);
					P_butane = P_flowrates(BUTANE);

% 					F_ethane = P_flowrates(ETHANE);
					% F_ethane = sum(P_flowrates(:)) + F_soln_ODE(i ,ETHANE);
					% ?? WHAT IS THE FLOWRATE INTO THE PLANT ??
					% F_ethane = F_ETHANE()
					F_ethane = F_ETHANE(select_1(i), select_2(i));

					if (~flowrates_valid(P_flowrates))
						disp("WARNING SOME FLOWATES MAY BE INVALID")
					end

					% Calculate the heat flux needed to keep reactor isothermal 
					heat_flux = 0;
					xi = get_xi(P_flowrates);
					F_steam = STEAM_TO_FEED_RATIO_MASS * F_ethane;
					heat_flux = heat_flux + heat_ethane(F_ethane, TEMP_ETHANE_FEED, TEMP_RXTR);
					heat_flux = heat_flux + heat_steam(F_steam, STEAM_CHOICE, PRESS_RXTR, TEMP_RXTR) ;
					heat_flux = heat_flux + heat_rxn(xi);

		% 			% Use the heat flux to calculate the fuel cost	
					[combusted_fuel_flow_rates, heat_flux_remaining] = fuel_combustion(heat_flux, P_flowrates);

					% Calculate how much natural gas you needed to combust
					F_natural_gas = natgas_combustion(heat_flux_remaining);

					% Determine how much of the product streams were combusted to keep the reactor isothermal	

					combusted_hydrogen = combusted_fuel_flow_rates(HYDROGEN);
					combusted_methane = combusted_fuel_flow_rates(METHANE);
					combusted_propane = combusted_fuel_flow_rates(PROPANE);
					combusted_butane = combusted_fuel_flow_rates(BUTANE);

		% 			% VALUE CREATED | Primary Products
					profit(i, 1) = profit(i, 1) + value_ethylene(P_ethylene);
					profit(i, 1) = profit(i, 1) + value_h2_chem(P_hydrogen - combusted_hydrogen);
					
					% VALUE CREATED | Non-combusted fuels 
					% The commented line can be removed or modified as per the context.
					% profit(i, 1) = profit(i, 1) + value_methane(P_methane - combusted_methane);
					
					profit(i, 1) = profit(i, 1) + value_propane(P_propane - combusted_propane);
					profit(i, 1) = profit(i, 1) + value_butane(P_butane - combusted_butane);  
					
					% COSTS INCURRED
					profit(i, 1) = profit(i, 1) - tax_C02(combusted_fuel_flow_rates, F_natural_gas);
					profit(i, 1) = profit(i, 1) - cost_steam(F_steam, COST_RATES_STEAM(STEAM_CHOICE, STEAM_COST_COL));
					profit(i, 1) = profit(i, 1) - value_ethane(F_ethane);
					profit(i, 1) = profit(i, 1) - cost_natural_gas_fuel(F_natural_gas);
					% WASTE STREAM HERE 

					if profit(i, 1) > 2 * 10^8
						fprintf("Profit = %f \n", profit(i, 1))
						% VALUE CREATED | Primary Products
						disp("value")
						value_ethylene(P_ethylene)
						value_h2_chem(P_hydrogen - combusted_hydrogen)
						
						% VALUE CREATED | Non-combusted fuels 
						% The commented line can be removed or modified as per the context.
						% profit(i, 1) = profit(i, 1) + value_methane(P_methane - combusted_methane);
						
						value_propane(P_propane - combusted_propane)
						value_butane(P_butane - combusted_butane) 
						disp("costs")
						% COSTS INCURRED
						tax_C02(combusted_fuel_flow_rates, F_natural_gas)
						cost_steam(F_steam, COST_RATES_STEAM(STEAM_CHOICE, STEAM_COST_COL))
						value_ethane(F_ethane)
						cost_natural_gas_fuel(F_natural_gas)
					end 
					% Sanity checking
					conserv_mass(i, 1) = F_ethane - sum(P_flowrates);
				end

				% PLOTTING_________________________________________________________________________
				

				


				col_names = {'V_rxtr [L] ', 'Hydrogen [kta]', 'Methane', ...
					'Ethylene', 'Propane', 'Butane','Ethane', 'conversion', ...
					'S1', 'S2', 'q0 [ L /s ]', 'Vol_plant [ L ]', 'q0 plant', 'cost reactor', 'profit', 'net profit', 'conserv mass'};
				soln_table = table( V_soln_ODE, F_soln_ODE(:, HYDROGEN), ...
							F_soln_ODE(:, METHANE), F_soln_ODE(:, ETHYLENE), ...
							F_soln_ODE(:, PROPANE), F_soln_ODE(:, BUTANE), ...
							F_soln_ODE(:, ETHANE), conversion,select_1, ...
							select_2,q0,V_plant,q0_plant,cost_rxt_vec,profit, profit - cost_rxt_vec, conserv_mass,'VariableNames',col_names)
	% 			soln_table.Properties.VariableNames = col_names;
	
	
	
	
				
	
				% Computer Selectivity vs conversion relationships 
	
				% Use Selectivity vs Conversion Relationships with lvl 2 & 3 balances 
				% to calculate the true feed flow rates into the reactor 
	
			end 
		end 
	end
end 



disp("The Script is done running ️")
% HELPER FUNCTIONS | PLOTTING______________________________________________

function z = plot_contour(x, y, z, options)
	global PSA_TOGGLE
	% Unpack options 
	x_label = options{1};
	y_label = options{2};
	plt_title = options{3};
	plt_saveName = options{4};
	
	if PSA_TOGGLE
    	stringValue = 'true';
	else
    	stringValue = 'false';
	end
	plt_title = plt_title + sprintf(" PSA %s ", stringValue);

	hold on 
	figure
    [C, h] = contourf(x, y, z); % Create filled contours
    clabel(C, h,  'FontSize', 10, 'Color', 'k', 'LabelSpacing', 200); % Customize label properties
	xlabel(x_label);
	ylabel(y_label);
	title(plt_title);
	saveas(gcf, plt_saveName);
	hold off
end 

function plot_3D(x, y, z, options)
	global PSA_TOGGLE

	% Unpack options 
    x_label = options{1};
    y_label = options{2};
    plt_title = options{3};
    plt_saveName = options{4};

	if PSA_TOGGLE
    	stringValue = 'true';
	else
    	stringValue = 'false';
	end
	plt_title = plt_title + sprintf(" PSA %s ", stringValue);

     % Create a new figure
    hold on; % Hold on to add multiple plot elements
    figure
	surf(x, y, z); % Create a 3D surface plot
    
    % Customizing the plot
    xlabel(x_label);
    ylabel(y_label);
    zlabel('Z Value'); % Add a label for the z-axis
    title(plt_title);
    colorbar; % Adds a color bar to indicate the scale of z values
%     shading interp; % Option for smoother color transition on the surface
    
    hold off; % Release the figure
    saveas(gcf, plt_saveName); % Save the figure to file
end


function plotFlowRatesForRow(row, flowRatesArray)
    % flowRatesArray is expected to be an array of matrices, where each matrix corresponds to a species' flow rates

    % Names of the gases for labeling purposes
    gasNames = {'Hydrogen', 'Methane', 'Ethylene', 'Propane', 'Butane', 'Ethane'};
    
    % Create a figure
    figure;
    hold on; % Hold on to plot all data on the same figure
    
    % Loop through each flow rate matrix in the array
    for i = 1:length(flowRatesArray)
        % Extract the specified row from the current matrix
        currentRow = flowRatesArray{i}(row, :);
        
        % Plot the current row with a marker
        plot(currentRow, '-o', 'DisplayName', gasNames{i});
    end
    
    % Adding plot features
    title(sprintf('Flow Rates for Row %d', row));
    xlabel('Selectivity 1 (S2 fixed)');
    ylabel('Flow Rate');
    legend('show');
    hold off; % Release the figure for other plots
end

% HELPER FUNCTIONS | HEAT ______________________________________________

function [combusted_fuel_flowrates, heatflux_left] = fuel_combustion(heat_flux, flowrates)
	global HYDROGEN METHANE ETHYLENE PROPANE BUTANE;
	global ENTHALPY_METHANE ENTHALPY_PROPANE ENTHALPY_BUTANE HEAT_CAPACITY_ETHANE;
	global MT_PER_KT G_PER_KT GJ_PER_KJ KJ_PER_GJ MOLMASS_METHANE KT_PER_G MOLMASS_BUTANE ...
			MOLMASS_PROPANE PSA_TOGGLE ENTHALPY_HYDROGEN MOLMASS_HYDROGEN

	% Note! : Longest Chain Hydrocarbons are cheapest to combust

	% initialize all values in the array to be zero 
	combusted_fuel_flowrates = flowrates * 0;

	% LOGIC : Goes through each heat source in order, returns if the heat flux supplied is sufficient.
	heatflux_left = heat_flux; 

	% (GJ / yr)           = (kt / yr)          * (g / kt) * (kJ / g)        * (GJ / kJ)
	Q_combust_all_hydrogen = flowrates(HYDROGEN) * G_PER_KT * ENTHALPY_HYDROGEN * GJ_PER_KJ;

	if (~PSA_TOGGLE)
		% Hydrogen
		if (heatflux_left > Q_combust_all_hydrogen)
			combusted_fuel_flowrates(HYDROGEN) = flowrates(HYDROGEN);
			heatflux_left = heatflux_left - Q_combust_all_hydrogen;
		else
			% (kt / yr)                       = ((GJ)                 ) * (KJ / GJ) *
			combusted_fuel_flowrates(HYDROGEN) = (heatflux_left) * KJ_PER_GJ * ...
				... % (mol / KJ)        * (g / mol)       * (kt / g)
				( 1 / ENTHALPY_HYDROGEN) * MOLMASS_HYDROGEN * KT_PER_G;
			heatflux_left = 0;
			return
		end
	end

	% (GJ / yr) 		  = (kt / yr)          * (g / kt) * (kJ / g)		* (GJ / kJ)
	Q_combust_all_methane = flowrates(METHANE) * G_PER_KT * ENTHALPY_METHANE * GJ_PER_KJ;
	
	% Methane
	if (heatflux_left > Q_combust_all_methane)
		combusted_fuel_flowrates(METHANE) = flowrates(METHANE);
		heatflux_left = heatflux_left - Q_combust_all_methane;
	else
		% (kt / yr) 					  = ((GJ)                 ) * (KJ / GJ) *
		combusted_fuel_flowrates(METHANE) = (heatflux_left) * KJ_PER_GJ * ...
			... % (mol / KJ) 		* (g / mol) 	  * (kt / g)
			( 1 / ENTHALPY_METHANE) * MOLMASS_METHANE * KT_PER_G;
		heatflux_left = 0;
		return
	end

	% (GJ / yr)           = (kt / yr)          * (g / kt) * (kJ / g)        * (GJ / kJ)
	Q_combust_all_propane = flowrates(PROPANE) * G_PER_KT * ENTHALPY_PROPANE * GJ_PER_KJ;

	% Propane
	if (heatflux_left > Q_combust_all_propane)
		combusted_fuel_flowrates(PROPANE) = flowrates(PROPANE);
		heatflux_left = heatflux_left - Q_combust_all_propane;
	else
		% (kt / yr)                        = ((GJ)                 ) * (KJ / GJ) *
		combusted_fuel_flowrates(PROPANE) = (heatflux_left) * KJ_PER_GJ * ...
			... % (mol / KJ)        * (g / mol)       * (kt / g)
			( 1 / ENTHALPY_PROPANE) * MOLMASS_PROPANE * KT_PER_G;
		heatflux_left = 0;
		return
	end

	% (GJ / yr)           = (kt / yr)          * (g / kt) * (kJ / g)        * (GJ / kJ)
	Q_combust_all_butane = flowrates(BUTANE) * G_PER_KT * ENTHALPY_BUTANE * GJ_PER_KJ;

	% Butane
	if (heatflux_left > Q_combust_all_butane)
		combusted_fuel_flowrates(BUTANE) = flowrates(BUTANE);
		heatflux_left = heatflux_left - Q_combust_all_butane;
	else
		% (kt / yr)                        = ((GJ)                 ) * (KJ / GJ) *
		combusted_fuel_flowrates(BUTANE) = (heatflux_left) * KJ_PER_GJ * ...
			... % (mol / KJ)        * (g / mol)       * (kt / g)
			( 1 / ENTHALPY_BUTANE) * MOLMASS_BUTANE * KT_PER_G;
		heatflux_left = 0;
		return
	end
end

%        GJ   =            (kta   ,  __         , bar      , C )
function heat = heat_steam(F_steam, STEAM_CHOICE, P_reactor, T_reactor)
	global COST_RATES_STEAM;
	global STEAM_PRESSURE_COL STEAM_TEMP_COL COST_RATES_STEAM G_PER_KT ...
			MOLMASS_WATER BAR_PER_PSIA C_TO_K HEAT_CAPACITY_WATER GJ_PER_KJ;

	P_steam = COST_RATES_STEAM(STEAM_CHOICE, STEAM_PRESSURE_COL); % [ psia ]
	T_steam = COST_RATES_STEAM(STEAM_CHOICE, STEAM_TEMP_COL);	  % [ C ]
	P_steam = P_steam * BAR_PER_PSIA;
	T_steam = T_steam + C_TO_K;
	T_reactor = T_reactor + C_TO_K;

	if (P_steam > P_reactor) % Adiabatic Expansion
		T_adibatic = (T_steam) * (P_reactor / P_steam);
		T_steam = T_adibatic;
	elseif (P_steam < P_reactor) % Compression
		W = compressor_work(T_reactor, P_steam, P_reactor);
		if ADD_COMPRESSOR_WORK_TO_STEAM_HEATFLUX
			heat = heat + W;
		end
		% I should add this to the heat flux probably ?? 
	end
	
	% KJ = kta     * (G / KT) * (mol / g)         * (KJ / MOL K)        * (K - K)
	heat = F_steam * G_PER_KT * (1/MOLMASS_WATER) * HEAT_CAPACITY_WATER * (T_reactor - T_steam);
	% GJ = KJ   * (KJ / GJ)
	heat = heat * GJ_PER_KJ;


	% Heat flux after temperture 



end

function T_f = adiabatic_temp(T_0, P_0, P_f)

	T_f = T_0 * ( P_0 / P_f);
end

function W = compressor_work(T, P_0, P_f)
	R = 8.314;		% [ J / mol K]

	W = - n * R * T * log(P_f / P_0);

end

% HELPER FUNCTIONS | TAXES______________________________________________

function cost = tax_C02(combusted_flowrates, F_natural_gas)
	global HYDROGEN METHANE ETHYLENE PROPANE BUTANE TAX_CO2_PER_MT;
	global MT_CO2_PER_KT_METHANE MT_CO2_PER_KT_PROPANE MT_CO2_PER_KT_BUTANE ...
	MT_CO2_PER_KT_NATURALGAS;

	% Calculate the cost per kt (in tax) of each combusted fuel
	methane = combusted_flowrates(METHANE);
	propane = combusted_flowrates(PROPANE);
	butane = combusted_flowrates(BUTANE);

	mt_c02 = 0;
	% kta  =  (MT)  + ( (kt fuel / yr) * (MT CO2 / KT FUEL) )
	mt_c02 = mt_c02 + methane * MT_CO2_PER_KT_METHANE;
	mt_c02 = mt_c02 + propane * MT_CO2_PER_KT_PROPANE;
	mt_c02 = mt_c02 + butane * MT_CO2_PER_KT_BUTANE;
	mt_c02 = mt_c02 + F_natural_gas * MT_CO2_PER_KT_NATURALGAS;

	cost = mt_c02 * TAX_CO2_PER_MT;
end

% HELPER FUNCTIONS | FUEL COSTS______________________________________________

function cost = cost_natural_gas_fuel(heat_flux_remaining)
	global VALUE_NATGAS_FUEL
	% $ / yr = (GJ) 		   * ($ / GJ)
	cost = heat_flux_remaining * VALUE_NATGAS_FUEL;
end 

% HELPER FUNCTIONS | FUEL FLOWRATES______________________________________________

function F_natural_gas = natgas_combustion(heat_flux_remaining)
	global KJ_PER_GJ ENTHALPY_NAT_GAS KT_PER_G MOLMASS_NATGAS;
	% output should be in kta, input is in GJ 

	%		kt			GJ				* (kJ / GJ) * (mol / kJ) *			(g / mol) *				(kt / g)
	F_natural_gas = heat_flux_remaining * KJ_PER_GJ * (1/ENTHALPY_NAT_GAS) * (MOLMASS_NATGAS) * KT_PER_G;

end

% FUNCTIONS | REACTOR ODE SYSTEM________________________________________________

function dFdV = reactionODEs(V, F, T, P, F_steam)
	global R_2 k1_f k1_r k2 k3 C_TO_K MOLMASS_METHANE MOLMASS_ETHANE MOLMASS_ETHYLENE ... 
		MOLMASS_PROPANE MOLMASS_HYDROGEN MOLMASS_BUTANE YR_PER_SEC G_PER_KT SEC_PER_YR KT_PER_G
	% INPUT UNITS 
	% V [ L ]
	% F [ kta ]
	% T [ Celcius ]
	% P [ bar ]

	% Change the input units so that evrything is consistent
	% P = P * ATM_PER_BAR;
	T = T + C_TO_K;
	
	% Product flow rate indicies 
	HYDROGEN = 1;
	METHANE = 2;
	ETHYLENE = 3;
	PROPANE = 4;
	BUTANE = 5;

	% Feed flow rate index
	ETHANE = 6;

	F_tot = sum(F) + F_steam;


	% Hydrogen = A
	dFAdV = (k1_f(T) * ( (F(ETHANE) * P) / (F_tot * R_2 * T) )   ) - ...
			(k1_r(T) * ( F(ETHYLENE) * F(HYDROGEN) * P^2) ) / (F_tot * R_2 * T)^2;
	
	% Methane = B
	dFBdV = (k2(T) * (F(ETHANE) * P)^2) / (F_tot * R_2 * T)^2;

	% Ethylene = C
	dFCdV = (k1_f(T) * (F(ETHANE) * P / (F_tot * R_2 * T))) - ...
			(k1_r(T) * (F(ETHYLENE) * F(HYDROGEN) * P^2) / (F_tot * R_2 * T)^2) - ...
			(k3(T) * (F(ETHANE) * F(ETHYLENE) * P^2) / (F_tot * R_2 * T)^2);

	% Propane = E
	dFEdV = k2(T) * (F(ETHANE) * P)^2 / (F_tot * R_2 * T)^2;

	% Butane = F
	dFFdV = (k3(T) * (F(ETHANE) * F(ETHYLENE) * P^2)) / (F_tot * R_2 * T)^2;

	% Ethane = D
	dFDdV = (-k1_f(T) * (F(ETHANE) * P / (F_tot * R_2 * T))) + ...
			(k1_r(T) * (F(ETHYLENE) * F(HYDROGEN) * P^2)/(F_tot * R_2 * T)^2) - ...
			(k2(T) * F(ETHANE)^2 * P^2 / (F_tot * R_2 * T)^2) - ...
			(k3(T) * F(ETHANE) * F(ETHYLENE) * P^2 / (F_tot * R_2 * T)^2);
	
	T = T - C_TO_K;

	dFdV = [dFAdV; dFBdV; dFCdV; dFEdV; dFFdV; dFDdV];
	
end

% function cost = cost_reactor(V_plant_input)
% 	global FT_PER_METER STEAM_TO_FEED_RATIO
% 	FT_PER_METER = 3.28084;
% 	% ??? WHAT ARE THE UNITS OF TIME
% 	%
% 	pi = 3.14159;
% 	D = 0.05;								% [m]
% 	V_plant_max = pi * (0.025)^2 * 20; 		%[m^3]
	
% 	% Reactors have a max length, so calculate the number of full size reactors
% 	% and add it to the cost of the one non-max length reactor

% 	cost = 0;

% 	% Find the Cost of the max-sized reactors
% 	num_of_additional_reactors = int64(V_plant_input / V_plant_max);	
% 	num_of_additional_reactors = double(num_of_additional_reactors);
	
% 	V_plant = V_plant_max;
% 	factor_1 = 3.25;
% 	factor_2 = (V_plant / (pi * (D/2)^2) * FT_PER_METER)^0.82;
% 	factor_3 = (101.9 * D * FT_PER_METER)^1.066;
% 	factor_4 = (1800 / 280);
% 	cost_max_reactor = factor_1 * factor_2 * factor_3 * factor_4; 
% 	cost = cost + num_of_additional_reactors * cost_max_reactor;
	

% 	% Find the cost of the non-max size reactor 
% 	V_plant = V_plant_input - V_plant_max * num_of_additional_reactors;
% 	if V_plant < 0
% 		V_plant = 0;
% 	end
% 	factor_1 = 3.25;
% 	factor_2 = (V_plant / (pi * (D/2)^2) * FT_PER_METER)^0.82;
% 	factor_3 = (101.9 * D * FT_PER_METER)^1.066;
% 	factor_4 = (1800 / 280);
% 	cost = cost + factor_1 * factor_2 * factor_3 * factor_4;


% end

function cost = cost_reactor(V_plant_input)
	global FT_PER_METER STEAM_TO_FEED_RATIO
	FT_PER_METER = 3.28084;
	% ??? WHAT ARE THE UNITS OF TIME
	%
	pi = 3.14159;
	D = 0.05;								% [m]
	V_plant_max = pi * (0.025)^2 * 20; 		%[m^3]
	
	% Reactors have a max length, so calculate the number of full size reactors
	% and add it to the cost of the one non-max length reactor

	cost = 0;

	% Find the Cost of the max-sized reactors
	num_of_additional_reactors = int64(V_plant_input / V_plant_max);	
	num_of_additional_reactors = double(num_of_additional_reactors);
	
	V_plant = V_plant_max;
	factor_1 = 4.18;
	factor_2 = (V_plant / (pi * (D/2)^2) * FT_PER_METER)^0.82;
	factor_3 = (101.9 * D * FT_PER_METER)^1.066;
	factor_4 = (1800 / 280);
	cost_max_reactor = factor_1 * factor_2 * factor_3 * factor_4; 
	cost = cost + num_of_additional_reactors * cost_max_reactor;
	

	% Find the cost of the non-max size reactor 
	V_plant = V_plant_input - V_plant_max * num_of_additional_reactors;
	if V_plant < 0
		V_plant = 0;
	end
	factor_1 = 4.18;
	factor_2 = (V_plant / (pi * (D/2)^2) * FT_PER_METER)^0.82;
	factor_3 = (101.9 * D * FT_PER_METER)^1.066;
	factor_4 = (1800 / 280);
	cost = cost + factor_1 * factor_2 * factor_3 * factor_4;


end

function cost = separation_cost(P_flowrates)
	clear; clc;
	R = 8.314; 		% [J / mol K ]
	T = 173.15; 	% [ K ] 
		% This is the average of boiling points non-wighted

	TEMP_SEPARATION = 173.15;

	%Using compositions from ASPEN [mol fraction ]
	z_ethane = 0.105622966291020;
	z_ethylene = 0.270187148320469;
	z_hydrogen = 0.310508998651457;
	z_methane = 1.25317745797863e-005;
	z_propane = 1.25317745797863e-005;
	z_butane = 4.03219013378326e-002;
	z_water = 0.273333921850062;

	% mol fractions AFTER separationS
	x_water = 1;
	x_ethane = 1;
	x_propane = 0.9997;
	x_butane = 0.0003;
	x_ethylene = 1;
	x_hydrogen = 1; 
	x_methane = 4.03293090303065e-004;

	% PSA pressure system (system not at current pressure)
	P_in = 2;
	P_H2 = 10;
	P_ME = 1;

	%Using flow rates from ASPEN [NOTE: FOR MATLAB USE THE VALUES FROM THE
	%SOLN_TABLE. WE USED THESE AS EXPECTED COSTS)

	% INPUT MY FLOWRATES INSTEAD OF USING THESE
	F_water = 237.1; %mol/s
		% Steam coming out of the reactor 
	F_LPG = 34.97; %mol/s
		% 
	F_ethylene = 234.3; %mol/s
	F_ethane = 91.6; %mol/s
	F_H2 = 0.8716; %mol/s
	F_ME = 2.69e-5; %mol/s

	% j / s = ...
	W_min_Sep_System = F_water*R*T*log(x_water/z_water) + ...
		F_LPG*R*T*log(x_propane/z_propane + x_butane/z_butane) + ...
		F_ethylene*R*T*log(x_ethylene/z_ethylene) + ...
		F_ethane*R*T*log(x_ethane/z_ethane) + ...
		R*T*(F_H2*log(P_H2/P_in)+ F_H2*log(x_hydrogen/z_hydrogen) + ...
		F_ME*log(x_methane/z_methane) + F_ME*log(P_ME/P_in)) ;

	OPEX_MAX_Sep_System = W_min_Sep_System/1e9 * 50 * 3 * 30.24e6
	OPEX_MIN_Sep_System = W_min_Sep_System/1e9 * 20 * 3 * 30.24e6

	CAPEX_MAX_Sep_System = 50 * W_min_Sep_System
	CAPEX_MIN_Sep_System = 30 * W_min_Sep_System

	TFCI_Max_Sep_System = 2.5 * CAPEX_MAX_Sep_System 
	TFCI_Min_Sep_System = 2.5 * CAPEX_MIN_Sep_System 
end

%        [$]  =                  ( kta   , C, bar)
function cost = cost_waste_stream(F_steam, T, P)
	global MOLMASS_WATER G_PER_KT YR_PER_SEC R_2 M3_PER_L
	% 15.15 is the flowrate from aspen 
	% replace the 15.15 with the m^3 / s of the flowrate of steam out of the reactor 
	% MAKE SURE ITS m^3 / s

	% This is the OUTLET temperature and pressure of the sep SYSTEM
	T_sep = T;
	P_sep = P;

	% q = 15.15/3600; %flow rate of waste water out of the plant m^3/s

	% mol/s = ( kt / yr) * (g / kt) * (mol / g)   * (yr / s)
	F_steam = F_steam * G_PER_KT * (MOLMASS_WATER) * YR_PER_SEC;

	% L / s =  (mol / s) * [ L bar / K mol ] * (Kelvin) / (bar)
	q = F_steam * R_2 * T / P

	% m^3 / s = (L / s) * (m^3 / L)
	q = q * M3_PER_L;

	a = 0.001 + 2e-4*q^(-0.6); %Source: Uldrich and Vasudevan
	b=0.1; %Source: Uldrich and Vasudevan
	CEPCI = 820; %Source: Lecture slides
	C_f = 3.0; %$/GJ

	%$/m^3 waste water
	cost_waste_water = a*CEPCI + b*C_f

	%$/yr = $/m^3 * m^3/hr * hr/yr
	% cpy_waste_water = cost_waste_water*15.15*8760 
	
	% m^3 / s = (m^3 / s) * (s / yr)
	q = q * SEC_PER_YR;
	% cpy_waste_water = cost_waste_water*(FLOWRATE WATER m^3 / yr)
	cost = cost_waste_water * q;

end













































