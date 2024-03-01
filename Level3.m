% WORKSPACE VARIABLES ARE NOT BEING CLEARED 
clc; close all; 

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
global MOLMASS_METHANE;

% USER NOTES____________________________________________________________________

% Note: The primary (high level) units of this script are ... 
% Mass			kta
% Energy		GJ
% Pressure 		Bar
% Temperature 	Celcius
% Moles			Gigamoles ???
% Value			Millions of Dollars ($ MM)

% [ __ ] THIS MEANS DIMENSIONLESS UNITS

% DESIGN PARAMETERS & USER INPUTS_______________________________________________

% Feed
STEAM_TO_FEED_RATIO = 1.0;		% [ __ ] 0.6 to 1.0

% Product
P_ETHYLENE_DES = 200;			% [	kta ]
	% Note! This design parameter's units are changed prior to the matrix def 

% Reactor Conditions
TEMP_RXTR = 800;				% [ C ] 775 to 825 C ???
PRESS_RXTR = 3;					% [ Bar ]
TEMP_ETHANE_FEED = 25;			% [ C ]
CONVERSION = 0.8;				% [ __ ]
USERINPUT_S1 = 0.60;				% [ __ ]
USERINPUT_S2 = 0.15; 			% [ __ ]

% Plotting : Tabs mean one input is dependent on another 
CONSOLE_OUTPUT_EFFECTIVE_VALUE_FUELS = true;
OUTPUT_LVL3_FLOWRATES_TO_CONSOLE = true;
	SANITY_CHECK_CALCULATIONS = true;
CALCULATE_ALL_SELECTIVITIES = true;
	PLOT_ECON_3D = true;
	PLOT_ECON_COUNTOUR = true;
%_______________________________________________________________________________
% DON'T TOUCH ANYTHING BELOW THIS LINE
%_______________________________________________________________________________

% CONSTANTS | PLOTTING_____________________________________________________

CONSOLE_SECTION_DIVIDER = ...
	"_____________________________________________________________________";
S1_MIN = 0.40;
S1_MAX = 1.00;
S1_POINTS = 100;
S2_MIN = 0.05;
S2_MAX = 0.20;
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

% Value
MMDOLLA_PER_DOLLA = 10^-6;	% [ $ MM / $]
DOLLA_PER_MMDOLLA = 10^6;	% [ $ / $ MM ]

% CONSTANTS | CHEMICAL_____________________________________________________

% Chemical | Molar Mass
MOLMASS_HYDROGEN = 1.008;				% [ g / mol ]
	% Source : ?? 
MOLMASS_METHANE = 16.04;				% [ g / mol ]
	% Source : ?? 
MOLMASS_CO2 = 44.01;					% [ g / mol ]
	% Source : ??
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

% Heat capacities 
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
ENTHALPY_METHANE = 890;					% [ kJ / mol ]	
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C74828&Mask=1
ENTHALPY_PROPANE = 2219.2;				% [ kJ / mol ]
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C74986&Mask=1
ENTHALPY_BUTANE = 2877.5;				% [ kJ / mol ]
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C106978&Mask=1
ENTHALPY_NAT_GAS = ENTHALPY_METHANE; 
	% Source : ???
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
% [ psia Temp[C] $/kg  kJ/kg ]
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
	STEAM_COST_ROW = 3;
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
TAX_CO2_PER_GJ_NATGAS = TAX_CO2_PER_GJ_METHANE; % ???

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

% FUNCTIONS | MATRIX (EXTENT OF RXN)_______________________________________

% (mol / yr)   = (kt / yr)        (g / kt) *  ( mol / g) 
P_ETHYLENE_DES = P_ETHYLENE_DES * G_PER_KT * (1 / MOLMASS_ETHYLENE);	

A = @(s1, s2)...
	[s1-1	,s1		,s1+1	;
     s2		,s2-1	,s2		;
     1		,0		,-1		];
b = [	0;	
		0;		
		P_ETHYLENE_DES	];

% FUNCTIONS | FLOWRATE_____________________________________________________

% (kta)	   =                        (mol / yr)
P_HYDROGEN = @(xi_1)				xi_1 * ...
	... % (g / mol)       * (kt / g) 
		MOLMASS_HYDROGEN * KT_PER_G;
P_METHANE = @(xi_2)					xi_2 * ...
			MOLMASS_METHANE * KT_PER_G;
P_ETHYLENE = @(xi_1, xi_3)			(xi_1 - xi_3) * ...
			MOLMASS_ETHYLENE * KT_PER_G;
P_PROPANE = @(xi_2)					xi_2 * ...
			MOLMASS_PROPANE * KT_PER_G;
P_BUTANE = @(xi_3)					xi_3 * ...
			MOLMASS_BUTANE * KT_PER_G;
F_ETHANE = @(xi_1, xi_2, xi_3)		(xi_1 + xi_2 + xi_3) * ...
			MOLMASS_ETHANE * KT_PER_G;

% FUNCTIONS | VALIDATION___________________________________________________

flowrates_valid = @( flowrates ) all(flowrates >= 0);

% FUNCTIONS | ECONOMICS____________________________________________________

% ($ / yr) =                    (kta) *   (MT / KT) * ($ / MT)
value_ethane = @(P_ethane) P_ethane * MT_PER_KT * VALUE_ETHANE;
value_ethylene = @(P_ethylene) P_ethylene * MT_PER_KT * VALUE_ETHYLENE;
value_h2_chem = @(P_h2_chem) P_h2_chem * MT_PER_KT * VALUE_HYDROGEN_CHEM;
value_methane = @(P_methane) P_methane * MT_PER_KT * VALUE_METHANE_FUEL;
value_propane = @(P_propane) P_propane * MT_PER_KT * VALUE_PROPANE_FUEL;
value_butane = @(P_butane) P_butane * MT_PER_KT * VALUE_BUTANE_FUEL;

% ($ / yr) = 							(kta) *  (kg / kt) * ($ / kg)
cost_steam = @(F_steam, steam_rate) F_steam * KG_PER_KT * steam_rate;

% FUNCTIONS | THEROMODYNAMICS______________________________________________
% (GJ / yr) =                          (kta) *  (g / KT)   * (mol gas/ g gas)    * (kJ / mol K)          * (GJ / KJ) * (K)
heat_ethane = @(F_ethane, T0, Tf) F_ethane * G_PER_KT * (1 / MOLMASS_ETHANE) * HEAT_CAPACITY_ETHANE * GJ_PER_KJ * (Tf - T0);

% (GJ / yr)  =      (mol / yr) * (kJ / mol) *   (GJ / kJ)
heat_rxn1 = @(xi_1) xi_1 * ENTHALPY_RXN_1 * GJ_PER_KJ;
heat_rxn2 = @(xi_2) xi_2 * ENTHALPY_RXN_2 * GJ_PER_KJ;
heat_rxn3 = @(xi_3) xi_3 * ENTHALPY_RXN_3 * GJ_PER_KJ;
heat_rxn = @(xi) heat_rxn1(xi(1)) + heat_rxn2(xi(2)) + heat_rxn3(xi(3)); 

% SCRIPT___________________________________________________________________

% Economics | Post-Tax Value of different fuel sources 
if (CONSOLE_OUTPUT_EFFECTIVE_VALUE_FUELS)
	disp(" [ $ / GJ ] ")
	EFFECTIVE_VALUE_METHANE_FUEL = VALUE_HYDROGEN_FUEL + TAX_CO2_PER_GJ_METHANE
	EFFECTIVE_VALUE_PROPANE_FUEL = VALUE_PROPANE_FUEL + TAX_CO2_PER_GJ_PROPANE
	EFFECTIVE_VALUE_BUTANE_FUEL = VALUE_BUTANE_FUEL + TAX_CO2_PER_GJ_BUTANE
	EFFECTIVE_VALUE_NAT_GAS_FUEL = VALUE_NATGAS_FUEL + TAX_CO2_PER_GJ_NATGAS
% 	EFFECTIVE_VALUE_NUM2_FUEL = VALUE_NATGAS_FUEL + TAX_CO2_PER_GJ_NUM2;
		% Not using number 2 fuel bc its too expensive 
end

if (OUTPUT_LVL3_FLOWRATES_TO_CONSOLE)
	xi = A(USERINPUT_S1, USERINPUT_S2) \ b;

	% Calculate the flow rates of each species (kta)
	P_hydrogen = P_HYDROGEN(xi(1));
	P_methane = P_METHANE(xi(2));
	P_ethylene = P_ETHYLENE(xi(1), xi(3));
	P_propane = P_PROPANE(xi(2));
	P_butane = P_BUTANE(xi(3));
	P_flowrates = [ P_hydrogen, P_methane, P_ethylene, P_propane, P_butane ];	

	disp(CONSOLE_SECTION_DIVIDER)
	if (flowrates_valid(P_flowrates))
		
		fprintf("Flowrates for the reactor given that s1 = %f, s2 = %f\n\n", ...
			USERINPUT_S1, USERINPUT_S2)

		disp(CONSOLE_SECTION_DIVIDER)
		disp("Level 2 Flowrates  in / out of the entire plant [ kt / yr ]")
		P_hydrogen = P_HYDROGEN(xi(1))
		P_methane = P_METHANE(xi(2))
		P_ethylene = P_ETHYLENE(xi(1), xi(3))
		P_propane = P_PROPANE(xi(2))
		P_butane = P_BUTANE(xi(3))

		disp("Fresh Feed Flowrate")
		F_ethane = F_ETHANE(xi(1), xi(2), xi(3))
		
		disp(CONSOLE_SECTION_DIVIDER)
		disp("Level 3 Flowrates [ kt / yr ] ")

		disp("Recycle Stream Flowrate")

		R_ethane = F_ethane * ((1-CONVERSION) / (CONVERSION))
		R_ethane = (P_ethylene/USERINPUT_S1) * ((1-CONVERSION)/CONVERSION)

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
			
			% Solve for extents of reaction
			xi = A(s1, s2) \ b;

			% Calculate the flow rates of each species (kta)
			P_hydrogen = P_HYDROGEN(xi(1));
			P_methane = P_METHANE(xi(2));
			P_ethylene = P_ETHYLENE(xi(1), xi(3));
			P_propane = P_PROPANE(xi(2));
			P_butane = P_BUTANE(xi(3));

			F_ethane = F_ETHANE(xi(1), xi(2), xi(3));

			P_flowrates = [ P_hydrogen, P_methane, P_ethylene, P_propane, P_butane ];
		
			if (flowrates_valid(P_flowrates))

				% Store for plotting (kta)
				ethylene_flowrates(i) = P_ETHYLENE(xi(1), xi(3));
				hydrogen_flowrates(i) = P_HYDROGEN(xi(1));
				methane_flowrates(i) = P_METHANE(xi(2));
				ethylene_flowrates(i) = P_ETHYLENE(xi(1), xi(3));
				propane_flowrates(i) = P_PROPANE(xi(2));
				butane_flowrates(i) = P_BUTANE(xi(3));
				ethane_flowrates(i) = F_ETHANE(xi(1), xi(2), xi(3));

				% Calculate the heat flux needed to keep reactor isothermal 
				heat_flux = 0;
				F_steam = STEAM_TO_FEED_RATIO * F_ethane;
				heat_flux = heat_flux + heat_ethane(P_ethylene, TEMP_ETHANE_FEED, TEMP_RXTR);
				heat_flux = heat_flux + heat_steam(F_steam, STEAM_50PSIA, PRESS_RXTR, TEMP_RXTR) ;
				heat_flux = heat_flux + heat_rxn(xi);

	% 			% Use the heat flux to calculate the fuel cost	
				[combusted_fuel_flow_rates, heat_flux_remaining] = fuel_combustion(heat_flux, P_flowrates);

				% Calculate how much natural gas you needed to combust
				F_natural_gas = natgas_combustion(heat_flux_remaining);

				% Determine how much of the product streams were combusted to keep the reactor isothermal	
				% Assume: no hydrogen is combusted
				combusted_methane = combusted_fuel_flow_rates(METHANE);
				combusted_propane = combusted_fuel_flow_rates(PROPANE);
				combusted_butane = combusted_fuel_flow_rates(BUTANE);

	% 			% VALUE CREATED | Primary Products
				profit(i) = profit(i) + value_ethylene(P_ethylene);
				profit(i) = profit(i) + value_h2_chem(P_hydrogen); % Assume no H2 combusted

				% VALUE CREATED | Non-combusted fuels 
				profit(i) = profit(i) + value_methane(P_methane - combusted_methane);
				profit(i) = profit(i) + value_propane(P_propane - combusted_propane);
				profit(i) = profit(i) + value_butane(P_butane - combusted_butane);	

				% COSTS INCURRED
				profit(i) = profit(i) - tax_C02(combusted_fuel_flow_rates, F_natural_gas);
	% 			profit(i) = profit(i) - cost_steam(F_steam, COST_RATES_STEAM(STEAM_COST_ROW,STEAM_50PSIA));
				profit(i) = profit(i) - value_ethane(F_ethane);
				profit(i) = profit(i) - cost_natural_gas_fuel(F_natural_gas);
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

disp("The Script is done running ⌛️")
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
    clabel(C, h,  'FontSize', 10, 'Color', 'k', 'LabelSpacing', 200); % Customize label properties
	xlabel(x_label);
	ylabel(y_label);
	title(plt_title);
	saveas(gcf, plt_saveName);
	hold off
end 

function plot_3D(x, y, z, options)
    % Unpack options 
    x_label = options{1};
    y_label = options{2};
    plt_title = options{3};
    plt_saveName = options{4};

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

function [combusted_fuel_flowrates, heatflux_left] = fuel_combustion(heat_flux, flowrates)
	global HYDROGEN METHANE ETHYLENE PROPANE BUTANE;
	global ENTHALPY_METHANE ENTHALPY_PROPANE ENTHALPY_BUTANE HEAT_CAPACITY_ETHANE;
	global MT_PER_KT G_PER_KT GJ_PER_KJ KJ_PER_GJ MOLMASS_METHANE KT_PER_G;

	% Note! : Longest Chain Hydrocarbons are cheapest to combust

	% initialize all values in the array to be zero 
	combusted_fuel_flowrates = flowrates * 0;

	% LOGIC : Goes through each heat source in order, returns if the heat flux supplied is sufficient.
	heatflux_left = heat_flux; 

	% (GJ / yr) 		  = (kt / yr)          * (g / kt) * (kJ / g)		* (GJ / kJ)
	Q_combust_all_methane = flowrates(METHANE) * G_PER_KT * ENTHALPY_METHANE * GJ_PER_KJ;
	% (GJ)   	  = (GJ)          - (GJ)
	heatflux_left = heatflux_left - Q_combust_all_methane;
	
	% Methane
	if (heatflux_left > 0)
		combusted_fuel_flowrates(METHANE) = flowrates(METHANE);
	else
		% (kt / yr) 					  = ((GJ)                 ) * (KJ / GJ) *
		combusted_fuel_flowrates(METHANE) = (Q_combust_all_methane) * KJ_PER_GJ * ...
			... % (mol / KJ) 		* (g / mol) 	  * (kt / g)
			( 1 / ENTHALPY_METHANE) * MOLMASS_METHANE * KT_PER_G;
		heatflux_left = 0;
		return
	end

    % Propane
	% (GJ / yr) 		  = (kt / yr)          * (g / kt) * (kJ / g)		* (GJ / kJ)
    Q_combust_all_propane = flowrates(PROPANE) * G_PER_KT * ENTHALPY_PROPANE * GJ_PER_KJ;
    heatflux_left = heatflux_left - Q_combust_all_propane;

    if (heatflux_left > 0)
        combusted_fuel_flowrates(PROPANE) = flowrates(PROPANE);
    else
        combusted_fuel_flowrates(PROPANE) = heat_flux / (ENTHALPY_PROPANE * GJ_PER_KJ);
        heatflux_left = 0;
        return
    end

    % Butane
    Q_combust_all_butane = flowrates(BUTANE) * ENTHALPY_BUTANE * GJ_PER_KJ;
    heatflux_left = heatflux_left - Q_combust_all_butane;

    if (heatflux_left > 0)
        combusted_fuel_flowrates(BUTANE) = flowrates(BUTANE);
    else
        combusted_fuel_flowrates(BUTANE) = (heatflux_left + Q_combust_all_butane) / (ENTHALPY_BUTANE * GJ_PER_KJ);
        heatflux_left = 0;
        return
    end

	%DEBUGGING 
	heatflux_from_combustion = heat_flux - heatflux_left;

end

function heat = heat_steam(F_steam, STEAM_50C, P_reactor, T_reactor)
	global COST_RATES_STEAM;
	global STEAM_PRESSURE_COL STEAM_TEMP_COL COST_RATES_STEAM G_PER_KT

	P_steam = COST_RATES_STEAM(STEAM_50C, STEAM_PRESSURE_COL);
	T_steam = COST_RATES_STEAM(STEAM_50C, STEAM_TEMP_COL);

	if (P_steam > P_reactor)
		T_adiabatic = T_steam * (P_reactor / P_steam);
		heat = F_steam * G_PER_KT; 
	end

	heat = 0; % REMOVE ME 

end

function cost = tax_C02(combusted_flowrates, F_natural_gas)
	global HYDROGEN METHANE ETHYLENE PROPANE BUTANE TAX_CO2_PER_MT;
	global MT_CO2_PER_KT_METHANE MT_CO2_PER_KT_PROPANE MT_CO2_PER_KT_BUTANE ...
	MT_CO2_PER_KT_NATURALGAS;

	cost = 0;
	% Calculate the cost per kt (in tax) of each combusted fuel
	combusted_flowrates;
	methane = combusted_flowrates(METHANE);
	propane = combusted_flowrates(PROPANE);
	butane = combusted_flowrates(BUTANE);

	mt_c02 = methane * MT_CO2_PER_KT_METHANE;
	mt_c02 = mt_c02 + propane * MT_CO2_PER_KT_PROPANE;
	mt_c02 = mt_c02 + butane * MT_CO2_PER_KT_BUTANE;
	mt_c02 = mt_c02 + F_natural_gas * MT_CO2_PER_KT_NATURALGAS;

	cost = mt_c02 * TAX_CO2_PER_MT;
	
	% Calculate the cost of the remaining natural gas C02 tax



end

function cost = cost_natural_gas_fuel(heat_flux_remaining)
	global VALUE_NATGAS_FUEL
	cost = heat_flux_remaining * VALUE_NATGAS_FUEL;
end 

function F_natural_gas = natgas_combustion(heat_flux_remaining)
	global KJ_PER_GJ ENTHALPY_NAT_GAS KT_PER_G MOLMASS_NATGAS;
	% output should be in kta, input is in GJ 

	%		kt			GJ				* (kJ / GJ) * (mol / kJ) *			(g / mol) *				(kt / g)
	F_natural_gas = heat_flux_remaining * KJ_PER_GJ * (1/ENTHALPY_NAT_GAS) * (MOLMASS_NATGAS) * KT_PER_G;

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











