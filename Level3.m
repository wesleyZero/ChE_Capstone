clc; clear; close all;
global S1_MIN S1_MAX S1_POINTS;
global S2_MIN S2_MAX S2_POINTS;
global INVALID_FLOWRATE;
global Fethyl_S1S2_plotOpt;
% global P_ETHYLENE;
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
global HYDROGEN METHANE ETHYLENE PROPANE BUTANE;
global ENTHALPY_METHANE ENTHALPY_PROPANE ENTHALPY_BUTANE HEAT_CAPACITY_ETHANE;



% DESIGN PARAMETERS_____________________________________________________________
STEAM_TO_FEED_RATIO = 0.6; %0.6 to 1.0

% Note: The primary units of this script are ... 
% Mass			kta
% Energy		GJ
% Pressure 		Bar
% Temperature 	Celcius

% CONSTANTS | HARD CODED________________________________________________________

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

% Chemical | Molar Mass
MOLMASS_PROPANE = 44.0956;				% [ g / mol ]
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C74986&Mask=1
MOLMASS_BUTANE = 58.1222;				% [ g / mol ]
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C106978&Mask=1
MOLMASS_ETHANE = 30.0690;				% [ g / mol ]
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C74840
MOLMASS_ETHYLENE = 28.0532;				% [ g / mol ];
	% Source = https://webbook.nist.gov/cgi/cbook.cgi?ID=74-85-1&Type=IR-SPEC&Index=QUANT-IR,20

% Economic | Fuel
VALUE_H2_FUEL = 3;			% [ $ / GJ ]
VALUE_CH4_FUEL = 3;			% [ $ / GJ ]
VALUE_C3H6_FUEL = 3;		% [ $ / GJ ]
VALUE_C4H8_FUEL = 3;		% [ $ / GJ ]
VALUE_NATGAS_FUEL = 3;		% [ $ / GJ ]
VALUE_NUM2OIL_FUEL = 4.5;	% [ $ / US Gallon ]


% Thermodynamics | Heats of Formation (at 25C)
HEAT_FORMATION_ETHANE = -83.8;			% [ kJ / mol  ] reference Temp = std
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C74840&Units=SI&Mask=1EFF
HEAT_FORMATION_METHANE = -74.87;		% [ kJ / mol  ] reference Temp = std
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C74828&Mask=1
HEAT_FORMATION_ETHYLENE = 52.47;		% [ kJ / mol ] reference Temp = std 
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C74851&Mask=1
HEAT_FORMATION_HYDROGEN = 0; 			% [ kJ / mol  ] reference Temp = std 
HEAT_FORMATION_PROPANE = -104.7;		% [ kJ / mol  ] reference Temp = std 
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C74986&Mask=1
HEAT_FORMATION_BUTANE = -125.6;			% [ kJ / mol  ] reference Temp = std 
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C106978&Mask=1

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

% Design params
P_ETHYLENE_DES = 200;			% [	kta ]
% F_ETHYLENE_MOLS = F_ETHYLENE * G_PER_KT * (1 / MOLMASS_ETHYLENE);

% CONSTANTS | FXNS OF CONSTANTS__________________________________________________

% Thermodynamics | Enthalpy of Reactions
ENTHALPY_RXN_1 = HEAT_FORMATION_HYDROGEN + HEAT_FORMATION_ETHYLENE ...
										- HEAT_FORMATION_ETHANE;
ENTHALPY_RXN_2 = HEAT_FORMATION_METHANE + HEAT_FORMATION_PROPANE ...
										- 2 * HEAT_FORMATION_ETHANE; 
ENTHALPY_RXN_3 = HEAT_FORMATION_ETHANE - HEAT_FORMATION_ETHANE ...
										- HEAT_FORMATION_ETHYLENE;

% Economics | Enviormental
TAX_CO2_PER_MT = 125;				% [ $ / MT ]
TAX_CO2_PER_GJ_METHANE = KJ_PER_GJ * (1 / ENTHALPY_METHANE) * CO2_TO_METHANE_COMBUSTION_STOICH * MT_PER_G * TAX_CO2_PER_MT;
TAX_CO2_PER_GJ_PROPANE = KJ_PER_GJ * (1 / ENTHALPY_PROPANE) * CO2_TO_PROPANE_COMBUSTION_STOICH * MT_PER_G * TAX_CO2_PER_MT;
TAX_CO2_PER_GJ_BUTANE = KJ_PER_GJ * (1 / ENTHALPY_BUTANE) * CO2_TO_BUTANE_COMBUSTION_STOICH * MT_PER_G * TAX_CO2_PER_MT;

% Economics | Post-Tax Value of different fuel sources
EFFECTIVE_VALUE_METHANE_FUEL = VALUE_H2_FUEL - TAX_CO2_PER_GJ_METHANE;
EFFECTIVE_VALUE_PROPANE_FUEL = VALUE_C3H6_FUEL - TAX_CO2_PER_GJ_PROPANE;
EFFECTIVE_VALUE_BUTANE_FUEL = VALUE_C4H8_FUEL - TAX_CO2_PER_GJ_BUTANE;

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

% A = @(s1, s2)...
% 	[s1-1	,s1		,s1+1;
%      s2		,s2-1	,s2;
%      1		,2		,1	];
% b = [0;		0;		P_ETHYLENE];

% writing the last equation to hold the Product instead of the feed
% constant
% A = @(s1, s2)...
% 	[s1-1	,s1		,s1+1;
%      s2		,s2-1	,s2;
%      1		,0		,-1	];
% b = [0;		0;		P_ETHYLENE];

% Isaiahs recommendation
A = @(s1, s2)...
	[s1-1	,s1		,s1+1;
     s2		,s2-1	,s2;
     1		,1		,1	];
b = [0;		0;		P_ETHYLENE_DES];

% FUNCTIONS | FLOWRATE_____________________________________________________

P_HYDROGEN = @(xi_1)			xi_1;
P_METHANE = @(xi_2)				xi_2;
P_ETHYLENE = @(xi_1, xi_3)		xi_1 - xi_3;
P_PROPANE = @(xi_2)				xi_2;
P_BUTANE = @(xi_3)				xi_3;

F_ETHANE = @(xi_1, xi_2, xi_3)	xi_1 + xi_2 + xi_3;

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
% kta * (g / kta) * (mol / g) * (kJ / mol) * (GJ / kJ) * K 
heat_ethane = @(F_ethane, T0, Tf) F_ethane * G_PER_KT * (1 / MOLMASS_ETHANE) * HEAT_CAPACITY_ETHANE * GJ_PER_KJ * (Tf - T0);

% Input: [ kta ] 	Output: [ GJ ]
heat_rxn1 = @(xi_1) xi_1 * ENTHALPY_RXN_1 ;
heat_rxn2 = @(xi_2) xi_2 * ENTHALPY_RXN_2;
heat_rxn3 = @(xi_3) xi_3 * ENTHALPY_RXN_3 ;
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

hydrogen_flowrates = (s1_mesh + s2_mesh) .* 0;
methane_flowrates = (s1_mesh + s2_mesh) .* 0;
ethylene_flowrates = (s1_mesh + s2_mesh) .* 0;
propane_flowrates = (s1_mesh + s2_mesh) .* 0;
butane_flowrates = (s1_mesh + s2_mesh) .* 0;
ethane_flowrates = (s1_mesh + s2_mesh) .* 0;

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
			hydrogen_flowrates(i) = P_HYDROGEN(xi(1));
			methane_flowrates(i) = P_METHANE(xi(2));
			ethylene_flowrates(i) = P_ETHYLENE(xi(1), xi(3));
			propane_flowrates(i) = P_PROPANE(xi(2));
			butane_flowrates(i) = P_BUTANE(xi(3));
			ethane_flowrates(i) = F_ETHANE(xi(1), xi(2), xi(3));

% 			% Calculate the heat flux needed to keep reactor isothermal 
			heat_flux = 0;
			F_steam = STEAM_TO_FEED_RATIO * F_ethane;
			heat_flux = heat_flux + heat_ethane(P_ethylene, T_ethane_feed, T_reactor);
% 			heat_flux = heat_flux + heat_steam(F_steam, STEAM_50C, P_reactor, T_reactor); 
			heat_flux = heat_flux + heat_rxn(xi);


			cost_of_heating_ethane = heat_ethane(P_ethylene, T_ethane_feed, T_reactor) * 3
			cost_of_heating_reaction = heat_rxn(xi) * 3
% 
% 			% Use the heat flux to calculate the fuel cost	
			[combusted_fuel_flow_rates, heat_flux_remaining] = fuel_combustion(heat_flux, flowrates);
			heat_flux_remaining
% 			cost_of_num2Fuel = heat_flux_remaining * 30 
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
			profit(i) = profit(i) - cost_natural_gas_fuel(heat_flux_remaining);
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

plot_3D(s1_mesh, s2_mesh, profit, PROFIT_S1S2_OPT);

methane_flowrates
propane_flowrates = propane_flowrates * 0;
% Prepare the array of flow rate matrices
flowRatesArray = {hydrogen_flowrates, methane_flowrates, ethylene_flowrates, propane_flowrates, butane_flowrates, ethane_flowrates};


% Call the function with the desired row
plotFlowRatesForRow(1, flowRatesArray); % To plot the first row across all matrices


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
	global MT_PER_KT G_PER_KT GJ_PER_KJ;

	combusted_fuel_flowrates = flowrates * 0;
	heatflux_left = heat_flux; 
	% Longest Chain Hydrocarbons are cheapest to combust
	Q_combust_all_methane = flowrates(METHANE) * ENTHALPY_METHANE * GJ_PER_KJ;
	heatflux_left = heatflux_left - Q_combust_all_methane;
	
	% Goes through each heat source in order, returns if the heat flux supplied is sufficient.
	if (heatflux_left > 0)
		combusted_fuel_flowrates(METHANE) = flowrates(METHANE);
	else
		combusted_fuel_flowrates(METHANE) = (heatflux_left + Q_combust_all_methane) * ( 1 / ENTHALPY_METHANE);
		heatflux_left = 0;
		return
	end

    % Start with Propane
    Q_combust_all_propane = flowrates(PROPANE) * ENTHALPY_PROPANE * GJ_PER_KJ;
    heatflux_left = heatflux_left - Q_combust_all_propane;

    if (heatflux_left > 0)
        combusted_fuel_flowrates(PROPANE) = flowrates(PROPANE);
    else
        combusted_fuel_flowrates(PROPANE) = heat_flux / (ENTHALPY_PROPANE * GJ_PER_KJ);
        heatflux_left = 0;
        return
    end

    % Then Butane
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

	heat = 0;

end

function cost = tax_C02(combusted_flowrates, heatflux_left)
	

end

function cost = cost_natural_gas_fuel(heat_flux_remaining)

	cost = heat_flux_remaining * 30;

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
    xlabel('Column Index');
    ylabel('Flow Rate');
    legend('show');
    hold off; % Release the figure for other plots
end











