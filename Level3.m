clc; clear; close all; 

% CONSTANTS________________________________________________________________

global S1_MIN S1_MAX S1_POINTS;
global S2_MIN S2_MAX S2_POINTS;
global INVALID_FLOWRATE;
global Fethyl_S1S2_plotOpt;
global FC2H6;
global M_C2H6;
global MT_PER_KT G_PER_KT GJ_PER_KJ;
global VALUE_ETHANE VALUE_ETHYLENE VALUE_H2_CHEM;
global COST_STEAM;
global VALUE_H2_FUEL VALUE_CH4_FUEL VALUE_C3H6_FUEL VALUE_C4H8_FUEL;
global VALUE_NATGAS_FUEL VALUE_NUM2OIL_FUEL;
global COST_CO2 COST_WASTESTREAM;
global ENTHALPY_PROPANE ENTHALPY_BUTANE;
global MOLMASS_PROPANE MOLMASS_BUTANE;
global PROFIT_S1S2_OPT;
global STEAM_COSTS;

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
FC2H6 = 200;			% [	kta ]

% Molar mass
M_C2H6 = 28.05;			% [ g / mol ]

% Unit conversions 
MT_PER_KT = 1000;		% [ kt / MT ]
G_PER_KT = 10^9;		% [ g / kt ]
GJ_PER_KJ = 10^-6;		% [ GJ / kJ ]

% Economic | Chemicals
VALUE_ETHANE = 200;		% [ $ / MT ]
VALUE_ETHYLENE = 900;	% [ $ / MT ]
VALUE_H2_CHEM = 1400;	% [ $ / MT ]

% Economic | Steam 
% [psia Temp[C] $/kg kJ/kg
COST_STEAM = [
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

% Economics | Enviormental
COST_CO2 = 125;				% [ $ / MT ]
COST_WASTESTREAM = NaN;		% Ulrich and Vasudevan

% Thermodynamics | Enthalpy of combustion of gas at standard conditions
ENTHALPY_PROPANE = 2219.2;		% [ kJ / mol ]
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C74986&Mask=1
ENTHALPY_BUTANE = 2877.5;		% [ kJ / mol ]
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C106978&Mask=1

% Chemical | Molar Mass
MOLMASS_PROPANE = 44.0956;		% [ g / mol ]
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C74986&Mask=1
MOLMASS_BUTANE = 58.1222;
	% Source : https://webbook.nist.gov/cgi/cbook.cgi?ID=C106978&Mask=1

% Chemical | Combustion



% SYSTEM OF EQUARTIONS (EXTENT OF RXN)_____________________________________

A = @(s1, s2)...
	[s1-1	,s1		,s1+1;
     s2		,s2-1	,s2;
     1		,2		,1	];
b = [0;		0;		FC2H6];

% FUNCTIONS | FLOWRATE_____________________________________________________

P_H2 = @(xi_1)			xi_1;
P_CH4 = @(xi_2)			xi_2;
P_C2H4 = @(xi_1, xi_3)	xi_1 - xi_3;
P_C3H8 = @(xi_2)		xi_2;
P_C4H10 = @(xi_3)		xi_3;

% FUNCTIONS | VALIDATION___________________________________________________

flowrates_valid = @( flowrates ) all(flowrates >= 0);

% FUNCTIONS | ECONOMICS____________________________________________________

% Inputs are in [ kta ] outputs are in [ $ ]
value_ethane = @(P_ethane) P_ethane * MT_PER_KT * VALUE_ETHANE;
value_ethylene = @(P_ethylene) P_ethylene * MT_PER_KT * VALUE_ETHYLENE;
value_h2_chem = @(P_h2_chem) P_h2_chem * MT_PER_KT * VALUE_H2_CHEM;

% SCRIPT___________________________________________________________________

% Iterates through each value of selectivities S1 and S2 to find the economic
% potential for different reaction conditions 
s1_domain = linspace(S1_MIN, S1_MAX, S1_POINTS);
s2_domain = linspace(S2_MIN, S2_MAX, S2_POINTS);
[s1_mesh, s2_mesh] = meshgrid(s1_domain, s2_domain);
ethylene_flowrates = (s1_mesh + s2_mesh) .* 0;
profit = (s1_mesh + s2_mesh) .* 0;	

i = 1;
for s1 = s1_domain
	for s2 = s2_domain
		
		% Solve for extents of reaction
		xi = A(s1, s2) \ b;

		% Calculate the flow rates of each species
		p_h2 = P_H2(xi(1));
		p_ch4 = P_CH4(xi(2));
		p_c2h4 = P_C2H4(xi(1), xi(3));
		p_c3h8 = P_C3H8(xi(2));
		p_c4h10 = P_C4H10(xi(3));
		
		flowrates = [ p_h2, p_ch4, p_c2h4, p_c3h8, p_c4h10 ];

		if (flowrates_valid(flowrates))
			ethylene_flowrates(i) = p_c2h4;

			% Value Created 
			disp("product h2 LPG")			
			value_ethylene(p_c2h4)
			value_h2_chem(p_h2)
			value_LPG(p_c3h8, p_c4h10)

			profit(i) = profit(i) + value_ethylene(p_c2h4);
			profit(i) = profit(i) + value_h2_chem(p_h2);
			profit(i) = profit(i) + value_LPG(p_c3h8, p_c4h10);

			% Costs incurred
			fuels = 1 : 10; % Placeholder value
			profit(i) = profit(i) - cost_C02(fuels);
			profit(i) = profit(i) - cost_steam(flowrates);
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

function value = value_LPG(P_propane, P_butane)
	% ASSUMPTION : VALUE OF LPG IS JUST SUM VALUE OF PROPANE + BUTANE BEING
	% COMBUSTED
	% kt * (g / kt) * (mol / g) * (kJ / mol) * (GJ / KJ) * ($ / GJ)
	global G_PER_KT MOLMASS_PROPANE ENTHALPY_PROPANE GJ_PER_KJ ...
		VALUE_C3H6_FUEL MOLMASS_BUTANE ENTHALPY_BUTANE VALUE_C4H8_FUEL;

	prop_val = P_propane * G_PER_KT * (1 / MOLMASS_PROPANE) * ...
							ENTHALPY_PROPANE * GJ_PER_KJ * VALUE_C3H6_FUEL;
	but_val = P_butane * G_PER_KT * (1 / MOLMASS_BUTANE) * ...
						ENTHALPY_BUTANE * GJ_PER_KJ * VALUE_C4H8_FUEL;
	value = prop_val + but_val;
end

function cost = cost_C02(fuels)
	p_h2 = fuels(1);				% [ kta ]
	p_methane = fuels(2);			% [ kta ]
	p_propane = fuels(3);			% [ kta ]
	p_butane = fuels(4);			% [ kta ]
	GJ_natural_gas = fuels(5);		% [ GJ ]
	gallons_2fuel_oil = fuels(6);	% [ gallons ]


	cost = 0; 
	
	
end

function cost = cost_steam(flowrates)
	cost = 0 
end
















