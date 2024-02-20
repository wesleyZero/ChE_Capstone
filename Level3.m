clc; clear; close all; 

% CONSTANTS________________________________________________________________

% Plotting 
S1_MIN = 0.05;
S1_MAX = 0.95;
S1_POINTS = 40;

S2_MIN = 0.05;
S2_MAX = 0.95;
S2_POINTS = 40;

INVALID_FLOW = 0;

cont_plt_opt = { ...
	'S_1 Selectivity', ...
	'S_2 Selectivity', ...
	'Ethylene Flowrate [kta]',...
	'P_ethylene_VS_S1_S2.jpg'};

% Design params
FC2H6 = 200;	% [	kta ]

% Molar mass
M_C2H6 = 28.05;	% [ g / mol ]

% Unit conversions 
metricTons_to_kiloton = 1000;

% Economic | Chemicals
value_ethane = 200;		% [ $ / MT ]
value_ethylene = 900;	% [ $ / MT ]

% Economic | Fuel
value_H2_fuel = 3;		% [ $ / GJ ]
value_CH4_fuel = 3;		% [ $ / GJ ]
value_C3H6_fuel = 3;	% [ $ / GJ ]
value_C4H8_fuel = 3;	% [ $ / GJ ]
value_NatGas_fuel = 3;	% [ $ / GJ ]
value_2FuelOil = 4.5;	% [ $ / US Gallon ]

% Economics | Enviormental




% SYSTEM OF EQUARTIONS (EXTENT OF RXN)_____________________________________

A = @(s1, s2)...
	[s1-1	,s1		,s1+1;
     s2		,s2-1	,s2;
     1		,2		,1	];
b = [0;		0;		FC2H6];

% FLOWRATE FUNCTIONS_______________________________________________________

P_H2 = @(xi_1)			xi_1;
P_CH4 = @(xi_2)			xi_2;
P_C2H4 = @(xi_1, xi_3)	xi_1 - xi_3;
P_C3H8 = @(xi_2)		xi_2;
P_C4H10 = @(xi_3)		xi_3;

% VALIDATION FUNCTIONS_____________________________________________________

flowrates_valid = @( flowrates ) ...
			all(flowrates >= 0);

% SCRIPT___________________________________________________________________

% Iterates through each value of selectivities S1 and S2 to find the economic
% potential for different reaction conditions 
s1_domain = linspace(S1_MIN, S1_MAX, S1_POINTS);
s2_domain = linspace(S2_MIN, S2_MAX, S2_POINTS);
[s1_mesh, s2_mesh] = meshgrid(s1_domain, s2_domain);
ethylene_flowrates = s1_mesh + s2_mesh;			% Just placeholder values

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
		else
			ethylene_flowrates(i) = INVALID_FLOW;
		end 
		i = i + 1;
	end 
end 

plot_contour(s1_mesh, s2_mesh, ethylene_flowrates, cont_plt_opt);



% HELPER FUNCTIONS | PLOTTING______________________________________________

function z = plot_contour(x, y, z, options)
	% Unpack options 
	x_label = options{1};
	y_label = options{2};
	plt_title = options{3};
	plt_saveName = options{4};

	hold on 
	contourf(x, y, z, "ShowText","on");
	xlabel(x_label);
	ylabel(y_label);
	title(plt_title);
	saveas(gcf, plt_saveName);
	hold off
end

function value = ethane_value(P_ethane)
	% inputs
	%		P_ethane [kta]
	% output
	%		value [$ USD]

	
end 
















