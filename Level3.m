clc; clear; close all; 




% CONSTANTS________________________________________________________________

% Plotting 
S1_MIN = 0.05;
S1_MAX = 0.95;
S1_POINTS = 4;

S2_MIN = 0.05;
S2_MAX = 0.95;
S2_POINTS = 4;

% Design params
FC2H6 = 200;	% [	kta ]

% Molar mass
M_C2H6 = 28.05;	% [ g / mol ]

% Unit conversions 
metricTons_to_kiloton = 1000;



% SYSTEM OF EQUARTIONS (EXTENT OF RXN)_____________________________________

A = @(s1, s2)...
	[s1-1	,s1		,s1+1;
     s2		,s2-1	,s2;
     1		,2		,1	];
b = [0;		0;		FC2H6];

% FLOWRATE FUNCTIONS_______________________________________________________

P_H2 = @(xi_1)			xi_1;
P_CH4 = @(xi_2)			xi_2;		% Verify with squad about this equation
P_C2H4 = @(xi_1, xi_3)	xi_1 - xi_3;
P_C3H8 = @(xi_2)		xi_2;
P_C4H10 = @(xi_3)		xi_3;

% VALIDITY EQUATIONS_______________________________________________________

% Makes sure the flowrates are physically possible
validSolution = @(p_h2, p_ch4, p_c2h4, p_c3h8, p_c4h10) ...
	

% SCRIPT___________________________________________________________________

% Iterates through each value of selectivities S1 and S2 to find the economic
% potential for different reaction conditions 


for s1 = linspace(S1_MIN, S1_MAX, S1_POINTS)
	for s2 = linspace(S2_MIN, S2_MAX, S2_POINTS)
		
		% Solve for extents of reaction
		xi = A(s1, s2) \ b; 
		disp("xi  is")
		xi 

		% Calculate the flow rates of each species
		p_h2 = P_H2(xi(1));
		p_ch4 = P_CH4(xi(2));
		p_c2h4 = P_C2H4(xi(1), xi(3));
		p_c3h8 = P_C3H8(xi(2));
		p_c4h10 = P_C4H10(xi(3));

	end 
end 