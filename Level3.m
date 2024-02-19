clc; clear; close all; 


% CONSTANTS________________________________________________________________

S1_MIN = 0.05;
S1_MAX = 0.95;
S1_POINTS = 50;

S2_MIN = 0.05;
S2_MAX = 0.95;
S2_POINTS = 50;


% SCRIPT___________________________________________________________________

% Code iterates through each value of S1 and S2 to find the economic
% potential for different reaction conditions 



for s1 = linspace(S1_MIN, S1_MAX, S1_POINTS)
	for s2 = linspace(S2_MIN, S2_MAX, S2_POINTS)
		fprintf("%f\n", s1)
	end 
end 