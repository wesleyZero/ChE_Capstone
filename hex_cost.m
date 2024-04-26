clc; clear; close all; 



% area is in ft^2 
A = [ 1.1 * 10^4, 3.4 * 10^5, 10^4 ];
cost = 0;
toggle_h2 = true;


for	a = A
	if a > 25000 
			if toggle_h2 
				hex_i_cost = (1800 / 280) * 101.3 * (25000^0.65) * 3.75;
			else
				hex_i_cost = (1800 / 280) * 101.3 * 25000^0.65;
			end
			cost = cost + hex_i_cost * floor(a / 25000);
			disp("ITS TOO BIG")
			fraction_of_a = (a/25000) -  floor(a / 25000);
			a = fraction_of_a * a;
	end
		
	if toggle_h2 
		hex_i_cost = (1800 / 280) * 101.3 * a^0.65 * 3.75;
	else
		hex_i_cost = (1800 / 280) * 101.3 * a^0.65;
	end
	cost = cost + hex_i_cost
end

disp("the total cost ")
cost

disp("________________________")
toggle_h2 = false;
for	a = A
	if a > 25000 
			if toggle_h2 
				hex_i_cost = (1800 / 280) * 101.3 * 25000^0.65 * 3.75;
			else
				hex_i_cost = (1800 / 280) * 101.3 * 25000^0.65;
			end
			cost = cost + hex_i_cost;
			disp("ITS TOO BIG: old and new a ")
% 			a
% 			a = ((a/25000) -  floor(a / 25000)) * a
	end
		
	if toggle_h2 
		hex_i_cost = (1800 / 280) * 101.3 * a^0.65 * 3.75;
	else
		hex_i_cost = (1800 / 280) * 101.3 * a^0.65;
	end
	cost = cost + hex_i_cost
end

disp("the total cost ")
cost