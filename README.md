# Capstone I 

This was the first capstone project, of a two-part series for my chemical engineering senior design class. I learned a lot on this first project and as a result, I cleaned up the code massively when doing the second project. All of the code was in a single script, as well as the the code being not very modular. Since the capstone II project has the same concepts, but with a better code base. Please look at my [Capstone II project](https://github.com/wesleyZero/capstone_II/tree/main) to evaluate my coding abilities. For this project, I would like to simply write what I learned to do better in the second project. 

## One cool thing I would like to point out on this project. 

It was really cool numerically approximating a system of differential equations to solve this problem. Just one highlight that makes this project different from the second capstone project. 

[reactionODEs](https://github.com/wesleyZero/ChE_Capstone/blob/newBranch/Level3.m#L1404-L1456)
```matlab
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
```



## What I learned (not) to do, and do in the second project 

1. I created a huge massive script, I (once again) vastly underestimated the scope of the project (how large it was going to become) so it became a massive mess. 

2. I put all of the files in the same directory (why???). I learned to put images in the /img folder and pdfs in the /pdf folder.

3. I did not use structures, and instead [used a TON of constants](https://github.com/wesleyZero/ChE_Capstone/blob/newBranch/Level3.m#L58-L370). At my internship at KARL STORZ I learned (from reading others code) that structures make things so much cleaner in MATLAB code. Why don't people just create classes and objects in MATLAB? I have no idea. The MATLAB-way seem to be to make a bunch of structures. 

4. Use [way too many lambda functions (anonymous functions)](https://github.com/wesleyZero/ChE_Capstone/blob/newBranch/Level3.m#L372-L444). I clearly like lambda functions and they are very nice for making some functions more readable. But the way I used them in this huge script wasn't making things easier to understand, it was making it harder. 

5. Did not make modular code. [Look at this disaster I created!!! omg](https://github.com/wesleyZero/ChE_Capstone/blob/newBranch/Level3.m#L640-L1110) This really should have been a collection of functions in a few different scripts. I kept digging a deeper grave instead of just putting parts of this monster into other functions. It's completely unreadable and impossible to maintain. I feel ashamed for creating it, but sometimes the best way to learn a lesson is to torture yourself by doing it the wrong way.... Never forget
