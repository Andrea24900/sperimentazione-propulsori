clear all
close all
clc

%% dati
D = 42e-3; %m
d = 9.94e-3; %m
qv_n = 96; %Nm^3/h
gamma = 1.4;
beta = d/D;
p2 = 4;   %bar
p_2= p2*101325; %Pa

q_m = qv_n * 1.29 /3600;

L1 = 25.4/42; %con D in mm
L2 = L1;
%interpolazione per trovare mu
mu_300=1.81e-5+(0.05e-5)/10*27; % mu a 300 K
ReD = 4*q_m/(pi*mu_300*D);

C = 0.5959 + 0.0312 * beta^2.1 - 0.184 * beta^8 + 0.0029 * beta^2.5 * (1e6/ReD)^0.75 + ...
    0.039 * L1 * beta^4/(1-beta^4) - 0.0337 * L2 * beta^3;

%calcolo di rho 2
R=287; %J/kg/K
T_2=300; %K
rho_2=p_2/(R*T_2);

%% iterazioni
C1 = (0.41 + 0.35*beta^4)/gamma;
C2 = C/sqrt(1-beta^4) * pi/4 * d^2 * sqrt(2*rho_2);
C3 = C2/sqrt(p_2);

epsilon = @(p_1) q_m-C3 * (1 - C1 +C1.*p_2./p_1) .* sqrt(-p_1.* p_2 + p_1.^2) ;

p1_sol = fzero(epsilon, 5e5);

delta_p=p1_sol-p_2;