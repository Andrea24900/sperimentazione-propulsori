clc
close all
clear

%% dati
D = 42; %mm
d = 9.94; %mm
qm_n = 96;
gamma = 1.4;
beta = d/D;
p2 = 4;   

k = 25*D/1e4;
qm = qm_n * 1.29 /3600;
mu = 1.8e-5;
l1 = 25.4/D;
l2 = l1;
Re = qm/(mu*D)

C = 0.5959 + 0.0312 * beta^2.1 - 0.184 * beta^8 + 0.0029 * beta^2.5 * (1e6/Re)^0.75 + ...
    0.09 * l1 * beta^4/(1-beta^4) - 0.0337 * l2 * beta^3;

%% iterazioni
C2 = 1 - (0.41 + 0.35*beta^4)/gamma;
C1 = C/sqrt(1-beta^4) * pi/4 * d^2;
C3 = qm/(C1/sqrt(2))

epsilon = @(p1) -C3 + (1 + C2) * sqrt(p1 * p2 - p1.^2) - C2*p2*sqrt(p2./p1 - 1)

p1_sol = fzero(epsilon, 5)