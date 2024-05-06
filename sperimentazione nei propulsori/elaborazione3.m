clear all
close all
clc

linewidth=1.2;
fontsize=13;

%% Definizione del tempo di campionamento

Ts=0.01; %[s]
fs=1/Ts;

%% Apertura dei files ed estrazione dei dati

short_vector=readmatrix("prova1600campioni.txt")';

long_vector=readmatrix("prova10000campioni.txt")';

load errori_stat.mat
load errori_sist.mat



%% Calcolo errore quantizzazione

numero_bit=12;

campo_inf=0; %[V]
campo_sup=10; %[V]

delta_V=campo_sup-campo_inf;

err_quant=delta_V/(2^(numero_bit+1));

%% Calcolo errore totale

err_tot_short=sqrt(err_quant^2+errore_sist^2+err_stat_short^2);

err_tot_long=sqrt(err_quant^2+errore_sist^2+err_stat_long^2);

%% Calcolo perdite irraggiamento




%% Calcolo errore totale