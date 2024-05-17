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

load errore_sist_short.mat

load errore_sist_long.mat



%% Calcolo errore quantizzazione

numero_bit=12;

campo_inf=0; %[V]
campo_sup=10; %[V]

delta_V=campo_sup-campo_inf;

err_quant_V=delta_V/(2^(numero_bit+1)); % in V

err_quant_mV=err_quant_V*1000; %in mV

% la relazione metrologica è del quarto ordine, quindi la propagazione
% degli errori dipende dal valore puntuale di tensione

%come valore di tensione di riferimento si sceglie il fondoscala in
%tensione della legge di calibrazione: questo è il più alto valore di
%tensione che compare nella legge metrologica e quindi produce la più alta
%derivata di sensibilità: l'errore di quantizzazione in temperatura risulta
%conservativo

tensione_FS_short=tensione_cal_short(end);
tensione_FS_long=tensione_cal_long(end);

% le due leggi sono del quarto ordine y=ax^4+bx^3+cx^2+dx+e
% la derivata è y'=4ax^3+3bx^2+2cx+d

derivata_legge_short=@(tens) 4*coefficienti_short(1).*tens.^3+3*coefficienti_short(2).*tens.^2+...
    2*coefficienti_short(3).*tens+coefficienti_short(4);

derivata_legge_long=@(tens) 4*coefficienti_long(1).*tens.^3+3*coefficienti_long(2).*tens.^2+...
    2*coefficienti_long(3).*tens+coefficienti_long(4);

err_quant_short=derivata_legge_short(tensione_FS_short)*err_quant_mV/100;
err_quant_long=derivata_legge_long(tensione_FS_long)*err_quant_mV/100;

%% Calcolo errore totale
err_tot_short=sqrt(err_quant_short^2+errore_sist_short^2+err_stat_short^2);

err_tot_long=sqrt(err_quant_long^2+errore_sist_long^2+err_stat_long^2);

%% Calcolo perdite irraggiamento





%% Calcolo errore totale