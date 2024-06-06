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

%% Individuazione massimi e minimi

max_short=max(short_vector);
max_long=max(long_vector);

min_short=min(short_vector);
min_long=min(long_vector);

fprintf("Max value of short dataset: %f \n",max_short)
fprintf("Max value of long dataset: %f \n",max_long)
fprintf("Min value of short dataset: %f \n",min_short)
fprintf("Min value of long dataset: %f \n",min_long)

%% Definizione classi

n_classes=10;

classes_short=linspace((min_short-0.5e-7),(max_short+0.5e-7),n_classes+1)';

classes_long=linspace((min_long-0.5e-7),(max_long+0.5e-7),n_classes+1)';

%%  Suddivisione dati nelle classi mediante istogramma

% hist_short=figure(1);
% h_short=histogram(short_vector,classes_short);
% xlabel("$T\, [^oC]$","Interpreter","latex","FontSize",fontsize)
% ylabel("Occorrenze","Interpreter","latex","FontSize",fontsize)
% %title("Serie corta",'interpreter','latex',"FontSize",fontsize)
% exportgraphics(hist_short,'istogrammashort.png','Resolution',600)
% hist_long=figure(2);
% h_long=histogram(long_vector,classes_long);
% xlabel("$T\, [^oC]$","Interpreter","latex","FontSize",fontsize)
% ylabel("Occorrenze","Interpreter","latex","FontSize",fontsize)
% %title("Serie lunga",'interpreter','latex',"FontSize",fontsize)
% exportgraphics(hist_long,'istogrammalong.png','Resolution',600)
%% Conteggio dati nelle classi

N_short=histcounts(short_vector,classes_short)';

N_long=histcounts(long_vector,classes_long)';

%% Calcolo frequenze relative

rel_freq_short=N_short/length(short_vector);

rel_freq_long=N_long/length(long_vector);

% rel_freq=figure(3);
% plot(rel_freq_short,'-*','LineWidth',linewidth);
% grid on
% hold on
% plot(rel_freq_long,'-o','LineWidth',linewidth);
% xlabel("$Classi$","Interpreter","latex","FontSize",fontsize)
% ylabel("$f$","Interpreter","latex","FontSize",fontsize)
% legend('Serie corta','Serie lunga','interpreter','latex','fontsize',fontsize)
% 
% %title("Frequenze relative per la serie corta",'interpreter','latex',"FontSize",fontsize)
% exportgraphics(rel_freq,'relboth.png','Resolution',600)
%% Calcolo frequenze cumulate

cum_freq_short=zeros(n_classes,1);
cum_freq_long=zeros(n_classes,1);

for i=1:n_classes
    cum_freq_short(i)=sum(rel_freq_short(1:i));
    cum_freq_long(i)=sum(rel_freq_long(1:i));
end

% cum_both=figure(5);
% plot(cum_freq_short,'*-','LineWidth',linewidth)
% grid on
% hold on
% plot(cum_freq_long,'o-','LineWidth',linewidth)
% xlabel("$Classi$","Interpreter","latex","FontSize",fontsize)
% ylabel("$F$","Interpreter","latex","FontSize",fontsize)
% %title("Frequenze cumulate per entrambe le serie",'interpreter','latex',"FontSize",fontsize)
% legend("Serie corta","Serie lunga",'interpreter','latex',"FontSize",fontsize,'location','southeast')
% exportgraphics(cum_both,'cumboth.png','Resolution',600)
cum_freq_rel_short=100*cum_freq_short/length(short_vector);
cum_freq_rel_long=100*cum_freq_long/length(long_vector);

%% Calcolo media, moda, mediana, skewness

mean_short=mean(short_vector)
mean_long=mean(long_vector)

median_short=median(short_vector)
median_long=median(long_vector)

[mode_short,freq_mode_short]=mode(short_vector);
[mode_long,freq_mode_long,count_long]=mode(long_vector);

skew_short=skewness_coefficient(short_vector)
skew_long=skewness_coefficient(long_vector)

%% Calcolo deviazione standard ed errore statistico

std_short=std(short_vector)
std_long=std(long_vector)

%per usare tinv bisogna prendere il livello di confidenza voluto, es. 95%,
%e aggiungere la metà del suo complemento a 1

t_95_short=tinv(0.975,length(short_vector)-1);
t_95_long=tinv(0.975,length(long_vector)-1);

sigmaTmean_short=std_short/sqrt(length(short_vector));
sigmaTmean_long=std_long/sqrt(length(long_vector));
err_stat_short=sigmaTmean_short*t_95_short;
err_stat_long=sigmaTmean_long*t_95_long;

save errori_stat.mat err_stat_short err_stat_long
%% Caratterizzazione temporale dei set

sampling_instants_short=[0:Ts:(length(short_vector)-1)*Ts];
sampling_instants_long=[0:Ts:(length(long_vector)-1)*Ts];

%% Plot dei due data set

time_short=figure(6);
plot(sampling_instants_short,short_vector,"LineWidth",linewidth);
grid on
xlabel("$t\, [s]$","Interpreter","latex","FontSize",fontsize)
ylabel("$T\, [^oC]$","Interpreter","latex","FontSize",fontsize)
%title("Time evolution of short data set","Interpreter","latex","FontSize",fontsize)
exportgraphics(time_short,'time_short.png','Resolution',600)
time_long=figure(7);
plot(sampling_instants_long,long_vector,"LineWidth",linewidth);
grid on
xlabel("$t\, [s]$","Interpreter","latex","FontSize",fontsize)
ylabel("$T\, [^oC]$","Interpreter","latex","FontSize",fontsize)
%
% title("Time evolution of long data set","Interpreter","latex","FontSize",fontsize)
exportgraphics(time_long,'time_long.png','Resolution',600)
%% Calcolo PSD dei dataset

% Dataset corto
% Note sulla finestratura e come viene usata in pwelch
% - l'argomento window scalare (M) indica di dividere i dataset in M
%   segmenti, su ciascuno viene poi usata una finestra di Hamming
% - l'argomento window vettoriale (lunghezza: L) divide i dataset in
%   segmenti di lunghezza L e su ciascuno applica la finestra specificata
%   in window mediante moltiplicazione del segmento per il vettore window

Number_of_windows = 10; %numero di finestre usate
L_window = floor(length(short_vector)/Number_of_windows);   % numero di valori nella finestra:
%   determinato dividendo la lunghezza del dataset per in numero di
%   finestre desiderate. Matlab usa 8 finestre di default. 
%   Aumentare il numero di finestre riduce la varianza e aumenta il bias.
%   Ridurre in numero di finestre riduce il bias e aumenta la varianza
window = hamming( L_window );   % hamming window with L_window values

n_overlap = floor(L_window/2);   % Overlap at 50% of window length: l'overlapping 
%   riduce il bias ma aumenta la correlazione

% Plot of the window
figure(8)
plot(window,"LineWidth",linewidth)
grid on
title('Window',"Interpreter","latex","FontSize",fontsize)

[ p_short, freq_psd ] = pwelch( short_vector, window, n_overlap, length(short_vector), fs, 'onesided' );

%      L'argomento 'onesided' calcola la PSD one sided e moltiplica la potenza a tutte le frequenze
    %      x2 al fine di conservare la potenza totale; ricordare opzione 'power' se mai servisse.
    %      Osservare che il fattore di scala è dovuto a
    %      G_uu(one-sided)=2*S_uu(two-sided)

% Plot:    
    freq_short=figure(9);
    plot( freq_psd, 10*log10(p_short), '-' ,"LineWidth",linewidth)
    % si usa 10*log10 invece che 20*log10 in quanto è insita l'operazione
    % di calcolo di sqrt(p_uu). Per riferimento vedere slide su "Power
    % Spectral Densities and Spectral Densities" di "Random+Processes".
    grid on;
    xlabel( 'Frequency [Hz]',"Interpreter","latex","FontSize",fontsize);
    ylabel( 'PSD [dB/Hz]' ,"Interpreter","latex","FontSize",fontsize);
 %   title('One-sided Power Spectral Density for short data set',...
 %       "Interpreter","latex","FontSize",fontsize)
 exportgraphics(freq_short,'freq_short.png','Resolution',600)


% Dataset lungo

Number_of_windows = 10; 
L_window = floor(length(long_vector)/Number_of_windows);   
window = hamming( L_window );   
n_overlap = floor(L_window/2);   



[ p_long, freq_psd ] = pwelch( long_vector, window, n_overlap, length(long_vector), fs, 'onesided' );

%      L'argomento 'onesided' calcola la PSD one sided e moltiplica la potenza a tutte le frequenze
    %      x2 al fine di conservare la potenza totale; ricordare opzione 'power' se mai servisse.
    %      Osservare che il fattore di scala è dovuto a
    %      G_uu(one-sided)=2*S_uu(two-sided)

% Plot:    
    freq_long=figure(10);
    plot( freq_psd, 10*log10(p_long), '-' ,"LineWidth",linewidth)
    grid on;
    xlabel( 'Frequency [Hz]' ,"Interpreter","latex","FontSize",fontsize)
    ylabel( 'PSD [dB/Hz]' ,"Interpreter","latex","FontSize",fontsize)
    exportgraphics(freq_long,'freq_long.png','Resolution',600)
    % title('One-sided Power Spectral Density for long data set',"Interpreter","latex","FontSize",fontsize)






















%% FUNZIONI----------------------------------------------------------------

function coeff=skewness_coefficient(x_vector)
    x_mean=mean(x_vector);
    x_std=std(x_vector);
    coeff=(1/length(x_vector)/x_std^3)*sum((x_vector(:)-x_mean).^3);
end


