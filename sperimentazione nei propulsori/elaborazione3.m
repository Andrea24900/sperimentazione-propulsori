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

tensione_FS_short=tensione_cal_short(1);
tensione_FS_long=tensione_cal_long(1);

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

%% Soluzione alternativa per quantizzazione
T_mean_short = mean(short_vector);
T_mean_long = mean(long_vector);

inverted_short=@(V) T_mean_short-(coefficienti_short(1).*V.^4+coefficienti_short(2).*V.^3+...
    coefficienti_short(3).*V.^2+coefficienti_short(4).*V+coefficienti_short(5));

inverted_long=@(V) T_mean_long-(coefficienti_long(1).*V.^4+coefficienti_long(2).*V.^3+...
    coefficienti_long(3).*V.^2+coefficienti_long(4).*V+coefficienti_long(5));

V_0_short=5.79;

V_0_long=6.29;

V_short=fzero(inverted_short,V_0_short);
V_long=fzero(inverted_long,V_0_long);

err_quant_short_alt=derivata_legge_short(V_short)*err_quant_mV/100;
err_quant_long_alt=derivata_legge_long(V_long)*err_quant_mV/100;

%% Calcolo errore totale alt
err_tot_short_alt=sqrt(err_quant_short_alt^2+errore_sist_short^2+err_stat_short^2);

err_tot_long_alt=sqrt(err_quant_long_alt^2+errore_sist_long^2+err_stat_long^2);

diff_short=(err_tot_short_alt-err_tot_short)/err_tot_short*100;
diff_long=(err_tot_long_alt-err_tot_long)/err_tot_long*100;
tensione_vect=linspace (3,8,100);


plot_confronto=figure(10);
yyaxis left
met_1=plot(tensione_FS_short,derivata_legge_short(tensione_FS_short),'o','color','green','LineWidth',linewidth,'MarkerSize',8);

grid on
hold on 


met_2=plot(V_short,derivata_legge_short(V_short),'square','LineWidth',linewidth,'Color','red','MarkerSize',8);

plot(tensione_vect,derivata_legge_short(tensione_vect),'-','lineWidth',linewidth,'MarkerSize',12); 
xlabel('$V$ [mV]','FontSize',fontsize,'Interpreter','latex')
ylabel('$dT/dV$ [$^o$C/mV]','FontSize',fontsize,'Interpreter','latex')
yyaxis right 
plot(tensione_vect,-inverted_short(tensione_vect)+T_mean_short,'lineWidth',linewidth)
ylabel('$T(V)$ [$^o$C]','FontSize',fontsize,'Interpreter','latex')
legend('$dT/dV$ 1$^o$ metodo','$dT/dV$ 2$^o$ metodo','interpreter','latex','fontsize',fontsize,'location','north')
exportgraphics(plot_confronto,'confronto.png','Resolution',600)

%% Calcolo perdite irraggiamento





%% Calcolo errore totale