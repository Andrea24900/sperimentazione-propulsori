clear all
close all
clc

linewidth=1.2;
fontsize=13;

load errori_stat.mat

% Max value of long dataset: 1449.917970  
% Min value of long dataset: 931.352290 

%% Tabella di calibrazione statica

tensione_cal=[0.786;1.791;2.431;2.784;3.158;3.551;3.963;4.395;4.844;5.311;5.793;...
    6.290;6.800;7.326;7.866;8.418;8.979;9.549;10.124;10.704;11.286]; %[mV]
temperatura_cal=[400;600;700;750;800;850;900;950;1000;1050;1100;1150;1200;1250;...
    1300;1350;1400;1450;1500;1550;1600];

%legge di calibrazione
cal_law=figure(1);
plot(tensione_cal,temperatura_cal,'-*','linewidth',linewidth)
grid on
xlabel('$V$ [mV]','Interpreter','latex','FontSize',fontsize)
ylabel('$T$ [$^o$C]','Interpreter','latex','FontSize',fontsize)
exportgraphics(cal_law,'calibration.png','Resolution',600)
ordine_max=4;


for j=1

    %% 5 tentativi di legge interpolante per dataset corto
    % Min value of short dataset: 953.745910->valore 950
    % Max value of short dataset: 1193.110960->1200
    tensione_cal_short=tensione_cal(9-j:12+j);
    temperatura_cal_short=temperatura_cal(9-j:12+j);
    primo_valore=temperatura_cal_short(1);
    ultimo_valore=temperatura_cal_short(end);
    % legge lineare 
    coefficienti_1_short=polyfit(tensione_cal_short,temperatura_cal_short,1);
    % legge quadratica
    coefficienti_2_short=polyfit(tensione_cal_short,temperatura_cal_short,2);
    % legge cubica
    coefficienti_3_short=polyfit(tensione_cal_short,temperatura_cal_short,3);
    % legge quartica
    coefficienti_4_short=polyfit(tensione_cal_short,temperatura_cal_short,4);
    % legge quinto ordine
    %coefficienti_5_short=polyfit(tensione_cal_short,temperatura_cal_short,5);
    %valutare le tensioni di calibrazione con le 4 leggi
    stima_T_short=zeros(length(temperatura_cal_short),ordine_max);
    
    stima_T_short(:,1)=polyval(coefficienti_1_short,tensione_cal_short);
    stima_T_short(:,2)=polyval(coefficienti_2_short,tensione_cal_short);
    stima_T_short(:,3)=polyval(coefficienti_3_short,tensione_cal_short);
    stima_T_short(:,4)=polyval(coefficienti_4_short,tensione_cal_short);
    %stima_T_short(:,5)=polyval(coefficienti_5_short,tensione_cal_short);
    
    errore_T_short=zeros(size(stima_T_short));
    errore_T_short_quadratico=zeros(size(stima_T_short));
    %valutare gli errori intesi come err=temperatura_cal-stima_T_n
    for i=1:ordine_max
    errore_T_short(:,i)=stima_T_short(:,i)-temperatura_cal_short;
    errore_T_short_quadratico(:,i)=errore_T_short(:,i).^2;

        stima_errore_int(i)=sqrt(sum(errore_T_short_quadratico(:,i))/(length(temperatura_cal_short)-(i+1)));
        %scelta di t95 per i 4 casi
    
        nu=length(temperatura_cal_short)-(i+1);
        t95(i)=tstudent(nu);
        errore_sistematico(i)=t95(i).*stima_errore_int(i);
    
    end
    
    [errore_sistematico_minimo,ordine_polinomio]=min(errore_sistematico);
    
    vettore_errori_sistematici(j,:)=[errore_sistematico_minimo,ordine_polinomio,primo_valore,ultimo_valore];
    
    % Calcolo di R^2 per tutti gli ordini
    
    for i=1:ordine_max
        for k=1:length(temperatura_cal_short)
            err(k)=stima_T_short(k,i)-mean(temperatura_cal_short);
            err_square(k)=err(k).^2;
            den(k)=temperatura_cal_short(k)-mean(temperatura_cal_short);
            den_square(k)=den(k).^2;
        end
        NUM=sum(err_square);
        DEN=sum(den_square);
        R_square=NUM/DEN;
        R_square_vector(i)=R_square;
    end
    Rsquare_matrix(j,:)=R_square_vector;
end
[migliore_soluzione,indice]=min(vettore_errori_sistematici(:,1));

ordine_migliore_soluzione=vettore_errori_sistematici(indice,2);


errore_sist_short=migliore_soluzione;

err_min_fig=figure(2);

sist=plot(errore_sistematico,'*','LineWidth',linewidth);
yyaxis left
grid on
hold on

best=plot(stima_errore_int,'x','LineWidth',1.5,'color',[0.4660 0.6740 0.1880],'MarkerSize',12);
set(gca,'XLim',[0 5])
set(gca,'xTick',0:5)
xlabel('Ordine del modello','Interpreter','latex','fontsize',fontsize)
ylabel('Errore $[^oC]$','Interpreter','latex','fontsize',fontsize)
yyaxis right
plot(t95,'o','LineWidth',linewidth)

ylabel('$t_{\textit{95}}$','Interpreter','latex','fontsize',fontsize)
legend([best sist],'Errore di regressione','Errore sistematico','interpreter','latex','fontsize',fontsize,'location','east')
exportgraphics(err_min_fig,'err_sist_short_test.png','Resolution',600)

%% Ricalcolo legge metrologica best case
j=indice;
tensione_cal_short=tensione_cal(9-j:12+j);
temperatura_cal_short=temperatura_cal(9-j:12+j);
coefficienti_short=polyfit(tensione_cal_short,temperatura_cal_short,ordine_migliore_soluzione);

%% Propagazione degli errori

err_tot_short_NOQUANT=sqrt(errore_sist_short^2+err_stat_short^2);

save errore_sist_short.mat errore_sist_short coefficienti_short...
    tensione_cal_short temperatura_cal_short err_tot_short_NOQUANT