clear all
close all
clc

linewidth=1.2;
fontsize=13;

load errori_stat.mat


%% Tabella di calibrazione statica

tensione_cal=[0.786;1.791;2.431;2.784;3.158;3.551;3.963;4.395;4.844;5.311;5.793;...
    6.290;6.800;7.326;7.866;8.418;8.979;9.549;10.124;10.704;11.286]; %[mV]
temperatura_cal=[400;600;700;750;800;850;900;950;1000;1050;1100;1150;1200;1250;...
    1300;1350;1400;1450;1500;1550;1600];

%legge di calibrazione
figure(1)
plot(tensione_cal,temperatura_cal,'-*','linewidth',linewidth)
grid on
xlabel('$V$ [mV]','Interpreter','latex','FontSize',fontsize)
ylabel('$T$ [$^o$C]','Interpreter','latex','FontSize',fontsize)

ordine_max=4;


for j=1:4

    %% 4 tentativi di legge interpolante per dataset lungo
    % Max value of long dataset: 1449.917970 -> 1450 -> 18 
    % Min value of long dataset: 931.352290 -> 900 -> 7

    tensione_cal_short=tensione_cal(8-j:17+j);
    temperatura_cal_short=temperatura_cal(8-j:17+j);
    % legge lineare 
    coefficienti_1_short=polyfit(tensione_cal_short,temperatura_cal_short,1);
    % legge quadratica
    coefficienti_2_short=polyfit(tensione_cal_short,temperatura_cal_short,2);
    % legge cubica
    coefficienti_3_short=polyfit(tensione_cal_short,temperatura_cal_short,3);
    % legge quartica
    coefficienti_4_short=polyfit(tensione_cal_short,temperatura_cal_short,4);
   
    %valutare le tensioni di calibrazione con le 4 leggi
    stima_T_short=zeros(length(temperatura_cal_short),ordine_max);
    
    stima_T_short(:,1)=polyval(coefficienti_1_short,tensione_cal_short);
    stima_T_short(:,2)=polyval(coefficienti_2_short,tensione_cal_short);
    stima_T_short(:,3)=polyval(coefficienti_3_short,tensione_cal_short);
    stima_T_short(:,4)=polyval(coefficienti_4_short,tensione_cal_short);
   
    
    errore_T_short=zeros(size(stima_T_short));
    errore_T_short_quadratico=zeros(size(stima_T_short));
    %valutare gli errori intesi come err=temperatura_cal-stima_T_n
    for i=1:ordine_max
    errore_T_short(:,i)=stima_T_short(:,i)-temperatura_cal_short;
    errore_T_short_quadratico(:,i)=errore_T_short(:,i).^2;
    %errore_T_short_quadratico_relFS(:,i)=errore_T_short_quadratico(:,i)./(max(temperatura_cal_short)^2)*100;
    end
    % figure(2)
    % plot(tensione_cal_short,errore_T_short_quadratico_relFS(:,1),'-*','LineWidth',linewidth)
    % grid on
    % hold on
    % for i=2:4
    % plot(tensione_cal_short,errore_T_short_quadratico_relFS(:,i),'-*','LineWidth',linewidth)
    % end
    % legend('Ordine 1','Ordine 2','Ordine 3','Ordine 4','interpreter','latex','fontsize',fontsize)
    % xlabel('$V$ [mV]','Interpreter','latex','FontSize',fontsize)
    % ylabel('$e^2_{\%FS}$','Interpreter','latex','FontSize',fontsize)
    
    for i=1:ordine_max
        stima_errore_int(i)=sqrt(sum(errore_T_short_quadratico(:,i))/(length(temperatura_cal_short)-(i+1)));
        %scelta di t95 per i 4 casi
    
        nu=length(temperatura_cal_short)-(i+1);
        t95(i)=tstudent(nu);
        errore_sistematico(i)=t95(i).*stima_errore_int(i);
    
    end
    
    [errore_sistematico_minimo,ordine_polinomio]=min(errore_sistematico);
    
    vettore_errori_sistematici(j,:)=[errore_sistematico_minimo,ordine_polinomio];
    
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


errore_sist_long=migliore_soluzione;

%% Propagazione degli errori

err_tot_short=sqrt(errore_sist_long^2+err_stat_long^2);

save errori_sist.mat errore_sist_long