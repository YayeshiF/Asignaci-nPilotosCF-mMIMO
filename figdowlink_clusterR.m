
close all;
clear;

%% Definimos simulación setup

%Numero de Monte-Carlo setups
nbrOfSetups = 250;%133

%Numero de channel realizations per setup
nbrOfRealizations = 100;

%Numbero de APs 
L = 100;

%Antenas por AP
N = 4;

%Numero de usuarios
K = [100]; %40

%Longitud bloque coherencia
tau_c = 200;

%Longitud de pilotos
tau_p = 10;

%Desviación estándar angular en el modelo de dispersión local (radianes)
ASD_varphi = deg2rad(15);  %angulo azimuth
ASD_theta = deg2rad(15);   %angulo de elevación

%% Parametros de propagación
%Potencia total de transmisión ascendente por UE (mW)
p = 100;

%Potencia total de transmisión descendente por AP (mW)
rho_tot = 200;

%Variables almacen resultados
SE_PMMSE_DCC = zeros(sum(K),nbrOfSetups); %P-MMSE(DCC) nuevo
SE_PMMSE_DCC_old = zeros(sum(K),nbrOfSetups);  %P-MMSE(DCC) old
%% Recorrer setups con configuración nueva
for n = 1:nbrOfSetups
    
    %Mostrar evolución simulación
    disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]);
    
    %Generar una configuración con UEs y APs en ubicaciones aleatorias y
    %criterio de igual repetición de pilotos
    [gainOverNoisedB,R,pilotIndex,D,D_small] = generateSetupUnicastRep(L,K,N,tau_p,1,n,ASD_varphi,ASD_theta);
    
    %Generar realizaciones de canal con estimaciones y matrices de correlación de error de estimación
    [Hhat,H,B,C] = canal_cluster(R,nbrOfRealizations,L,K,N,tau_p,pilotIndex,p);
    

    %% Cell-Free Massive MIMO con DCC

    gainOverNoise = db2pow(gainOverNoisedB); %calculo SNR en dBs
    
    %Calcular la asignación de potencia en (6.36) para la precodificación distribuida
    rho_dist = zeros(L,sum(K));
    
    for l = 1:L
        
        %Extraer qué UEs son utilizadas por el AP l
        servedUEs = find(D(l,:)==1);
        
        %Calcular el denominador en (6.36)
        normalizationAPl = sum(sqrt(gainOverNoise(l,servedUEs)));
        
        for ind = 1:length(servedUEs)
            
            rho_dist(l,servedUEs(ind)) = rho_tot*sqrt(gainOverNoise(l,servedUEs(ind)))/normalizationAPl;
            

        end
        
    end
    
    %Calcular SE para el caso de P-MMSE (DCC) nuevo
        [SE_P_MMSE] ...
            = functionComputeSE_downlink_P_MMSE(Hhat,H,D,C,tau_c,tau_p,nbrOfRealizations,N,sum(K),L,p,gainOverNoisedB,rho_tot);
    
    %Guardar valores de SE
    SE_PMMSE_DCC(:,n) = SE_P_MMSE;
    
    
    %Limpiar los valores de las matrices
    clear Hhat H B C R;
    
end
%% Recorrer setups con configuración antigua
for n = 1:nbrOfSetups
    
    %Mostrar evolución simulación
    disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]);
    
    %Generar una configuración con UEs y APs en ubicaciones aleatorias
    [gainOverNoisedB,R,pilotIndex,D,D_small] = generateSetupUnicast(L,K,N,tau_p,1,n,ASD_varphi,ASD_theta);
    
    %Generar realizaciones de canal con estimaciones y matrices de correlación de error de estimación
    [Hhat,H,B,C] = canal_cluster(R,nbrOfRealizations,L,K,N,tau_p,pilotIndex,p);
    
    %% Cell-Free Massive MIMO con DCC
    
    %Calcular la asignación de potencia en (6.36) para la precodificación distribuida
    rho_dist = zeros(L,sum(K));
    gainOverNoise = db2pow(gainOverNoisedB);
   
    
    for l = 1:L
        
        %Extraer qué UEs son utilizadas por el AP l
        servedUEs = find(D(l,:)==1);
        
        %Calcular el denominador en (6.36)
        normalizationAPl = sum(sqrt(gainOverNoise(l,servedUEs)));
        
        for ind = 1:length(servedUEs)
            
            rho_dist(l,servedUEs(ind)) = rho_tot*sqrt(gainOverNoise(l,servedUEs(ind)))/normalizationAPl;
            

        end
        
    end
    
    %Calcular SE para el caso de P-MMSE (DCC) antiguo
        [SE_P_MMSE] ...
            = functionComputeSE_downlink_P_MMSE(Hhat,H,D,C,tau_c,tau_p,nbrOfRealizations,N,sum(K),L,p,gainOverNoisedB,rho_tot);
    
    %Guardar valores de SE
    SE_PMMSE_DCC_old(:,n) = SE_P_MMSE;
    
    %Limpiar los valores de las matrices
    clear Hhat H B C R;
    
end


%% Graficar resultados de la simulación
% Graficar figura SE
figure;
title 'Comparación de SE por usuario'
hold on; box on;
set(gca,'fontsize',16);

plot(sort(SE_PMMSE_DCC(:)),linspace(0,1,sum(K)*nbrOfSetups),'k:','LineWidth',2); %configuración nueva
plot(sort(SE_PMMSE_DCC_old(:)),linspace(0,1,sum(K)*nbrOfSetups),'r:','LineWidth',2);%configuración old

xlabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
ylabel('CDF','Interpreter','Latex');
legend({'Algoritmo mejorado','Algoritmo base'},'Interpreter','Latex','Location','SouthEast');
xlim([0 12]);
grid on;
grid minor;

figure;
title 'Comparación de SE total en la distribución'
hold on; box on;
set(gca,'fontsize',16);

plot(sort(sum(SE_PMMSE_DCC)),linspace(0,1,nbrOfSetups),'k:','LineWidth',2); %configuración nueva
plot(sort(sum(SE_PMMSE_DCC_old)),linspace(0,1,nbrOfSetups),'r:','LineWidth',2);%configuración old

xlabel('Sum spectral efficiency [bit/s/Hz]','Interpreter','Latex');
ylabel('CDF','Interpreter','Latex');
legend({'Algoritmo actualizado','Algoritmo base'},'Interpreter','Latex','Location','SouthEast');
ylim([0 1]);
grid on;

