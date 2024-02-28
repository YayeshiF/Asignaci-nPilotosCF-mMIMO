function [gainOverNoisedB,R,pilotIndex,D,APpositions,UEpositions,distances, masterAPs] = generateSetupUnicastRep(L,K,N,tau_p,nbrOfSetups,seed,ASD_varphi,ASD_theta)

%% Define simulation setup

%Establecer el número de semilla si se especifica un valor diferente de cero
if (nargin>5)&&(seed>0)
    rng(seed)
end

%Tamaño área cobertura(as a square with wrap-around)
squareLength = 1000; %meter

%Tamaño del cluster
C = length(K); %numero de clusters
squareLengthCluster = 150*ones(1,C); %tamaño en metros

%Ancho de banda (Hz)
B = 20e6;

%Factor de ruido (in dB)
noiseFigure = 7;

%Calcular la potencia de ruido (in dBm)
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;

%Parámetros de pérdida de trayectoria (5.42)
alpha = 36.7;
constantTerm = -30.5;

%Desviación estándar del desvanecimiento de sombra (5.43)
sigma_sf = 4;

%Distancia de decorrelación del desvanecimiento de sombra (5.43)
decorr = 9;

%Diferencia altura AP y UE (in meters)
distanceVertical = 10;

%Definir espacio antena (in number of wavelengths)
antennaSpacing = 1/2; %longitud de onda

% Número de veces que se repite cada piloto
rep = sum(K) / tau_p; 

%Variables almacenan resultados
gainOverNoisedB = zeros(L,sum(K),nbrOfSetups);
R = zeros(N,N,L,sum(K),nbrOfSetups);
distances = zeros(L,sum(K),nbrOfSetups);
pilotIndex = zeros(sum(K),nbrOfSetups);
D = zeros(L,sum(K),nbrOfSetups);
masterAPs = zeros(sum(K),1); %indice de MasterAP para cada UE

%% Recorrer setups
for n = 1:nbrOfSetups   
    %Ubicaciones aleatorias de AP con distribución uniforme
    APpositions = (rand(L,1) + 1i*rand(L,1)) * squareLength;   
    %Para calcular ubicaciones UEs
    UEpositions = zeros(sum(K),1);   
    %Calcular ubicaciones alternativas de AP con wrap around
    wrapHorizontal = repmat([-squareLength 0 squareLength],[3 1]);
    wrapVertical = wrapHorizontal';
    wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
    APpositionsWrapped = repmat(APpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[L 1]);   
    %Guardar matriz de correlación de Shadowing
    shadowCorrMatrix = sigma_sf^2*ones(sum(K),sum(K));
    shadowAPrealizations = zeros(sum(K),L); 
    
    %Añadir UEs en clusters
    %Generar centro ubicación de clusters
    Clusterpositions = (rand(C,1) + 1i*rand(C,1)) * squareLength;
    % Permutar todos los usuarios
    allUsers = randperm(sum(K));

    % Inicializar el índice para recorrer el vector permutado de usuarios
    userIndex = 1;

    for c = 1:C
        cluster_users = allUsers(userIndex:userIndex + K(c) - 1);

        for k = 1:K(c)
            % Generar ubicación de UE aleatoria en el área
            UEpositions(cluster_users(k)) = Clusterpositions(c) + (rand(1, 1) + 1i * rand(1, 1)) * squareLengthCluster(c);
        end

        % Actualizar el índice de usuarios
        userIndex = userIndex + K(c);
    end


    UEpositions = UEpositions.';
    UEpositions = UEpositions(:);
    plot(UEpositions, 'b*') %graficar posiciones de UEs/clusters

    for k = 1:sum(K)
        %Calcular distancias asumiendo que los AP están a 10 m por encima de los UEs
        [distanceAPstoUE,whichpos] = min(abs(APpositionsWrapped - repmat(UEpositions(k),size(APpositionsWrapped))),[],2);
        distances(:,k,n) = sqrt(distanceVertical^2+distanceAPstoUE.^2);           
        %Si no es el primer UE
        if k-1>0                
            %Calcular distancias desde el nuevo UE prospectivo hasta todos los demás UEs
            shortestDistances = zeros(k-1,1);                
            for i = 1:k-1
                shortestDistances(i) = min(abs(UEpositions(k) - UEpositions(i) + wrapLocations));
            end                
            %Calcular la media condicional y la desviación estándar necesarias para
            % obtener las nuevas realizaciones de desvanecimiento de sombra, cuando los anteriores
            % realizaciones de desvanecimiento de sombra de los UEs ya se han generado.
            %This computation is based on Theorem 10.2 in "Fundamentals of
            %Statistical Signal Processing: Estimation Theory" by S. Kay
            newcolumn = sigma_sf^2*2.^(-shortestDistances/decorr);
            term1 = newcolumn'/shadowCorrMatrix(1:k-1,1:k-1);
            meanvalues = term1*shadowAPrealizations(1:k-1,:);
            stdvalue = sqrt(sigma_sf^2 - term1*newcolumn);                
        else %Si es el primer UE             
            %Agregar el UE y comenzar a almacenar los valores de
            %correlación de desvanecimiento de shadowing
            meanvalues = 0;
            stdvalue = sigma_sf;
            newcolumn = [];                
        end            
        %Generar las realizaciones de desvanecimiento de sombra
        shadowing = meanvalues + stdvalue*randn(1,L);
        %Calcular la ganancia del canal dividida por la potencia de ruido
        gainOverNoisedB(:,k,n) = constantTerm - alpha*log10(distances(:,k,n)) + shadowing' - noiseVariancedBm;                                  
        %Actualizar la matriz de correlación de desvanecimiento de sombra y almacenar las realizaciones
        shadowCorrMatrix(1:k-1,k) = newcolumn;
        shadowCorrMatrix(k,1:k-1) = newcolumn';
        shadowAPrealizations(k,:) = shadowing;  
        %Recorrer los APs
        for l = 1:L            
            %CCalcular el ángulo nominal entre el UE k y el AP l
            angletoUE_varphi = angle(UEpositions(k)-APpositionsWrapped(l,whichpos(l))); %angulo azimuth
            angletoUE_theta = asin(distanceVertical/distances(l,k,n));  %angulo de elevación
            %Generar la matriz de correlación espacial utilizando el modelo de dispersión local
            %(2.18) y distribución angular gaussiana escalando las matrices normalizadas con la ganancia del canal
            if nargin>6
                R(:,:,l,k,n) = db2pow(gainOverNoisedB(l,k,n))*functionRlocalscattering(N,angletoUE_varphi,angletoUE_theta,ASD_varphi,ASD_theta,antennaSpacing);
            else
                R(:,:,l,k,n) = db2pow(gainOverNoisedB(l,k,n))*eye(N);  %If angular standard deviations are not specified, set i.i.d. fading
            end
        end
    end
    for k = 1:sum(K)
        %Determinar el masterAP para cada UE buscando el AP con
        %condición de canal
        [~,master] = max(gainOverNoisedB(:,k,n));
        D(master,k,n) = 1;
        masterAPs(k) = master;        
       %Asignar pilotos ortogonales a los primeros tau_p UEs
        %Algorithm 4.1
        if k <= tau_p
            
            pilotIndex(k,n) = k;
            
        else %Además UEs
            %Calcular la interferencia acumulada en cada UE para piloto
            %respecto su masterAP
            %Algorithm 4.1
            pilotinterference = zeros(tau_p,1);
            
            for t = 1:tau_p
                
                pilotinterference(t) = sum(db2pow(gainOverNoisedB(master,pilotIndex(1:k-1,n)==t,n)));
                
            end

            % Asegurarse que ningun piloto es asignado más de rep veces
            count = 0;
             while true
                % Encuentra el piloto con la menor interferencia
                [~, bestpilot] = min(pilotinterference);
            
                % Verifica si el piloto se asigna más de rep veces
                if sum(pilotIndex(1:k-1,n) == bestpilot) >= rep
                    pilotinterference(bestpilot) = Inf;  % Excluir este piloto de la consideración
                    count = count + 1;
                
                % Si se han probado todos los pilotos, romper el bucle
                    if count == tau_p
                        break;
                    end
                else
                    % Asigna el piloto si no se ha asignado más de rep veces
                    pilotIndex(k, n) = bestpilot;
                    break;
                end
            end
        end
    end

    % Cada AP sirve al UE con la mejor condición de canal en cada uno de
    % los pilotos en la configuración cell-free
    for l = 1:L       
        for t = 1:tau_p          
            pilotUEs = find(t==pilotIndex(:,n));
            [~,UEindex] = max(gainOverNoisedB(l,pilotUEs,n));
            D(l,pilotUEs(UEindex),n) = 1;           
        end        
    end
end  

    
    


