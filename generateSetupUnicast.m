function [gainOverNoisedB,R,pilotIndex,D,APpositions,UEpositions,distances, masterAPs] = generateSetupUnicast(L,K,N,tau_p,nbrOfSetups,seed,ASD_varphi,ASD_theta)

%% Define simulation setup

%Set the seed number if it is specified other than zero
if (nargin>5)&&(seed>0)
    rng(seed)
end

%Size of the coverage area (as a square with wrap-around)
squareLength = 1000; %meter

%Size of the clusters of users
C = length(K); %number of clusters of users
squareLengthCluster = 150*ones(1,C); %meter

%Communication bandwidth (Hz)
B = 20e6;

%Noise figure (in dB)
noiseFigure = 7;

%Compute noise power (in dBm)
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;

%Pathloss parameters for the model in (5.42)
alpha = 36.7;
constantTerm = -30.5;

%Standard deviation of the shadow fading in (5.43)
sigma_sf = 4;

%Decorrelation distance of the shadow fading in (5.43)
decorr = 9;

%Height difference between an AP and a UE (in meters)
distanceVertical = 10;

%Define the antenna spacing (in number of wavelengths)
antennaSpacing = 1/2; %Half wavelength distance

%Prepare to save results
gainOverNoisedB = zeros(L,sum(K),nbrOfSetups);
R = zeros(N,N,L,sum(K),nbrOfSetups);
distances = zeros(L,sum(K),nbrOfSetups);
pilotIndex = zeros(sum(K),nbrOfSetups);
D = zeros(L,sum(K),nbrOfSetups);
masterAPs = zeros(sum(K),1); %the indices of master AP of each UE k 

%% Go through all setups
for n = 1:nbrOfSetups   
    %Random AP locations with uniform distribution
    APpositions = (rand(L,1) + 1i*rand(L,1)) * squareLength;   
    %Prepare to compute UE locations
    UEpositions = zeros(sum(K),1);   
    %Compute alternative AP locations by using wrap around
    wrapHorizontal = repmat([-squareLength 0 squareLength],[3 1]);
    wrapVertical = wrapHorizontal';
    wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
    APpositionsWrapped = repmat(APpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[L 1]);   
    %Prepare to store shadowing correlation matrix
    shadowCorrMatrix = sigma_sf^2*ones(sum(K),sum(K));
    shadowAPrealizations = zeros(sum(K),L); 
    %Add UEs in clusters
    %Generate center of users' cluster location
    Clusterpositions = (rand(C,1) + 1i*rand(C,1)) * squareLength;
    % Permutar todos los usuarios
    allUsers = randperm(sum(K));

    % Inicializar el índice para recorrer el vector permutado de usuarios
    userIndex = 1;

    for c = 1:C
        cluster_users = allUsers(userIndex:userIndex + K(c) - 1);

        for k = 1:K(c)
            % Generate a random UE location in the area
            UEpositions(cluster_users(k)) = Clusterpositions(c) + (rand(1, 1) + 1i * rand(1, 1)) * squareLengthCluster(c);
        end

        % Actualizar el índice de usuarios
        userIndex = userIndex + K(c);
    end


    UEpositions = UEpositions.';
    UEpositions = UEpositions(:);
    plot(UEpositions, 'b*') %graficar posiciones de UEs/clusters

    for k = 1:sum(K)
        %Compute distances assuming that the APs are 10 m above the UEs
        [distanceAPstoUE,whichpos] = min(abs(APpositionsWrapped - repmat(UEpositions(k),size(APpositionsWrapped))),[],2);
        distances(:,k,n) = sqrt(distanceVertical^2+distanceAPstoUE.^2);           
        %If this is not the first UE
        if k-1>0                
            %Compute distances from the new prospective UE to all other UEs
            shortestDistances = zeros(k-1,1);                
            for i = 1:k-1
                shortestDistances(i) = min(abs(UEpositions(k) - UEpositions(i) + wrapLocations));
            end                
            %Compute conditional mean and standard deviation necessary to
            %obtain the new shadow fading realizations, when the previous
            %UEs' shadow fading realization have already been generated.
            %This computation is based on Theorem 10.2 in "Fundamentals of
            %Statistical Signal Processing: Estimation Theory" by S. Kay
            newcolumn = sigma_sf^2*2.^(-shortestDistances/decorr);
            term1 = newcolumn'/shadowCorrMatrix(1:k-1,1:k-1);
            meanvalues = term1*shadowAPrealizations(1:k-1,:);
            stdvalue = sqrt(sigma_sf^2 - term1*newcolumn);                
        else %If this is the first UE                
            %Add the UE and begin to store shadow fading correlation values
            meanvalues = 0;
            stdvalue = sigma_sf;
            newcolumn = [];                
        end            
        %Generate the shadow fading realizations
        shadowing = meanvalues + stdvalue*randn(1,L); 
        %Compute the channel gain divided by noise power
        gainOverNoisedB(:,k,n) = constantTerm - alpha*log10(distances(:,k,n)) + shadowing' - noiseVariancedBm;                                  
        %Update shadowing correlation matrix and store realizations
        shadowCorrMatrix(1:k-1,k) = newcolumn;
        shadowCorrMatrix(k,1:k-1) = newcolumn';
        shadowAPrealizations(k,:) = shadowing;  
        %Go through all APs
        for l = 1:L            
            %Compute nominal angle between UE k and AP l
            angletoUE_varphi = angle(UEpositions(k)-APpositionsWrapped(l,whichpos(l))); %azimuth angle
            angletoUE_theta = asin(distanceVertical/distances(l,k,n));  %elevation angle
            %Generate spatial correlation matrix using the local
            %scattering model in (2.18) and Gaussian angular distribution
            %by scaling the normalized matrices with the channel gain
            if nargin>6
                R(:,:,l,k,n) = db2pow(gainOverNoisedB(l,k,n))*functionRlocalscattering(N,angletoUE_varphi,angletoUE_theta,ASD_varphi,ASD_theta,antennaSpacing);
            else
                R(:,:,l,k,n) = db2pow(gainOverNoisedB(l,k,n))*eye(N);  %If angular standard deviations are not specified, set i.i.d. fading
            end
        end
    end
    for k = 1:sum(K)
        %Determine the master AP for UE k by looking for AP with best
        %channel condition
        [~,master] = max(gainOverNoisedB(:,k,n));
        D(master,k,n) = 1;
        masterAPs(k) = master;        
       %Assign orthogonal pilots to the first tau_p UEs according to
        %Algorithm 4.1
        if k <= tau_p
            
            pilotIndex(k,n) = k;
            
        else %Assign pilot for remaining UEs
            
            %Compute received power to the master AP from each pilot
            %according to Algorithm 4.1
            pilotinterference = zeros(tau_p,1);
            
            for t = 1:tau_p
                
                pilotinterference(t) = sum(db2pow(gainOverNoisedB(master,pilotIndex(1:k-1,n)==t,n)));
                
            end
            
            %Find the pilot with the least receiver power according to
            %Algorithm 4.1
            [~,bestpilot] = min(pilotinterference);
            pilotIndex(k,n) = bestpilot;
            
        end
    end
    %Each AP serves the UE with the strongest channel condition on each of
    %the pilots in the cell-free setup
    for l = 1:L       
        for t = 1:tau_p          
            pilotUEs = find(t==pilotIndex(:,n));
            [~,UEindex] = max(gainOverNoisedB(l,pilotUEs,n));
            D(l,pilotUEs(UEindex),n) = 1;           
        end        
    end
end  

    
    


