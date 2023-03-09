function [gainOverNoisedB,R,pilotIndex,D,D_small,APpositions,UEpositions,distances] = generateSetup(L,K,N,tau_p,nbrOfSetups,seed,APselection,alpha,pilotAssignment,ASD_varphi,ASD_theta)
%This function generates realizations of the simulation setup described in
%Section 5.3.
%
%INPUT:
%L               = Number of APs per setup
%K               = Number of UEs in the network
%N               = Number of antennas per AP
%tau_p           = Number of orthogonal pilots
%nbrOfSetups     = Number of setups with random UE and AP locations
%seed            = Seed number of pseudorandom number generator
%ASD_varphi      = Angular standard deviation in the local scattering model
%                  for the azimuth angle (in radians)
%ASD_theta       = Angular standard deviation in the local scattering model
%                  for the elevation angle (in radians)
%
%OUTPUT:
%gainOverNoisedB = Matrix with dimension L x K x nbrOfSetups where
%                  element (l,k,n) is the channel gain (normalized by the
%                  noise variance) between AP l and UE k in setup n
%R               = Matrix with dimension N x N x L x K x nbrOfSetups
%                  where (:,:,l,k,n) is the spatial correlation matrix
%                  between AP l and UE k in setup n, normalized by noise
%pilotIndex      = Matrix with dimension K x nbrOfSetups containing the
%                  pilots assigned to the UEs
%D               = DCC matrix with dimension L x K x nbrOfSetups where (l,k,n)
%                  is one if AP l serves UE k in setup n and zero otherwise
%                  for cell-free setup
%D_small         = DCC matrix with dimension L x K x nbrOfSetups where (l,k,n)
%                  is one if AP l serves UE k in setup n and zero otherwise
%                  for small-cell setup
%APpositions     = Vector of length L with the AP locations, where the real
%                  part is the horizontal position and the imaginary part
%                  is the vertical position
%UEpositions     = Vector of length K with UE positions, measured in the
%                  same way as APpositions
%distances       = Matrix with same dimension as gainOverNoisedB containing
%                  the distances in meter between APs and UEs
%
%This Matlab function was developed to generate simulation results to:
%
%Ozlem Tugfe Demir, Emil Bjornson and Luca Sanguinetti (2021),
%"Foundations of User-Centric Cell-Free Massive MIMO",
%Foundations and Trends in Signal Processing: Vol. 14: No. 3-4,
%pp 162-472. DOI: 10.1561/2000000109
%
%This is version 1.0 (Last edited: 2021-01-31)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%monograph as described above.


%% Define simulation setup

%Set the seed number if it is specified other than zero
if (nargin>5)&&(seed>0)
    rng(seed)
end

%Size of the coverage area (as a square with wrap-around)
squareLength = 1000; %meter

%Communication bandwidth (Hz)
B = 100e6;%20e6;

%Frequency center (GHz)
fc_GHz = 3.5;%2

%Noise figure (in dB)
noiseFigure = 8;%7;

%Compute noise power (in dBm)
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;

%Pathloss parameters for the model in (5.42)
% alpha = 36.7;
constantTerm = -30.5;

%Standard deviation of the shadow fading in (5.43)
sigma_sf = 4;

%Decorrelation distance of the shadow fading in (5.43)
decorr = 9;

%Height difference between an AP and a UE (in meters)
distanceVertical = 10;
HeightUE_m = 1.65;

%Define the antenna spacing (in number of wavelengths)
antennaSpacing = 1/2; %Half wavelength distance


%Prepare to save results
gainOverNoisedB = zeros(L,K,nbrOfSetups);
R = zeros(N,N,L,K,nbrOfSetups);
distances = zeros(L,K,nbrOfSetups);
pilotIndex = zeros(K,nbrOfSetups);
D = zeros(L,K,nbrOfSetups);
D_small = zeros(L,K,nbrOfSetups);

masterAPs = zeros(K,1); %the indices of master AP of each UE k


%% Go through all setups
for n = 1:nbrOfSetups
    
    %Random AP locations with uniform distribution
    %     APpositions = (rand(L,1) + 1i*rand(L,1)) * squareLength;
    
    %Random AP locations with Poisson distribution
    [APpositions,~,~] = MF_AF_GS_generateSetupAPs(squareLength,squareLength,distanceVertical+HeightUE_m,L,1,'Poisson');
    
    
    %Prepare to compute UE locations
    UEpositions = zeros(K,1);
    
    
    %Compute alternative AP locations by using wrap around
    wrapHorizontal = repmat([-squareLength 0 squareLength],[3 1]);
    wrapVertical = wrapHorizontal';
    wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
    APpositionsWrapped = repmat(APpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[L 1]);
    
    %Prepare to store shadowing correlation matrix
    shadowCorrMatrix = sigma_sf^2*ones(K,K);
    shadowAPrealizations = zeros(K,L);
    
    
    
    %Add UEs
    for k = 1:K
        
        %Generate a random UE location in the area
        UEposition = (rand(1,1) + 1i*rand(1,1)) * squareLength;
        
        %Compute distances assuming that the APs are 10 m above the UEs
        [distanceAPstoUE,whichpos] = min(abs(APpositionsWrapped - repmat(UEposition,size(APpositionsWrapped))),[],2);
        distances(:,k,n) = sqrt(distanceVertical^2+distanceAPstoUE.^2);
        
        %If this is not the first UE
        if k-1>0
            
            %Compute distances from the new prospective UE to all other UEs
            shortestDistances = zeros(k-1,1);
            
            for i = 1:k-1
                shortestDistances(i) = min(abs(UEposition - UEpositions(i) + wrapLocations));
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
        
        %Compute the channel gain divided by noise power, for Fc_GHz = 2 defined in TS 36.814, Table B.1.2.1-1
        %         gainOverNoisedB(:,k,n) = constantTerm - alpha*log10(distances(:,k,n)) + shadowing' - noiseVariancedBm;
        
        %Compute the channel gain divided by noise power, for any Fc_GHz defined in TR 38.901, Table 7.4.1-1
        gainOverNoisedB(:,k,n) = AF_LargeScaleGainOverNoise('Umi',fc_GHz,noiseVariancedBm,distances(:,k,n),distanceVertical+HeightUE_m,HeightUE_m);
        gainOverNoisedB(:,k,n) = gainOverNoisedB(:,k,n) + shadowing';
        
        
        %Update shadowing correlation matrix and store realizations
        shadowCorrMatrix(1:k-1,k) = newcolumn;
        shadowCorrMatrix(k,1:k-1) = newcolumn';
        shadowAPrealizations(k,:) = shadowing;
        
        %Store the UE position
        UEpositions(k) = UEposition;
        
        
        %Determine the master AP for UE k by looking for AP with best
        %channel condition
        [~,master] = max(gainOverNoisedB(:,k,n));
        D(master,k,n) = 1;
        masterAPs(k) = master;
        
        if  strcmp(pilotAssignment,'nonJoint_minContamination') || strcmp(pilotAssignment,'joint_minContamination')
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
        elseif strcmp(pilotAssignment,'nonJoint_random') || strcmp(pilotAssignment,'joint_random')
            tau_p = min(tau_p,K);
            pilotIndex(:,n) = [1:tau_p randi(tau_p,1,K-tau_p)]';
        end
        
        
        
        %Go through all APs
        for l = 1:L
            
            %Compute nominal angle between UE k and AP l
            angletoUE_varphi = angle(UEpositions(k)-APpositionsWrapped(l,whichpos(l))); %azimuth angle
            angletoUE_theta = asin(distanceVertical/distances(l,k,n));  %elevation angle
            %Generate spatial correlation matrix using the local
            %scattering model in (2.18) and Gaussian angular distribution
            %by scaling the normalized matrices with the channel gain
            if nargin>9
                R(:,:,l,k,n) = db2pow(gainOverNoisedB(l,k,n))*functionRlocalscattering(N,angletoUE_varphi,angletoUE_theta,ASD_varphi,ASD_theta,antennaSpacing);
            else
                R(:,:,l,k,n) = db2pow(gainOverNoisedB(l,k,n))*eye(N);  %If angular standard deviations are not specified, set i.i.d. fading
            end
        end
        
    end
    
    
    %Each AP serves the UE with the strongest channel condition on each of
    %the pilots in the cell-free setup
    if strcmp(APselection,'BookEmil') || strcmp(APselection,'ClusterControlCmax')
        for l = 1:L
            
            for t = 1:tau_p
                
                pilotUEs = find(t==pilotIndex(:,n));
                [~,UEindex] = max(gainOverNoisedB(l,pilotUEs,n));
                D(l,pilotUEs(UEindex),n) = 1;
                
            end
            
        end

    end
    
    %Determine the AP serving each UE in the small-cell setup according to
    %(5.47) by considering only the APs from the set M_k for UE k, i.e.,
    %where D(:,k,n) is one.
    for k=1:K
        
        tempmat = -inf*ones(L,1);
        tempmat(D(:,k,n)==1,1) = gainOverNoisedB(D(:,k,n)==1,k,n);
        [~,servingAP] = max(tempmat);
        D_small(servingAP,k,n) = 1;
        
    end
    
end

if strcmp(APselection,'CanonicalCF')
    D = ones(L,K);
            
elseif strcmp(APselection,'ClusterControlCmax')
            
    E = D;
    U_max = tau_p; 
%     alpha = 0.75;
    Cmax = floor((alpha*U_max*L)/K);
    clusterSize_UEs = sum(D,1);
    overloadedUEs_id = find(clusterSize_UEs>Cmax);
    nbrOfAPsToExclude = sum(D(:,overloadedUEs_id),1) - Cmax;
    gainOverNoise_Lin = db2pow(gainOverNoisedB);
    
    for indexUE = 1:length(overloadedUEs_id)
        servingAPs_id = find(D(:,overloadedUEs_id(indexUE)) == 1);
        [~, APsToExclude_id] = sort(gainOverNoise_Lin(servingAPs_id,overloadedUEs_id(indexUE)),'ascend');
        
        APsToExclude_id = servingAPs_id(APsToExclude_id);
        APsToExclude_id = APsToExclude_id(1:nbrOfAPsToExclude(indexUE));
        D(APsToExclude_id, overloadedUEs_id(indexUE)) = 0;
    end
    
    A = 0;
    
end

end

% Authors: Marx Freitas
%          André Fernandes
%          Gilvan Borges

% Generates the AP distribution in the coverage area. The APS can be distributed in conventional
% geometries or randomly

%INPUT:
%startAPPosX_m  = the length in meters of the X-axis of the APs positions in the coverage area
%startAPPosY_m  = the length in meters of the Y-axis of the APs positions in the coverage area
%HeightAP_m     = the height of APs in meters
%nbrAPsPerArray = the number of APs per array
%nbrArrays      = the number of arrays utilized to form the geometry of the APs distributions
%APsCoverage    = the type of coverage selected in the main code

%OUTPUT:
%APPositionsXYZ_m = matrix with dimension L x 3 is the 3D position of the APs in the coverage area
%nbrOfAPs         = the number of APs used in each AP distribution (APsCoverage)

function [APpositions,APPositionsXYZ_m,nbrOfAPs] = MF_AF_GS_generateSetupAPs(startAPPosX_m,startAPPosY_m,HeightAP_m,nbrAPsPerArray,nbrArrays,APsCoverage)

switch APsCoverage
    
    %=====================================================================%
    % For a single array of APs
    %=====================================================================%
    
    % Author: Marx Freitas
    case 'Line'
        
        % Set the APs locations at XYZ plane
        APPositionsXYZ_m = [ linspace(-startAPPosX_m, startAPPosX_m, nbrAPsPerArray)' zeros(nbrAPsPerArray, 1)...
            HeightAP_m * ones(nbrAPsPerArray, 1) ];
        
        % Computes the number of APs
        nbrOfAPs = nbrAPsPerArray;
        
        %=====================================================================%
        % For randomic APs, considering a gaussian distribution
        %=====================================================================%
        
        % Author: Marx Freitas
    case 'Randomic'
        
        % Set the APs locations at XYZ plane
        APPositionsXYZ_m = [ startAPPosX_m * randn(nbrAPsPerArray*nbrArrays, 1)...
            startAPPosY_m * randn(nbrAPsPerArray*nbrArrays, 1) HeightAP_m * ones(nbrAPsPerArray*nbrArrays, 1) ];
        
        % Computes the number of APs
        nbrOfAPs = nbrAPsPerArray*nbrArrays;
        
        %=====================================================================%
        % For randomic APs, considering a uniform distribution
        %=====================================================================%
        
        % Author: Marx Freitas
    case 'RandomicUniform'
        
        % Set the APs locations at XYZ plane
        APPositionsXYZ_m = [ -startAPPosX_m+(2*startAPPosX_m)*rand(nbrAPsPerArray*nbrArrays, 1)...
            -startAPPosY_m+(2*startAPPosY_m)*rand(nbrAPsPerArray*nbrArrays, 1) ...
            HeightAP_m * ones(nbrAPsPerArray*nbrArrays, 1) ];
        
        % Computes the number of APs
        nbrOfAPs = nbrAPsPerArray*nbrArrays;
        
        %=====================================================================%
        % For randomic APs, considering a Poison process
        %=====================================================================%
        
        % Author: Gilvan Soares
    case 'Poisson'
        
        % Computes the number of APs
        nbrOfAPs = nbrAPsPerArray*nbrArrays;
        
        % Set the APs locations at XYZ plane
        APPositionsXYZ_m = GS_apsDistribution(nbrOfAPs,[startAPPosX_m; startAPPosY_m], HeightAP_m, 'regularity', 1);
        APPositionsXYZ_m(:,1) = (2.*APPositionsXYZ_m(:,1)) - startAPPosX_m;
        APPositionsXYZ_m(:,2) = (2.*APPositionsXYZ_m(:,2)) - startAPPosY_m;
        
        %=====================================================================%
        % For a Matrix of APs
        %=====================================================================%
        
        % Author: Marx Freitas
    case 'Matrix'
        
        % Creates the first columns of the matrix
        APPositionsXYZ_m = repmat( [zeros(nbrAPsPerArray,1) linspace(-startAPPosY_m, startAPPosY_m, nbrAPsPerArray)' ...
            HeightAP_m.*ones(nbrAPsPerArray, 1)],[nbrArrays 1] );
        
        auxAPPositionsXYZ_m = linspace(-startAPPosX_m, startAPPosX_m, nbrArrays);
        
        % In order to create the matrix, copy the first collum at different locations of the XYZ plane
        for indexAP_XYZ = 1: length(auxAPPositionsXYZ_m)
            % Set the APs locations at XYZ plane
            APPositionsXYZ_m( (indexAP_XYZ-1)*nbrAPsPerArray+1:indexAP_XYZ*nbrAPsPerArray,1 ) = auxAPPositionsXYZ_m(indexAP_XYZ);
        end
        
        % Computes the number of APs
        nbrOfAPs = nbrAPsPerArray*nbrArrays;
        
        %=====================================================================%
        % For a "Diagonal X" of APs
        %=====================================================================%
        
        % Author: Marx Freitas
    case 'Xdiagonal'
        
        % Computes the number of APs
        nbrOfAPs = 2*nbrAPsPerArray;
        
        % Set the APs locations at XYZ plane
        APPositionsXYZ_m = [zeros(nbrOfAPs,1) zeros(nbrOfAPs,1) HeightAP_m.*ones(nbrOfAPs, 1)];
        APPositionsXYZ_m(1:nbrAPsPerArray,1) = linspace(startAPPosX_m,-startAPPosX_m,nbrAPsPerArray)';
        APPositionsXYZ_m(nbrAPsPerArray+1:nbrOfAPs,1) = linspace(-startAPPosX_m,startAPPosX_m,nbrAPsPerArray)';
        APPositionsXYZ_m(1:nbrAPsPerArray,2) = linspace(startAPPosY_m,-startAPPosY_m,nbrAPsPerArray)';
        APPositionsXYZ_m(nbrAPsPerArray+1:nbrOfAPs,2) = APPositionsXYZ_m(1:nbrAPsPerArray,2);
        
        %=====================================================================%
        % For APs serving a square area
        %=====================================================================%
        % Authors: Marx Freitas and André Fernandes
        
        % Author: André Fernandes
        % Edited by Marx Freitas
    case 'Rectangle_two_array'
        
        % Calculates the number of APs
        nbrOfAPs = 2*(nbrAPsPerArray - 1);
        
        % Genetares the base of the positions of APs. The Z coordinates are constant for all the dimensions below for array 1(oriental and south sides)
        Array1=linspace(0,2*(startAPPosY_m+startAPPosX_m), nbrAPsPerArray)';
        APPositionsXYZarr1_m =zeros(length(Array1),3);
        
        for indexArray1=1:length(Array1)
            if Array1(indexArray1)<=2*startAPPosY_m
                APPositionsXYZarr1_m(indexArray1,:)=[startAPPosX_m,startAPPosY_m-Array1(indexArray1),HeightAP_m];
            else
                APPositionsXYZarr1_m(indexArray1,:)=[startAPPosX_m-(Array1(indexArray1)-2*(startAPPosY_m)),-startAPPosY_m,HeightAP_m];
            end
        end
        
        % For array 2(north and ocidental sides)
        APPositionsXYZarr2_m(:,:) = -APPositionsXYZarr1_m(:,:);
        APPositionsXYZarr2_m(end,:) = []; APPositionsXYZarr2_m(1,:) = []; % MF: In order to avoid AP overlap.
        
        % Set the APs locations at XYZ plane
        APPositionsXYZ_m = vertcat(APPositionsXYZarr1_m,APPositionsXYZarr2_m);
        
        clear Array1 Array2 APPositionsXYZarr1_m APPositionsXYZarr2_m
        
        % Author: Marx Freitas
    case 'Rectangle'
        
        % Place the first array of APs in the region of the XY plane, where Y is negative, and X vary from -Xo to Xo
        APPositionsXYZdim1_m = [ linspace(startAPPosX_m, -startAPPosX_m, nbrAPsPerArray)'...
            zeros(nbrAPsPerArray, 1) HeightAP_m*ones(nbrAPsPerArray, 1) ];
        APPositionsXYZdim1_m = [ APPositionsXYZdim1_m(:,1) APPositionsXYZdim1_m(:,2)-startAPPosY_m APPositionsXYZdim1_m(:,3) ];
        
        % Place the second array of APs in the region of the XY plane, where Y is positive, and X vary from -Xo to Xo
        APPositionsXYZdim2_m = [APPositionsXYZdim1_m(:,1) APPositionsXYZdim1_m(:,2)+2*startAPPosY_m APPositionsXYZdim1_m(:,3)];
        
        % Place the third array of APs in the region of the XY plane, where X is negative, and Y vary from -Yo to Yo
        APPositionsXYZdim3_m = [ -startAPPosX_m.*ones(nbrAPsPerArray, 1) linspace(-startAPPosY_m, startAPPosY_m, nbrAPsPerArray)'...
            HeightAP_m*ones(nbrAPsPerArray, 1) ];
        APPositionsXYZdim3_m = APPositionsXYZdim3_m(2:end-1,1:end,1:end);
        
        % Place the fourth array of APs in the region of the XY plane, where X is positive, and Y vary from -Yo to Yo
        APPositionsXYZdim4_m = [ startAPPosX_m.*ones(nbrAPsPerArray, 1) linspace(+startAPPosY_m, -startAPPosY_m, nbrAPsPerArray)'...
            HeightAP_m*ones(nbrAPsPerArray, 1) ];
        APPositionsXYZdim4_m = APPositionsXYZdim4_m(2:end-1,1:end,1:end);
        
        % Set the APs locations at XYZ plane
        APPositionsXYZ_m = vertcat(APPositionsXYZdim1_m,APPositionsXYZdim2_m, APPositionsXYZdim3_m,APPositionsXYZdim4_m);
        
        % Compute the number of APs
        nbrOfAPs = 4*(nbrAPsPerArray-1);
        
        clear auxAPPosition APPositionsXYZdim1_m APPositionsXYZdim2_m...
            APPositionsXYZdim3_m APPositionsXYZdim4_m
        
    case 'RectangleXdiagonal'
        
        % Computes the rectangle distributions of the APs. Estimates the four dimensions of the rectangle
        APPositionsXYZdim1_m = [ linspace(startAPPosX_m, -startAPPosX_m, nbrAPsPerArray)'...
            zeros(nbrAPsPerArray, 1) HeightAP_m*ones(nbrAPsPerArray, 1) ];
        APPositionsXYZdim1_m = [ APPositionsXYZdim1_m(:,1) APPositionsXYZdim1_m(:,2)-startAPPosY_m APPositionsXYZdim1_m(:,3) ];
        
        APPositionsXYZdim2_m = APPositionsXYZdim1_m; APPositionsXYZdim2_m(:,2) = APPositionsXYZdim2_m(:,2)+2*startAPPosY_m;
        
        APPositionsXYZdim3_m = [ -startAPPosX_m.*ones(nbrAPsPerArray, 1)...
            linspace(-startAPPosY_m, +startAPPosY_m, nbrAPsPerArray)' APPositionsXYZdim1_m(:,3) ];
        APPositionsXYZdim3_m = APPositionsXYZdim3_m(2:end-1,1:end,1:end);
        
        APPositionsXYZdim4_m = [ startAPPosX_m.*ones(nbrAPsPerArray, 1)...
            linspace(+startAPPosY_m, -startAPPosY_m, nbrAPsPerArray)' APPositionsXYZdim1_m(:,3) ];
        APPositionsXYZdim4_m = APPositionsXYZdim4_m(2:end-1,1:end,1:end);
        
        APPositionsRecXYZ_m = vertcat(APPositionsXYZdim1_m,APPositionsXYZdim2_m,APPositionsXYZdim3_m,APPositionsXYZdim4_m);
        
        % Computes Xdiagonal distributions of the APs
        numberOfAPsXdiag = 4*(nbrAPsPerArray);
        
        APPositionsXYZXdiag_m = [zeros(numberOfAPsXdiag,1) zeros(numberOfAPsXdiag,1) HeightAP_m.*ones(numberOfAPsXdiag, 1)];
        APPositionsXYZXdiag_m(1:numberOfAPsXdiag/2,1) = linspace(startAPPosX_m,-startAPPosX_m,numberOfAPsXdiag/2)';
        APPositionsXYZXdiag_m(numberOfAPsXdiag/2+1:numberOfAPsXdiag,1) = linspace(-startAPPosX_m,startAPPosX_m,numberOfAPsXdiag/2)';
        APPositionsXYZXdiag_m(1:numberOfAPsXdiag/2,2) = linspace(startAPPosY_m,-startAPPosY_m,numberOfAPsXdiag/2)';
        APPositionsXYZXdiag_m(numberOfAPsXdiag/2+1:numberOfAPsXdiag,2) = APPositionsXYZXdiag_m(1:numberOfAPsXdiag/2,2);
        
        % Find APs overlap in the code processing
        indexXdiag = find(abs(APPositionsXYZXdiag_m(:,1)) == startAPPosX_m);
        
        % In order to avoid AP overlap, exclude the repeated positions.
        APPositionsXYZXdiag_m(indexXdiag,:) = [];
        
        APPositionsXYZ_m = vertcat(APPositionsRecXYZ_m,APPositionsXYZXdiag_m);
        
        % Computes the number of APs
        nbrOfAPs = 8*(nbrAPsPerArray-1);
        
        clear APPositionsXYZdim1_m APPositionsXYZdim2_m APPositionsXYZdim3_m APPositionsXYZdim4_m...
            APPositionsRecXYZ_m APPositionsXYZXdiag_m indexXdiag numberOfAPsXdiag
        
end

% Adjusting the APs locations at XYZ plane. With this, the coverage area will be located only in the plane (x>0,y>0)
APPositionsXYZ_m(:,1) = (APPositionsXYZ_m(:,1) + startAPPosX_m)/2;
APPositionsXYZ_m(:,2) = (APPositionsXYZ_m(:,2) + startAPPosY_m)/2;
APpositions = APPositionsXYZ_m(:,1) + 1i*APPositionsXYZ_m(:,2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Creates access point (AP) positions randomly on any 3D plane delimited
% by a user-defined parallelogram region. The distribution of APs can vary
% gradually from a pure random process (Binomial Point Process) to another
% with more regularity (Hard Core Point Process), according to the user.
%
%
% SYNTAX:
%   apsPositions = aps_distribution(numberOfAps)
%   apsPositions = aps_distribution(numberOfAps, dimensions)
%   apsPositions = aps_distribution(numberOfAps, dimensions, startPosition)
%   apsPositions = aps_distribution(---, name, value)
%   apsPositions = aps_distribution(---, name1, value1, name2, value2)
%
%
% INPUT:
%   numberOfAps:
%      Define the number of APs will be generated.
%
%   dimensions:
%      Matrix of dimensions 2-by-3, in which the lines represent
%      three-dimensional vectors. These two vectors define the plane and
%      the parallelogram region in which the APs will be generated. Each
%      vector represents a side of the parallelogram whose vertex lies at
%      the origin of the coordinate system. The norms of these vectors
%      represent the lengths of the sides of the parallelogram. For
%      convenience, the user can enter a 2-by-2 matrix instead of 2-by-3,
%      in which case a third null coordinate will automatically be added
%      to the two two-dimensional vectors. The user can also enter a
%      2-by-1 column in the form [width ; height], in this case the matrix
%      elements represent the width and height of a rectangle in the xy
%      plane, that is, it is equivalent to a 2-by-3 matrix in the form
%      [width 0 0; 0 length 0]. When not specified, it assumes the default
%      value [1 0 0 ; 0 1 0], that is, a square with side 1 in the xy
%      plane.
%
%   startPosition:
%      Vector of dimensions 1-by-3, representing the initial position of
%      the simulation parallelogram region. For convenience, the user can
%      enter a 1-by-2 vector in the form [x0 y0], which is equivalent to
%      [x0 y0 0]. The user can also enter a scalar z0, in this case it is
%      equivalent to entering [0 0 z0]. When not specified, it assumes the
%      default value [0 0 0].
%
%   name, value:
%      Input argument pair that specifies some settings. Can assume:
%      ---------------------------------------------------------------
%      |      name      |                   value                    |
%      ---------------------------------------------------------------
%      |  'regularity'  |      Real value in the range [0, 1].       |
%      |                                                             |
%      |   'attempts'   |        Non-negative integer value.         |
%      ---------------------------------------------------------------
%      The name 'regularity' defines the degree of regularity in the
%      distribution of the generated APs. Being 0 a distribution with
%      minimal regularity (Binomial Point Process) and 1 with maximum
%      regularity achieved by the algorithm (Hard Core Point Process). The
%      closer to 1, the longer the processing time. When not specified,
%      the default is 1/2.
%      The name 'attempts' defines the maximum number of attempts used by
%      the algorithm to satisfy the regularity condition in each AP. When
%      not specified, the default is 30. Particular attention should be
%      paid to increasing the maximum number of attempts, with the risk of
%      exponential growth in processing time.
%
%
% OUTPUT:
%   apsPositions:
%      Matrix of dimension numberOfAps-by-3, where the rows are the
%      three-dimensional coordinates of each generated AP.
%
%
% EXEMPLES:
% % Exemple 1 - Regularly distribute 400 APs on the roof of a 400m x 200m
% % shed and 6m ceiling height:
%      numberOfAps = 400;
%      shedDimensions = [400; 200];
%      heightOfAps = 12;
%      apsPositions = aps_distribution(numberOfAps, ...
%          shedDimensions, heightOfAps, 'regularity', 1);
%
%      figure
%      plot3(apsPositions(:, 1), apsPositions(:, 2), ...
%          apsPositions(:, 3), '*')
%      grid
%
% % Exemple 2 - Distribute 100 APs on a 5m square panel on a side,
% % inclined 45 degrees from the floor:
%      numberOfAps = 100;
%      side= 5;
%      apsPositions = aps_distribution(numberOfAps, ...
%          [side 0 0; 0 -side side], [0 side 0]);
%
%      figure
%      plot3(apsPositions(:, 1), apsPositions(:, 2), ...
%          apsPositions(:, 3), '*')
%      grid
%
% % Exemple 3 - Distribute 1000 APs on a parallelogram randomly defined in
% % space:
%      apsPositions = aps_distribution(1000, [rand(2, 3)]);
%
%      figure
%      plot3(apsPositions(:, 1), apsPositions(:, 2), ...
%          apsPositions(:, 3), '*')
%      grid
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Created by:    Gilvan Soares Borges
% Contact:       gilvan.borges@ifpa.edu.br
% Last modified: 07/10/2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function apsPositions3D = GS_apsDistribution(numberOfAps, varargin)

[width, height, linearTransformMatrix, startPosition, regularity, attempts] = ...
    check_input_arguments(numberOfAps, varargin); % Auxiliary function: Formats input arguments and avoids possible input errors.


%% DEFINING THE DEFAULT PARTITION CELL:

du = sqrt(width*height/numberOfAps);

if (width >= du) && (height >= du) % default case.
    dx = du;
    dy = du;
else % particular cases.
    if (width  < du)
        dx = width;
        dy = height/numberOfAps;
    end
    if (height < du)
        dx = width /numberOfAps;
        dy = height;
    end
end


%% GENERATING 2D AP POSITIONS WITH A UNIFORM BPP (BINOMIAL POINT PROCESS):

xPartition = 0:dx:width;
yPartition = 0:dy:height;

nXPartition = length(xPartition);
nYPartition = length(yPartition);
initialNumberAps = nXPartition*nYPartition;

[xCellPositions, yCellPositions] = meshgrid(xPartition, yPartition);

xApsPositions = xCellPositions + dx*rand(nYPartition, nXPartition);
yApsPositions = yCellPositions + dy*rand(nYPartition, nXPartition);

apsPositions2D = [xApsPositions(:) yApsPositions(:)];


%% INCREASING FIRST-ORDER STATIONARITY BASED ON HCPP (HARD CORE POINT PROCESS):

minDist0 = regularity*du;

listNeighbors = list_neighbors(nXPartition, nYPartition); % Auxiliary function: Presents the neighbors of each AP.

j = 0;
pairsCloseNeighbors = zeros(4*initialNumberAps - 3*(nXPartition + nYPartition) + 2, 2); % Pre-allocation with the maximum possible number of access points nearby.
for i = 1:length(listNeighbors) % Identify nearby APs.
    neighborsOfI = listNeighbors{i}(listNeighbors{i} > i);
    neighborsNearI = neighborsOfI(  dist( apsPositions2D(neighborsOfI, :), apsPositions2D(i, :).' ) < minDist0  );
    nNeighborsNearI = length(neighborsNearI);
    
    pairsCloseNeighbors(j+1:j+nNeighborsNearI, 1) = i;
    pairsCloseNeighbors(j+1:j+nNeighborsNearI, 2) = neighborsNearI;
    
    j = j + nNeighborsNearI;
end
pairsCloseNeighbors(j+1:end, :) = [];
shiftedAps = best_shifted_aps(pairsCloseNeighbors); % Auxiliary function: Among the nearby APs, those that will be repositioned are chosen.

for i = shiftedAps % Reposition the selected nearby APs.
    counter = 0;
    lessWorstDist = min(  dist( apsPositions2D(listNeighbors{i}, :), apsPositions2D(i, :).' )  );
    
    while (counter < attempts) && (lessWorstDist <= minDist0) % Reposition while you are near other APs or if the attempt limit is reached.
        newApPosition = [(xCellPositions(i) + dx*rand())  ; (yCellPositions(i) + dy*rand())];
        
        minDist = min(  dist( apsPositions2D(listNeighbors{i}, :), newApPosition )  );
        
        if (minDist > lessWorstDist)
            lessWorstDist = minDist;
            apsPositions2D(i, :) = newApPosition.';
        end
        
        counter = counter + 1;
    end
end


%% ADJUSTING THE NUMBER OF APS:

invalidAps = [... % APs outside the rectangular region.
    nYPartition * find(  apsPositions2D((1:nXPartition-1).*nYPartition, 2) > height  ) ; ...
    initialNumberAps - find(  apsPositions2D(initialNumberAps-(1:nYPartition-1), 1) > width  ) ; ...
    initialNumberAps + find( (apsPositions2D(initialNumberAps, 1) > width) || (apsPositions2D(initialNumberAps, 2) > height) )  -  1 ...
    ];

initialValidNumberAps = initialNumberAps - length(invalidAps);
missingAps = numberOfAps - initialValidNumberAps;

apsPositions2D = [apsPositions2D ; zeros(missingAps, 2)];
for i = 1:(missingAps) % Adding missing APs.
    lessWorstDist = 0; counter = 0;
    
    while (lessWorstDist <= minDist0) && (counter < attempts)
        ind = randi(initialNumberAps);
        newApPosition = [(xCellPositions(ind) + dx*rand()) ; (yCellPositions(ind) + dy*rand())];
        
        if (  ( newApPosition(1) <= width ) && ( newApPosition(2) <= height )  )
            minDist = min(  dist( [apsPositions2D(listNeighbors{ind}, :) ; apsPositions2D(ind, :)], newApPosition )  );
            
            if (minDist > lessWorstDist)
                lessWorstDist = minDist;
                apsPositions2D(initialNumberAps + i, :) = newApPosition.';
            end
        end
        
        counter = counter + 1;
    end
end

apsPositions2D(invalidAps, :) = [];                                          % Excluding APs outside the rectangular region.
apsPositions2D(randi(initialValidNumberAps, -missingAps, 1), :) = [];        % Excluding surplus APs.


%% LINEAR TRANSFORMATION FROM 2D TO 3D:

apsPositions3D = startPosition + apsPositions2D*linearTransformMatrix;


%% FOR QUICK TEST:

% figure
% plot(apsPositions2D(:, 1), apsPositions2D(:, 2), '*')
% hold on
% plot(xCellPositions, yCellPositions, 'r')
% plot(xCellPositions.', yCellPositions.', 'r')
% axis([0, width, 0, height])
%
% figure
% plot3(apsPositions3D(:, 1), apsPositions3D(:, 2), apsPositions3D(:, 3), '*')
% grid


end











%% ======================== AUXILIARY FUNCTIONS 1 ========================

% Formats input arguments and avoids possible input errors.
function [width, height, linearTransformMatrix, startPosition, regularity, attempts] = check_input_arguments(numberOfAps, varargin)

%% DEFAULT VALUES:

dimensions = [1 0 0; 0 1 0];  booleanDimensionsDefault = true;
startPosition = [0 0 0];      booleanP0Default         = true;
regularity = 1/2; % With "regularity = 1/2", only 15% of APs need to be repositioned. With "regularity = 1", 60% of APs need to be repositioned.
attempts   = 30;  % With up to 30 attempts, 100% of nearby APs are guaranteed to be repositioned correctly, as long as "regularity = 1/2".


%% ALLOCATING INPUT ARGUMENTS:

names = []; values = [];
numberExtraInput = length(varargin{1});

if (numberExtraInput == 1) % aps_distribution(numberOfAps, dimensions)
    dimensions = varargin{1}{1}; booleanDimensionsDefault = false;
    
elseif (numberExtraInput == 2)
    if ~ischar(varargin{1}{1}) % aps_distribution(numberOfAps, dimensions, startPosition)
        dimensions       = varargin{1}{1}; booleanDimensionsDefault = false;
        startPosition = varargin{1}{2};    booleanP0Default         = false;
        
    else % aps_distribution(numberOfAps, Name, Value)
        names{1}  = varargin{1}{1};
        values{1} = varargin{1}{2};
        
    end
    
elseif (numberExtraInput == 3) % aps_distribution(numberOfAps, dimensions, Name, Value)
    dimensions = varargin{1}{1}; booleanDimensionsDefault = false;
    names{1}   = varargin{1}{2};
    values{1}  = varargin{1}{3};
    
elseif (numberExtraInput == 4)
    if ~ischar(varargin{1}{1}) % aps_distribution(numberOfAps, dimensions, startPosition, Name, Value)
        dimensions       = varargin{1}{1}; booleanDimensionsDefault = false;
        startPosition = varargin{1}{2};    booleanP0Default         = false;
        names{1}   = varargin{1}{3};
        values{1}  = varargin{1}{4};
        
    else % aps_distribution(numberOfAps, Name1, Value1, Name2, Value2)
        names{1}  = varargin{1}{1};
        values{1} = varargin{1}{2};
        names{2}  = varargin{1}{3};
        values{2} = varargin{1}{4};
    end
    
elseif (numberExtraInput == 5) % aps_distribution(numberOfAps, dimensions, Name1, Value1, Name2, Value2)
    dimensions = varargin{1}{1}; booleanDimensionsDefault = false;
    names{1}   = varargin{1}{2};
    values{1}  = varargin{1}{3};
    names{2}   = varargin{1}{4};
    values{2}  = varargin{1}{5};
    
elseif (numberExtraInput == 6) % aps_distribution(numberOfAps, dimensions, startPosition, Name1, Value1, Name2, Value2)
    dimensions       = varargin{1}{1}; booleanDimensionsDefault = false;
    startPosition = varargin{1}{2};    booleanP0Default         = false;
    names{1}         = varargin{1}{3};
    values{1}        = varargin{1}{4};
    names{2}         = varargin{1}{5};
    values{2}        = varargin{1}{6};
    
elseif (numberExtraInput > 6)
    error('The number of input arguments is invalid.')
end


%% CHECKING THE INPUT ARGUMENTS "numberOfAps":

if ( ~isnumeric(numberOfAps) || ~isreal(numberOfAps) || ~isscalar(numberOfAps) ...
        || (numberOfAps <= 0) || ((numberOfAps - fix(numberOfAps)) ~= 0) )
    error('The input argument "numberOfAps" is invalid.')
end


%% CHECKING THE INPUT ARGUMENT "dimensions":

if ~booleanDimensionsDefault
    if (   ~isnumeric(dimensions) || ~isreal(dimensions) || ...
            ~(size(dimensions, 1) == 2) || ~(size(dimensions, 2) < 4)   )
        error('The input argument "dimensions" is invalid.')
    end
    
    if (size(dimensions, 2) == 1)
        if ( (dimensions(1) > 0) || (dimensions(2) > 0) )
            dimensions = [dimensions(1) 0  0; 0 dimensions(2) 0];
        else
            error('The width and height dimensions of the parallelogram region of the simulation should be positive.')
        end
        
    elseif (size(dimensions, 2) == 2)
        dimensions = [dimensions [0;0]];
    end
    
    if (rank(dimensions) < 2)
        error('The sides of the simulation parallelogram region cannot be collinear.')
    end
end

width  = norm(dimensions(1, :));
height = norm(dimensions(2, :));

linearTransformMatrix = [dimensions(1, :)./width ; dimensions(2, :)./height];


%% CHECKING THE INPUT ARGUMENT "startPosition":

if ~booleanP0Default
    if (  ~isnumeric(startPosition) || ~isreal(startPosition) ...
            || ~isvector(startPosition) || (length(startPosition) > 3)  )
        error('The input argument "startPosition" is invalid.')
    end
    
    if (length(startPosition) == 1); startPosition = [0 0 startPosition]; end
    if (length(startPosition) == 2); startPosition(3) = 0;                end
end


%% CHECKING THE INPUT ARGUMENT PAIR "name" and "value":

for i = 1:length(names)
    switch names{i}
        case 'regularity'
            regularity = values{i};
            
            if (  ~isnumeric(regularity) || ~isreal(regularity) || ~isscalar(regularity)  )
                error('The input argument "regularity" is invalid.')
            end
            if (regularity < 0)
                regularity = 0;
                warning('The "regularity" input argument cannot take values ??outside the range of zero to one. Value zero was assumed.')
            end
            if (regularity > 1)
                regularity = 1;
                warning('The "regularity" input argument cannot take values ??outside the range of zero to one. Value one was assumed.')
            end
            
        case 'attempts'
            attempts = values{i};
            if ( ~isnumeric(attempts) || ~isreal(attempts) || ~isscalar(attempts) ...
                    || (attempts < 0) || ((attempts - fix(attempts)) ~= 0) )
                error('The input argument "attempts" is invalid.')
            end
            
        otherwise
            error([ 'The input argument "' names{i} '" is invalid.' ])
    end
end


end





%% ======================== AUXILIARY FUNCTIONS 2 ========================

% For each AP, it determines its neighbors.
function listNeighbors = list_neighbors(nXPartition, nYPartition)

numberAps = nXPartition*nYPartition;
listNeighbors = cell(numberAps, 1);

for k = 1:numberAps
    if (k == 1)
        listNeighbors{k} = [(k+1) ; ...
            (k+nYPartition) ; ...
            (k+nYPartition+1)    ];
        
    elseif(k == nYPartition)
        listNeighbors{k} = [(k-1); ...
            (k+nYPartition-1); ...
            (k+nYPartition)       ];
        
    elseif (k == numberAps-nYPartition+1)
        listNeighbors{k} = [(k-nYPartition); ...
            (k-nYPartition+1); ...
            (k+1)                    ];
        
    elseif (k == numberAps)
        listNeighbors{k} = [(k-nYPartition-1); ...
            (k-nYPartition); ...
            (k-1)                    ];
        
    elseif ( (k > 1) && (k < nYPartition) )
        listNeighbors{k} = [(k-1); ...
            (k+1); ...
            (k+nYPartition-1); ...
            (k+nYPartition); ...
            (k+nYPartition+1)     ];
        
    elseif (mod(k-1, nYPartition) == 0)
        listNeighbors{k} = [(k-nYPartition); ...
            (k-nYPartition+1);  ...
            (k+1); ...
            (k+nYPartition); ...
            (k+nYPartition+1)     ];
        
    elseif (mod(k, nYPartition) == 0)
        listNeighbors{k} = [(k-nYPartition-1); ...
            (k-nYPartition); ...
            (k-1); ...
            (k+nYPartition-1); ...
            (k+nYPartition)       ];
        
    elseif ( (k > numberAps-nYPartition+1) && (k < numberAps) )
        listNeighbors{k} = [(k-nYPartition-1); ...
            (k-nYPartition); ...
            (k-nYPartition+1); ...
            (k-1); ...
            (k+1)                 ];
        
    else
        listNeighbors{k} = [(k-nYPartition-1); ...
            (k-nYPartition); ...
            (k-nYPartition+1); ...
            (k-1); ...
            (k+1); ...
            (k+nYPartition-1); ...
            (k+nYPartition); ...
            (k+nYPartition+1)     ];
    end
    
end


end





%% ======================== AUXILIARY FUNCTIONS 3 ========================

% Among the nearby access points, it is chosen which are the most suitable
% for relocating positions. The basic idea is to prioritize the indexes
% that are most repeated.
function shiftedAps = best_shifted_aps(pairsCloseNeighbors)

if isempty(pairsCloseNeighbors)
    shiftedAps = [];
    return
end

frequency =  histcounts(pairsCloseNeighbors(:), (0:max(max(pairsCloseNeighbors))) + 1/2).';
pairsCloseNeighbors = sortrows([ max(frequency(pairsCloseNeighbors(:, 1)), frequency(pairsCloseNeighbors(:, 2))) pairsCloseNeighbors ], 1, 'descend');
pairsCloseNeighbors(:, 1) = [];

nPairsCloseNeighbors = size(pairsCloseNeighbors, 1);
shiftedAps = zeros(1, nPairsCloseNeighbors);
for i = 1:nPairsCloseNeighbors
    Ap1Priority = frequency(pairsCloseNeighbors(i, 1));
    Ap2Priority = frequency(pairsCloseNeighbors(i, 2));
    
    if (Ap1Priority == Ap2Priority)
        if (Ap1Priority == 1)
            shiftedAps(i) = pairsCloseNeighbors(i, randi(2));
        else
            balance1 = length(find(shiftedAps(1:i-1) == pairsCloseNeighbors(i, 1))) - length(find(pairsCloseNeighbors(1:i-1, :) == pairsCloseNeighbors(i, 1)));
            balance2 = length(find(shiftedAps(1:i-1) == pairsCloseNeighbors(i, 2))) - length(find(pairsCloseNeighbors(1:i-1, :) == pairsCloseNeighbors(i, 2)));
            
            if (balance1 > balance2)
                shiftedAps(i) = pairsCloseNeighbors(i, 1);
            else
                shiftedAps(i) = pairsCloseNeighbors(i, 2);
            end
            
        end
        
    elseif (Ap1Priority > Ap2Priority)
        shiftedAps(i) = pairsCloseNeighbors(i, 1);
    else
        shiftedAps(i) = pairsCloseNeighbors(i, 2);
    end
    
end

shiftedAps = unique(shiftedAps);

end

function gainOverNoisedB = AF_LargeScaleGainOverNoise(model,fc_GHz,noiseVariancedBm,distancesUEAP_m,h_AP_m,h_UE_m)
heigth_diference=h_AP_m-h_UE_m;
if strcmp(model,'Umi')==1
    %outdoor street canyon
    
    %variables declaration
    distance2D=sqrt(distancesUEAP_m.^2-heigth_diference^2);
    selector_LOS_NLOS=rand(size(distancesUEAP_m));
    PrLOS=zeros(size(distancesUEAP_m));
    PL=zeros(size(distancesUEAP_m));
    d_bp_line=4*heigth_diference*fc_GHz*10^9/(300*10^6);
    
    %set LOS probabilities for all cases on model
    PrLOS(distance2D<=18)=1;
    if sum(PrLOS)~=length(distance2D)
        %PrLOS(distance2D>18)=ones(size(distance2D))*18./distance2D+exp(-distance2D(distance2D>18)/36).*(1-(ones(size(distance2D))*18./distance2D));
        PrLOS(distance2D>18)=ones(size(distance2D(distance2D>18)))*18./distance2D(distance2D>18)+exp(-distance2D(distance2D>18)/36).*(1-(ones(size(distance2D(distance2D>18)))*18./distance2D(distance2D>18)));
    end
    %select 1 for LOS and 0 for NLOS
    selector_LOS_NLOS=(PrLOS>selector_LOS_NLOS);
    
    %Path loss calculations
    PL((selector_LOS_NLOS==1)==(distance2D<=d_bp_line))=32.4+21*log10(distancesUEAP_m((selector_LOS_NLOS==1)==(distance2D<=d_bp_line)))+20*log10(fc_GHz);
    PL((selector_LOS_NLOS==1)==(distance2D>d_bp_line))=32.4+40*log10(distancesUEAP_m((selector_LOS_NLOS==1)==(distance2D>d_bp_line)))+20*log10(fc_GHz)...
        -9.5*log10(d_bp_line^2+heigth_diference^2);
    %auxiliar variables for ease the calculation of NLOS path Loss
    aux=zeros(size(distancesUEAP_m));
    aux((selector_LOS_NLOS==0)==(distance2D<=d_bp_line))=32.4+21*log10(distancesUEAP_m((selector_LOS_NLOS==0)==(distance2D<=d_bp_line)))+20*log10(fc_GHz);
    aux((selector_LOS_NLOS==0)==(distance2D>d_bp_line))=32.4+40*log10(distancesUEAP_m((selector_LOS_NLOS==0)==(distance2D>d_bp_line)))+20*log10(fc_GHz)...
        -9.5*log10(d_bp_line^2+heigth_diference^2);
    aux2=zeros(size(distancesUEAP_m));
    aux2((selector_LOS_NLOS==0))=22.4+35.3*log10(distancesUEAP_m((selector_LOS_NLOS==0)))+21.3*log10(fc_GHz)...
        -0.3*(h_UE_m-1.5);
    PL((selector_LOS_NLOS==0))=max(aux(selector_LOS_NLOS==0),aux2(selector_LOS_NLOS==0));
    
    %Gain over Noise calculation
    gainOverNoisedB=-PL-noiseVariancedBm;
elseif strcmp(model,'InH-Mixed')==1
    %Used for Offices with divisories
    
    %variables declaration
    distance2D=sqrt(distancesUEAP_m.^2-heigth_diference^2);
    selector_LOS_NLOS=rand(size(distancesUEAP_m));
    PrLOS=zeros(size(distancesUEAP_m));
    PL=zeros(size(distancesUEAP_m));
    
    %set LOS probabilities for all cases on model
    PrLOS(distance2D<=1.2)=1;
    PrLOS((distance2D<=6.5)==(distance2D>1.2))=exp(-(distance2D((distance2D<=6.5)==(distance2D>1.2))-1.2)/4.7);
    PrLOS(distance2D>6.5)=exp(-(distance2D(distance2D>6.5)-6.5)/32.6)*0.32;
    
    %select 1 for LOS and 0 for NLOS
    selector_LOS_NLOS=(PrLOS>selector_LOS_NLOS);
    
    %Path loss calculations
    PL(selector_LOS_NLOS==1)=32.4+17.3*log10(distancesUEAP_m(selector_LOS_NLOS==1))+20*log10(fc_GHz);
    PL(selector_LOS_NLOS==0)=max(17.3+38.3*log10(distancesUEAP_m(selector_LOS_NLOS==0))+24.9*log10(fc_GHz)...
        ,32.4+17.3*log10(distancesUEAP_m(selector_LOS_NLOS==0))+20*log10(fc_GHz));
    
    %Gain over Noise calculation
    gainOverNoisedB=-PL-noiseVariancedBm;
elseif strcmp(model,'InH-Open')==1
    %Used for open Offices
    
    %variables declaration
    distance2D=sqrt(distancesUEAP_m.^2-heigth_diference^2);
    selector_LOS_NLOS=rand(size(distancesUEAP_m));
    PrLOS=zeros(size(distancesUEAP_m));
    PL=zeros(size(distancesUEAP_m));
    
    %set LOS probabilities for all cases on model
    PrLOS(distance2D<=5)=1;
    PrLOS((distance2D<=49)==(distance2D>5))=exp(-(distance2D((distance2D<=49)==(distance2D>5))-5)/70.8);
    PrLOS(distance2D>49)=exp(-(distance2D(distance2D>49)-49)/211.17)*0.54;
    
    %select 1 for LOS and 0 for NLOS
    selector_LOS_NLOS=(PrLOS>selector_LOS_NLOS);
    
    %Path loss calculations
    PL(selector_LOS_NLOS==1)=32.4+17.3*log10(distancesUEAP_m(selector_LOS_NLOS==1))+20*log10(fc_GHz);
    PL(selector_LOS_NLOS==0)=max(17.3+38.3*log10(distancesUEAP_m(selector_LOS_NLOS==0))+24.9*log10(fc_GHz)...
        ,32.4+17.3*log10(distancesUEAP_m(selector_LOS_NLOS==0))+20*log10(fc_GHz));
    
    %Gain over Noise calculation
    gainOverNoisedB=-PL-noiseVariancedBm;
elseif strcmp(model,'InF-SL')==1
    %Used for Big machineries composed of regular metallic surfaces.
    %For example: several mixed production areas with open spaces and storage/commissioning areas
    
    %variables declaration
    distance2D=sqrt(distancesUEAP_m.^2-heigth_diference^2);
    d_clutter_m=10;
    clutter_density=0.2;
    selector_LOS_NLOS=rand(size(distancesUEAP_m));
    
    %set LOS probabilities for all cases on model
    PrLOS=exp(-distance2D/(-d_clutter_m/log(1-clutter_density)));
    
    %select 1 for LOS and 0 for NLOS
    selector_LOS_NLOS=(PrLOS>selector_LOS_NLOS);
    
    %Path loss calculations
    PL(selector_LOS_NLOS==1)=31.84+25.5*log10(distancesUEAP_m(selector_LOS_NLOS==1))+19*log10(fc_GHz);
    PL(selector_LOS_NLOS==0)=max(33+25.5*log10(distancesUEAP_m(selector_LOS_NLOS==0))+20*log10(fc_GHz)...
        ,31.84+25.5*log10(distancesUEAP_m(selector_LOS_NLOS==0))+19*log10(fc_GHz));
    
    %Gain over Noise calculation
    gainOverNoisedB=-PL-noiseVariancedBm;
elseif strcmp(model,'InF-DL')==1
    %Used for Small to medium metallic machinery and objects with irregular structure.
    %For example: assembly and production lines surrounded by mixed small-sized machineries.
    
    %variables declaration
    d_clutter_m=2;
    distance2D=sqrt(distancesUEAP_m.^2-heigth_diference^2);
    clutter_density=0.7;
    selector_LOS_NLOS=rand(size(distancesUEAP_m));
    
    %set LOS probabilities for all cases on model
    PrLOS=exp(-distance2D/(-d_clutter_m/log(1-clutter_density)));
    
    %select 1 for LOS and 0 for NLOS
    selector_LOS_NLOS=(PrLOS>selector_LOS_NLOS);
    
    %Path loss calculations
    PL(selector_LOS_NLOS==1)=31.84+25.5*log10(distancesUEAP_m(selector_LOS_NLOS==1))+19*log10(fc_GHz);
    PL(selector_LOS_NLOS==0)=max(18.6+35.7*log10(distancesUEAP_m(selector_LOS_NLOS==0))+20*log10(fc_GHz)...
        ,31.84+25.5*log10(distancesUEAP_m(selector_LOS_NLOS==0))+19*log10(fc_GHz));
    
    %Gain over Noise calculation
    gainOverNoisedB=-PL-noiseVariancedBm;
else
    
end
end



