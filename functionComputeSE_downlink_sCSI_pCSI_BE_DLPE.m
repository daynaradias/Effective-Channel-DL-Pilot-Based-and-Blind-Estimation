function  [SE_MR, SE_LP_MMSE, SE_P_MMSE, SE_P_RZF, ...
    NMSE_MR, NMSE_LP_MMSE, NMSE_P_MMSE, NMSE_P_RZF,...
    SE_MR_pCSI, SE_LP_MMSE_pCSI, SE_P_MMSE_pCSI, SE_P_RZF_pCSI,...
    SE_MR_BE,SE_LP_MMSE_BE,SE_P_MMSE_BE,SE_P_RZF_BE,...
    NMSE_MR_BE,NMSE_LP_MMSE_BE,NMSE_P_MMSE_BE,NMSE_P_RZF_BE,...
    SE_MR_BE_InfTau,SE_LP_MMSE_BE_InfTau,SE_P_MMSE_BE_InfTau,SE_P_RZF_BE_InfTau,...
    NMSE_MR_BE_InfTau,NMSE_LP_MMSE_BE_InfTau,NMSE_P_MMSE_BE_InfTau,NMSE_P_RZF_BE_InfTau,...
    SE_MR_DLPE,SE_LP_MMSE_DLPE,SE_P_MMSE_DLPE,SE_P_RZF_DLPE,...
    NMSE_MR_DLPE,NMSE_LP_MMSE_DLPE,NMSE_P_MMSE_DLPE,NMSE_P_RZF_DLPE]...
    = functionComputeSE_downlink_sCSI_pCSI_BE_DLPE(capacityBound,tau_c,tau_up,pilotIndexUL,tau_dp,nbrOfRealizations,K,pilotAssignmentMethod,...
    CHandFP_MR,CHandFP_LP_MMSE,CHandFP_P_MMSE,CHandFP_P_RZF,...
    interUserGains_MR,interUserGains_LP_MMSE,interUserGains_P_MMSE,interUserGains_P_RZF)
%% DL pilot estimation

%Compute the prelog factor assuming only downlink data transmission
prelogFactor_DLPE = (1 - (tau_up + tau_dp)/tau_c);


if strcmp(pilotAssignmentMethod,'nonJoint_minContamination')
    %Min pilot contamination DL pilot reuse
    [~,W_MR] = pilotAssignment(K,tau_dp,CHandFP_MR);
    [~,W_LP_MMSE] = pilotAssignment(K,tau_dp,CHandFP_LP_MMSE);
    [~,W_P_MMSE] = pilotAssignment(K,tau_dp,CHandFP_P_MMSE);
    [~,W_P_RZF] = pilotAssignment(K,tau_dp,CHandFP_P_RZF);
elseif strcmp(pilotAssignmentMethod,'nonJoint_random')
    pilotIndexDL = [1:tau_dp randi(tau_dp,1,K-tau_dp)]';
    W = zeros(K,K);
    for i =1:K
        for j =1:K
            W (i,j) = isequal (pilotIndexDL(i),  pilotIndexDL(j));
        end
    end
    W_MR = W;
    W_LP_MMSE = W;
    W_P_MMSE = W;
    W_P_RZF = W;
elseif strcmp(pilotAssignmentMethod,'joint_minContamination')
    %Min pilot contamination DL pilot reuse
    [~,W_MR] = jointULandDLpilotAssignment(K,tau_dp,CHandFP_MR,pilotIndexUL);%
    [~,W_LP_MMSE] = jointULandDLpilotAssignment(K,tau_dp,CHandFP_LP_MMSE,pilotIndexUL);%
    [~,W_P_MMSE] = jointULandDLpilotAssignment(K,tau_dp,CHandFP_P_MMSE,pilotIndexUL);%
    [~,W_P_RZF] = jointULandDLpilotAssignment(K,tau_dp,CHandFP_P_RZF,pilotIndexUL);%
elseif strcmp(pilotAssignmentMethod,'joint_random')
    
    pilotIndexDL = [1:tau_dp randi(tau_dp,1,K-tau_dp)]';
    
    for i = 1:K
        for j=1:K
            if i ~= j
                aux =  (isequal(pilotIndexUL(i),pilotIndexUL(j)) + ...
                    isequal(pilotIndexDL(i),pilotIndexDL(j))) == 2;
                
                if aux == 1
                    
                    pilotIndexDL(j) = randi(tau_dp,1);
                    
                end
            end
        end
    end
        
    W = zeros(K,K);
    for i =1:K
        for j =1:K
            W (i,j) = isequal (pilotIndexDL(i),  pilotIndexDL(j));
        end
    end
    W_MR = W;
    W_LP_MMSE = W;
    W_P_MMSE = W;
    W_P_RZF = W;
end

[pilotEstimatedEffectiveChannel_MR,NMSE_MR_DLPE] = DLbeamformingTraining(K,nbrOfRealizations,tau_dp,W_MR,interUserGains_MR);
[pilotEstimatedEffectiveChannel_LP_MMSE,NMSE_LP_MMSE_DLPE] = DLbeamformingTraining(K,nbrOfRealizations,tau_dp,W_LP_MMSE,interUserGains_LP_MMSE);
[pilotEstimatedEffectiveChannel_P_MMSE,NMSE_P_MMSE_DLPE] = DLbeamformingTraining(K,nbrOfRealizations,tau_dp,W_P_MMSE,interUserGains_P_MMSE);
[pilotEstimatedEffectiveChannel_P_RZF,NMSE_P_RZF_DLPE] = DLbeamformingTraining(K,nbrOfRealizations,tau_dp,W_P_RZF,interUserGains_P_RZF);

switch capacityBound
    
    case 'Ngo2017'
        % (Ngo, 2017) capacity lower bound
        SE_MR_DLPE = capacity_bound_ngo(K,nbrOfRealizations,pilotEstimatedEffectiveChannel_MR,interUserGains_MR);        
        SE_LP_MMSE_DLPE = capacity_bound_ngo(K,nbrOfRealizations,pilotEstimatedEffectiveChannel_LP_MMSE,interUserGains_LP_MMSE);
        SE_P_MMSE_DLPE = capacity_bound_ngo(K,nbrOfRealizations,pilotEstimatedEffectiveChannel_P_MMSE,interUserGains_P_MMSE);
        SE_P_RZF_DLPE = capacity_bound_ngo(K,nbrOfRealizations,pilotEstimatedEffectiveChannel_P_RZF,interUserGains_P_RZF);

    case 'Ngo2013'
        % (Ngo, 2013) capacity lower bound
        [SE_MR_DLPE,~] = capacity_bound_ngo_2013(K,nbrOfRealizations,pilotEstimatedEffectiveChannel_MR,interUserGains_MR); 
        [SE_LP_MMSE_DLPE,~] = capacity_bound_ngo_2013(K,nbrOfRealizations,pilotEstimatedEffectiveChannel_LP_MMSE,interUserGains_LP_MMSE);        
        [SE_P_MMSE_DLPE,~] = capacity_bound_ngo_2013(K,nbrOfRealizations,pilotEstimatedEffectiveChannel_P_MMSE,interUserGains_P_MMSE);
        [SE_P_RZF_DLPE,~] = capacity_bound_ngo_2013(K,nbrOfRealizations,pilotEstimatedEffectiveChannel_P_RZF,interUserGains_P_RZF);
    
    case 'Ngo2013_2'
        % (Ngo, 2013) capacity lower bound 2        
        SE_MR_DLPE = capacity_lower_bound_ngo(K,nbrOfRealizations,pilotEstimatedEffectiveChannel_MR,interUserGains_MR);       
        SE_LP_MMSE_DLPE = capacity_lower_bound_ngo(K,nbrOfRealizations,pilotEstimatedEffectiveChannel_LP_MMSE,interUserGains_LP_MMSE);
        SE_P_MMSE_DLPE = capacity_lower_bound_ngo(K,nbrOfRealizations,pilotEstimatedEffectiveChannel_P_MMSE,interUserGains_P_MMSE);
        SE_P_RZF_DLPE = capacity_bound_ngo(K,nbrOfRealizations,pilotEstimatedEffectiveChannel_P_RZF,interUserGains_P_RZF);

    case 'Hardening'
        
        SE_MR_DLPE = zeros(K,1);
        SE_LP_MMSE_DLPE = zeros(K,1);
        SE_P_MMSE_DLPE = zeros(K,1);
        SE_P_RZF_DLPE = zeros(K,1);
        
        for k=1:K
            SE_MR_DLPE(k) = mean(log2( 1 + (abs(pilotEstimatedEffectiveChannel_MR(k,:)).^2) ./ (1+ var( reshape (sum(interUserGains_MR(k,:,:),2), [1,nbrOfRealizations]) - pilotEstimatedEffectiveChannel_MR(k,:) ,1,2) )));
            SE_LP_MMSE_DLPE(k) = mean(log2( 1 + (abs(pilotEstimatedEffectiveChannel_LP_MMSE(k,:)).^2) ./ (1+ var( reshape (sum(interUserGains_LP_MMSE(k,:,:),2), [1,nbrOfRealizations]) - pilotEstimatedEffectiveChannel_LP_MMSE(k,:) ,1,2) )));
            SE_P_MMSE_DLPE(k) = mean(log2( 1 + (abs(pilotEstimatedEffectiveChannel_P_MMSE(k,:)).^2) ./ (1+ var( reshape (sum(interUserGains_P_MMSE(k,:,:),2), [1,nbrOfRealizations]) - pilotEstimatedEffectiveChannel_P_MMSE(k,:) ,1,2) )));
            SE_P_RZF_DLPE(k) = mean(log2( 1 + (abs(pilotEstimatedEffectiveChannel_P_RZF(k,:)).^2) ./ (1+ var( reshape (sum(interUserGains_P_RZF(k,:,:),2), [1,nbrOfRealizations]) - pilotEstimatedEffectiveChannel_P_RZF(k,:) ,1,2) )));
        end
             
    case 'Interdonato2019'
        % (Interdonato, 2019) capacity bound
        [SE_MR_DLPE,~] = capacity_bound_interdonato(K,pilotEstimatedEffectiveChannel_MR,interUserGains_MR);
        [SE_LP_MMSE_DLPE,~] = capacity_bound_interdonato(K,pilotEstimatedEffectiveChannel_LP_MMSE,interUserGains_LP_MMSE);       
        [SE_P_MMSE_DLPE,~] = capacity_bound_interdonato(K,pilotEstimatedEffectiveChannel_P_MMSE,interUserGains_P_MMSE);
        [SE_P_RZF_DLPE,~] = capacity_bound_interdonato(K,pilotEstimatedEffectiveChannel_P_RZF,interUserGains_P_RZF);

    case 'UnF'
        % Use and then forget capacity bound
        [SE_MR_DLPE,~] = use_and_then_forget_bound(K,pilotEstimatedEffectiveChannel_MR,interUserGains_MR); 
        [SE_LP_MMSE_DLPE,~] = use_and_then_forget_bound(K,pilotEstimatedEffectiveChannel_LP_MMSE,interUserGains_LP_MMSE);
        [SE_P_MMSE_DLPE,~] = use_and_then_forget_bound(K,pilotEstimatedEffectiveChannel_P_MMSE,interUserGains_P_MMSE);
        [SE_P_RZF_DLPE,~] = use_and_then_forget_bound(K,pilotEstimatedEffectiveChannel_P_RZF,interUserGains_P_RZF);
end

%% Blind estimation        
%Compute the prelog factor assuming only downlink data transmission
prelogFactor = (1-tau_up/tau_c);
tau_d = tau_c - tau_up;

[blind_estimated_effective_channel_MR,NMSE_MR_BE,blind_estimated_effective_channel_MR_InfTau,NMSE_MR_BE_InfTau]= blind_estimation(tau_d,K,interUserGains_MR,nbrOfRealizations);
[blind_estimated_effective_channel_LP_MMSE,NMSE_LP_MMSE_BE,blind_estimated_effective_channel_LP_MMSE_InfTau,NMSE_LP_MMSE_BE_InfTau]= blind_estimation(tau_d,K,interUserGains_LP_MMSE,nbrOfRealizations);
[blind_estimated_effective_channel_P_MMSE,NMSE_P_MMSE_BE,blind_estimated_effective_channel_P_MMSE_InfTau,NMSE_P_MMSE_BE_InfTau]= blind_estimation(tau_d,K,interUserGains_P_MMSE,nbrOfRealizations);
[blind_estimated_effective_channel_P_RZF,NMSE_P_RZF_BE,blind_estimated_effective_channel_P_RZF_InfTau,NMSE_P_RZF_BE_InfTau]= blind_estimation(tau_d,K,interUserGains_P_RZF,nbrOfRealizations);


switch capacityBound
    
    case 'Ngo2017'
        % (Ngo, 2017) capacity lower bound
        SE_MR_BE = capacity_bound_ngo(K,nbrOfRealizations,blind_estimated_effective_channel_MR,interUserGains_MR);
        SE_LP_MMSE_BE = capacity_bound_ngo(K,nbrOfRealizations,blind_estimated_effective_channel_LP_MMSE,interUserGains_LP_MMSE);        
        SE_P_MMSE_BE = capacity_bound_ngo(K,nbrOfRealizations,blind_estimated_effective_channel_P_MMSE,interUserGains_P_MMSE);
        SE_P_RZF_BE = capacity_bound_ngo(K,nbrOfRealizations,blind_estimated_effective_channel_P_RZF,interUserGains_P_RZF);

    case  'Ngo2013'
        % (Ngo, 2013) capacity lower bound
        [SE_MR_BE,~] = capacity_bound_ngo_2013(K,nbrOfRealizations,blind_estimated_effective_channel_MR,interUserGains_MR);
        [SE_LP_MMSE_BE,~] = capacity_bound_ngo_2013(K,nbrOfRealizations,blind_estimated_effective_channel_LP_MMSE,interUserGains_LP_MMSE);       
        [SE_P_MMSE_BE,~] = capacity_bound_ngo_2013(K,nbrOfRealizations,blind_estimated_effective_channel_P_MMSE,interUserGains_P_MMSE);
        [SE_P_RZF_BE,~] = capacity_bound_ngo_2013(K,nbrOfRealizations,blind_estimated_effective_channel_P_RZF,interUserGains_P_RZF);
    
    case  'Ngo2013_2'
        % (Ngo, 2013) capacity lower bound  2
        SE_MR_BE = capacity_lower_bound_ngo(K,nbrOfRealizations,blind_estimated_effective_channel_MR,interUserGains_MR);        
        SE_LP_MMSE_BE = capacity_lower_bound_ngo(K,nbrOfRealizations,blind_estimated_effective_channel_LP_MMSE,interUserGains_LP_MMSE);
        SE_P_MMSE_BE = capacity_lower_bound_ngo(K,nbrOfRealizations,blind_estimated_effective_channel_P_MMSE,interUserGains_P_MMSE);
        SE_P_RZF_BE = capacity_lower_bound_ngo(K,nbrOfRealizations,blind_estimated_effective_channel_P_RZF,interUserGains_P_RZF);

    case 'Hardening'
        
        SE_MR_BE = zeros(K,1);
        SE_LP_MMSE_BE = zeros(K,1);
        SE_P_MMSE_BE = zeros(K,1);
        SE_P_RZF_BE = zeros(K,1);
        
        for k=1:K
            SE_MR_BE(k) =mean(log2( 1 + (abs(blind_estimated_effective_channel_MR(k,:)).^2) ./ (1+ var( reshape (sum(interUserGains_MR(k,:,:),2), [1,nbrOfRealizations]) - blind_estimated_effective_channel_MR(k,:) ,1,2) )));
            SE_LP_MMSE_BE(k) = mean(log2( 1 + (abs(blind_estimated_effective_channel_LP_MMSE(k,:)).^2) ./ (1+ var( reshape (sum(interUserGains_LP_MMSE(k,:,:),2), [1,nbrOfRealizations]) - blind_estimated_effective_channel_LP_MMSE(k,:) ,1,2) )));
            SE_P_MMSE_BE(k) = mean(log2( 1 + (abs(blind_estimated_effective_channel_P_MMSE(k,:)).^2) ./ (1+ var( reshape (sum(interUserGains_P_MMSE(k,:,:),2), [1,nbrOfRealizations]) - blind_estimated_effective_channel_P_MMSE(k,:) ,1,2) )));
            SE_P_RZF_BE(k) = mean(log2( 1 + (abs(blind_estimated_effective_channel_P_RZF(k,:)).^2) ./ (1+ var( reshape (sum(interUserGains_P_RZF(k,:,:),2), [1,nbrOfRealizations]) - blind_estimated_effective_channel_P_RZF(k,:) ,1,2) )));

        end
        
    case 'Interdonato2019'
        % (Interdonato, 2019) capacity bound    
        [SE_MR_BE,~] = capacity_bound_interdonato(K,blind_estimated_effective_channel_MR,interUserGains_MR);     
        [SE_LP_MMSE_BE,~] = capacity_bound_interdonato(K,blind_estimated_effective_channel_LP_MMSE,interUserGains_LP_MMSE);     
        [SE_P_MMSE_BE,~] = capacity_bound_interdonato(K,blind_estimated_effective_channel_P_MMSE,interUserGains_P_MMSE);
        [SE_P_RZF_BE,~] = capacity_bound_interdonato(K,blind_estimated_effective_channel_P_RZF,interUserGains_P_RZF);

    case 'UnF'
        %Use and then forget capacity bound
        [SE_MR_BE,~] = use_and_then_forget_bound(K,blind_estimated_effective_channel_MR,interUserGains_MR);       
        [SE_LP_MMSE_BE,~] = use_and_then_forget_bound(K,blind_estimated_effective_channel_LP_MMSE,interUserGains_LP_MMSE);       
        [SE_P_MMSE_BE,~] = use_and_then_forget_bound(K,blind_estimated_effective_channel_P_MMSE,interUserGains_P_MMSE);
        [SE_P_RZF_BE,~] = use_and_then_forget_bound(K,blind_estimated_effective_channel_P_RZF,interUserGains_P_RZF);      
end

switch capacityBound
    
    case 'Ngo2017'
        % (Ngo, 2017) capacity lower bound
        SE_MR_BE_InfTau = capacity_bound_ngo(K,nbrOfRealizations,blind_estimated_effective_channel_MR_InfTau,interUserGains_MR);
        SE_LP_MMSE_BE_InfTau = capacity_bound_ngo(K,nbrOfRealizations,blind_estimated_effective_channel_LP_MMSE_InfTau,interUserGains_LP_MMSE);        
        SE_P_MMSE_BE_InfTau = capacity_bound_ngo(K,nbrOfRealizations,blind_estimated_effective_channel_P_MMSE_InfTau,interUserGains_P_MMSE);
        SE_P_RZF_BE_InfTau = capacity_bound_ngo(K,nbrOfRealizations,blind_estimated_effective_channel_P_RZF_InfTau,interUserGains_P_RZF);

    case  'Ngo2013'
        % (Ngo, 2013) capacity lower bound
        [SE_MR_BE_InfTau,~] = capacity_bound_ngo_2013(K,nbrOfRealizations,blind_estimated_effective_channel_MR_InfTau,interUserGains_MR);
        [SE_LP_MMSE_BE_InfTau,~] = capacity_bound_ngo_2013(K,nbrOfRealizations,blind_estimated_effective_channel_LP_MMSE_InfTau,interUserGains_LP_MMSE);       
        [SE_P_MMSE_BE_InfTau,~] = capacity_bound_ngo_2013(K,nbrOfRealizations,blind_estimated_effective_channel_P_MMSE_InfTau,interUserGains_P_MMSE);
        [SE_P_RZF_BE_InfTau,~] = capacity_bound_ngo_2013(K,nbrOfRealizations,blind_estimated_effective_channel_P_RZF_InfTau,interUserGains_P_RZF);
    
    case  'Ngo2013_2'
        % (Ngo, 2013) capacity lower bound  2
        SE_MR_BE_InfTau = capacity_lower_bound_ngo(K,nbrOfRealizations,blind_estimated_effective_channel_MR_InfTau,interUserGains_MR);        
        SE_LP_MMSE_BE_InfTau = capacity_lower_bound_ngo(K,nbrOfRealizations,blind_estimated_effective_channel_LP_MMSE_InfTau,interUserGains_LP_MMSE);
        SE_P_MMSE_BE_InfTau = capacity_lower_bound_ngo(K,nbrOfRealizations,blind_estimated_effective_channel_P_MMSE_InfTau,interUserGains_P_MMSE);
        SE_P_RZF_BE_InfTau = capacity_lower_bound_ngo(K,nbrOfRealizations,blind_estimated_effective_channel_P_RZF_InfTau,interUserGains_P_RZF);

    case 'Hardening'
        
        SE_MR_BE_InfTau = zeros(K,1);
        SE_LP_MMSE_BE_InfTau = zeros(K,1);
        SE_P_MMSE_BE_InfTau = zeros(K,1);
        SE_P_RZF_BE_InfTau = zeros(K,1);
        
        for k=1:K
            SE_MR_BE_InfTau(k) =mean(log2( 1 + (abs(blind_estimated_effective_channel_MR_InfTau(k,:)).^2) ./ (1+ var( reshape (sum(interUserGains_MR(k,:,:),2), [1,nbrOfRealizations]) - blind_estimated_effective_channel_MR_InfTau(k,:) ,1,2) )));
            SE_LP_MMSE_BE_InfTau(k) = mean(log2( 1 + (abs(blind_estimated_effective_channel_LP_MMSE_InfTau(k,:)).^2) ./ (1+ var( reshape (sum(interUserGains_LP_MMSE(k,:,:),2), [1,nbrOfRealizations]) - blind_estimated_effective_channel_LP_MMSE_InfTau(k,:) ,1,2) )));
            SE_P_MMSE_BE_InfTau(k) = mean(log2( 1 + (abs(blind_estimated_effective_channel_P_MMSE_InfTau(k,:)).^2) ./ (1+ var( reshape (sum(interUserGains_P_MMSE(k,:,:),2), [1,nbrOfRealizations]) - blind_estimated_effective_channel_P_MMSE_InfTau(k,:) ,1,2) )));
            SE_P_RZF_BE_InfTau(k) = mean(log2( 1 + (abs(blind_estimated_effective_channel_P_RZF_InfTau(k,:)).^2) ./ (1+ var( reshape (sum(interUserGains_P_RZF(k,:,:),2), [1,nbrOfRealizations]) - blind_estimated_effective_channel_P_RZF_InfTau(k,:) ,1,2) )));

        end
        
    case 'Interdonato2019'
        % (Interdonato, 2019) capacity bound    
        [SE_MR_BE_InfTau,~] = capacity_bound_interdonato(K,blind_estimated_effective_channel_MR_InfTau,interUserGains_MR);     
        [SE_LP_MMSE_BE_InfTau,~] = capacity_bound_interdonato(K,blind_estimated_effective_channel_LP_MMSE_InfTau,interUserGains_LP_MMSE);     
        [SE_P_MMSE_BE_InfTau,~] = capacity_bound_interdonato(K,blind_estimated_effective_channel_P_MMSE_InfTau,interUserGains_P_MMSE);
        [SE_P_RZF_BE_InfTau,~] = capacity_bound_interdonato(K,blind_estimated_effective_channel_P_RZF_InfTau,interUserGains_P_RZF);

    case 'UnF'
        %Use and then forget capacity bound
        [SE_MR_BE_InfTau,~] = use_and_then_forget_bound(K,blind_estimated_effective_channel_MR_InfTau,interUserGains_MR);       
        [SE_LP_MMSE_BE_InfTau,~] = use_and_then_forget_bound(K,blind_estimated_effective_channel_LP_MMSE_InfTau,interUserGains_LP_MMSE);       
        [SE_P_MMSE_BE_InfTau,~] = use_and_then_forget_bound(K,blind_estimated_effective_channel_P_MMSE_InfTau,interUserGains_P_MMSE);
        [SE_P_RZF_BE_InfTau,~] = use_and_then_forget_bound(K,blind_estimated_effective_channel_P_RZF_InfTau,interUserGains_P_RZF);      
end
%% sCSI and pCSI
%for sCSI
estimated_effective_channel_MR = diag(mean(interUserGains_MR,3)).*ones(K,nbrOfRealizations);
estimated_effective_channel_LP_MMSE = diag(mean(interUserGains_LP_MMSE,3)).*ones(K,nbrOfRealizations);
estimated_effective_channel_P_MMSE = diag(mean(interUserGains_P_MMSE,3)).*ones(K,nbrOfRealizations);
estimated_effective_channel_P_RZF = diag(mean(interUserGains_P_RZF,3)).*ones(K,nbrOfRealizations);

NMSE_MR = diag(mean(abs(mean(interUserGains_MR,3).*ones(K,K,nbrOfRealizations) - interUserGains_MR).^2,3)./mean(abs(interUserGains_MR).^2,3));
NMSE_LP_MMSE = diag(mean(abs(mean(interUserGains_LP_MMSE,3).*ones(K,K,nbrOfRealizations) - interUserGains_LP_MMSE).^2,3)./mean(abs(interUserGains_LP_MMSE).^2,3));
NMSE_P_MMSE = diag(mean(abs(mean(interUserGains_P_MMSE,3).*ones(K,K,nbrOfRealizations) - interUserGains_P_MMSE).^2,3)./mean(abs(interUserGains_P_MMSE).^2,3));
NMSE_P_RZF = diag(mean(abs(mean(interUserGains_P_RZF,3).*ones(K,K,nbrOfRealizations) - interUserGains_P_RZF).^2,3)./mean(abs(interUserGains_P_RZF).^2,3));

%for pCSI
perfectEstimatedEffectiveChannel_MR = zeros(K,nbrOfRealizations);
perfectEstimatedEffectiveChannel_LP_MMSE = zeros(K,nbrOfRealizations);
perfectEstimatedEffectiveChannel_P_MMSE  = zeros(K,nbrOfRealizations);
perfectEstimatedEffectiveChannel_P_RZF  = zeros(K,nbrOfRealizations);

for k = 1:K
    perfectEstimatedEffectiveChannel_MR(k,:) = interUserGains_MR(k,k,:);
    perfectEstimatedEffectiveChannel_LP_MMSE(k,:) = interUserGains_LP_MMSE(k,k,:);
    perfectEstimatedEffectiveChannel_P_MMSE(k,:) = interUserGains_P_MMSE(k,k,:);
    perfectEstimatedEffectiveChannel_P_RZF(k,:) = interUserGains_P_RZF(k,k,:);
end

switch capacityBound
    
    case 'Ngo2017'
        % (Ngo, 2017) capacity lower bound
        SE_MR = capacity_bound_ngo(K,nbrOfRealizations,estimated_effective_channel_MR,interUserGains_MR);
        SE_LP_MMSE = capacity_bound_ngo(K,nbrOfRealizations,estimated_effective_channel_LP_MMSE,interUserGains_LP_MMSE);      
        SE_P_MMSE = capacity_bound_ngo(K,nbrOfRealizations,estimated_effective_channel_P_MMSE,interUserGains_P_MMSE);
        SE_P_RZF = capacity_bound_ngo(K,nbrOfRealizations,estimated_effective_channel_P_RZF,interUserGains_P_RZF);

        SE_MR_pCSI = capacity_bound_ngo(K,nbrOfRealizations,perfectEstimatedEffectiveChannel_MR,interUserGains_MR);
        SE_LP_MMSE_pCSI = capacity_bound_ngo(K,nbrOfRealizations,perfectEstimatedEffectiveChannel_LP_MMSE,interUserGains_LP_MMSE);      
        SE_P_MMSE_pCSI = capacity_bound_ngo(K,nbrOfRealizations,perfectEstimatedEffectiveChannel_P_MMSE,interUserGains_P_MMSE);
        SE_P_RZF_pCSI = capacity_bound_ngo(K,nbrOfRealizations,perfectEstimatedEffectiveChannel_P_RZF,interUserGains_P_RZF);

    case  'Ngo2013'
        % (Ngo, 2013) capacity lower bound
        SE_MR = capacity_bound_ngo_2013(K,nbrOfRealizations,estimated_effective_channel_MR,interUserGains_MR);        
        SE_LP_MMSE = capacity_bound_ngo_2013(K,nbrOfRealizations,estimated_effective_channel_LP_MMSE,interUserGains_LP_MMSE);      
        SE_P_MMSE = capacity_bound_ngo_2013(K,nbrOfRealizations,estimated_effective_channel_P_MMSE,interUserGains_P_MMSE);
        SE_P_RZF = capacity_bound_ngo_2013(K,nbrOfRealizations,estimated_effective_channel_P_RZF,interUserGains_P_RZF);

        %Compute SE with MR precoding, assuming blind estimation at the UE
        SE_MR_pCSI = capacity_bound_ngo_2013(K,nbrOfRealizations,perfectEstimatedEffectiveChannel_MR,interUserGains_MR);       
        SE_LP_MMSE_pCSI = capacity_bound_ngo_2013(K,nbrOfRealizations,perfectEstimatedEffectiveChannel_LP_MMSE,interUserGains_LP_MMSE);              
        SE_P_MMSE_pCSI = capacity_bound_ngo_2013(K,nbrOfRealizations,perfectEstimatedEffectiveChannel_P_MMSE,interUserGains_P_MMSE);
        SE_P_RZF_pCSI = capacity_bound_ngo_2013(K,nbrOfRealizations,perfectEstimatedEffectiveChannel_P_RZF,interUserGains_P_RZF);

    case  'Ngo2013_2'
        % (Ngo, 2013) capacity lower bound  2      
        %Compute SE with MR precoding, assuming blind estimation at the UE
        SE_MR = capacity_lower_bound_ngo(K,nbrOfRealizations,estimated_effective_channel_MR,interUserGains_MR);       
        SE_LP_MMSE = capacity_lower_bound_ngo(K,nbrOfRealizations,estimated_effective_channel_LP_MMSE,interUserGains_LP_MMSE);              
        SE_P_MMSE = capacity_lower_bound_ngo(K,nbrOfRealizations,estimated_effective_channel_P_MMSE,interUserGains_P_MMSE);
        SE_P_RZF = capacity_lower_bound_ngo(K,nbrOfRealizations,estimated_effective_channel_P_RZF,interUserGains_P_RZF);

        %Compute SE with MR precoding, assuming blind estimation at the UE
        SE_MR_pCSI = capacity_lower_bound_ngo(K,nbrOfRealizations,perfectEstimatedEffectiveChannel_MR,interUserGains_MR);        
        SE_LP_MMSE_pCSI = capacity_lower_bound_ngo(K,nbrOfRealizations,perfectEstimatedEffectiveChannel_LP_MMSE,interUserGains_LP_MMSE);              
        SE_P_MMSE_pCSI = capacity_lower_bound_ngo(K,nbrOfRealizations,perfectEstimatedEffectiveChannel_P_MMSE,interUserGains_P_MMSE);
        SE_P_RZF_pCSI = capacity_lower_bound_ngo(K,nbrOfRealizations,perfectEstimatedEffectiveChannel_P_RZF,interUserGains_P_RZF);

    case 'Hardening'
        
        SE_MR = zeros(K,1);
        SE_LP_MMSE = zeros(K,1);
        SE_P_MMSE = zeros(K,1);
        SE_P_RZF = zeros(K,1);

        for k=1:K
            SE_MR(k) =mean(log2( 1 + (abs(estimated_effective_channel_MR(k,:)).^2) ./ (1+ var( reshape (sum(interUserGains_MR(k,:,:),2), [1,nbrOfRealizations]) - estimated_effective_channel_MR(k,:) ,1,2) )));
            SE_LP_MMSE(k) = mean(log2( 1 + (abs(estimated_effective_channel_LP_MMSE(k,:)).^2) ./ (1+ var( reshape (sum(interUserGains_LP_MMSE(k,:,:),2), [1,nbrOfRealizations]) - estimated_effective_channel_LP_MMSE(k,:) ,1,2) )));
            SE_P_MMSE(k) = mean(log2( 1 + (abs(estimated_effective_channel_P_MMSE(k,:)).^2) ./ (1+ var( reshape (sum(interUserGains_P_MMSE(k,:,:),2), [1,nbrOfRealizations]) - estimated_effective_channel_P_MMSE(k,:) ,1,2) )));
            SE_P_RZF(k) = mean(log2( 1 + (abs(estimated_effective_channel_P_RZF(k,:)).^2) ./ (1+ var( reshape (sum(interUserGains_P_RZF(k,:,:),2), [1,nbrOfRealizations]) - estimated_effective_channel_P_RZF(k,:) ,1,2) )));
 
        end
        
        SE_MR_pCSI = zeros(K,1);
        SE_LP_MMSE_pCSI = zeros(K,1);
        SE_P_MMSE_pCSI = zeros(K,1);
        SE_P_RZF_pCSI = zeros(K,1);

        for k=1:K
            SE_MR_pCSI(k) =mean(log2( 1 + (abs(perfectEstimatedEffectiveChannel_MR(k,:)).^2) ./ (1+ var( reshape (sum(interUserGains_MR(k,:,:),2), [1,nbrOfRealizations]) - perfectEstimatedEffectiveChannel_MR(k,:) ,1,2) )));
            SE_LP_MMSE_pCSI(k) = mean(log2( 1 + (abs(perfectEstimatedEffectiveChannel_LP_MMSE(k,:)).^2) ./ (1+ var( reshape (sum(interUserGains_LP_MMSE(k,:,:),2), [1,nbrOfRealizations]) - perfectEstimatedEffectiveChannel_LP_MMSE(k,:) ,1,2) )));
            SE_P_MMSE_pCSI(k) = mean(log2( 1 + (abs(perfectEstimatedEffectiveChannel_P_MMSE(k,:)).^2) ./ (1+ var( reshape (sum(interUserGains_P_MMSE(k,:,:),2), [1,nbrOfRealizations]) - perfectEstimatedEffectiveChannel_P_MMSE(k,:) ,1,2) )));
            SE_P_RZF_pCSI(k) = mean(log2( 1 + (abs(perfectEstimatedEffectiveChannel_P_RZF(k,:)).^2) ./ (1+ var( reshape (sum(interUserGains_P_RZF(k,:,:),2), [1,nbrOfRealizations]) - perfectEstimatedEffectiveChannel_P_RZF(k,:) ,1,2) )));

        end
        
    case 'Interdonato2019'
        % (Interdonato, 2019) capacity bound             
        [SE_MR,~] = capacity_bound_interdonato(K,estimated_effective_channel_MR,interUserGains_MR);        
        [SE_LP_MMSE,~] = capacity_bound_interdonato(K,estimated_effective_channel_LP_MMSE,interUserGains_LP_MMSE);        
        [SE_P_MMSE,~] = capacity_bound_interdonato(K,estimated_effective_channel_P_MMSE,interUserGains_P_MMSE);
        [SE_P_RZF,~] = capacity_bound_interdonato(K,estimated_effective_channel_P_RZF,interUserGains_P_RZF);
              
        [SE_MR_pCSI,~] = capacity_bound_interdonato(K,perfectEstimatedEffectiveChannel_MR,interUserGains_MR);
        [SE_LP_MMSE_pCSI,~] = capacity_bound_interdonato(K,perfectEstimatedEffectiveChannel_LP_MMSE,interUserGains_LP_MMSE);        
        [SE_P_MMSE_pCSI,~] = capacity_bound_interdonato(K,perfectEstimatedEffectiveChannel_P_MMSE,interUserGains_P_MMSE);
        [SE_P_RZF_pCSI,~] = capacity_bound_interdonato(K,perfectEstimatedEffectiveChannel_P_RZF,interUserGains_P_RZF);
    
    case 'UnF'
        %Use and then forget capacity bound
        %Compute SE with MR precoding, assuming blind estimation at the UE
        [SE_MR,~] = use_and_then_forget_bound(K,estimated_effective_channel_MR,interUserGains_MR);        
        [SE_LP_MMSE,~] = use_and_then_forget_bound(K,estimated_effective_channel_LP_MMSE,interUserGains_LP_MMSE);        
        [SE_P_MMSE,~] = use_and_then_forget_bound(K,estimated_effective_channel_P_MMSE,interUserGains_P_MMSE);
        [SE_P_RZF,~] = use_and_then_forget_bound(K,estimated_effective_channel_P_RZF,interUserGains_P_RZF);

         %Compute SE with MR precoding, assuming blind estimation at the UE
        [SE_MR_pCSI,~] = use_and_then_forget_bound(K,perfectEstimatedEffectiveChannel_MR,interUserGains_MR);        
        [SE_LP_MMSE_pCSI,~] = use_and_then_forget_bound(K,perfectEstimatedEffectiveChannel_LP_MMSE,interUserGains_LP_MMSE);                
        [SE_P_MMSE_pCSI,~] = use_and_then_forget_bound(K,perfectEstimatedEffectiveChannel_P_MMSE,interUserGains_P_MMSE);
        [SE_P_RZF_pCSI,~] = use_and_then_forget_bound(K,perfectEstimatedEffectiveChannel_P_RZF,interUserGains_P_RZF);

end
%% SE with prelog 
SE_MR = prelogFactor*SE_MR;
SE_LP_MMSE = prelogFactor*SE_LP_MMSE;
SE_P_MMSE = prelogFactor*SE_P_MMSE;
SE_P_RZF = prelogFactor*SE_P_RZF;

SE_MR_BE = prelogFactor*SE_MR_BE;
SE_LP_MMSE_BE = prelogFactor*SE_LP_MMSE_BE;
SE_P_MMSE_BE = prelogFactor*SE_P_MMSE_BE;
SE_P_RZF_BE = prelogFactor*SE_P_RZF_BE;

SE_MR_BE_InfTau = prelogFactor*SE_MR_BE_InfTau;
SE_LP_MMSE_BE_InfTau = prelogFactor*SE_LP_MMSE_BE_InfTau;
SE_P_MMSE_BE_InfTau = prelogFactor*SE_P_MMSE_BE_InfTau;
SE_P_RZF_BE_InfTau = prelogFactor*SE_P_RZF_BE_InfTau;

SE_MR_pCSI = prelogFactor*SE_MR_pCSI;
SE_LP_MMSE_pCSI = prelogFactor*SE_LP_MMSE_pCSI;
SE_P_MMSE_pCSI = prelogFactor*SE_P_MMSE_pCSI;
SE_P_RZF_pCSI = prelogFactor*SE_P_RZF_pCSI;
 
SE_MR_DLPE = prelogFactor_DLPE*SE_MR_DLPE;
SE_LP_MMSE_DLPE = prelogFactor_DLPE*SE_LP_MMSE_DLPE;
SE_P_MMSE_DLPE = prelogFactor_DLPE*SE_P_MMSE_DLPE;
SE_P_RZF_DLPE = prelogFactor_DLPE*SE_P_RZF_DLPE;

end

%% Capacity bound functions

% Calculates a lower bound capacity based on the equation (30) from the article 'Downlink Training in Cell-Free Massive MIMO:
% A Blessing in Disguise.'
%
% SYNTAX:
%   spectralEfficiency = capacity_bound_interdonato(numberOfUEs,estimatedEffectiveChannel,effectiveChannel)
%
% INPUT:
%   numberOfUEs: number of users
%   estimatedEffectiveChannel: Matrix of dimention numberOfUEs-by-nbrOfRealizations containig the 
%      estimated Effective Channel of UE k (α^^kk).       
%   effectiveChannel: Matrix of dimention numberOfUEs-by-numberOfUEs-by-nbrOfRealizations  
%      containig all effective channels from UE k to UE k' (αkk').
%
% OUTPUT:
%   spectralEfficiency:
%     Vector of dimension numberOfUEs-by-1, containing the spectral efficiency in bits/sec/Hz of each UE.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Created by:    Daynara Dias Soyza
% Contact:       daynara@ufpa.br
% Last modified: 08/19/2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spectralEfficiency,SINR] = capacity_bound_interdonato(numberOfUEs,estimatedEffectiveChannel,effectiveChannel_all)

%Prepare to save simulation results
SINR = zeros(numberOfUEs,1);
spectralEfficiency = zeros(numberOfUEs,1);

effectiveChannel = zeros(size(estimatedEffectiveChannel));
for aux1 = 1:numberOfUEs
effectiveChannel(aux1,:) = effectiveChannel_all(aux1,aux1,:);
end

for i = 1:numberOfUEs
   
    SINR(i) = ( mean ( ( abs( estimatedEffectiveChannel(i,:) ) ) .^2) ) ./ ...
        (mean( ( abs( ((effectiveChannel(i,:))) - estimatedEffectiveChannel(i,:) ) ).^2 )  + sum(mean(  abs ( effectiveChannel_all(i,[1:i-1 i+1:end],:) ) .^2 ,3),2) + 1 );

    spectralEfficiency(i) = log2 (1 + SINR(i));
end

end

% Calculates a lower bound capacity based on the equation (44) from the 
% article 'Ngo, H. Q., G. Larsson, E., “No downlink pilots are needed in 
% Massive MIMO�?. IEEE Transactions on Communications, 2017.'
%
% SYNTAX:
%   spectralEfficiency = capacity_bound_ngo(numberOfUEs,totalDLpower,powerCoefficients,nbrOfRealizations,x,y,z)
%
% INPUT:
%   numberOfUEs: number of users
%   totalDLpower: total downlink power
%   powerCoefficients: Vector of dimension 1-by-numberOfUEs containing the 
%                      power control coefficients.
%   nbrOfRealizations: number of channel realizations
%   x: Matrix of dimention numberOfUEs-by-nbrOfRealizations containig the 
%      estimated Effective Channel of UE k (α^^kk(1)).       
%   y: Matrix of dimention numberOfUEs-by-nbrOfRealizations containig the
%      effective channel of UE k (αkk).
%   z: Matrix of dimention numberOfUEs-by-numberOfUEs-by-nbrOfRealizations  
%      containig all effective channels from 
%      UE k to UE k' (|αkk'|).
%
% OUTPUT:
%   spectralEfficiency:
%     Vector of dimension numberOfUEs-by-1, containing the spectral 
%         efficiency in bits/sec/Hz of each UE.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Created by:    Daynara Dias Soyza
% Contact:       daynara@ufpa.br
% Last modified: 08/19/2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [spectralEfficiency] = capacity_bound_ngo(numberOfUEs,nbrOfRealizations,x,z)

y = zeros(numberOfUEs,nbrOfRealizations);
for aux = 1:numberOfUEs
y(numberOfUEs,:) = z (numberOfUEs,numberOfUEs,:);
end

%Size of the numberOfSamples by numberOfSamples grid over which the density 
%is computed, numberOfSamples has to be a power of 2
numberOfSamples = 2^7;%2^10;%nbrOfRealizations;%
x = abs(x);
y = abs(y);
z = (abs(z)).^2;
epslon = 1e-20;
%Prepare to save simulation results
espectationOfAllEffectiveChannelsSquaredGivenTheirEstimation = zeros(numberOfSamples,numberOfUEs);
SINR = zeros(numberOfSamples,numberOfUEs);
spectralEfficiency = zeros(numberOfUEs,1);


for i = 1:numberOfUEs
    
    %Probability density function (PDF)
    xRange = 0;% max(x(i,:)) - min(x(i,:));
    xMin = min(x(i,:)) - xRange/10;
    xMax = max(x(i,:)) + xRange/10;
    [~,xPdf,X,~] = kde(x(i,:), numberOfSamples, xMin, xMax);
    deltaX = X(1,2) -  X(1,1);%(X(1,2) -  X(1,1)) * ones(numberOfSamples,1);
    xPdf(xPdf<epslon) = epslon;
    
    %Joint probability density function (PDF)   
    yxData = [transpose(abs(y(i,:))) transpose(x(i,:))];
    yxRange = 0;%max(yxData,[],1) - min(yxData,[],1);
    yxMax = max(yxData,[],1) + yxRange/10;
    yxMin = min(yxData,[],1) - yxRange/10;
    [~, yxPdf, Y, ~] = kde2d(yxData, numberOfSamples, yxMin, yxMax);
    deltaY = Y(1,2) -  Y(1,1);%(Y(:,2) -  Y(:,1)) * ones(1,numberOfSamples);
    
    %Eq. (47)
    espectationOfEffectiveChannelGivenEstimation = sum( Y .* deltaY .* yxPdf , 2 )./xPdf;

    for j =1:numberOfUEs
        
        %Joint probability density function (PDF)
        zxData = [squeeze(z(i,j,:)) transpose(x(i,:))];
        zxRange = 0;%max(zxData,[],1) - min(zxData,[],1);
        zxMax = max(zxData,[],1) + zxRange/10;
        zxMin = min(zxData,[],1) - zxRange/10;
        [~, zxPdf, Z, ~]=kde2d(zxData, numberOfSamples, zxMin, zxMax); 
        deltaZ = Z(1,2) -  Z(1,1);%(Z(:,2) -  Z(:,1)) * ones(1,numberOfSamples);
        
        %Eq.(48)
        espectationOfAllEffectiveChannelsSquaredGivenTheirEstimation(:,j) =  (sum(Z .* deltaZ .* zxPdf , 2 ))./xPdf;
       
    end    
    %Eq.(46)
    SINR(:,i) = ( (  ( abs ( espectationOfEffectiveChannelGivenEstimation ) ) .^2 ) ./...
        (1 +  sum( espectationOfAllEffectiveChannelsSquaredGivenTheirEstimation ,2)  -  ( abs ( espectationOfEffectiveChannelGivenEstimation ) ) .^2 ) );
    
    SINR(SINR<0) = 0;
    
    spectralEfficiency(i) = sum( xPdf .* deltaX .* ( log2( 1 + SINR(:,i) ) ) );

end

end


% Calculates a lower bound capacity based on the equation (X) from the article 'XX.'
%
% SYNTAX:
%   spectralEfficiency = capacity_bound_ngo(numberOfUEs,totalDLpower,powerCoefficients,estimatedEffectiveChannels,effectiveChannels)
%
% INPUT:
%   numberOfUEs: number of users
%   totalDLpower: total downlink power
%   powerCoefficients: Vector of dimension 1-by-numberOfUEs containing the
%                      power control coefficients.
%   estimatedEffectiveChannels: Matrix of dimention numberOfUEs-by-by-nbrOfRealizations  containig the
%      estimated Effective Channel of UE k (α^^kk).
%   effectiveChannels: Matrix of dimention numberOfUEs-by-by-nbrOfRealizations
%      containig effective channels from UE k to UE k (αkk).
%
% OUTPUT:
%   spectralEfficiency:
%     Vector of dimension numberOfUEs-by-1, containing the spectral efficiency in bits/sec/Hz of each UE.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Created by:    Daynara Dias Soyza
% Contact:       daynara@ufpa.br
% Last modified: 08/19/2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spectralEfficiency,SINR] = capacity_bound_ngo_2013(numberOfUEs,nbrOfRealizations,estimatedEffectiveChannels,effectiveChannels)

%Prepare to save simulation results
SINR = zeros(numberOfUEs,nbrOfRealizations);
spectralEfficiency = zeros(numberOfUEs,1);


effectiveChannels2 = zeros(numberOfUEs,nbrOfRealizations);
for aux =1:numberOfUEs
effectiveChannels2(aux,:) = effectiveChannels(aux,aux,:);
end
for i = 1:numberOfUEs
       
    SINR(i,:) =  ( ( abs( estimatedEffectiveChannels(i,:) ) ) .^2)  ./ ...
        ( sum(mean( abs( effectiveChannels2 - estimatedEffectiveChannels ).^2 , 2))  + sum( (abs (estimatedEffectiveChannels([1:i-1 i+1:end],:))).^2) + 1 );
    
    spectralEfficiency(i) = mean (log2 (1 +  SINR(i,:) ) ,2);
    
end
SINR = mean(SINR,2);
end


% Calculates a lower bound capacity based on the equation (12) from the article 'Massive MU-MIMO Downlink TDD Systems with
% Linear Precoding and Downlink Pilots'
%
% SYNTAX:
%   spectralEfficiency = capacity_bound_interdonato(numberOfUEs,nbrOfRealizations,estimatedEffectiveChannels,effectiveChannels)
%
% INPUT:
%   numberOfUEs: number of users
%   estimatedEffectiveChannels: Matrix of dimention numberOfUEs-by-nbrOfRealizations containig the 
%      estimated Effective Channel of UE k (α^^kk).       
%   effectiveChannels: Matrix of dimention numberOfUEs-by-numberOfUEs-by-nbrOfRealizations  
%      containig all effective channels from UE k to UE k' (αkk').
%
% OUTPUT:
%   spectralEfficiency:
%     Vector of dimension numberOfUEs-by-1, containing the spectral efficiency in bits/sec/Hz of each UE.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Created by:    Daynara Dias Soyza
% Contact:       daynara@ufpa.br
% Last modified: 08/19/2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spectralEfficiency = capacity_lower_bound_ngo(numberOfUEs,nbrOfRealizations,estimatedEffectiveChannels,effectiveChannels)

%Prepare to save simulation results
SINR = zeros(numberOfUEs,nbrOfRealizations);
spectralEfficiency = zeros(numberOfUEs,1);

for i = 1:numberOfUEs

    SINR(i,:) = (   ( ( abs( estimatedEffectiveChannels(i,:) ) ) .^2) ) ./ ...
        ( mean( ( abs( transpose(squeeze(effectiveChannels(i,i,:))) - estimatedEffectiveChannels(i,:) ) ).^2 )  +   transpose(squeeze(sum( ( abs ( effectiveChannels(i,[1:i-1 i+1:end],:) ) ).^2 ,2))) + 1 );

    spectralEfficiency(i) = mean(log2 (1 + SINR(i,:)));
    
end

end

% Calculates a lower bound capacity based on the equation (59) from the article 'Ngo, H. Q., G. Larsson, E., “No downlink pilots are needed in Massive
% MIMO�?. IEEE Transactions on Communications, 2017.'
%
% SYNTAX:
%   spectralEfficiency = use_and_then_forget_bound_ngo(numberOfUEs,estimatedEffectiveChannel,effectiveChannel)
%
% INPUT:
%   numberOfUEs: number of users
%   estimatedEffectiveChannel: Matrix of dimention numberOfUEs-by-nbrOfRealizations containig the 
%      estimated Effective Channel of UE k (α^^kk).       
%   effectiveChannel: Matrix of dimention numberOfUEs-by-numberOfUEs-by-nbrOfRealizations  
%      containig all effective channels from UE k to UE k' (αkk').
%
% OUTPUT:
%   spectralEfficiency:
%     Vector of dimension numberOfUEs-by-1, containing the spectral efficiency in bits/sec/Hz of each UE.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Created by:    Daynara Dias Soyza
% Contact:       daynara@ufpa.br
% Last modified: 08/19/2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [spectralEfficiency,SINR] = use_and_then_forget_bound(numberOfUEs,estimatedEffectiveChannel,effectiveChannel_all)

%Prepare to save simulation results
SINR = zeros(numberOfUEs,1);
spectralEfficiency = zeros(numberOfUEs,1);
temp = zeros(numberOfUEs,1);

[~,nbrOfRealizations] = size(estimatedEffectiveChannel);
effectiveChannel = zeros(size(estimatedEffectiveChannel));
for aux1 = 1:numberOfUEs
effectiveChannel(aux1,:) = effectiveChannel_all(aux1,aux1,:);
end

for i = 1:numberOfUEs

    for j =1:numberOfUEs
        
        if j ~= i
            aux2 = zeros(1,nbrOfRealizations);
            aux2(1,:) = effectiveChannel_all(i,j,:);
            temp(j) =  mean( (abs ( aux2 ./ estimatedEffectiveChannel(i,:) ) ).^2 ) ;
            
        end
        
    end

    SINR(i) = ( ( abs( mean( ((effectiveChannel(i,:))) ./ estimatedEffectiveChannel(i,: ) ) ) ) .^2 ) ./ ...
        ( var( ((effectiveChannel(i,:))) ./ estimatedEffectiveChannel(i,: ) )  +  sum(temp) ...
        + mean( 1./( abs ( estimatedEffectiveChannel(i,:) ) ).^2  ) );

    spectralEfficiency(i) = log2 (1 + SINR(i));
    
    temp = zeros(numberOfUEs,1);
 
end


end

