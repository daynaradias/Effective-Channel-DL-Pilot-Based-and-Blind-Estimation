% tic
%Empty workspace and close figures
% close all;clear;
%% Define simulation setup

%Number of APs
% L = 100;

%Number of antennas per AP
% N = 1;

%Number of UEs in the network
% K = 20;

%Length of coherence block
% tau_c = 200;

%Length of UL and DL pilot sequences
% tau_up = 10;
% tau_dp = 10;

% Options for capacityBound: 'Interdonato2019', 'UnF','Ngo2013','Ngo2013_2','Ngo2017','Hardening'
% capacityBound      = 'UnF';

%Total downlink transmit power per AP (mW)
% rho_tot = 200;

% figures and .mat name
% name = [capacityBound '_L' num2str(L) '_N' num2str(N) '_K' num2str(K) '_rhoD' num2str(rho_tot)...
%     '_tauC' num2str(tau_c) '_tauUP' num2str(tau_up) '_tauDP' num2str(tau_dp)];

function mainFunction(L,N,K,tau_c,tau_up,tau_dp,capacityBound, rho_tot, pilotAssignmentMethod, APselection , alpha, name )

%Number of Monte-Carlo setups
nbrOfSetups = 100;

%Number of channel realizations per setup
nbrOfRealizations = 500;

%Angular standard deviation in the local scattering model (in radians)
ASD = deg2rad(15);


%% Propagation parameters

%Total uplink transmit power per UE (mW)
p = 100;

%% Prepare to save simulation results

%Prepare to save simulation results (average SEs)
SE_P_MMSE_tot   = zeros(K,nbrOfSetups);
SE_P_RZF_tot    = zeros(K,nbrOfSetups);
SE_LP_MMSE_tot  = zeros(K,nbrOfSetups);
SE_MR_tot       = zeros(K,nbrOfSetups);

Gen_SE_P_MMSE_tot  = zeros(K,nbrOfSetups);
Gen_SE_P_RZF_tot   = zeros(K,nbrOfSetups);
Gen_SE_LP_MMSE_tot = zeros(K,nbrOfSetups);
Gen_SE_MR_tot      = zeros(K,nbrOfSetups);

%Prepare to store the hardening
CH_MR      = zeros(K,nbrOfSetups);
CH_LP_MMSE = zeros(K,nbrOfSetups);
CH_P_MMSE  = zeros(K,nbrOfSetups);
CH_P_RZF   = zeros(K,nbrOfSetups);

%Prepare to store the favorable propagation
FP_MR      = zeros(K*K - K,nbrOfSetups);
FP_LP_MMSE = zeros(K*K - K,nbrOfSetups);
FP_P_MMSE  = zeros(K*K - K,nbrOfSetups);
FP_P_RZF   = zeros(K*K - K,nbrOfSetups);

%Prepare to SE and NMSE of different estimation methods
SE_MR_sCSI_tot       = zeros(K,nbrOfSetups);
SE_LP_MMSE_sCSI_tot  = zeros(K,nbrOfSetups);
SE_P_MMSE_sCSI_tot   = zeros(K,nbrOfSetups);
SE_P_RZF_sCSI_tot   = zeros(K,nbrOfSetups);

SE_MR_pCSI_tot       = zeros(K,nbrOfSetups);
SE_LP_MMSE_pCSI_tot  = zeros(K,nbrOfSetups);
SE_P_MMSE_pCSI_tot   = zeros(K,nbrOfSetups);
SE_P_RZF_pCSI_tot   = zeros(K,nbrOfSetups);

SE_MR_BE_tot         = zeros(K,nbrOfSetups);
SE_LP_MMSE_BE_tot    = zeros(K,nbrOfSetups);
SE_P_MMSE_BE_tot     = zeros(K,nbrOfSetups);
SE_P_RZF_BE_tot     = zeros(K,nbrOfSetups);

SE_MR_BE_InfTau_tot       = zeros(K,nbrOfSetups);
SE_LP_MMSE_BE_InfTau_tot  = zeros(K,nbrOfSetups);
SE_P_MMSE_BE_InfTau_tot   = zeros(K,nbrOfSetups);
SE_P_RZF_BE_InfTau_tot   = zeros(K,nbrOfSetups);

SE_MR_DLPE_tot       = zeros(K,nbrOfSetups);
SE_LP_MMSE_DLPE_tot  = zeros(K,nbrOfSetups);
SE_P_MMSE_DLPE_tot   = zeros(K,nbrOfSetups);
SE_P_RZF_DLPE_tot   = zeros(K,nbrOfSetups);

NMSE_MR_tot          = zeros(K,nbrOfSetups);
NMSE_LP_MMSE_tot     = zeros(K,nbrOfSetups);
NMSE_P_MMSE_tot      = zeros(K,nbrOfSetups);
NMSE_P_RZF_tot      = zeros(K,nbrOfSetups);

NMSE_MR_BE_tot       = zeros(K,nbrOfSetups);
NMSE_LP_MMSE_BE_tot  = zeros(K,nbrOfSetups);
NMSE_P_MMSE_BE_tot   = zeros(K,nbrOfSetups);
NMSE_P_RZF_BE_tot   = zeros(K,nbrOfSetups);

NMSE_MR_BE_InfTau_tot       = zeros(K,nbrOfSetups);
NMSE_LP_MMSE_BE_InfTau_tot  = zeros(K,nbrOfSetups);
NMSE_P_MMSE_BE_InfTau_tot   = zeros(K,nbrOfSetups);
NMSE_P_RZF_BE_InfTau_tot   = zeros(K,nbrOfSetups);

NMSE_MR_DLPE_tot     = zeros(K,nbrOfSetups);
NMSE_LP_MMSE_DLPE_tot= zeros(K,nbrOfSetups);
NMSE_P_MMSE_DLPE_tot = zeros(K,nbrOfSetups);
NMSE_P_RZF_DLPE_tot = zeros(K,nbrOfSetups);

Ee_MR_tot = zeros(1,nbrOfSetups);
Ee_LP_MMSE_tot = zeros(1,nbrOfSetups);
Ee_P_MMSE_tot = zeros(1,nbrOfSetups);
Ee_P_RZF_tot = zeros(1,nbrOfSetups);

Gen_Ee_MR_tot = zeros(1,nbrOfSetups);
Gen_Ee_LP_MMSE_tot = zeros(1,nbrOfSetups);
Gen_Ee_P_MMSE_tot = zeros(1,nbrOfSetups);
Gen_Ee_P_RZF_tot = zeros(1,nbrOfSetups);

Ee_MR_BE_tot = zeros(1,nbrOfSetups);
Ee_LP_MMSE_BE_tot = zeros(1,nbrOfSetups);
Ee_P_MMSE_BE_tot = zeros(1,nbrOfSetups);
Ee_P_RZF_BE_tot = zeros(1,nbrOfSetups);

Ee_MR_DLPE_tot = zeros(1,nbrOfSetups);
Ee_LP_MMSE_DLPE_tot = zeros(1,nbrOfSetups);
Ee_P_MMSE_DLPE_tot = zeros(1,nbrOfSetups);
Ee_P_RZF_DLPE_tot = zeros(1,nbrOfSetups);

%Prepare to store DCC matrices
D_tot = zeros(L,K,nbrOfSetups);
% APselection =
% 'BookEmil';%'BookEmil';'CanonicalCF';'ClusterControlCmax'
% pilotAssignmentMethod = 'joint_minContamination';%'joint_random';'joint_minContamination';'nonJoint_random';'nonJoint_minContamination';

%% Go through all setups
parfor n = 1:nbrOfSetups
    
    %Display simulation progress
    disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]);
       
    %% Generate one setup with UEs and APs at random locations
    [gainOverNoisedB,R,pilotIndex,D,D_small] = generateSetup(L,K,N,tau_up,1,0,APselection,alpha,pilotAssignmentMethod,ASD,ASD);
    %Save the DCC matrix
    D_tot(:,:,n) = D;
    
    %% Generate channel realizations with estimates and estimation
    %error correlation matrices
    [Hhat,H,B,C] = functionChannelEstimates(R,nbrOfRealizations,L,K,N,tau_up,pilotIndex,p);
    
    %% Downlink
    % Compute the power allocation in (6.36) for distributed precoding
    rho_dist = zeros(L,K);
    
    gainOverNoise = db2pow(gainOverNoisedB);
    
    for l = 1:L
        
        %Extract which UEs are served by AP l
        servedUEs = find(D(l,:)==1);
        
        %Compute denominator in (6.36)
        normalizationAPl = sum(sqrt(gainOverNoise(l,servedUEs)));
        
        for ind = 1:length(servedUEs)
            
            rho_dist(l,servedUEs(ind)) = rho_tot*sqrt(gainOverNoise(l,servedUEs(ind)))/normalizationAPl;
            
        end
        
    end
    
    
    %Compute DL SE using precoders and results in Section 6 for centralized
    %and distributed uplink operations for DCC
    
    %Compute SEs for DCC case using hardening bound 
    % statistical CSI and perfect CSI
    [SE_MMSE, SE_P_MMSE, SE_P_RZF,  ...
        SE_L_MMSE, SE_LP_MMSE, SE_MR, ...
        Gen_SE_P_MMSE, Gen_SE_P_RZF, Gen_SE_LP_MMSE, Gen_SE_MR, ...
        interUserGains_MMSE, interUserGains_P_MMSE, interUserGains_P_RZF,...
        interUserGains_L_MMSE, interUserGains_LP_MMSE, interUserGains_MR,cont_MR,...
        w_MR, w_LP_MMSE, w_P_MMSE, w_PRZF] ...
        = functionComputeInterUserGains(Hhat,H,D,B,C,tau_c,tau_up,nbrOfRealizations,N,K,L,p,R,pilotIndex,rho_dist,gainOverNoisedB,rho_tot);
    
    % channel hardening and favorable propagation
    CHandFP_MR = real((var(interUserGains_MR,1,3))./(diag(mean(interUserGains_MR,3)).^2).*ones(K,K));
    CHandFP_LP_MMSE = real((var(interUserGains_LP_MMSE,1,3))./(diag(mean(interUserGains_LP_MMSE,3)).^2).*ones(K,K));
    CHandFP_P_MMSE = real((var(interUserGains_P_MMSE,1,3))./(diag(mean(interUserGains_P_MMSE,3)).^2).*ones(K,K));
    CHandFP_P_RZF = real((var(interUserGains_P_RZF,1,3))./(diag(mean(interUserGains_P_RZF,3)).^2).*ones(K,K));
    
    % Blind estimation and DL pilot-based estimation
    [SE_MR_sCSI, SE_LP_MMSE_sCSI, SE_P_MMSE_sCSI, SE_P_RZF_sCSI, ...
        NMSE_MR, NMSE_LP_MMSE, NMSE_P_MMSE, NMSE_P_RZF,...
        SE_MR_pCSI, SE_LP_MMSE_pCSI, SE_P_MMSE_pCSI, SE_P_RZF_pCSI,...
        SE_MR_BE,SE_LP_MMSE_BE,SE_P_MMSE_BE,SE_P_RZF_BE,...
        NMSE_MR_BE,NMSE_LP_MMSE_BE,NMSE_P_MMSE_BE,NMSE_P_RZF_BE,...
        SE_MR_BE_InfTau,SE_LP_MMSE_BE_InfTau,SE_P_MMSE_BE_InfTau,SE_P_RZF_BE_InfTau,...
        NMSE_MR_BE_InfTau,NMSE_LP_MMSE_BE_InfTau,NMSE_P_MMSE_BE_InfTau,NMSE_P_RZF_BE_InfTau,...
        SE_MR_DLPE,SE_LP_MMSE_DLPE,SE_P_MMSE_DLPE,SE_P_RZF_DLPE,...
        NMSE_MR_DLPE,NMSE_LP_MMSE_DLPE,NMSE_P_MMSE_DLPE,NMSE_P_RZF_DLPE]...
            = functionComputeSE_downlink_sCSI_pCSI_BE_DLPE...
            (capacityBound,tau_c,tau_up,pilotIndex,tau_dp,nbrOfRealizations,K,pilotAssignmentMethod,...
            CHandFP_MR,CHandFP_LP_MMSE,CHandFP_P_MMSE,CHandFP_P_RZF,...
            interUserGains_MR,interUserGains_LP_MMSE,interUserGains_P_MMSE,interUserGains_P_RZF);

        % energy efficiency
        [Ee_MR , Ee_LP_MMSE, Ee_P_MMSE, Ee_P_RZF,...
            Gen_Ee_MR , Gen_Ee_LP_MMSE, Gen_Ee_P_MMSE, Gen_Ee_P_RZF,...
            Ee_MR_BE , Ee_LP_MMSE_BE, Ee_P_MMSE_BE, Ee_P_RZF_BE,...
            Ee_MR_DLPE , Ee_LP_MMSE_DLPE, Ee_P_MMSE_DLPE, Ee_P_RZF_DLPE] = ...
            EnergyEfficiency (w_MR,w_LP_MMSE,w_P_MMSE,w_PRZF,...
            SE_MR , SE_LP_MMSE, SE_P_MMSE, SE_P_RZF,...
            Gen_SE_MR , Gen_SE_LP_MMSE, Gen_SE_P_MMSE, Gen_SE_P_RZF,...
            SE_MR_BE, SE_LP_MMSE_BE,SE_P_MMSE_BE,SE_P_RZF_BE,...
            SE_MR_DLPE, SE_LP_MMSE_DLPE, SE_P_MMSE_DLPE, SE_P_RZF_DLPE,...
            L,N,D,nbrOfRealizations);
        
    %Save results
    SE_P_MMSE_tot(:,n)      = SE_P_MMSE;
    SE_P_RZF_tot(:,n)       = SE_P_RZF;
    SE_LP_MMSE_tot(:,n)     = SE_LP_MMSE;
    SE_MR_tot(:,n)          = SE_MR;
      
    Gen_SE_P_MMSE_tot(:,n)  = Gen_SE_P_MMSE;
    Gen_SE_P_RZF_tot(:,n)   = Gen_SE_P_RZF;
    Gen_SE_LP_MMSE_tot(:,n) = Gen_SE_LP_MMSE;
    Gen_SE_MR_tot(:,n)      = Gen_SE_MR;
    
    
    CH_MR(:,n)              = diag(CHandFP_MR);
    CH_LP_MMSE(:,n)         = diag(CHandFP_LP_MMSE);
    CH_P_MMSE(:,n)          = diag(CHandFP_P_MMSE);
    CH_P_RZF(:,n)           = diag(CHandFP_P_RZF);
    
    idx = eye(K,K);
    FP_MR(:,n)              = CHandFP_MR(~idx);
    FP_LP_MMSE(:,n)         = CHandFP_LP_MMSE(~idx);
    FP_P_MMSE(:,n)          = CHandFP_P_MMSE(~idx);
    FP_P_RZF(:,n)           = CHandFP_P_RZF(~idx);
    
    SE_MR_sCSI_tot(:,n)        = SE_MR_sCSI;
    SE_LP_MMSE_sCSI_tot(:,n)   = SE_LP_MMSE_sCSI;
    SE_P_MMSE_sCSI_tot(:,n)    = SE_P_MMSE_sCSI;
    SE_P_RZF_sCSI_tot(:,n)    = SE_P_RZF_sCSI;
    
    SE_MR_pCSI_tot(:,n)        = SE_MR_pCSI;
    SE_LP_MMSE_pCSI_tot(:,n)   = SE_LP_MMSE_pCSI;
    SE_P_MMSE_pCSI_tot(:,n)    = SE_P_MMSE_pCSI;
    SE_P_RZF_pCSI_tot(:,n)    = SE_P_RZF_pCSI;
       
    SE_MR_BE_tot(:,n)          = SE_MR_BE;
    SE_LP_MMSE_BE_tot(:,n)     = SE_LP_MMSE_BE;
    SE_P_MMSE_BE_tot(:,n)      = SE_P_MMSE_BE;
    SE_P_RZF_BE_tot(:,n)      = SE_P_RZF_BE;
    
    SE_MR_BE_InfTau_tot(:,n)       = SE_MR_BE_InfTau;
    SE_LP_MMSE_BE_InfTau_tot(:,n)  = SE_LP_MMSE_BE_InfTau;
    SE_P_MMSE_BE_InfTau_tot(:,n)   = SE_P_MMSE_BE_InfTau;
    SE_P_RZF_BE_InfTau_tot(:,n)   = SE_P_RZF_BE_InfTau;
    
    SE_MR_DLPE_tot(:,n)        = SE_MR_DLPE;
    SE_LP_MMSE_DLPE_tot(:,n)   = SE_LP_MMSE_DLPE;
    SE_P_MMSE_DLPE_tot(:,n)    = SE_P_MMSE_DLPE;
    SE_P_RZF_DLPE_tot(:,n)    = SE_P_RZF_DLPE;
    
    NMSE_MR_tot (:,n)          = NMSE_MR;
    NMSE_LP_MMSE_tot(:,n)      = NMSE_LP_MMSE;
    NMSE_P_MMSE_tot(:,n)       = NMSE_P_MMSE;
    NMSE_P_RZF_tot(:,n)       = NMSE_P_RZF;
    
    NMSE_MR_BE_tot (:,n)       = NMSE_MR_BE;
    NMSE_LP_MMSE_BE_tot(:,n)   = NMSE_LP_MMSE_BE;
    NMSE_P_MMSE_BE_tot(:,n)    = NMSE_P_MMSE_BE;
    NMSE_P_RZF_BE_tot(:,n)    = NMSE_P_RZF_BE;
    
    NMSE_MR_BE_InfTau_tot(:,n)       = NMSE_MR_BE_InfTau;
    NMSE_LP_MMSE_BE_InfTau_tot(:,n)  = NMSE_LP_MMSE_BE_InfTau;
    NMSE_P_MMSE_BE_InfTau_tot(:,n)   = NMSE_P_MMSE_BE_InfTau;
    NMSE_P_RZF_BE_InfTau_tot(:,n)   = NMSE_P_RZF_BE_InfTau;
    
    NMSE_MR_DLPE_tot (:,n)     = NMSE_MR_DLPE;
    NMSE_LP_MMSE_DLPE_tot (:,n)= NMSE_LP_MMSE_DLPE;
    NMSE_P_MMSE_DLPE_tot (:,n) = NMSE_P_MMSE_DLPE;
    NMSE_P_RZF_DLPE_tot (:,n) = NMSE_P_RZF_DLPE;
    
    
    Ee_MR_tot(:,n) = Ee_MR;
    Ee_LP_MMSE_tot(:,n) = Ee_LP_MMSE;
    Ee_P_MMSE_tot(:,n) = Ee_P_MMSE;
    Ee_P_RZF_tot(:,n) = Ee_P_RZF;
    
    Gen_Ee_MR_tot(:,n) = Gen_Ee_MR;
    Gen_Ee_LP_MMSE_tot(:,n) = Gen_Ee_LP_MMSE;
    Gen_Ee_P_MMSE_tot(:,n) = Gen_Ee_P_MMSE;
    Gen_Ee_P_RZF_tot(:,n) = Gen_Ee_P_RZF;
    
    Ee_MR_BE_tot(:,n) = Ee_MR_BE;
    Ee_LP_MMSE_BE_tot(:,n) = Ee_LP_MMSE_BE;
    Ee_P_MMSE_BE_tot(:,n) = Ee_P_MMSE_BE;
    Ee_P_RZF_BE_tot(:,n) = Ee_P_RZF_BE;
    
    Ee_MR_DLPE_tot(:,n) = Ee_MR_DLPE;
    Ee_LP_MMSE_DLPE_tot(:,n) = Ee_LP_MMSE_DLPE;
    Ee_P_MMSE_DLPE_tot(:,n) = Ee_P_MMSE_DLPE;
    Ee_P_RZF_DLPE_tot(:,n) = Ee_P_RZF_DLPE;
        
    %Remove large matrices at the end of analyzing this setup
%     clear Hhat H B C R;
    Hhat = []; H = []; B = []; C = []; R = [];
    
    
    
end
% % toc
% % Plot spectral efficiency results
% figure;
% hold on
% box on;
% plot(sort((SE_MR_tot(:))),linspace(0,1,K*nbrOfSetups),'k','LineWidth',2);
% plot(sort((SE_LP_MMSE_tot(:))),linspace(0,1,K*nbrOfSetups),'b','LineWidth',2);
% plot(sort((SE_P_MMSE_tot(:))),linspace(0,1,K*nbrOfSetups),'r','LineWidth',2);
% plot(sort((SE_P_RZF_tot(:))),linspace(0,1,K*nbrOfSetups),'g','LineWidth',2);
% xlabel('SE(bit/s/Hz)');
% ylabel('CDF','Interpreter','Latex');
% legend({'MR','LP-MMSE','P-MMSE','P-RZF'},'Interpreter','tex','Location','SouthEast');
% saveas(gcf,[name 'SE'],'fig')
% 
% figure;
% hold on
% box on;
% plot(sort((Gen_SE_MR_tot(:))),linspace(0,1,K*nbrOfSetups),'k','LineWidth',2);
% plot(sort((Gen_SE_LP_MMSE_tot(:))),linspace(0,1,K*nbrOfSetups),'b','LineWidth',2);
% plot(sort((Gen_SE_P_MMSE_tot(:))),linspace(0,1,K*nbrOfSetups),'r','LineWidth',2);
% plot(sort((Gen_SE_P_RZF_tot(:))),linspace(0,1,K*nbrOfSetups),'g','LineWidth',2);
% xlabel('SE(bit/s/Hz)');
% ylabel('CDF','Interpreter','Latex');
% legend({'MR(pCSI)','LP-MMSE(pCSI)','P-MMSE(pCSI)','P-RZF(pCSI)'},'Interpreter','tex','Location','SouthEast');
% saveas(gcf,[name 'Gen_SE'],'fig')
% 
% % Plot Channel hardening results
% figure;
% hold on
% box on;
% plot(sort((CH_MR(:))),linspace(0,1,K*nbrOfSetups),'k','LineWidth',2);
% plot(sort((CH_LP_MMSE(:))),linspace(0,1,K*nbrOfSetups),'b','LineWidth',2);
% plot(sort((CH_P_MMSE(:))),linspace(0,1,K*nbrOfSetups),'r','LineWidth',2);
% plot(sort((CH_P_RZF(:))),linspace(0,1,K*nbrOfSetups),'g','LineWidth',2);
% xlabel('Channel hardening');
% ylabel('CDF','Interpreter','Latex');
% legend({'MR','LP-MMSE','P-MMSE','P-RZF'},'Interpreter','tex','Location','SouthEast');
% saveas(gcf,[name 'ChD'],'fig')
% 
% % Plot favorable propagation results
% figure;
% hold on
% box on;
% plot(sort((FP_MR(:))),linspace(0,1,(K^2 -K)*nbrOfSetups),'k','LineWidth',2);
% plot(sort((FP_LP_MMSE(:))),linspace(0,1,(K^2 -K)*nbrOfSetups),'b','LineWidth',2);
% plot(sort((FP_P_MMSE(:))),linspace(0,1,(K^2 -K)*nbrOfSetups),'r','LineWidth',2);
% plot(sort((FP_P_RZF(:))),linspace(0,1,(K^2 -K)*nbrOfSetups),'g','LineWidth',2);
% xlabel('Favorable propagation');
% ylabel('CDF','Interpreter','Latex');
% legend({'MR','LP-MMSE','P-MMSE','P-RZF'},'Interpreter','tex','Location','SouthEast');
% saveas(gcf,[name 'FP'],'fig')
% 
% 
% %% Plot simulation results
% figure; hold on; box on;
% plot(sort(SE_MR_sCSI_tot(:)),linspace(0,1,K*nbrOfSetups),'k-','LineWidth',2);
% plot(sort(SE_MR_BE_tot(:)),linspace(0,1,K*nbrOfSetups),'k--','LineWidth',2);
% plot(sort(SE_MR_DLPE_tot(:)),linspace(0,1,K*nbrOfSetups),'k:','LineWidth',2);
% plot(sort(SE_MR_pCSI_tot(:)),linspace(0,1,K*nbrOfSetups),'m','LineWidth',2);
% plot(sort(SE_MR_BE_InfTau_tot(:)),linspace(0,1,K*nbrOfSetups),'c--','LineWidth',2);
% xlabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
% ylabel('CDF','Interpreter','Latex');
% legend({'MR (sCSI)',...
%     'MR (BE)',...
%     'MR (DLPE)',...
%     'MR (pCSI)','MR (BE: \tau_c = \infty)'},'Interpreter','Latex','Location','SouthEast');
% % xlim([0 14]);
% saveas(gcf,[name 'SE_MR'],'fig')
% 
% figure; hold on; box on;
% plot(sort(SE_LP_MMSE_sCSI_tot(:)),linspace(0,1,K*nbrOfSetups),'b-','LineWidth',2);
% plot(sort(SE_LP_MMSE_BE_tot(:)),linspace(0,1,K*nbrOfSetups),'b--','LineWidth',2);
% plot(sort(SE_LP_MMSE_DLPE_tot(:)),linspace(0,1,K*nbrOfSetups),'b:','LineWidth',2);
% plot(sort(SE_LP_MMSE_pCSI_tot(:)),linspace(0,1,K*nbrOfSetups),'m','LineWidth',2);
% plot(sort(SE_LP_MMSE_BE_InfTau_tot(:)),linspace(0,1,K*nbrOfSetups),'c--','LineWidth',2);
% xlabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
% ylabel('CDF','Interpreter','Latex');
% legend({'LP-MMSE (sCSI)',...
%     'LP-MMSE (BE)',...
%     'LP-MMSE (DLPE)',...
%     'LP-MMSE (pCSI)','MR (BE: \tau_c = \infty)'},'Interpreter','Latex','Location','SouthEast');
% % xlim([0 14]);
% saveas(gcf,[name 'SE_LP_MMSE'],'fig')
% 
% figure; hold on; box on;
% plot(sort(SE_P_MMSE_sCSI_tot(:)),linspace(0,1,K*nbrOfSetups),'r-','LineWidth',2);
% plot(sort(SE_P_MMSE_BE_tot(:)),linspace(0,1,K*nbrOfSetups),'r--','LineWidth',2);
% plot(sort(SE_P_MMSE_DLPE_tot(:)),linspace(0,1,K*nbrOfSetups),'r:','LineWidth',2);
% plot(sort(SE_P_MMSE_pCSI_tot(:)),linspace(0,1,K*nbrOfSetups),'m','LineWidth',2);
% plot(sort(SE_P_MMSE_BE_InfTau_tot(:)),linspace(0,1,K*nbrOfSetups),'c--','LineWidth',2);
% xlabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
% ylabel('CDF','Interpreter','Latex');
% legend({'P-MMSE (sCSI)',...
%     'P-MMSE (BE)',...
%     'P-MMSE (DLPE)',...
%     'P-MMSE (pCSI)','MR (BE: \tau_c = \infty)'},'Interpreter','Latex','Location','SouthEast');
% % xlim([0 14]);
% saveas(gcf,[name 'SE_P_MMSE'],'fig')
% 
% figure; hold on; box on;
% plot(sort(SE_P_RZF_sCSI_tot(:)),linspace(0,1,K*nbrOfSetups),'g-','LineWidth',2);
% plot(sort(SE_P_RZF_BE_tot(:)),linspace(0,1,K*nbrOfSetups),'g--','LineWidth',2);
% plot(sort(SE_P_RZF_DLPE_tot(:)),linspace(0,1,K*nbrOfSetups),'g:','LineWidth',2);
% plot(sort(SE_P_RZF_pCSI_tot(:)),linspace(0,1,K*nbrOfSetups),'m','LineWidth',2);
% plot(sort(SE_P_RZF_BE_InfTau_tot(:)),linspace(0,1,K*nbrOfSetups),'c--','LineWidth',2);
% xlabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
% ylabel('CDF','Interpreter','Latex');
% legend({'P-RZF (sCSI)',...
%     'P-RZF(BE)',...
%     'P-RZF (DLPE)',...
%     'P-RZF (pCSI)','MR (BE: \tau_c = \infty)'},'Interpreter','Latex','Location','SouthEast');
% % xlim([0 14]);
% saveas(gcf,[name 'SE_P_RZF'],'fig')
% 
% figure; 
% semilogx(sort(NMSE_MR_tot(:)),linspace(0,1,K*nbrOfSetups),'k-','LineWidth',2);
% hold on; box on;
% semilogx(sort(NMSE_MR_BE_tot(:)),linspace(0,1,K*nbrOfSetups),'k--','LineWidth',2);
% semilogx(sort(NMSE_MR_DLPE_tot(:)),linspace(0,1,K*nbrOfSetups),'k:','LineWidth',2);
% semilogx(sort(NMSE_MR_BE_InfTau_tot(:)),linspace(0,1,K*nbrOfSetups),'c--','LineWidth',2);
% xlabel('Normalized Mean Square Error','Interpreter','Latex');
% ylabel('CDF','Interpreter','Latex');
% legend({'MR (sCSI)',...
%     'MR (BE)',...
%     'MR (DLPE)','MR (BE: \tau_c = \infty)'},'Interpreter','Latex','Location','SouthEast');
% saveas(gcf,[name 'NMSE_MR'],'fig')
% 
% figure;
% semilogx(sort(NMSE_LP_MMSE_tot(:)),linspace(0,1,K*nbrOfSetups),'b-','LineWidth',2);
% hold on; box on;
% semilogx(sort(NMSE_LP_MMSE_BE_tot(:)),linspace(0,1,K*nbrOfSetups),'b--','LineWidth',2);
% semilogx(sort(NMSE_LP_MMSE_DLPE_tot(:)),linspace(0,1,K*nbrOfSetups),'b:','LineWidth',2);
% semilogx(sort(NMSE_LP_MMSE_BE_InfTau_tot(:)),linspace(0,1,K*nbrOfSetups),'c--','LineWidth',2);
% xlabel('Normalized Mean Square Error','Interpreter','Latex');
% ylabel('CDF','Interpreter','Latex');
% legend({'LP-MMSE (sCSI)',...
%     'LP-MMSE (BE)',...
%     'LP-MMSE (DLPE)','MR (BE: \tau_c = \infty)'},'Interpreter','Latex','Location','SouthEast');
% saveas(gcf,[name 'NMSE_LP_MMSE'],'fig')
% 
% figure; 
% semilogx(sort(NMSE_P_MMSE_tot(:)),linspace(0,1,K*nbrOfSetups),'r-','LineWidth',2);
% hold on; box on;
% semilogx(sort(NMSE_P_MMSE_BE_tot(:)),linspace(0,1,K*nbrOfSetups),'r--','LineWidth',2);
% semilogx(sort(NMSE_P_MMSE_DLPE_tot(:)),linspace(0,1,K*nbrOfSetups),'r:','LineWidth',2);
% semilogx(sort(NMSE_P_MMSE_BE_InfTau_tot(:)),linspace(0,1,K*nbrOfSetups),'c--','LineWidth',2);
% xlabel('Normalized Mean Square Error','Interpreter','Latex');
% ylabel('CDF','Interpreter','Latex');
% legend({'P-MMSE (sCSI)',...
%     'P-MMSE (BE)',...
%     'P-MMSE (DLPE)','MR (BE: \tau_c = \infty)'},'Interpreter','Latex','Location','SouthEast');
% saveas(gcf,[name 'NMSE_P_MMSE'],'fig')
% 
% figure; 
% semilogx(sort(NMSE_P_RZF_tot(:)),linspace(0,1,K*nbrOfSetups),'g-','LineWidth',2);
% hold on; box on;
% semilogx(sort(NMSE_P_RZF_BE_tot(:)),linspace(0,1,K*nbrOfSetups),'g--','LineWidth',2);
% semilogx(sort(NMSE_P_RZF_DLPE_tot(:)),linspace(0,1,K*nbrOfSetups),'g:','LineWidth',2);
% semilogx(sort(NMSE_P_RZF_BE_InfTau_tot(:)),linspace(0,1,K*nbrOfSetups),'c--','LineWidth',2);
% xlabel('Normalized Mean Square Error','Interpreter','Latex');
% ylabel('CDF','Interpreter','Latex');
% legend({'P-RZF (sCSI)',...
%     'P-RZF (BE)',...
%     'P-RZF (DLPE)','MR (BE: \tau_c = \infty)'},'Interpreter','Latex','Location','SouthEast');
% saveas(gcf,[name 'NMSE_P_RZF'],'fig')
save(name)
end
