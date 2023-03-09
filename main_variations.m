clear
%Setup padrão:
K=5;% number of UEs
tau_up = 10;%number of UL pilots
tau_dp = 10;%number of DL pilots
tau_c = 200;%coherence interval length
L=10;%100;%number of APs
N=4;%number of Antennas per AP
rho_tot = 200;%DL transmit power in mW
capacityBound = 'UnF';%'Interdonato2019'; %
APselection = 'BookEmil';%'BookEmil';'CanonicalCF';'ClusterControlCmax';'ClusterControlCmax'
alpha = 0.75;
pilotAssignmentMethod = 'nonJoint_minContamination';%'joint_random';'joint_minContamination';'nonJoint_random';'nonJoint_minContamination';

%% 1) vary pilot assignment method
varpilotAssignmentMethod = [{'joint_random'},{'nonJoint_minContamination'},{'nonJoint_random'}];%{'joint_minContamination'},
for i = 1:length(varpilotAssignmentMethod)
    pilotAssignmentMethod1 = char(varpilotAssignmentMethod(i));
    name1 = [pilotAssignmentMethod1 '_' APselection '_' capacityBound '_L' num2str(L) '_N' num2str(N) '_K' num2str(K) '_rhoD' num2str(rho_tot)...
    '_tauC' num2str(tau_c) '_tauUP' num2str(tau_up) '_tauDP' num2str(tau_dp)];
    mainFunction(L,N,K,tau_c,tau_up,tau_dp,capacityBound, rho_tot, pilotAssignmentMethod1, APselection, alpha, name1 )
    close all;
end
%% 2) vary clustering stringency  parameter alpha
varAlpha = [0.25,0.5, 0.75];%0.25, 
for i = 1:length(varAlpha)
    APselection1 = 'ClusterControlCmax';
    alpha1 = varAlpha(i);
    name1 = [pilotAssignmentMethod '_' APselection1 '_alpha' num2str(100*alpha1) '_' capacityBound '_L' num2str(L) '_N' num2str(N) '_K' num2str(K) '_rhoD' num2str(rho_tot)...
    '_tauC' num2str(tau_c) '_tauUP' num2str(tau_up) '_tauDP' num2str(tau_dp)];
    mainFunction(L,N,K,tau_c,tau_up,tau_dp,capacityBound, rho_tot,pilotAssignmentMethod, APselection1, alpha1, name1 )
    close all;
end
%% 3) vary AP selection
varAPselection = [{'ClusterControlCmax'},{'CanonicalCF'}];%
for i = 1:length(varAPselection)
    APselection1 = char(varAPselection(i));
    name1 = [pilotAssignmentMethod '_' APselection1 '_alpha' num2str(100*alpha) '_' capacityBound '_L' num2str(L) '_N' num2str(N) '_K' num2str(K) '_rhoD' num2str(rho_tot)...
    '_tauC' num2str(tau_c) '_tauUP' num2str(tau_up) '_tauDP' num2str(tau_dp)];
    mainFunction(L,N,K,tau_c,tau_up,tau_dp,capacityBound, rho_tot,pilotAssignmentMethod, APselection1, alpha, name1 )
    close all;
end
%% 4) vary DL power of APs
power1 =[1, 25, 50, 200, 400, 600, 800, 1000, 5000];%[25, 50, 200, 400, 600, 800, 1000]; %
for i = 1:length(power1)
    rho_tot1 = power1(i);
    name1 = [capacityBound '_L' num2str(L) '_N' num2str(N) '_K' num2str(K) '_rhoD' num2str(rho_tot1)...
    '_tauC' num2str(tau_c) '_tauUP' num2str(tau_up) '_tauDP' num2str(tau_dp)];
    mainFunction(L,N,K,tau_c,tau_up,tau_dp,capacityBound, rho_tot1, name1 )
    close all;
end

%% 5)vary coherence interval
coherency2 = [800,1000 ,2000 ,3000];%[25,50, 100, 300, 400, 500,600,700,800,1000 ,2000 ,3000];%[];%%[50, 100, 200, 300, 400, 500];
for j = 1:length(coherency2)
    tau_c2 = coherency2(j);
    name2 = [capacityBound '_L' num2str(L) '_N' num2str(N) '_K' num2str(K) '_rhoD' num2str(rho_tot)...
    '_tauC' num2str(tau_c2) '_tauUP' num2str(tau_up) '_tauDP' num2str(tau_dp)];
    mainFunction(L,N,K,tau_c2,tau_up,tau_dp,capacityBound, rho_tot, name2 )
    close all;
end

%% 6) number of antennas per AP
nbrOfAntennasPerAP3 = [1, 2, 3, 5];%[1, 2, 3, 4, 5];
for k = 1:length(nbrOfAntennasPerAP3)
    N3 = nbrOfAntennasPerAP3(k);
    name3 = [capacityBound '_L' num2str(L) '_N' num2str(N3) '_K' num2str(K) '_rhoD' num2str(rho_tot)...
    '_tauC' num2str(tau_c) '_tauUP' num2str(tau_up) '_tauDP' num2str(tau_dp)];
    mainFunction(L,N3,K,tau_c,tau_up,tau_dp,capacityBound, rho_tot, name3 )
    close all;
end

%% 7) vary DL pilot length
DLpilots7 = [1,5,15,20];%[5,10,15,20];%25,30,35,40];%[5,10,15,20,25,30,35,40];
for o = 1:length(DLpilots7)
    tau_dp7 = DLpilots7(o);
    name7 = [capacityBound '_L' num2str(L) '_N' num2str(N) '_K' num2str(K) '_rhoD' num2str(rho_tot)...
    '_tauC' num2str(tau_c) '_tauUP' num2str(tau_up) '_tauDP' num2str(tau_dp7)];
    mainFunction(L,N,K,tau_c,tau_up,tau_dp7,capacityBound, rho_tot, name7 )
    close all
end
%% 8) vary UL pilot length
ULpilots6 = [5,15,20];%
for n = 1:length(ULpilots6)
    tau_up6 = ULpilots6(n);
    name6 = [capacityBound '_L' num2str(L) '_N' num2str(N) '_K' num2str(K) '_rhoD' num2str(rho_tot)...
    '_tauC' num2str(tau_c) '_tauUP' num2str(tau_up6) '_tauDP' num2str(tau_dp)];
    mainFunction(L,N,K,tau_c,tau_up6,tau_dp,capacityBound, rho_tot, name6 )
    close all
end

