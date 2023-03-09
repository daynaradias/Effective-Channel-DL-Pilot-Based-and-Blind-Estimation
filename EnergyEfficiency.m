function [Ee_MR , Ee_LP_MMSE, Ee_P_MMSE, Ee_P_RZF,...
    Gen_Ee_MR , Gen_Ee_LP_MMSE, Gen_Ee_P_MMSE, Gen_Ee_P_RZF,...
    Ee_MR_BE , Ee_LP_MMSE_BE, Ee_P_MMSE_BE, Ee_P_RZF_BE,...
    Ee_MR_DLPE , Ee_LP_MMSE_DLPE, Ee_P_MMSE_DLPE, Ee_P_RZF_DLPE] ...
    = EnergyEfficiency (w_MR,w_LP_MMSE,w_P_MMSE,w_PRZF,...
            SE_MR , SE_LP_MMSE, SE_P_MMSE, SE_P_RZF,...
            Gen_SE_MR , Gen_SE_LP_MMSE, Gen_SE_P_MMSE, Gen_SE_P_RZF,...
            SE_MR_BE, SE_LP_MMSE_BE,SE_P_MMSE_BE,SE_P_RZF_BE,...
            SE_MR_DLPE, SE_LP_MMSE_DLPE, SE_P_MMSE_DLPE, SE_P_RZF_DLPE,...
            L,N,D,nbrOfRealizations)
        
w_MR      = abs(w_MR).^2;
w_LP_MMSE = abs(w_LP_MMSE).^2;
w_P_MMSE  = abs(w_P_MMSE).^2 ;
w_PRZF    = abs(w_PRZF).^2;
%w_L_MMSE = abs(w_L_MMSE).^2;
%w_MMSE = abs(w_L_MMSE).^2;

mean_w_MR      = sum(w_MR,3)/nbrOfRealizations;
mean_w_LP_MMSE = sum(w_LP_MMSE,3)/nbrOfRealizations;
mean_w_P_MMSE  = sum(w_P_MMSE,3)/nbrOfRealizations;
mean_w_PRZF    = sum(w_PRZF,3)/nbrOfRealizations;
%mean_w_L_MMSE = sum(w_LP_MMSE,3)/nbrOfRealizations;
%mean_w_MMSE = sum(w_LP_MMSE,3)/nbrOfRealizations;

power_w_MR      = sum(mean_w_MR,2);
power_w_LP_MMSE = sum(mean_w_LP_MMSE,2);
power_w_P_MMSE  = sum(mean_w_P_MMSE,2);
power_w_PRZF    = sum(mean_w_PRZF,2);
%power_w_L_MMSE = sum(mean_w_LP_MMSE,2);
%power_w_MMSE = sum(mean_w_LP_MMSE,2);

BW_Hz = 100e6;
%Noise figure (in dB)
noiseFigure = 8;%7;
%Compute noise power (in dBm)
noiseVariance_dBm = -174 + 10*log10(BW_Hz) + noiseFigure;

Ee_MR      = MF_EnergyEfficiency (BW_Hz,D,SE_MR,power_w_MR,L,N,noiseVariance_dBm);
Ee_LP_MMSE = MF_EnergyEfficiency (BW_Hz,D,SE_LP_MMSE,power_w_LP_MMSE,L,N,noiseVariance_dBm);
Ee_P_MMSE  = MF_EnergyEfficiency (BW_Hz,D,SE_P_MMSE,power_w_P_MMSE,L,N,noiseVariance_dBm);
Ee_P_RZF   = MF_EnergyEfficiency (BW_Hz,D,SE_P_RZF,power_w_PRZF,L,N,noiseVariance_dBm);

Gen_Ee_MR      = MF_EnergyEfficiency (BW_Hz,D,Gen_SE_MR,power_w_MR,L,N,noiseVariance_dBm);
Gen_Ee_LP_MMSE = MF_EnergyEfficiency (BW_Hz,D,Gen_SE_LP_MMSE,power_w_LP_MMSE,L,N,noiseVariance_dBm);
Gen_Ee_P_MMSE  = MF_EnergyEfficiency (BW_Hz,D,Gen_SE_P_MMSE,power_w_P_MMSE,L,N,noiseVariance_dBm);
Gen_Ee_P_RZF   = MF_EnergyEfficiency (BW_Hz,D,Gen_SE_P_RZF,power_w_PRZF,L,N,noiseVariance_dBm);

Ee_MR_BE      = MF_EnergyEfficiency (BW_Hz,D,SE_MR_BE,power_w_MR,L,N,noiseVariance_dBm);
Ee_LP_MMSE_BE = MF_EnergyEfficiency (BW_Hz,D,SE_LP_MMSE_BE,power_w_LP_MMSE,L,N,noiseVariance_dBm);
Ee_P_MMSE_BE  = MF_EnergyEfficiency (BW_Hz,D,SE_P_MMSE_BE,power_w_P_MMSE,L,N,noiseVariance_dBm);
Ee_P_RZF_BE   = MF_EnergyEfficiency (BW_Hz,D,SE_P_RZF_BE,power_w_PRZF,L,N,noiseVariance_dBm);

Ee_MR_DLPE      = MF_EnergyEfficiency (BW_Hz,D,SE_MR_DLPE,power_w_MR,L,N,noiseVariance_dBm);
Ee_LP_MMSE_DLPE = MF_EnergyEfficiency (BW_Hz,D,SE_LP_MMSE_DLPE,power_w_LP_MMSE,L,N,noiseVariance_dBm);
Ee_P_MMSE_DLPE  = MF_EnergyEfficiency (BW_Hz,D,SE_P_MMSE_DLPE,power_w_P_MMSE,L,N,noiseVariance_dBm);
Ee_P_RZF_DLPE   = MF_EnergyEfficiency (BW_Hz,D,SE_P_RZF_DLPE,power_w_PRZF,L,N,noiseVariance_dBm);

end

% Ptotal = sum(P_m) + sum(Pbh_m), for l varying from 1 to L

% ("P_m") is the power consumption at the m-th AP due to the amplifier and the circuit power consumption part
% (including the power consumption of the transceiver chains and the power consumed for signal processing)
% ("Pbh_m") is the power consumed by the backhaul link connecting the CPU and the m-th AP.

% P_m = (1/alpha_m)*power_w_P_MMSE + N*Ptc_m; ou P_m = (1/alpha_m)*E[Xl] + N*Ptc_m;
% ("alpha_m") is such that: 0 < alpha_m ? 1 is the power amplifier efficiency
% ("Ptc_m") is the internal power required to run the circuit components (e.g. converters, mixers, and filters)
% related to each antenna of the m-th AP


% Pbh_m = P0_m + (B*Se*Pbt_m), where:
% ("P0_m") is a fixed power consumption of each backhaul (traffic-independentpower) which may depend on the
% distances between the APs and the CPU and the system topology
% ("Pbt_m") is the traffic-dependent power (in Watt per bit/s)
% ("B") is the system bandwidth
% ("Se") is the sum spectral efficiency of each user "k"
% ("Ee") is the energy efficiency

% Note: The backhaul is also used to transfer the power allocation coefficients, synchronization signals, etc.
% This is done once per large-scale fading realization which stays constant for many coherence intervals.
% Therefore, they neglect the power consumed by this processing.


function Ee = MF_EnergyEfficiency (BW_Hz,D,Se_precoding,power_w_precoding,nbrOfAPs,N,noiseVariance_dBm)

%% Input Parameters (from table II)
alpha_m = 0.4;
Ptc_m   = 0.2; % in W
Pbt_m   = 0.25/10^9; % 0.25 W per Gb/s, but the unit is in W/Gb/s
P0_m    = 0.825;

alpha_m = alpha_m*ones(nbrOfAPs,1);
Ptc_m   = Ptc_m.*ones(nbrOfAPs,1);
Pbt_m   = Pbt_m.*ones(nbrOfAPs,1);
Se_m    = zeros(nbrOfAPs,1); % SE of the users served by the m-th AP
aux = zeros(nbrOfAPs,1);

%% Computing Eqs. (18), (39) - backhaul power consumption
noiseVariance_w = 10.^((noiseVariance_dBm./10)-3); % Converting from dBm to w

for indexAP = 1:nbrOfAPs
  aux(indexAP,1) = sum(power_w_precoding((indexAP-1)*N+1:indexAP*N));
end

power_w_precoding = aux;
power   = noiseVariance_w.*(power_w_precoding./1e3); % Converting from mW to W
%power   = (power_w_precoding./1e3); % Converting from mW to W

P_m = (1./alpha_m).*power + N.*Ptc_m; % Eq. (18)

%If the AP is off it is not consuming energy
APs_off_id = sum(D,2) == 0;
P_m(APs_off_id) = eps;

for indexAP = 1:nbrOfAPs
    ServedUEs_id = D(indexAP,:) == 1;
    Se_m(indexAP) = sum(Se_precoding(ServedUEs_id));
end


% Computing Eq. (19)
P0_m = P0_m.*ones(nbrOfAPs,1); % In watts
Pbh_m = P0_m + (BW_Hz.*Se_m.*Pbt_m); % Eq.(19)
Ptotal = sum(P_m) + sum(Pbh_m); % Eq.(17)
Se = sum(Se_precoding);
Ee = (BW_Hz.*Se)/Ptotal; % Eq.(21)

% Falta definir (Se), (B)
end