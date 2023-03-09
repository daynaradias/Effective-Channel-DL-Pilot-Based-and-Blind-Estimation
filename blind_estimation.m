% Blind estimation method, based on "Ngo, H. Q., G. Larsson, E., “No downlink pilots are needed in Massive
% MIMO�?. IEEE Transactions on Communications, 2017."
% Author: Daynara Dias Soyza
% Date: 06/16/2020
function [estimated_effective_channel_UnCorr,NMSE_blind,estimated_effective_channel_InfTau,NMSE_blind_InfTau] = blind_estimation(tau_d,K,effective_channel,nbrOfRealizations)
% tic

% symbol1= zeros(1,K,1,tau_d);
% for j=1:K
%     symbol1(1,j,1,:) = transpose(nrSymbolModulate(randi([0 1],tau_d,1),'BPSK'));
% end

BPSK = sqrt(2)/2*([1 + 1i*1; -1 - 1i*1]);
symbol = zeros(1,K,1,tau_d); 
for j=1:K
    for t =1:tau_d
        symbol(1,j,1,t) = BPSK(randi(2));
    end
end


%E{|αkk|}
mean_effective_channel = (mean(abs(effective_channel),3));

% Algorithm 1 Blind Downlink Channel Estimation Method
estimated_effective_channel_UnCorr = zeros(K,nbrOfRealizations);
NMSE_blind = zeros(K,1);

estimated_effective_channel_InfTau = zeros(K,nbrOfRealizations);
NMSE_blind_InfTau = zeros(K,1);

for i = 1:K
    
    meanSquaredEffectiveChannel = mean( sum( abs( effective_channel(i,[1:i-1 i+1:end],:) ).^2 ,2) , 3);

    Wdown = sqrt(0.5)*(randn(1,1,nbrOfRealizations,tau_d) + 1i*randn(1,1,nbrOfRealizations,tau_d));
    
    receiveSignal = effective_channel(i,i,:) .* symbol(1,i,1,:) ...
        + sum(effective_channel(i,[1:i-1 i+1:end],:).* symbol(1,[1:i-1 i+1:end],1,:) ,2) + Wdown;
       
    epsilonUnCorr = (sum(abs(receiveSignal(:,:,:,2:end)).^2,4)/(tau_d-1));
   
    ind = epsilonUnCorr > (1 + meanSquaredEffectiveChannel);
    
    estimated_effective_channel_UnCorr(i,ind) = sqrt( ( epsilonUnCorr(1,1,ind) - 1 - meanSquaredEffectiveChannel ) );
    
    estimated_effective_channel_UnCorr(i,ind==0) = mean_effective_channel(i,i);
    
    NMSE_blind(i) = (mean( (abs(permute(estimated_effective_channel_UnCorr(i,:),[3 1 2]) - effective_channel(i,i,:)) ).^2 ,3) ) ./ ( mean( abs( effective_channel(i,i,:) ).^2, 3 ) );    

%%  For tau_d -> inf
    epsilonInfTau = sum ( abs(effective_channel(i,:,:)).^2, 2 ) + 1;
    
    indInfTau = epsilonInfTau > (1 + meanSquaredEffectiveChannel);
    
    estimated_effective_channel_InfTau(i,indInfTau) = sqrt( ( epsilonInfTau(1,1,indInfTau) - 1 - meanSquaredEffectiveChannel ) );
    
    estimated_effective_channel_InfTau(i,indInfTau==0) = mean_effective_channel(i,i);
    
    NMSE_blind_InfTau(i) = (mean( (abs(permute(estimated_effective_channel_InfTau(i,:),[3 1 2]) - effective_channel(i,i,:)) ).^2 ,3) ) ./ ( mean( abs( effective_channel(i,i,:) ).^2, 3 ) );    

end

end