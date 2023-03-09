

function  [pilotEstimatedEffectiveChannel,NMSE] = DLbeamformingTraining(nbrOfUEs,nbrOfRealizations,tau_DL_pilot,W,effectiveChannels)

pilotEstimatedEffectiveChannel = zeros(nbrOfUEs,nbrOfRealizations);

NMSE = zeros(nbrOfUEs,1);
    
for k = 1:nbrOfUEs
    
    noise_DL_pilot = sqrt(0.5)*(randn(1,1,nbrOfRealizations) + 1i*randn(1,1,nbrOfRealizations));
    
    receivedDLpilotSignal =  sqrt(tau_DL_pilot) *  effectiveChannels (k,k,:) + ...
        sqrt(tau_DL_pilot)* sum( effectiveChannels(k,[1:k-1 k+1:end],:) .* W (k,[1:k-1 k+1:end]), 2)  ...
        + noise_DL_pilot;
    
    A(1,:) =  effectiveChannels(k,k,:);
    B(1,:) = receivedDLpilotSignal;
    covariance1 = cov( A , B);
    
    pilotEstimatedEffectiveChannel1 = mean (effectiveChannels (k,k,:), 3) ...
        +  (covariance1 (1,2) / covariance1 (2,2) ).* ( receivedDLpilotSignal - mean(receivedDLpilotSignal , 3) );
    
    NMSE(k) = (mean( (abs( pilotEstimatedEffectiveChannel1 - effectiveChannels(k,k,:) ) ).^2 ,3) ) ./ ( mean( abs( effectiveChannels(k,k,:) ).^2, 3 ) );
    
    pilotEstimatedEffectiveChannel(k,:) = pilotEstimatedEffectiveChannel1;
      
end

end