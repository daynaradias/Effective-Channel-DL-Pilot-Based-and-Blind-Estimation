function [pilotIndexDL,W] = pilotAssignment(nbrOfUEs,tau_dp,CHandFP)

pilotIndexDL = zeros(nbrOfUEs,1);

for k = 1:nbrOfUEs
    if k <= tau_dp
        pilotIndexDL(k) = k;
    else
        %Compute received power from to the master AP from each pilot
        pilotinterference = zeros(tau_dp,1);
        
        for t = 1:tau_dp
            pilotinterference(t) = sum( CHandFP( k, pilotIndexDL(1:k-1)==t) );
        end
        
        %Find the pilot with the least receiver power
        [~,bestpilot] = min(pilotinterference);
        pilotIndexDL(k) = bestpilot;
    end    
end
% pilotIndexDL = [1:tau_dp randi(tau_dp,1,nbrOfUEs-tau_dp)]';

W = zeros(nbrOfUEs,nbrOfUEs);
for k =1:nbrOfUEs
    for j =1:nbrOfUEs
        W (k,j) = isequal (pilotIndexDL(k),  pilotIndexDL(j));
    end
end


% ChD = diag(CHandFP);
% UEpriority = QoS.*(w*DS + (1-w)*ChD);
% MAX = 1; MIN = 0;
% UEpriority = (MAX - MIN)* ((UEpriority - min(UEpriority))./(max(UEpriority) - min(UEpriority))) + MIN;
% % [~,indexUEpriority] = sort(UEpriority,1,'descend');
% %UEs with lower priority first
% [~,indexUEpriority] = sort(UEpriority,1,'ascend');

end