%method1
function [pilotIndexDL,W] = jointULandDLpilotAssignment(nbrOfUEs,tau_dp,CHandFP,pilotIndexUL)

pilotIndexDL = zeros(nbrOfUEs,1);
DLpilotinterference = zeros(tau_dp,nbrOfUEs);
for k = 1:nbrOfUEs
    if k <= tau_dp
        pilotIndexDL(k) = k;
    else
        %Compute received power from to the master AP from each pilot
        pilotinterference = zeros(tau_dp,1);
        
        for t = 1:tau_dp
            pilotinterference(t) = sum( CHandFP( k, pilotIndexDL(1:k-1)==t) );
        end
        DLpilotinterference(:,k) = pilotinterference;
        %Find the pilot with the least receiver power
        [~,bestpilot] = min(pilotinterference);
        pilotIndexDL(k) = bestpilot;
    end    
end


% W_up = zeros(nbrOfUEs,nbrOfUEs);
% W_dp = zeros(nbrOfUEs,nbrOfUEs);
% W1 = zeros(nbrOfUEs,nbrOfUEs);
% for m = 1:nbrOfUEs
%     for n=1:nbrOfUEs
%         if m ~= n 
%         W_up(m,n) = isequal(pilotIndexUL(m),pilotIndexUL(n));
%         W_dp(m,n) = isequal(pilotIndexDL(m),pilotIndexDL(n));
%         W1(m,n) =  (W_up(m,n)+ W_dp(m,n)) == 2;
%         end
%     end
% end
% 
% solve1 = sum(W1(:));

for i = 1:nbrOfUEs
    for j=1:nbrOfUEs
        if i ~= j
            aux =  (isequal(pilotIndexUL(i),pilotIndexUL(j)) + ...
                isequal(pilotIndexDL(i),pilotIndexDL(j))) == 2;
            
            if aux == 1
                DLpilotinterference(pilotIndexDL(j),j) = max(DLpilotinterference(:,j));
                
                [~,bestpilot] = min(DLpilotinterference(:,j));
                
                pilotIndexDL(j) = bestpilot;
                
            end
        end
    end
end

% W_up = zeros(nbrOfUEs,nbrOfUEs);
% W_dp = zeros(nbrOfUEs,nbrOfUEs);
% W2 = zeros(nbrOfUEs,nbrOfUEs);
% for m = 1:nbrOfUEs
%     for n=1:nbrOfUEs
%         if m ~= n 
%         W_up(m,n) = isequal(pilotIndexUL(m),pilotIndexUL(n));
%         W_dp(m,n) = isequal(pilotIndexDL(m),pilotIndexDL(n));
%         W2(m,n) =  (W_up(m,n)+ W_dp(m,n)) == 2;
%         end
%     end
% end
% 
% solve2 = sum(W2(:));

W = zeros(nbrOfUEs,nbrOfUEs);
for i =1:nbrOfUEs
    for j =1:nbrOfUEs
        W (i,j) = isequal (pilotIndexDL(i),  pilotIndexDL(j));
    end
end
