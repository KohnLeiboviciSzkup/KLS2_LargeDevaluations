function [measure diff_M] = KLS2_measure(s,r)              
    
%% Computational objects
numZ = length(r.z_grid);
numA = length(r.a_grid);           
  
%Productivity
pi = r.z_pi;
P = r.z_P;

Mnew = zeros(numA,numZ);
Mnew(1,:) = pi'; 

ap_ind = r.ap_ind;

iter = 0;
diff_M = 1;

while diff_M>s.eps_sm

    Mold = Mnew;
    Mnew = zeros(numA,numZ);
    
    for i=1:numA %Old asset state
        for j=1:numZ %Old productivity state
            
            Mnew(ap_ind(i,j),:) = Mnew(ap_ind(i,j),:) + Mold(i,j)*P(j,:);
            
        end
    end

    diff_M = norm(Mnew(:)-Mold(:))/norm(1+Mold(:));
    
%     if mod(iter,1000)==0 %True if iter is divisible by 1000
%         disp(['diff_M: ' num2str(diff_M)]);        
%         if iter==5000
%             break;
%         end
%     end    
    
    iter = iter + 1;
    
end

measure = Mnew;  

end
