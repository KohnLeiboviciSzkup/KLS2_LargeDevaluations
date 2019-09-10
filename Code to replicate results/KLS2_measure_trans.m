function M = KLS2_measure_trans(r,N,apt_ind,M0) 

    numZ = length(r.z_grid);
    numA = length(r.a_grid);

    M = zeros(numA,numZ,N);
    M(:,:,1) = M0;
    M(:,:,2) = M0; %Unexpected shocks happen at start of period 2, after period 1 decisions are made and states of period 2 determined
    Mnew = M0;

    for n=2:N-1

        Mold = Mnew;
        Mnew = zeros(numA,numZ);

        for i=1:numA %Old asset state
            for j=1:numZ %Old productivity state

                Mnew(apt_ind(i,j,n),:) = Mnew(apt_ind(i,j,n),:) + Mold(i,j)*r.z_P(j,:);

            end
        end

        M(:,:,n+1) = Mnew;

    end

            
end
