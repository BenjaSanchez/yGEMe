%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LCC = localClustCoeff(M)
% Computes the local clustering coefficient, the i-th position shows the
% fraction of existing connections over all possible connections between
% the neighborhood of metabolite i (not including i itself).
% 
% INPUT:    Connectivity matrix (M)
% OUTPUT:   Local Clustering Coefficient (LCC)
% 
% Benjamín Sánchez. Last edited: 2014-12-03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function LCC = localClustCoeff(M)

m   = length(M);
LCC = zeros(m,1);
for i = 1:m
    %Find neighbors and compute maximum possible links between them:
    neighbors = find(M(i,:));
    K         = length(neighbors);
    Nmax      = K*(K-1)/2;
    N         = 0;
    
    %Count number of existinting links between neihborhood:
    for j = neighbors
        for k = neighbors
            if j > k && M(j,k) == 1
                N = N+1;
            end
        end
    end
    
    %Calculate local clustering coefficient:
    if K == 1 || K == 0
        LCC(i) = 1;
    else
        LCC(i) = N/Nmax;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%