%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GCC = globalClustCoeff(M)
% Computes the global clustering coefficient, which indicates the fraction
% of closed triplets among all triplets, both closed and open.
% 
% INPUT: Connectivity matrix (M)
% OUTPUT: Global Clustering Coefficient (GCC)
% 
% Benjamín Sánchez. Last edited: 2014-11-20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GCC = globalClustCoeff(M)

m   = length(M);
C_a = 0;    %Triplet connection (any, open or closed)
C_c = 0;    %Triplet closed connection

for i = 1:m
    for j = 1:m
        for k = 1:m
            if i~=j && j~=k && k~=i
                if M(i,j)+M(j,k)+M(k,i) == 2
                    %Open triplet:
                    C_a = C_a+1;
                elseif M(i,j)+M(j,k)+M(k,i) == 3
                    %Closed triplet:
                    C_a = C_a+1;
                    C_c = C_c+1;
                end
            end
        end
    end
    disp(['Computing GCC: ' num2str(i) '/' num2str(m) ' nodes ready']);
end

GCC = C_c/C_a;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%