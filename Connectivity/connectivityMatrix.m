%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M = conectivity_matrix(S)
% Calculates the connectivity matrix, in which 2 metabolites i and j are
% conected (and therefore the position (i,j) and (j,i) are = 1) if they are
% both present in at least one shared reaction.
% 
% INPUT: Stochiometric matrix (S)
% OUTPUT: Connectivity matrix (M)
% 
% Benjamín Sánchez. Last edited: 2014-11-20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M = connectivityMatrix(S)

[m,n] = size(S);
M     = zeros(m);

for i = 1:n
    %Find non-zero elements in S for reaction i:
    mets_in_rxn = find(S(:,i))';
    %Connect all metabolites that are present in reaction i:
    for j = mets_in_rxn
        for k = mets_in_rxn
            if j ~= k
                M(j,k) = 1;
            end
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%