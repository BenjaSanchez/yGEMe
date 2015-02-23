%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UM = unreach_mets(model,SPM)
% Finds the % of unreachable metabolites (from the extracellular media).
%
% INPUT:    a GEM model (in COBRA or RAVEN format) and SPM (shortest paths
%           matrix, in which the (i,j) position shows the number of steps
%           for getting from i to j)
% OUTPUT:   The % of unreachable metabolites
%
% Benjamín J. Sánchez. Last edited: 2015-02-19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function UM = unreach_mets(model,SPM)

[m,n] = size(model.S);
UM    = ones(m,1);

%Find extracelluar metabollites:
extra_pos = zeros(m,1);
for i = 1:n
    rxn_coeffs = model.S(:,i);
    rxn_mets   = find(rxn_coeffs);
    if length(rxn_mets) == 1 && sum(rxn_coeffs) == -1
        extra_pos(rxn_mets) = 1;
    end
end

%Try to find a path from any extracellular metabolite to each metabolite:
for i = 1:m
    for j = find(extra_pos)'
        if i == j || SPM(j,i) ~= 0
            UM(i) = 0;
        end
    end
end

UM = sum(UM)/m*100;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%