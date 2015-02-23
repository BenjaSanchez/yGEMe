%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [mets_u,rxns_u,comps] = uniqueParts(model)
% Finds the unique metabolites, unique reactions and compartments of a GEM
% model.
%
% INPUT: a GEM model (in COBRA or RAVEN format)
% OUTPUTS: a list of unique metabolites (mets_u), unique reactions (rxns_u)
%          and compartments (comps).
% 
% Benjamín J. Sánchez. Last edited: 2014-12-01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mets_u,rxns_u,comps] = uniqueParts(model)

%Start from the complete lists:
mets_u = model.metNames;
rxns_u = model.rxnNames;
comps  = model.compNames;

%Reactions: First form unique S-matrix, in which stochiometric coefficients
%will be assigned to only one of the non-unique metabolites (the last one).
S     = full(model.S);
[M,N] = size(S);
S_u   = zeros(M,N);
for i = 1:N
    %For each rxn, obtain a list with the involved metabolites and the
    %stochiometric coefficients:
    rxn_col   = S(:,i);
    rxn_pos   = find(rxn_col);
    rxn_mets  = mets_u(rxn_pos);
    rxn_coefs = rxn_col(rxn_pos);
    
    for j = 1:length(rxn_mets)
        %For each metabolite in the list, find the last position in which
        %it is repeated in the model. Therefore, the rxn will involve only
        %the last of each non-unique metabolites:
        for k = 1:M
            if strcmp(rxn_mets(j),mets_u(k))
                K = k;
            end
        end
        S_u(K,i) = rxn_coefs(j);
    end
    disp(['Computing unique matrix: ' num2str(i) '/' num2str(N) ' rows ready']);
end

%Now, find in S_u equal columns, and assign those rxns as repeated:
for i = 1:N-1
    for j = i:N
        dist = S_u(:,i)-S_u(:,j);
        if sum(dist.^2) == 0
            rxns_u(j) = rxns_u(i);
        end
    end
end

%Delete repeated names:
mets_u = delete_repeated(mets_u);
rxns_u = delete_repeated(rxns_u);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%