%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ep_int = epistasis(dataset)
% Studies epistatic interactions, comparing single KO and double KO data
% (works for both model and experimental results).
%
% INPUT:    dataset including SDc (single KO data in complex media) and
%           DDc (double KO data in complex media).
% OUTPUT:   epistatic interactions (the first 2 rows are the combination of
%           genes, and the third has the epistatic score.  
%
% Benjamín J. Sánchez. Last edited: 2014-11-26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ep_int = epistasis(dataset)

DDc    = dataset.DDc;   %Double Gene deletions
ep_int = DDc(:,1:3);    %Epistatic Interactions
m      = length(DDc);

disp('Calculating experimental epistasis...')
for i = 1:m
    ep_int{i,3} = DDc{i,5} - DDc{i,3}*DDc{i,4};
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%