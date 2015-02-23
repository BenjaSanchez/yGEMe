%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = standardizeModel(model,toolbox)
% Standardize yeast GEM models for avoiding errors in comparison.
%
% INPUT:    The raw SBML model from the corresponding publication, and the
%           toolbox name (RAVEN or COBRA).
% OUTPUT:   The standardized model. Main corrections applied are:
%           *Removal of any compartment reference in MetNames
%           *Adition of the compNames vector (if not present)
%           *Adition of the metComps vector (if not present)
%
% Benjamín J. Sánchez. Last edited: 2014-11-12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = standardizeModel(model,toolbox)

mets = model.metNames;
M    = length(mets);

%Adition of the compNames and metComps vectors (if not present) & removal
%of any compartment reference in MetNames:
if strcmp(toolbox,'COBRA')
    comps = cell(M,1);
    for i = 1:M
        name      = mets{i};
        m         = length(name);
        %Split the met name in 2: The compartment and the metabolite name.
        %Afterwards, assign each to the corresponding list.
        k         = strfind(name,'[');
        n         = length(k);
        mets{i} = name(1:k(n)-2);
        comps{i}  = name(k(n)+1:m-1);
    end
    model.metNames  = mets;
    
    %Create both compNames (list with the copartment names) and metComps
    %(compartment in which each metabolite is, using indexing of
    %compNames):
    comps_u         = delete_repeated(comps);
    met_comps       = zeros(M,1);
    for i = 1:M
        for j = 1:length(comps_u)
            if strcmp(comps{i},comps_u{j})
                met_comps(i) = j;
            end
        end
    end
    model.compNames = comps_u;
    model.metComps  = met_comps;
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%