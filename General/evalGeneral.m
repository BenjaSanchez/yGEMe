%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% general = evalGeneral(model,toolbox)
% Counts all relevant features of a GEM model.
%
% INPUT:    a GEM model (in COBRA or RAVEN format) and the toolbox name.
% OUTPUT:   a vector with the following information:
%           #rxns
%           #mets
%           #genes
%           #compartments
%           #unique rxns
%           #unique mets
%           #GPR associations/#genes
%           % EC codes in reactions
%           % ChEBI identifiers in metabolites
%           % ChEBI identifiers in metabolites
%
% Benjamín J. Sánchez. Last edited: 2015-02-23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function general = evalGeneral(model,toolbox)

general    = zeros(10,1);

%Count rxns, mets and genes:
general(1) = length(model.rxns);
general(2) = length(model.mets);
general(3) = length(model.genes);

%Count comps, unique mets and unique rxns:
[mets_u,rxns_u,comps] = uniqueParts(model);
general(4) = length(comps);
general(5) = length(rxns_u);
general(6) = length(mets_u);

%Count GPR associations:
general(7) = sum(sum(full(model.rxnGeneMat)))/length(model.genes);

%Count % EC numbers in reactions:
if strcmp(toolbox,'COBRA')
    general(8) = length(delete_empty(model.rxnECNumbers))/length(model.rxns)*100;
elseif strcmp(toolbox,'RAVEN')
    general(8) = length(delete_empty(model.eccodes))/length(model.rxns)*100;
end

%Count % ChEBI identifiers in metabolites:
if strcmp(toolbox,'COBRA')
    general(9) = length(delete_empty(model.metChEBIID))/length(model.mets)*100;
elseif strcmp(toolbox,'RAVEN')
    general(9) = length(delete_empty(model.metMiriams))/length(model.mets)*100;
end

%Count % InChI identifiers in metabolites:
if strcmp(toolbox,'COBRA')
    general(10) = length(delete_empty(model.metInChIString))/length(model.mets)*100;
elseif strcmp(toolbox,'RAVEN')
    general(10) = length(delete_empty(model.inchis))/length(model.mets)*100;
end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%