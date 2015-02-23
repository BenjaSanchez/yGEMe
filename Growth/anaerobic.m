%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = anaerobic(model)
% Transforms the model to anaerobic (OBS: MODEL SPECIFIC)
%
% INPUT:    a GEM model (in COBRA or RAVEN format)
% OUTPUT:   The modified model
%
% Benjamín J. Sánchez. Last edited: 2014-12-08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = anaerobic(model)

if strcmp(model.description,'yeast_4.05.xml')
    %No modifications suggested in Dobson et al. 2010
    rxns2block = {'O2 transport'};
    rxns2open  = {};

elseif strcmp(model.description,'yeast_5.01_model.xml')
    %Modifications performed as suggested by Heavner et al. 2012
    rxns2block = {'oxygen exchange'
                  'lipid pseudoreaction'};
    rxns2open  = {'lipid pseudoreaction [no 14-demethyllanosterol, no ergosta-5,7,22,24(28)-tetraen-3beta-ol]'
                  'ergosterol exchange'
                  'lanosterol exchange'
                  'zymosterol exchange'
                  'phosphatidate exchange'};
    
elseif strcmp(model.description,'yeast_6.06_cobra.xml') || strcmp(model.description,'yeast_7.11_cobra.xml') 
    %Modifications performed as suggested by Heavner et al. 2013
    rxns2block = {'oxygen exchange'};
    rxns2open  = {'ergosterol exchange'
                  'lanosterol exchange'
                  'zymosterol exchange'
                  'episterol exchange'
                  'fecosterol exchange'
                  'ergosta-5,7,22,24(28)-tetraen-3beta-ol exchange'};
end

%Block corresponding rxns:
for i = 1:length(rxns2block)
    rxn_pos           = find(strcmp(rxns2block(i),model.rxnNames))';
    model.lb(rxn_pos) = 0;
    model.ub(rxn_pos) = 0;
end

%Open corresponding rxns:
for i = 1:length(rxns2open)
    rxn_pos           = find(strcmp(rxns2open(i),model.rxnNames))';
    model.lb(rxn_pos) = -1000;
    model.ub(rxn_pos) = +1000;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%