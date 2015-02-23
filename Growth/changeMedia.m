%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = changeMedia(model,media)
% Changes the media composition.
%
% INPUT:    a GEM model (in COBRA or RAVEN format) and the media type
%           ('minimal' or 'complex')
% OUTPUT:   The modified model for the corresponding media
%
% Benjamín J. Sánchez. Last edited: 2014-12-05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = changeMedia(model,media)

%LB parameter (manually optimized):
b = 0.1;

%Common constraints:
model = allowUptake(model,'O2',1000,media);
model = allowUptake(model,'oxygen',1000,media);
model = allowUptake(model,'sulphate',1000,media);
model = allowUptake(model,'sulfate',1000,media);
model = allowUptake(model,'phosphate',1000,media);

if strcmp(media,'nitrogen_lim')
    model = allowUptake(model,'glucose',1000,media);
    model = allowUptake(model,'ammonium',1,media);
    model = allowUptake(model,'NH3',1,media);
    
else
    model = allowUptake(model,'glucose',10,media);
    model = allowUptake(model,'ammonium',1000,media);
    model = allowUptake(model,'NH3',1000,media);
    
    if strcmp(media,'complex')
        %Aminoacids:
        model = allowUptake(model,'alanine',b,media);
        model = allowUptake(model,'cysteine',b,media);
        model = allowUptake(model,'aspartate',b,media);
        model = allowUptake(model,'glutamate',b,media);
        model = allowUptake(model,'phenylalanine',b,media);
        model = allowUptake(model,'glycine',b,media);
        model = allowUptake(model,'histidine',b,media);
        model = allowUptake(model,'isoleucine',b,media);
        model = allowUptake(model,'lysine',b,media);
        model = allowUptake(model,'leucine',b,media);
        model = allowUptake(model,'methionine',b,media);
        model = allowUptake(model,'asparagine',b,media);
        model = allowUptake(model,'proline',b,media);
        model = allowUptake(model,'glutamine',b,media);
        model = allowUptake(model,'arginine',b,media);
        model = allowUptake(model,'serine',b,media);
        model = allowUptake(model,'threonine',b,media);
        model = allowUptake(model,'valine',b,media);
        model = allowUptake(model,'tryptophan',b,media);
        model = allowUptake(model,'tyrosine',b,media);
        %Nucleotides:
        model = allowUptake(model,'adenine',b,media);
        model = allowUptake(model,'thymine',b,media);
        model = allowUptake(model,'guanine',b,media);
        model = allowUptake(model,'cytosine',b,media);
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%