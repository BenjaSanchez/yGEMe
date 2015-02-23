%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model_SD = singleDel(model,exp_SD,toolbox,media)
% Performs single gene KO analysis to all genes of the model.
%
% INPUT:    a GEM model (in COBRA or RAVEN format), the experimental SD
%           analysis, the toolbox name and the media type (minimal or
%           complex)
% OUTPUT:   model_SD, a mx2 cell with the gene names in the first column
%           and the fitness score in the second column (mu(KO)/mu(WT))
%
% Benjamín J. Sánchez. Last edited: 2014-12-02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model_SD = singleDel(model,exp_SD,toolbox,media)

model    = changeMedia(model,media);
m        = length(exp_SD);
model_SD = cell(m,2);

if strcmp(toolbox,'COBRA')
    WT_sol  = optimizeCbModel(model);
    for i = 1:m
        KO_model      = deleteModelGenes(model,exp_SD{i,1});
        KO_sol        = optimizeCbModel(KO_model);
        model_SD{i,1} = exp_SD{i,1};
        model_SD{i,2} = abs(KO_sol.f/WT_sol.f);
        disp(['Calculating SD analysis in ' media ' media: ' num2str(i) '/' num2str(m) ' genes ready.'])
    end
elseif strcmp(toolbox,'RAVEN')
    WT_sol  = solveLP(model);
    for i = 1:m
        KO_model      = removeGenes(model,exp_SD{i,1},false,true);
        KO_sol        = solveLP(KO_model);
        model_SD{i,1} = exp_SD{i,1};
        model_SD{i,2} = abs(KO_sol.f/WT_sol.f);
        disp(['Calculating SD analysis in ' media ' media: ' num2str(i) '/' num2str(m) ' genes ready.'])
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%