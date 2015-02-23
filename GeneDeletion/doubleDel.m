%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model_DD = doubleDel(model,exp_DD,toolbox,media)
% Performs double gene KO analysis to all gene cpombinations of the model.
%
% INPUT:    a GEM model (in COBRA or RAVEN format), the experimental DD
%           analysis, the toolbox name and the media type (minimal or
%           complex)
% OUTPUT:   model_DD, a mx5 cell with the model gene names in the 2 first
%           columns, the corresponding SD fitness in the next 2 columns,
%           and the DD fitness in the last column
%
% Benjamín J. Sánchez. Last edited: 2014-12-02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model_DD = doubleDel(model,exp_DD,toolbox,media)

model    = changeMedia(model,media);
m        = length(exp_DD);
model_DD = cell(m,5);

if strcmp(toolbox,'COBRA')
    WT_sol  = optimizeCbModel(model);
    for i = 1:m
        sKO1_model    = deleteModelGenes(model,exp_DD{i,1});
        sKO2_model    = deleteModelGenes(model,exp_DD{i,2});
        dKO_model     = deleteModelGenes(sKO1_model,exp_DD{i,2});
        sKO1_sol      = optimizeCbModel(sKO1_model);
        sKO2_sol      = optimizeCbModel(sKO2_model);
        dKO_sol       = optimizeCbModel(dKO_model);
        model_DD{i,1} = exp_DD{i,1};
        model_DD{i,2} = exp_DD{i,2};
        model_DD{i,3} = abs(sKO1_sol.f/WT_sol.f);
        model_DD{i,4} = abs(sKO2_sol.f/WT_sol.f);
        model_DD{i,5} = abs(dKO_sol.f/WT_sol.f);
        disp(['Calculating DD analysis in ' media ' media: ' num2str(i) '/' num2str(m) ' combinations ready.'])
    end
elseif strcmp(toolbox,'RAVEN')
    WT_sol  = solveLP(model);
    for i = 1:m
        sKO1_model = removeGenes(model,exp_DD{i,1},false,true);
        sKO2_model = removeGenes(model,exp_DD{i,2},false,true);
        try
            dKO_model = removeGenes(sKO1_model,exp_DD{i,2},false,true);
        catch
            dKO_model = sKO1_model;
        end
        sKO1_sol      = solveLP(sKO1_model);
        sKO2_sol      = solveLP(sKO2_model);
        dKO_sol       = solveLP(dKO_model);
        model_DD{i,1} = exp_DD{i,1};
        model_DD{i,2} = exp_DD{i,2};
        model_DD{i,3} = abs(sKO1_sol.f/WT_sol.f);
        model_DD{i,4} = abs(sKO2_sol.f/WT_sol.f);
        model_DD{i,5} = abs(dKO_sol.f/WT_sol.f);
        disp(['Calculating DD analysis in ' media ' media: ' num2str(i) '/' num2str(m) ' combinations ready.'])
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%