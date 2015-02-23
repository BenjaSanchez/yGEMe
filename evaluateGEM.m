%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% results = evaluateGEM(file_name,toolbox)
% Performs a complete evaluation of the yeast GEM.
%
% INPUT:    The .mat file name containing the model structure and the
%           toolbox to be used for performing the optimizations
%           (COBRA or RAVEN).
% OUTPUT:   Results structure with the following fields:
%           *name: GEM name
%           *toolbox: Toolbox name
%           *model: GEM
%           *general: General results
%           *connectivity: Connectivity results
%           *growth: Growth results
%           *geneDeletion: Gene deletion results
%           For a detailed explanation of each field go to the
%           corresponding function.
%
% Benjamín J. Sánchez. Last edited: 2015-02-23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = evaluateGEM(file_name,toolbox)

results.name    = file_name;
results.toolbox = toolbox;

addpath General
addpath Connectivity
addpath Growth
addpath GeneDeletion

if strcmp(toolbox,'COBRA')
    initCobraToolbox
end

model = load(file_name,'model');
model = model.model;

%Standardize model (metabolite naming convention, create compartment
%vectors, etc):
model = standardizeModel(model,toolbox);

results.model = model;

%Count basic model features:
results.general = evalGeneral(model,toolbox);

%Check unreachable metabolites, blocked reactions and calculate other
%connectivity measurements:
results.connectivity = evalConnectivity(model,toolbox);

%Evaluate kinetic predictions of model using experimental data:
results.growth = evalGrowth(model,toolbox);

%Evaluate genetic predictions of model using experimental data:
results.geneDeletion = evalGeneDeletion(model,toolbox);

rmpath General
rmpath Connectivity
rmpath Growth
rmpath GeneDeletion

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%