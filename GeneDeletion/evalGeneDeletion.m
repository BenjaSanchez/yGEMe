%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% results = evalGeneDeletion(model,toolbox)
% Studies genetic predictions of the model.
%
% INPUT:    a GEM model (in COBRA or RAVEN format) and the toolbox name
% OUTPUT:   a structure with the following fields:
%           *exp_KO   (experimental KO information)
%           *model_KO (model KO results)
%           *metrics  (all the predictions scores and metrics)
%
% Benjamín J. Sánchez. Last edited: 2015-02-23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = evalGeneDeletion(model,toolbox)

%Load relevant experimental data:
exp_KO = load_exp(model);

%Calculate experimental epistatic interactions:
exp_KO.epi = epistasis(exp_KO);

%Single Deletion Analysis in Minimum Media:
model_KO.SDm = singleDel(model,exp_KO.SDm,toolbox,'carbon_lim');

%Single Deletion Analysis in Complex Media:
model_KO.SDc = singleDel(model,exp_KO.SDc,toolbox,'complex');

%Double Deletion Analysis in Complex Media:
model_KO.DDc = doubleDel(model,exp_KO.DDc,toolbox,'complex');

%Calculate epistatic interactions in model:
model_KO.epi = epistasis(model_KO);

%Calculate all metrics associated to the predictions:
metrics.SDm = checkPredictions(exp_KO.SDm,model_KO.SDm);
metrics.SDc = checkPredictions(exp_KO.SDc,model_KO.SDc);
metrics.DDc = checkPredictions(exp_KO.DDc,model_KO.DDc);
metrics.epi = checkPredictions(exp_KO.epi,model_KO.epi);

%Save results:
results.exp_KO   = exp_KO;
results.model_KO = model_KO;
results.metrics  = metrics;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%