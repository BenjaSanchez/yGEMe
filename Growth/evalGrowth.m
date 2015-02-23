%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% results = evalGrowth(model,toolbox)
% Studies kinetic predictions of the model.
%
% INPUT:    a GEM model (in COBRA or RAVEN format) and the toolbox name
% OUTPUT:   *results.data: a cell file with 5 columns
%               *Media ('carbon_lim' or 'nitrogen_lim')
%               *O2 presence ('aerobic' or 'anaerobic')
%               *Substrate uptake (mmol/gDWh)
%               *Experimental measure of growth rate [1/h]
%               *Experimental error of measurement [1/h]
%               *Model prediction of growth rate [1/h]
%           *results.score: Squared differences sum between model
%            predictions and experimental data
%
% Benjamín J. Sánchez. Last edited: 2015-02-23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = evalGrowth(model,toolbox)

%Load exp data:
cd Growth/exp_data
[exp_data,~]   = xlsread('exp_data.xlsx',1,'C2:E100');
[~,conditions] = xlsread('exp_data.xlsx',1,'A2:B100');
results.data   = cell(length(exp_data),7);
cd ../..

for i = 1:length(exp_data)
    results.data{i,1} = conditions{i,1};
    results.data{i,2} = conditions{i,2};
    results.data{i,3} = exp_data(i,1);
    results.data{i,4} = exp_data(i,2);
    results.data{i,5} = exp_data(i,3);
    
    %Change media:
    model_i = changeMedia(model,conditions{i,1});
    if strcmp(conditions{i,2},'anaerobic')
        model_i = anaerobic(model_i);
    end
    
    %Solve FBA:
    if strcmp(toolbox,'COBRA')
        sol = optimizeCbModel(model_i);
    elseif strcmp(toolbox,'RAVEN')
        sol = solveLP(model_i);
    end
    
    results.data{i,6} = sol.x(find_exc(model,conditions{i,1}));
    results.data{i,7} = abs(sol.f);
    
end

%Calculate difference score:
diff          = cell2mat(results.data(:,4)) - cell2mat(results.data(:,7));
results.score = sum(diff.^2);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pos = find_exc(model,condition)
%Finds the position of an active exchange rxn.

if strcmp(condition,'nitrogen_lim')
    met_name = 'ammonium';
elseif strcmp(condition,'carbon_lim')
    met_name = 'glucose';
end

%Find all metabolites named as the query (usually more than one will appear
%because they have the same name in all compartments). The prefixes
%'alpha-D-', 'D-' and 'L-' will also be tried (useful for sugars and
%aminoacids):
met_pos = find(strcmp(model.metNames,['alpha-D-' met_name]))';
if isempty(met_pos)
    met_pos = find(strcmp(model.metNames,['D-' met_name]))';
    if isempty(met_pos)
        met_pos = find(strcmp(model.metNames,['L-' met_name]))';
        if isempty(met_pos)
            met_pos = find(strcmp(model.metNames,met_name))';
        end
    end
end

%For all hits find an exchange rxn with LB<0 if uptake or LB>0 if
%excretion. If found, return pos:
for i = met_pos
    met_rxns = find(model.S(i,:));
    for j = met_rxns
        if length(find(model.S(:,j))) == 1
            if sum(model.S(:,j)) == -1 && model.lb(j) < 0
                pos = j;
            elseif sum(model.S(:,j)) == 1 && model.ub(j) > 0
                pos = j;
            end
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%