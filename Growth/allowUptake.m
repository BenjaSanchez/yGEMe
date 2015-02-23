%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = allowUptake(model,met_name,b,media)
% Finds the exchange reaction of the corresponding metabolite and assigns
% the appropiate lower bound.
%
% INPUT:    a GEM model (in COBRA or RAVEN format), the metabolite name,
%           the value for the exchange constraint and the media type
%           ('minimal' or 'complex')
% OUTPUT:   The modified model
%
% Benjamín J. Sánchez. Last edited: 2014-12-05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = allowUptake(model,met_name,b,media)

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

%For all hits find an exchange rxn. If found, allow uptake:
for i = met_pos
    met_rxns = find(model.S(i,:));
    success = false;
    for j = met_rxns
        if length(find(model.S(:,j))) == 1
            %Assign to the first exchange rxn (in-take or out-take) the
            %corresponding bound (lower or upper, respectively):
            if sum(model.S(:,j)) == -1 && ~success
                model.lb(j) = -b;
                model.ub(j) = 1000;
                success = true;
                disp([met_name ' was correctly allowed for uptake in ' media ' media.'])
            elseif sum(model.S(:,j)) == 1 && ~success
                model.lb(j) = -1000;
                model.ub(j) = b;
                success = true;
                disp([met_name ' was correctly allowed for uptake in ' media ' media.'])
            
            %Block any other exchange rxn for the same metabolite:
            elseif sum(model.S(:,j)) == -1 && success
                model.lb(j) = 0;
                model.ub(j) = 0;
            elseif sum(model.S(:,j)) == 1 && success
                model.lb(j) = 0;
                model.ub(j) = 0;
            end
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%