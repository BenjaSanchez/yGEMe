%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ZR = zero_rxns(model,toolbox)
% Finds the % of reactions that cannot carry flux (if all uptakes are
% allowed), using flux variability analysis (FVA).
%
% INPUT:    a GEM model (in COBRA or RAVEN format) and the toolbox name.
% OUTPUT:   The % of reactions that cannot carry flux
%
% Benjamín J. Sánchez. Last edited: 2014-11-28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ZR = zero_rxns(model,toolbox)

[~,n] = size(model.S);
ZR    = ones(n,1);

%Allow uptake for all exchange rxns:
for i = 1:n
    rxn_coeffs = model.S(:,i);
    if length(find(rxn_coeffs)) == 1
        model.lb(i) = -1000;
        model.ub(i) = +1000;
    end
end

%Perform FVA to each reaction:
for i = 1:n
    model.c = zeros(n,1);
    
    %Maximize & minimize the reaction:
    model.c(i) = 1;
    if strcmp(toolbox,'COBRA')
        FBA_max    = optimizeCbModel(model);
        model.c(i) = -1;
        FBA_min    = optimizeCbModel(model);
    elseif strcmp(toolbox,'RAVEN')
        FBA_max    = solveLP(model);
        model.c(i) = -1;
        FBA_min    = solveLP(model);
    end
    
    %Check if the rxn can carry flux:
    if ~isempty(FBA_max.f)
        if abs(FBA_max.f) > 1e-6
            ZR(i) = 0;
        end
    end
    if ~isempty(FBA_min.f)
        if abs(FBA_min.f) > 1e-6
            ZR(i) = 0;
        end
    end
    disp(['Performing FVA: ' num2str(i) '/' num2str(n) ' rxns ready.'])
end

ZR = sum(ZR)/n*100;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%