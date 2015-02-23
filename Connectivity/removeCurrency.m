%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mncm = removeCurrency(model,Mfull)
% Removes currency metabolites from connectivity matrix.
%
% INPUT:    a GEM model (in COBRA or RAVEN format) and the full
%           connectivity matrix.
% OUTPUT:   Connectivity matrix with no currency metabolites: Water,
%           proton, CO2, O2, phosphate3-, diphosphate4-, ammonium, ATP,
%           ADP, AMP, NAD+, NADH, NADP+, NADPH.
%
% Benjamín J. Sánchez. Last edited: 2014-12-09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M = removeCurrency(model,M)

currency_mets = {'H2O','H+','carbon dioxide','oxygen','phosphate', ...
                 'diphosphate','ammonium','ATP','ADP','AMP','NAD', ...
                 'NAD(+)','NADH','NADP(+)','NADPH'};
             
del_pos = zeros(1,length(M));
for i = 1:length(currency_mets)
    %Find all rows/columns of M associated with a currency metabolite (in
    %any compartment):
    met_pos = strcmp(model.metNames,currency_mets{i})';
    if sum(met_pos) == 0
        disp(['WARNING: could not remove ' currency_mets{i} ' from M.'])
    else
        del_pos = del_pos + met_pos;
    end
end

%Remove empty rows and/or columns:
for i = 1:length(M)
    column = M(:,i);
    row    = M(i,:);
    if sum(column.^2) == 0 || sum(row.^2) == 0
        met_pos    = zeros(1,length(M));
        met_pos(i) = 1;
        del_pos    = del_pos + met_pos;
    end
end

%Remove rows & columns from the connectivity matrix:
M(find(del_pos),:) = [];
M(:,find(del_pos)) = [];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%