%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% results = evalConnectivity(model)
% Performs connectivity tests to the model.
%
% INPUT:    a GEM model (in COBRA or RAVEN format) and the toolbox name.
% OUTPUT:   a structure with the following fields:
%           *Mfull -> Full Connectivity matrix analysis
%           *Mncm  -> Connectivity matrix analysis with no currency
%                     metabolites (Herrgard et al. 2008)
%           *UM    -> % unreachable metabolites, 1x1
%           *ZR    -> % reactions that cannot carry flux, 1x1
%           The first 2 fields have the following subfields:
%           *M   (connectivity matrix, mxm)
%           *GCC (global clustering coefficient, 1x1)
%           *LCC (local clustering coefficient, mx1)
%           *ACC (average clustering coefficient, 1x1)
%           *ND  (node degree, mx1)
%           *AND (average node degree, 1x1)
%           *SPM (shortest path matrix, mxm)
%           *SPD (shortest paths diversity, mxm)
%           *SPA (shortest paths alternatives, mxm)
%           *CPL (characteristic path length, 1x1)
%           *D   (network diameter, 1x1)
%           *APD (average path diversity, 1x1)
%           *HPD (highest path diversity, 1x1)
%           *BC  (betweenness centrality, mx1)
%           *ABC (average betweenness centrality, 1x1)
%           *HBC (highest betweenness centrality, 1x1)
%
% Benjamín J. Sánchez. Last edited: 2014-12-04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = evalConnectivity(model,toolbox)

%Compute connectivity matrix:
Mfull = connectivityMatrix(full(model.S));

%Remove currency metabolites from connectivity matrix to form the second
%connectivity matrix:
Mncm = removeCurrency(model,Mfull);

%Perform all analysis for both connectivity matrices:
for i = 1:2
    if i == 1
        M = Mfull;
    else
        M = Mncm;
    end
    
    m     = length(M);
    res.M = M;
    
    %Compute global, local and average clustering coefficients:
    res.GCC = globalClustCoeff(M);
    res.LCC = localClustCoeff(M);
    res.ACC = mean(res.LCC);

    %The i-th node degree is simply how much links the i-th metabolite has:
    res.ND  = sum(M)';
    res.AND = mean(res.ND);

    %Obtain all shortest paths between all pair of nodes in the network:
    [SPM,SPD,SPA] = shortestPaths(M);
    res.SPM       = SPM;
    res.SPD       = SPD;
    res.SPA       = SPA;

    %Compute characteristic path length, diameter and average path diversity:
    res.CPL = sum(sum(SPM))/(m*(m-1));
    res.D   = max(max(SPM));
    res.APD = sum(sum(SPD))/(m*(m-1)/2);
    res.HPD = max(max(SPD));

    %Calculate beyweenness centrality for each node, and find the average and
    %highest values:
    res.BC  = betweennessCentral(SPA);
    res.ABC = mean(res.BC);
    res.HBC = max(res.BC);

    if i == 1
        results.Mfull = res;
    else
        results.Mncm = res;
    end
end

%Find % unreachable metabolites:
results.UM = unreach_mets(model,results.Mfull.SPM);

%Find % reactions that cannot carry flux:
results.ZR = zero_rxns(model,toolbox);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%