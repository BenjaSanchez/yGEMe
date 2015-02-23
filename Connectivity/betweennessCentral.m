%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BC = betweennessCentral(SPA)
% Computes the betweenness centrality of each node of a network, the i-th
% position indicates the fraction of times node i is present in the
% shortest paths between nodes j and k, averaged among all possible pairs
% (j,k).
% 
% INPUT:    Shortest Paths Alternatives Array (SPA)
% OUTPUT:   Betweenness Centrality (BC)
% 
% Benjamín Sánchez. Last edited: 2014-11-20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function BC = betweennessCentral(SPA)

m  = length(SPA);
BC = zeros(m,1);
for i = 1:m
    BCi = zeros(m);
    B   = 0;
    for j = 1:m-1
        for k = j+1:m
            SPAjk = SPA{j,k};
            if i ~= j && i~= k && ~isempty(SPAjk)
                %For each (j,k) pair, search all paths between j and k. If
                %i comes up, store it in paths_i:
                paths_i = 0;
                [R,S]   = size(SPAjk);
                for r = 1:R
                    for s = 1:S
                        if SPAjk(r,s) == i
                            paths_i = paths_i + 1;
                        end
                    end
                end
                %Calculate the betweenness centrality of i in path j->k:
                BCi(j,k) = paths_i/R;
                B        = B+1;
            end
        end        
    end
    %Calculate BC for i, averaging all BCi's:
    BC(i) = sum(sum(BCi))/B;
    disp(['Calculating Betweenness Centrality: Node ' num2str(i) ' ready.']);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%