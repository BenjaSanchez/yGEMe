%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [d,p] = dijsktra(M,start)
% Uses the Dijsktra algorithm to find the shortest path from a starting
% point to every other point in a graph.
% 
% INPUT:  M (matrix representing the graph), and start (the starting point)
% OUTPUT: d (shortest distance to each point) and p (previous points, all
%            possibilities)
% 
% Benjamín Sánchez. Last edited: 2014-11-20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [d,p] = dijsktra(M,start)

m = length(M);
S = zeros(m,1);     %Visited nodes
d = inf(1,m);       %Distances from start
p = zeros(1,m);     %Previous
A = ones(1,m);      %Alternatives

d(start)   = 0;
p(1,start) = start;

for i = 1:m
    
    %Find closest:
    dmin = inf;
    for j = 1:m
        if d(j) < dmin && isempty(find(S==j))
            S(i) = j;
            dmin = d(j);
        end
    end
    
    if S(i) ~= 0
        %Find neighbors:
        nb = find(M(S(i),:));
        %Take out the ones in S:
        for j = 1:length(nb)
            if ~isempty(find(S==nb(j)));
                nb(j) = 0;
            end
        end

        %Update distances and previous:
        for j = 1:length(nb)
            k = nb(j);
            if k > 0
                new_d = d(S(i)) + M(S(i),k);
                
                %If a shorter path exists, replace it as the only one so
                %far:
                if new_d < d(k)
                    d(k)   = new_d;
                    p(1,k) = S(i);
                    if A(k) > 1
                        p(2:A(k),k) = zeros(2:A(k),1);
                        A(k)   = 1;
                    end
                %If a equally short path exists, then add the alternative:
                elseif new_d == d(k)
                    A(k)      = A(k)+1;
                    p(A(k),k) = S(i);
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%