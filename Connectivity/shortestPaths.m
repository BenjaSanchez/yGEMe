%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [SPM,SPD,SPA] = shortestPaths(M)
% Computes the shortest path matrix and the corresponding paths.
% 
% INPUT:    *Connectivity matrix (M)
% OUTPUTS:  *Shortest Path Matrix (SPM): The (i,j) position contains the
%            shortest distance between node i and node j.
%           *Shortest Path Diversity (SPD): The (i,j) position contains the
%            number of alternatives there are for going from i to j in
%            SPM(i,j) steps.
%           *Shortest Path Alternatives (SPA): The {i,j} position contains
%            all SPD(i,j) paths, as a matrix in which the rows are the
%            different paths from i to j.
% 
% Benjamín Sánchez. Last edited: 2014-12-09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [SPM,SPD,SPA] = shortestPaths(M)

%Perform the dijsktra algorithm for all nodes:
m   = length(M);
SPM = zeros(m,m);
P   = cell(m,1);
for i = 1:m
    [SPM(i,:),P{i}] = dijsktra(M,i);
    for j = 1:m
        if SPM(i,j) == inf
            SPM(i,j) = 0;
        end
    end
    disp(['computing SPM: ' num2str(i) '/' num2str(m) ' nodes ready'])
end

%Find all shortest paths between each pair of nodes. Store them in SPA, and
%store the number of shortest paths in SPD:
SPA = cell(m);
SPD = zeros(m);
for i = 1:m
    for j = i+1:m
        assignin('base','paths',[]);
        find_path(i,j,P{i})
        SPA{i,j}     = evalin('base','paths');
        [SPD(i,j),~] = size(SPA{i,j});
    end
    disp(['Generating shortest paths: Path ' num2str(i) ' ready.']);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function find_path(start,new_path,Pi)
%Recursive function that generates all shortest paths between 2 nodes.

first = new_path(1);
if first == start
    %Shortest path completed. Store it in a workspace variable:
    paths = [evalin('base','paths');new_path];
    assignin('base','paths',paths);
else
    %Add a previous step to the path. Iterate the function for all
    %possibilities:
    prev = Pi(:,first);
    for i = 1:length(prev)
        if Pi(i,first) > 0
            find_path(start,[Pi(i,first) new_path],Pi)
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%