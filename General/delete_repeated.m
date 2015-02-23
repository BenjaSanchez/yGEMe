%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% list = delete_repeated(list)
% Recieves a list, removes all repeated elements and returns the list with
% only non repeated elements.
%
% NOTE: Must be a string or cell array
%
% Benjamín J. Sánchez. Last edited: 2014-11-12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function list = delete_repeated(list)

N          = length(list);
delete_pos = zeros(1,N);
for i = 1:N-1
    for j = i+1:N
        if strcmp(list{i},list{j})
            delete_pos(j) = 1;
        end
    end
end
list(find(delete_pos)) = [];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%