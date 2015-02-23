%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% list = delete_empty(list)
% Recieves a list, removes all empty elements and returns the list with only
% non empty elements. If the list has more than 1 column, the decision will
% be made based on the last column.
%
% INPUT:    Complete list
% OUTPUT:   List without empty spaces
%
% Benjamín J. Sánchez. Last edited: 2014-11-24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function list = delete_empty(list)

[N,M]      = size(list);
delete_pos = zeros(1,N);
for i = 1:N
    if isempty(list{i,M})
        delete_pos(i) = 1;
    end
end
list(find(delete_pos)) = [];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%