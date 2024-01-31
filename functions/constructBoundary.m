function [B1, B2] = constructBoundary(V, E, T)

%% Construct Boundary Operator B1


B1 = zeros(size(V,2), size(E,1));
for j = 1:size(V,2)
    for k = 1: size(E, 1)
        if V(j) == E(k, 1)
            B1(j,k) = -1; % node j is the head of edge j
        elseif V(j) == E(k, 2)
            B1(j,k) = 1;
        end
    end
end



%% Construct Boundary Operator B2
B2 = zeros(size(E,1), size(T,1));
for m = 1:size(E,1)
    for n = 1: size(T, 1)
        if ismember(E(m, :), T(n, :))
            [~, indiciesOfMissingElements] = find(ismember(T(n, :), E(m, :))==0);
            temp = T(n, :); temp(indiciesOfMissingElements) = [];
            v1 = find(ismember(T(n, :), E(m,1)) == 1);
            v2 = find(ismember(T(n, :), E(m,2)) == 1);
            if mod(v2-v1,3) == 1
                B2(m,n) = 1;
            else 
                B2(m,n) = -1;
            end
        end
    end
end


end
