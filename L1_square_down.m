function L1_down_square = L1_square_down(L1_down, l, L1_down_square)

ind = zeros(1, 2*size(l,2));
ind(1:2:size(ind,2)) = l*2-1;
ind(2:2:size(ind,2)) = l*2;

L1_down_square(ind, ind) = L1_down + L1_down_square(ind, ind);

end

