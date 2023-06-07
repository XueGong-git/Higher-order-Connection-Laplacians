function [L1_up_square, L1_down_square] = L1_square(L1_up, L1_down, l, L1_up_square, L1_down_square)

L1_up_square(l(1:3), l(1:3)) = L1_up + L1_up_square(l(1:3), l(1:3));
L1_up_square(l(3:5), l(3:5)) = L1_up + L1_up_square(l(3:5), l(3:5));

L1_down_square(l(1:3), l(1:3)) = L1_down;
L1_down_square(l(3:5), l(3:5)) = L1_down;

end

