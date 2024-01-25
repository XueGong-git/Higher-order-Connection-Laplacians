function [L1_up_square, L1_down_square] = L1_square_b(L1_up, L1_down, l, L1_up_square, L1_down_square)

l1=[l(2), l(3),l(1)];
l2=[l(4),l(3),l(5)];

L1_up_square(l1, l1) = L1_up + L1_up_square(l1, l1);
L1_up_square(l2, l2) = L1_up + L1_up_square(l2, l2);

L1_down_square(l1, l1) = L1_down;
L1_down_square(l2, l2) = L1_down;

end

