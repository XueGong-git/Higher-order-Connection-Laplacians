function L1_up_square = L1_square_up(L1_up, l, L1_up_square)

l1=[l(2), l(3),l(1)];
l2=[l(4),l(3),l(5)];

L1_up_square(l1, l1) = L1_up + L1_up_square(l1, l1);
L1_up_square(l2, l2) = L1_up + L1_up_square(l2, l2);

end

