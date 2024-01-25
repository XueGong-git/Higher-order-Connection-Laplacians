                      function L1_up_square = L1_square_up_type2(L1_up1, L1_up2, l, L1_up_square)

l1=[2*l(2)-1, 2*l(2), 2*l(3)-1, 2*l(3), 2*l(1)-1, 2*l(1)];

l2=[2*l(4)-1, 2*l(4), 2*l(3)-1, 2*l(3), 2*l(5)-1, 2*l(5)];

L1_up_square(l1, l1) = L1_up1 + L1_up_square(l1, l1);
L1_up_square(l2, l2) = L1_up2 + L1_up_square(l2, l2);

%L1_down_square(l1, l1) = L1_down1;
%L1_down_square(l2, l2) = L1_down2;

end

