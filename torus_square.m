% Triangle 4->1, 4->2, 4->3, 1->2, 2->3, 3->1 and 1->2->3
clear
%%% B1
% edges {1, 2}, {1, 3}, {2, 3}, {1, 4}, {2, 4}, {3, 4}
l = [1 10 19 4 13; 
    2 13 20 5 16;
    3 16 21 6 10;
    4 11 22 7 14;
    5 14 23 8 17;
    6 17 24 9 11;
    7 12 25 1 15;
    8 15 26 2 18;
    9 18 27 3 12;];
i1 = 1i;

L1_up_square = zeros(27, 27); L1_down_square = zeros(27, 27);

max_step = 96;
x_grid = (2*pi/max_step) * linspace(1,max_step, max_step);

for n=1:max_step

delta=2*pi*n/max_step;

L1_up=[1, -exp(i1*delta), exp(-i1*delta);
        -exp(-i1*delta), 1, -exp(i1*delta); 
        exp(i1*delta), -exp(-i1*delta),1];
    
L1_down=[2,exp(-i1*delta),-exp(i1*delta); 
           exp(i1*delta),2,exp(-i1*delta); 
           -exp(-i1*delta),exp(i1*delta),2];
       
for row = 1:size(l,1)
    [L1_up_square, L1_down_square] = L1_square(L1_up, L1_down, l(row,:), L1_up_square, L1_down_square);
    
end


Lup =  L1_up_square; Ldown = L1_down_square;
anti_commutator(n) = sum(sum((Lup*Ldown-Ldown*Lup).^2)); % squared frobenius norm of anti-commutator
[V, D] = eigs(Lup, 27, 'smallestreal');
V_up(n,:,:) = reshape(V, [1, 27, 27]); lambda_up(n,:) = diag(D);
[V, D] = eigs(Ldown, 27, 'smallestreal');
V_down(n,:,:) = reshape(V, [1, 27, 27]); lambda_down(n,:) = diag(D);
[V, D] = eigs(Ldown+Lup, 27, 'smallestreal');
V_L(n,:,:) = reshape(V, [1, 27, 27]); lambda_L(n,:) = diag(D);
hodge_norm1(n) = norm(Lup*Ldown,"fro"); % Frobenius norm of Lup*Ldown
hodge_norm2(n) = norm(Ldown*Lup,"fro"); % calculate   

end


figure
subplot(2,2,1);
plot(x_grid, lambda_up(:,1))
hold on 
for ind = 2:27
    plot(x_grid, lambda_up(:,ind))
end
%plot(x_grid, 1+2*cos(x_grid), '-.')
hold off
title('Spectrum of L_{[1]}^{M, up}');
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
% Add a label to the x-axis
xlabel('\delta')

subplot(2,2,2);
plot(x_grid, lambda_down(:,1))
hold on 
for ind = 2:27
    plot(x_grid, lambda_down(:,ind))
end
plot(x_grid, 2+2*cos(x_grid-pi/3), '-.')
hold off
title('Spectrum of L_{[1]}^{M, down}');
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
% Add a label to the x-axis
xlabel('\delta')

subplot(2,2,3);
plot(x_grid, lambda_L(:,1))
hold on 
for ind = 2:27
    plot(x_grid, lambda_L(:,ind))
end
plot(x_grid, 3+2*sqrt(3)*sin(x_grid), '-.')
hold off
title('Spectrum of L_{[1]}^{M}');
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
% Add a label to the x-axis
xlabel('\delta')
%set(gca, 'FontSize', 14);
saveas(gcf, 'plots\case1.eps', 'epsc');

subplot(2,2,4);
plot(x_grid, anti_commutator)
title('Commutator');
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
% Add a label to the x-axis
xlabel('\delta')
ylim([-1 1])
saveas(gcf, 'plots\case1.eps', 'epsc');


figure
plot(x_grid, hodge_norm1)
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
title('Frobenius Norm of Lup*Ldown');

figure
plot(x_grid, hodge_norm2)
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
title('Frobenius Norm of Ldown*Lup');


