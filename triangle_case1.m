% Triangle 1->2, 2->3, 3->1 and 1->2->3
clear
%%% B1
% edges {1, 2}, {1, 3}, {2, 3}
B1 = [-1   -1   0 ; 
       1  0   -1 ; 
       0   1   1 ];

%%% B2
B2 = [ 1 ;
       -1 ;
       1 ];

%% Undirected Version

% Hodge Laplacian
L0u = B1 * B1'; % 0-up Laplacian
L1d = B1' * B1; % 1-up Laplacian, n_edge*n_edge matrix
L1u = B2 * B2'; % 1-up Laplacian, n_edge*n_edge matrix
L2d = B2' * B2; % 2-down Laplacian, n_triangle*n_triangle matrix

% The undirected infered adjacency matrix
A0u = diag(diag(L0u)) - L0u; % this is the same as the standard adjacency matrix for standard graph
D1u = diag(diag(L1u));
A1u = D1u - L1u; % infered adjacency between edges
D1d = diag(diag(L1d));
A1d = D1d - L1d; % infered adjacency between edges
i1 = 1*i;

max_step = 96;
% Triangle 1->2, 2->3, 3->1 and 1->2->3
x_grid = (2*pi/max_step) * linspace(1,max_step, max_step);


for n=1:max_step
    
    delta=2*pi*n/max_step;
    Lup=[1, -exp(i1*delta), exp(-i1*delta);
        -exp(-i1*delta), 1, -exp(i1*delta); 
        exp(i1*delta), -exp(-i1*delta),1];
    
    Ldown=[2,exp(-i1*delta),-exp(i1*delta); 
           exp(i1*delta),2,exp(-i1*delta); 
           -exp(-i1*delta),exp(i1*delta),2];
    
    anti_commutator(n) = sum(sum((Lup*Ldown-Ldown*Lup).^2)); % squared frobenius norm of anti-commutator
    [V, D] = eigs(Lup, 3, 'smallestreal');
    V_up(n,:,:) = reshape(V, [1, 3, 3]); lambda_up(n,:) = diag(D);
    [V, D] = eigs(Ldown, 3, 'smallestreal');
    V_down(n,:,:) = reshape(V, [1, 3, 3]); lambda_down(n,:) = diag(D);
    [V, D] = eigs(Ldown+Lup, 3, 'smallestreal');
    V_L(n,:,:) = reshape(V, [1, 3, 3]); lambda_L(n,:) = diag(D);
    hodge_norm1(n) = norm(Lup*Ldown,"fro"); % Frobenius norm of Lup*Ldown
    hodge_norm2(n) = norm(Ldown*Lup,"fro"); % calculate  

end


figure
subplot(2,2,1);
plot(x_grid, lambda_up(:,1))
hold on 
plot(x_grid, lambda_up(:,2))
plot(x_grid, lambda_up(:,3))
plot(x_grid, 1+2*cos(x_grid), '-.')
hold off
title('Spectrum of L_{[1]}^{M, up}');
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
% Add a label to the x-axis
xlabel('\delta')

subplot(2,2,2);
plot(x_grid, lambda_down(:,1))
hold on 
plot(x_grid, lambda_down(:,2))
plot(x_grid, lambda_down(:,3))
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
plot(x_grid, lambda_L(:,2))
plot(x_grid, lambda_L(:,3))
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
