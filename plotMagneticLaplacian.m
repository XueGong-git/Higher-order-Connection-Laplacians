function [B1, B2] = plotMagneticLaplacian(V, E, T, max_step)


%nNode =  size(V,2); 
nEdge = size(E,1);
%nTriag = size(T,1);
[B1, B2] = constructBoundary(V, E, T);



% Hodge Laplacian
%L0u = B1 * B1'; % 0-up Laplacian
L1d = B1' * B1; % 1-up Laplacian, n_edge*n_edge matrix
L1u = B2 * B2'; % 1-up Laplacian, n_edge*n_edge matrix
%L2d = B2' * B2; % 2-down Laplacian, n_triangle*n_triangle matrix


% solve the eigenvalues

% 1-down Laplacian
%[V1d,D1d] = eig(L1d); % find clusters of eddges by common nodes 

% 1-up Laplacian
%[V1u,D1u] = eig(L1u); % zero eigenvalues find connected components of edges through shared triangles

% 1 Laplacian
%[V1,D1] = eig(L1u+L1d); % zero eigenvalues find connected components of edges through shared triangles

% 2-down Laplacian
%[V2d,D2d] = eig(L2d); % find clusters of triangles by common edges 




%% Magnetic Laplacian
%max_step = 96;
x_grid = (2*pi/max_step) * linspace(1,max_step, max_step);
commutator = zeros(1, max_step); hodge_norm1 = zeros(1, max_step); hodge_norm2 = zeros(1, max_step); 
V_up=zeros(max_step, nEdge*2, nEdge*2); V_down=zeros(max_step, nEdge*2, nEdge*2); V_L = zeros(max_step, nEdge*2, nEdge*2);
lambda_up = zeros(max_step, nEdge*2); lambda_down = zeros(max_step, nEdge*2); 
lambda_L = zeros(max_step, nEdge*2);

for iter=1:max_step
    
delta=2*pi*iter/max_step;
[Lup, Ldown] = constructMagneticLaplacian(V, E, T, delta, L1u, L1d);
commutator(iter) = sum(sum((Lup*Ldown-Ldown*Lup).^2)); % squared frobenius norm of anti-commutator

[Vec, D] = eigs(Lup, nEdge*2, 'smallestreal');
V_up(iter,:,:) = reshape(Vec, [1, nEdge*2, nEdge*2]); 
lambda_up(iter,:) = diag(D);

[Vec, D] = eigs(Ldown, nEdge*2, 'smallestreal');
V_down(iter,:,:) = reshape(Vec, [1, nEdge*2, nEdge*2]);
lambda_down(iter,:) = diag(D);

[Vec, D] = eigs(Ldown+Lup, nEdge*2, 'smallestreal');
V_L(iter,:,:) = reshape(Vec, [1, nEdge*2, nEdge*2]); 
lambda_L(iter,:) = diag(D);

hodge_norm1(iter) = norm(Lup * Ldown,"fro"); % Frobenius norm of Lup*Ldown
hodge_norm2(iter) = norm(Ldown * Lup,"fro"); % calculate  

end


%}



figure
subplot(2,2,1);
plot(x_grid, lambda_up(:,1))
hold on 
for ind = 2:nEdge*2
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
for ind = 2:nEdge*2
    plot(x_grid, lambda_down(:,ind))
end
%plot(x_grid, 2+2*cos(x_grid-pi/3), '-.')
hold off
title('Spectrum of L_{[1]}^{M, down}');
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
% Add a label to the x-axis
xlabel('\delta')

subplot(2,2,3);
plot(x_grid, lambda_L(:,1))
hold on 
for ind = 2:nEdge*2
    plot(x_grid, lambda_L(:,ind))
end
%plot(x_grid, 3+2*sqrt(3)*sin(x_grid), '-.')
hold off
title('Spectrum of L_{[1]}^{M}');
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
% Add a label to the x-axis
xlabel('\delta')
%set(gca, 'FontSize', 14);
saveas(gcf, 'plots\case1.eps', 'epsc');

subplot(2,2,4);
plot(x_grid, commutator)
title('Commutator');
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
% Add a label to the x-axis
xlabel('\delta')
%ylim([-1 1])
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

end