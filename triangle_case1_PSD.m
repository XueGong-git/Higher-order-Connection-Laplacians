% Triangle 1->2, 2->3, 3->1 and 1->2->3
clear
close all
%%% B1
% edges {1, 2}, {2, 3}, {3, 1}
B1 = [-1   0   1 ; 
       1  -1  0 ; 
       0   1   -1 ];

%%% B2
B2 = [ 1 ;
       1 ;
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
D1u_mod = sum(abs(A1u), 2); % modified 1-up degree matrix

D1d = diag(diag(L1d));
A1d = D1d - L1d; % infered adjacency between edges
D1d_mod = sum(abs(A1d), 2); % modified 1-down degree matrix

i1 = 1i;

max_step = 96;
% Triangle 1->2, 2->3, 3->1 and 1->2->3
x_grid = (2*pi/max_step) * linspace(1,max_step, max_step);
%delta = 2*pi/4;

for n=1:max_step
    
   delta=2*pi*n/max_step;
    Lup=[2, -exp(i1*delta), exp(-i1*delta);
        -exp(-i1*delta), 2, -exp(i1*delta); 
        exp(i1*delta), -exp(-i1*delta),2];
    
    Ldown=[2,exp(-i1*delta),-exp(i1*delta); 
           exp(i1*delta),2,exp(-i1*delta); 
           -exp(-i1*delta),exp(i1*delta),2];
    
    %
    Lup = kron(Lup, eye(2)); Ldown = kron(Ldown, eye(2));
    anti_commutator(n) = sum(sum((Lup*Ldown-Ldown*Lup).^2)); % squared frobenius norm of anti-commutator
    
    sizeMat = size(Lup,1);
    
    [V, D] = eigs(Lup, sizeMat, 'smallestreal');
    V_up(n,:,:) = reshape(V, [1, sizeMat, sizeMat]); lambda_up(n,:) = diag(D);
    
    [V, D] = eigs(Ldown, sizeMat, 'smallestreal');
    V_down(n,:,:) = reshape(V, [1, sizeMat, sizeMat]); lambda_down(n,:) = diag(D);
    
    [V, D] = eigs(Ldown+Lup, sizeMat, 'smallestreal');
    V_L(n,:,:) = reshape(V, [1, sizeMat, sizeMat]); lambda_L(n,:) = diag(D);
    hodge_norm1(n) = norm(Lup*Ldown,"fro"); % Frobenius norm of Lup*Ldown
    hodge_norm2(n) = norm(Ldown*Lup,"fro"); % calculate  

end


figure
subplot(2,2,1);
plot(x_grid, lambda_up(:,1))
hold on 
for ind = 2:sizeMat
    plot(x_grid, lambda_up(:,ind))
end
plot(x_grid, 2+2*cos(x_grid), '-.')
hold off
title('Spectrum of $L_{1}^{c, up}$', 'Interpreter', 'latex');
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
ylim([-0.1 4.1])

% Add a label to the x-axis
xlabel('\delta')

subplot(2,2,2);
plot(x_grid, lambda_down(:,1))
hold on 
for ind = 2:sizeMat
    plot(x_grid, lambda_down(:,ind))
end
plot(x_grid, 2+2*cos(x_grid-pi/3), '-.')
hold off
title('Spectrum of $L_{1}^{c, down}$', 'Interpreter', 'latex');
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
ylim([-0.1 4.1])
% Add a label to the x-axis
xlabel('\delta')

subplot(2,2,3);
plot(x_grid, lambda_L(:,1))
hold on 
for ind = 2:sizeMat
    plot(x_grid, lambda_L(:,ind))
end
plot(x_grid, 4+2*sqrt(3)*sin(x_grid), '-.')
hold off
title('Spectrum of $L_{1}^{c}$', 'Interpreter', 'latex');
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
% Add a label to the x-axis
xlabel('\delta')
%set(gca, 'FontSize', 14);
saveas(gcf, 'plots\case1.eps', 'epsc');

subplot(2,2,4);
plot(x_grid, anti_commutator)
title('Commutator', 'Interpreter', 'latex');
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
% Add a label to the x-axis
xlabel('\delta')
ylim([-1 1])

lines = findall(gcf, 'Type', 'line');
% Loop through each line and set color to blue
for i = 1:numel(lines)
    lines(i).Color =  [0, 0.4470, 0.7410];
end

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

