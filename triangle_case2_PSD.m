% Triangle 1->2, 2->3, 3->1 and 3->2->1
clear
close all
%%% B1
% edges {1, 2}, {1, 3}, {2, 3}
B1 = [-1   0   1 ; 
       1  -1   0 ; 
       0   1   -1 ];

%%% B2 
B2 = [ -1 ;
       -1 ;
       -1 ];

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

%%
% Triangle slight frustration 1->2, 2->3, 3->1 and 3->2->1
Lup=zeros(6,6);
Ldown=zeros(6,6);

i1 = 1i;
for n=1:100
    
    delta=(2*pi/100)*n;
    Lup(1,:)=[2, 0,                0, -exp(-i1*delta),   0, exp(i1*delta)];
    Lup(2,:)=[0, 2,                -exp(-i1*delta), 0,   exp(i1*delta), 0];
    
    Lup(3,:)=[0, -exp(i1*delta),    2,0,                0,-exp(-i1*delta)];
    Lup(4,:)=[-exp(i1*delta), 0,    0,2,                -exp(-i1*delta),0];
    
    Lup(5,:)=[0, exp(-i1*delta),  0,-exp(i1*delta),   2,0];
    Lup(6,:)=[exp(-i1*delta), 0,  -exp(i1*delta),0,   0,2];
    
    
    Ldown(1,:)=[2,0,               exp(-i1*delta),0,   -exp(i1*delta),0];
    Ldown(2,:)=[0,2,               0,exp(-i1*delta),   0,-exp(i1*delta)];
    
    Ldown(3,:)=[exp(i1*delta),0,  2,0,                 exp(-i1*delta),0];
    Ldown(4,:)=[0,exp(i1*delta),  0,2,                 0,exp(-i1*delta)];
    
    Ldown(5,:)=[-exp(-i1*delta),0,  exp(i1*delta),0,     2,0];
    Ldown(6,:)=[0,-exp(-i1*delta),  0,exp(i1*delta),     0,2];
    
    anti_commutator(n)=sum(sum((Lup*Ldown-Ldown*Lup).^2));
    lambda_up(n,:)=sort(eig((Lup)));
    lambda_down(n,:)=sort(eig((Ldown)));
    lambda_L(n,:)=sort(eig((Ldown+Lup)));
    hodge_norm1(n) = norm(Lup*Ldown,"fro"); % calculate  
    hodge_norm2(n) = norm(Ldown*Lup,"fro"); % calculate  

end
x_grid = (2*pi/100) * linspace(1,100);

figure

subplot(2,2,1);
plot(x_grid, lambda_up(:,1))
hold on 
plot(x_grid, lambda_up(:,2))
plot(x_grid, lambda_up(:,3))
plot(x_grid, lambda_up(:,4))
plot(x_grid, lambda_up(:,5))
plot(x_grid, lambda_up(:,6))
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
plot(x_grid, lambda_down(:,2))
plot(x_grid, lambda_down(:,3))
plot(x_grid, lambda_down(:,4))
plot(x_grid, lambda_down(:,5))
plot(x_grid, lambda_down(:,6))
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
plot(x_grid, lambda_L(:,2))
plot(x_grid, lambda_L(:,3))
plot(x_grid, lambda_L(:,4))
plot(x_grid, lambda_L(:,5))
plot(x_grid, lambda_L(:,6))
%plot(0.7227, 0, '*')
%plot(x_grid, 3+4cos()
hold off
title('Spectrum of $L_{1}^{c}$', 'Interpreter', 'latex');
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
ylim([-0.1 8.1])
% Add a label to the x-axis
xlabel('\delta')

subplot(2,2,4);
plot(x_grid, anti_commutator)
title('Anticommutator', 'Interpreter', 'latex');
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

%set(gca, 'FontSize', 14);
saveas(gcf, 'plots\case2.eps', 'epsc');

figure
plot(x_grid, hodge_norm1)
title('Frobenius Norm of Lup*Ldown');

figure
plot(x_grid, hodge_norm2)
title('Frobenius Norm of Ldown*Lup');
