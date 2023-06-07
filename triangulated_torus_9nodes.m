clear
%  9 nodes
%  27 edges 
%  18 triangle 
% create boundary operators
%%% B1


% Node list
V = [1 2 3 4 5 6 7 8 9];

% Edge list
E = [2 1; 1 3; 1 4; 5 1; 7 1; 8 1;
     3 2; 4 2; 2 6; 7 2; 9 2;
     3 5; 6 3; 3 8; 9 3;
     4 5; 6 4; 4 8; 9 4;
     5 9; 5 6; 7 5;
     6 7; 8 6;
     7 8; 9 7;
     8 9];
        
% Triangle list
T = [1 4 8; 1 3 8; 4 5 9; 4 8 9; 5 1 3; 5 9 3;
    3 8 6; 3 2 6; 8 9 7; 8 6 7; 9 3 2; 9 7 2; 
    2 6 4; 2 1 4; 6 7 5; 6 4 5; 7 2 1; 7 5 1;];


nNode =  size(V,2); nEdge = size(E,1); nTriag = size(T,1);
%% Construct Boundary Operator B1
B1 = zeros(size(V,2), size(E,1));
for j = 1:size(V,2)
    for k = 1: size(E, 1)
        if V(j) == E(k, 1)
            B1(j,k) = -1; % node j is the head of edge j
        elseif V(j) == E(k, 2)
            B1(j,k) = 1;
        end
    end
end



%% Construct Boundary Operator B2
B2 = zeros(size(E,1), size(T,1));
for m = 1:size(E,1)
    for n = 1: size(T, 1)
        if ismember(E(m, :), T(n, :))
            [dummy, indiciesOfMissingElements] = find(ismember(T(n, :), E(m, :))==0);
            temp = T(n, :); temp(indiciesOfMissingElements) = [];
            if temp == E(m, :)
                B2(m,n) = 1;
            else 
                B2(m,n) = -1;
            end
        end
    end
end


%% 



% Hodge Laplacian
L0u = B1 * B1'; % 0-up Laplacian
L1d = B1' * B1; % 1-up Laplacian, n_edge*n_edge matrix
L1u = B2 * B2'; % 1-up Laplacian, n_edge*n_edge matrix
L2d = B2' * B2; % 2-down Laplacian, n_triangle*n_triangle matrix

% The undirected infered adjacency matrix
A0_up = diag(diag(L0u)) - L0u; % this is the same as the standard adjacency matrix for standard graph

D1_up = diag(diag(L1u));
A1_up = diag(diag(L1u)) - L1u; % infered adjacency between edges

D1_down = diag(diag(L1d));
A1_down = diag(diag(L1d)) - L1d; % this is the same as the standard adjacency matrix for standard graph

% solve the eigenvalues
% 0-up Laplacian
[V0u,D0u] = eig(L0u); % same as normal Laplacian

% 1-down Laplacian
[V1d,D1d] = eig(L1d); % find clusters of eddges by common nodes 

% 1-up Laplacian
[V1u,D1u] = eig(L1u); % zero eigenvalues find connected components of edges through shared triangles

% 1 Laplacian
[V1,D1] = eig(L1u+L1d); % zero eigenvalues find connected components of edges through shared triangles

% 2-down Laplacian
[V2d,D2d] = eig(L2d); % find clusters of triangles by common edges 




%% Directed Version
max_step = 96;
x_grid = (2*pi/max_step) * linspace(1,max_step, max_step);

for iter=1:max_step
    delta=2*pi*iter/max_step;
% Up Laplacian
% Calculate rotation matrix
T_up = zeros(size(E,1), size(E,1));
Theta_up = zeros(size(E,1), size(E,1));
P_up = zeros(4*size(E,1), 4*size(E,1));
sx = [0 1; 1 0]; sy = [0 -1i; 1i 0]; sz = [1 0; 0 -1];

for n = 1:size(T,1) % triangle
    % find edges
    edgeInd = find(ismember(E(:,2), T(n, :)) == 1 & ismember(E(:,1), T(n, :)) == 1);
    for x = 1:3 % edge
        e1 = edgeInd(x);
        for y = 1:3
            e2 = edgeInd(y);
            if e1 ~= e2
                v1 = find(ismember(T(n, :), E(e1,1)) == 1);
                v2 = find(ismember(T(n, :), E(e1,2)) == 1);
                if E(e1,2) == E(e2,1) % tail of e1 == head of e2
                    % find index of common nodes in triangle
                    
                    if  mod(v2-v1,3) == 1 % case 1
                        Theta_up(e1, e2) = - delta; % ijk; jki; kij
                        P_up(2*e1-1: 2*e1, 2*e2-1: 2*e2) = eye(2);
                        T_up(2*e1-1: 2*e1, 2*e2-1: 2*e2) = exp(-1i*delta)*eye(2);
                    else % case 4
                        Theta_up(e1, e2) = delta; 
                        P_up(2*e1-1: 2*e1, 2*e2-1: 2*e2) = sx;
                        T_up(2*e1-1: 2*e1, 2*e2-1: 2*e2) = exp(1i*delta)*sx;

                    end
                    
                elseif E(e1,1) == E(e2,2) % head of e1 == tail of e2
                    if  mod(v2-v1,3) == 2 % case 3: kji; jik; ikj
                        Theta_up(e1, e2) =  -delta; 
                        P_up(2*e1-1: 2*e1, 2*e2-1: 2*e2) = eye(2);
                        T_up(2*e1-1: 2*e1, 2*e2-1: 2*e2) = exp(-1i*delta)*eye(2);

                    else % case 3 : ijk; jki; kij
                        Theta_up(e1, e2) = delta; 
                        P_up(2*e1-1: 2*e1, 2*e2-1: 2*e2) = sx;
                        T_up(2*e1-1: 2*e1, 2*e2-1: 2*e2) = exp(1i*delta)*sx;

                    end
                    
                elseif E(e1,1) == E(e2,1) % head of e1 == head of e2
                    if  mod(v2-v1,3) == 1 % case 5 ijk; jki; kij
                        Theta_up(e1, e2) = - delta; % ijk; jki; kij
                        P_up(2*e1-1: 2*e1, 2*e2-1: 2*e2) = sy;
                        T_up(2*e1-1: 2*e1, 2*e2-1: 2*e2) = exp(-1i*delta)*sy;

                    else % case 6
                        Theta_up(e1, e2) = delta; 
                        P_up(2*e1-1: 2*e1, 2*e2-1: 2*e2) = sy;
                        T_up(2*e1-1: 2*e1, 2*e2-1: 2*e2) = exp(1i*delta)*sy;

                    end
                elseif E(e1,2) == E(e2,2) % tail of e1 == tail of e2
                    if  mod(v2-v1,3) == 1 % case 7 ijk; jki; kij
                        Theta_up(e1, e2) = - delta; % ijk; jki; kij
                        P_up(2*e1-1: 2*e1, 2*e2-1: 2*e2) = sz;
                        T_up(2*e1-1: 2*e1, 2*e2-1: 2*e2) = exp(-1i*delta)*sz;

                    else % case 8
                        Theta_up(e1, e2) = delta; 
                        P_up(2*e1-1: 2*e1, 2*e2-1: 2*e2) = sz;
                        T_up(2*e1-1: 2*e1, 2*e2-1: 2*e2) = exp(1i*delta)*sz;

                    end
                end
            end              
        end
    end
end

Lup = kron(D1_up, eye(2)) - T_up.* kron(A1_up, ones(2));




%% Down Laplacian

for n = 1:nNode % triangle
    % find edges
    [edgeInd c] = find(ismember(E(), V(n)) == 1 );
    for x = 1:size(edgeInd,1) % edge
        e1 = edgeInd(x);
        for y = 1:size(edgeInd,1)
            e2 = edgeInd(y);
            if e1 ~= e2
                v1 = find(ismember(E(e1,:), V(n)) == 1);
                v2 = find(ismember(E(e2,:), V(n)) == 1);
                if v1 == 1 & v2 == 1 % j-i, j-k
                        Theta_down(e1, e2) = 0; 
                        P_down(2*e1-1: 2*e1, 2*e2-1: 2*e2) = sy;
                        T_down(2*e1-1: 2*e1, 2*e2-1: 2*e2) = sy;
                elseif v1 == 1 & v2 == 2  % case 2
                        Theta_down(e1, e2) = -delta; 
                        P_down(2*e1-1: 2*e1, 2*e2-1: 2*e2) = eye(2);
                        T_down(2*e1-1: 2*e1, 2*e2-1: 2*e2) = exp(-1i*delta)*eye(2);
                elseif v1 == 2 & v2 == 1  % case 1
                        Theta_down(e1, e2) = delta; 
                        P_down(2*e1-1: 2*e1, 2*e2-1: 2*e2) = eye(2);
                        T_down(2*e1-1: 2*e1, 2*e2-1: 2*e2) = exp(1i*delta)*eye(2);
                elseif v1 == 2 & v2 == 2  % case 4
                        Theta_down(e1, e2) = 0; 
                        P_down(2*e1-1: 2*e1, 2*e2-1: 2*e2) = sz;
                        T_down(2*e1-1: 2*e1, 2*e2-1: 2*e2) = sz;

                end
            end
                    
               

        end
    end
end              

Ldown = kron(D1_down, eye(2)) - T_down.* kron(A1_down, ones(2));


anti_commutator(iter) = sum(sum((Lup*Ldown-Ldown*Lup).^2)); % squared frobenius norm of anti-commutator
[Vec, D] = eigs(Lup, 54, 'smallestreal');
V_up(n,:,:) = reshape(Vec, [1, 54, 54]); 
lambda_up(iter,:) = diag(D);
[Vec, D] = eigs(Ldown, 54, 'smallestreal');
V_down(iter,:,:) = reshape(Vec, [1, 54, 54]); lambda_down(iter,:) = diag(D);
[Vec, D] = eigs(Ldown+Lup, 54, 'smallestreal');
V_L(iter,:,:) = reshape(Vec, [1, 54, 54]); lambda_L(iter,:) = diag(D);
hodge_norm1(iter) = norm(Lup*Ldown,"fro"); % Frobenius norm of Lup*Ldown
hodge_norm2(iter) = norm(Ldown*Lup,"fro"); % calculate  

end


%}



figure
subplot(2,2,1);
plot(x_grid, lambda_up(:,1))
hold on 
for ind = 2:54
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
for ind = 2:54
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
for ind = 2:54
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
plot(x_grid, anti_commutator)
title('Anticommutator');
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
