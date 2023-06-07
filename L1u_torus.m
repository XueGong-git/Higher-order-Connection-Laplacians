%  7 nodes
%  21 edges 
%  14 triangle 
% create boundary operators
%%% B1
B1 = zeros(7, 21);
B1(1, 1:6) = 1;
B1(2, 7:11) = 1; B1(2, 1) = -1;
B1(3, 12:15) = 1; B1(3, 2) = -1;  B1(3, 7) = -1;  
B1(4, 16:18) = 1; B1(4, 3) = -1;  B1(4, 8) = -1;  B1(4, 12) = -1;  
B1(5, 19:20) = 1; B1(5, 4) = -1;  B1(5, 9) = -1;  B1(5, 13) = -1;  B1(5, 16) = -1;  
B1(6, 21) = 1; B1(6, 5) = -1;  B1(6, 10) = -1;  B1(6, 14) = -1;  B1(6, 17) = -1;  B1(6, 19) = -1;  
B1(7, 21) = -1; B1(7, 6) = -1;  B1(7, 11) = -1;  B1(7, 15) = -1;  B1(7, 18) = -1;  B1(7, 20) = -1;  


%%% B2
B2 = zeros(21, 14);
B2(1, 1:2) = 1;
B2(2, 3:4) = 1;
B2(3, [1, 3]) = -1; 
B2(4, 5:6) = 1;
B2(5, [2, 5]) = -1;
B2(6, [4, 6]) = -1;
B2(7, [7,8]) = -1;
B2(8, [1,9]) = -1;
B2(9, [7,9]) = -1;
B2(10, [2,11]) = 1;
B2(11, [8,11]) = -1;
B2(12, [3,12]) = 1;
B2(13, [7,13]) = 1;
B2(14, [12,13]) = 1;
B2(15, [4,8]) = 1;
B2(16, [9,14]) = [1 -1];
B2(17, [10, 12]) = 1;
B2(18, [10,14]) = -1;
B2(19, [5,13]) = 1;
B2(20, [6,14]) = 1;
B2(21, [10,11]) = 1;

%% Undirected Version


% Hodge Laplacian
L0u = B1 * B1'; % 0-up Laplacian
L1d = B1' * B1; % 1-up Laplacian, n_edge*n_edge matrix
L1u = B2 * B2'; % 1-up Laplacian, n_edge*n_edge matrix
L2d = B2' * B2; % 2-down Laplacian, n_triangle*n_triangle matrix

% The undirected infered adjacency matrix
A0_up = diag(diag(L0u)) - L0u; % this is the same as the standard adjacency matrix for standard graph
D1_up = diag(diag(L1u));
A1_up = diag(diag(L1u)) - L1u; % infered adjacency between edges

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
