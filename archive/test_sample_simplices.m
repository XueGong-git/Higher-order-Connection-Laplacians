%%% 6 nodes, 8 edges, 3 triangles 
% create boundary operators
%%% B1
B1 = [-1 0 -1 0 0 0 0 0; 1 -1 0 -1 0 0 0 0; 0 1 1 0 -1 0 0 0; 0 0 0 1 1 -1 -1 0; 0 0 0 0 0 1 0 -1; 0 0 0 0 0 0 1 1];

%%% B2
B2 = [1 0 0 ; 1 1 0; -1 0 0; 0 -1 0; 0 1 0; 0 0 1; 0 0 -1; 0 0 1];

% Hodge Laplacian
L0u = B1 * B1';
L1d = B1' * B1; %  n_edge*n_edge matrix

L1u = B2 * B2'; % n_edge*n_edge matrix
L2d = B2' * B2; % n_triangle*n_triangle matrix



%solve the eigenvalues
[V0u,D0u] = eig(L0u); % same as normal Laplacian, node 1,2,3 on the left, node 4 in the middle, node 5,6 on the right
[V1u,D1u] = eig(L1u); % zero eigenvalues find connected components of edges through shared triangles
[V1d,D1d] = eig(L1d); % find clusters of eddges by common nodes 
[V2d,D2d] = eig(L2d); % find clusters of triangles by common triangles 