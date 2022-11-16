%%% Test Magnetic Laplacian 
g = 1/3;
%%% 6 nodes, 8 edges, 3 triangles 
% create boundary operators
%%% B1
B1 = [-1 0 -1 0 0 0 0 0; 1 -1 0 -1 0 0 0 0; 0 1 1 0 -1 0 0 0; 0 0 0 1 1 -1 -1 0; 0 0 0 0 0 1 0 -1; 0 0 0 0 0 0 1 1];

%%% B2
B2 = [1 0 0 ; 1 1 0; -1 0 0; 0 -1 0; 0 1 0; 0 0 1; 0 0 -1; 0 0 1];

%%% Omega - assume direction is the same as the alignment

Omega0 = 0.5*[0 1 1 0 0 0; 1 0 1 1 0 0; 1 1 0 1 0 0; 0 1 1 0 0 1; 0 0 0 0 0 1; 0 0 0 1 1 0];
Omega1 = 0.5*[0 1 1 0 0 0 0 0; 1 0 1 1 1 0 0 0; 1 1 0 0 0 0 0 0; 0 1 0 0 1 0 0 0; 0 1 0 1 0 0 0 0; 0 0 0 0 0 0 1 1; 0 0 0 0 0 1 0 1; 0 0 0 0 0 1 1 0];

%%% Phi - assume direction is the same as the alignment

Phi0 = 2*Omega0;
Phi1 = 2*Omega1;

% Hodge Laplacian
L0u = B1 * B1'; % n_node*n_node matrix
L1d = B1' * B1; % n_edge*n_edge matrix

L1u = B2 * B2'; % n_edge*n_edge matrix
L2d = B2' * B2; % n_triangle*n_triangle matrix


%solve the eigenvalues
[V_0,D_0, Phi_0] = magHodge(Omega0,g); % this should be the same as the Magnetic Laplacian
[V_1,D_1, Phi_1] = magHodge(Omega1,g); % find directed clusters of edges according to triangle