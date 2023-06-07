%  3 nodes {1}, {2}, {3},
%  3 edges {1,2}, {1,3}, {2,3}
%  1 triangle {1,2,3}
% create boundary operators
%%% B1
B1 = [-1 -1 0; 1 0 -1; 0 1 1];

%%% B2
B2 = [1; -1; 1];

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

%% Directed version {1,2,3}, {1,2}, {2,3}, {3,1}
g = 1/6;

%% New definition of L1Mp with Pauli matrices
delta = 2*pi*g;
One = ones(2);
I2 = eye(2);
Im = exp(-delta*i)*I2;
Ip = exp(delta*i)*I2;
O = zeros(2,2);

I2 = eye(2); Im = exp(-delta*i)*I2; Ip = exp(delta*i)*I2;
Sx = [0, 1; 1, 0]; Sxm = exp(-delta*i)*Sx; Sxp = exp(delta*i)*Sx;
Sy = [0, -1i; 1i, 0]; Sym = exp(-delta*i)*Sx; Syp = exp(delta*i)*Sx;
Sz = [1, 0; 0, -1]; Szm = exp(-delta*i)*Sx; Szp = exp(delta*i)*Sx;

T =[I2 Syp Im; Sym I2 Szp; Ip Szm I2];
Lm1u = kron(D1_up, I2) - T .* kron(A1_up, One); % matrix product is not Hermitian

% 1-up Laplacian
[Vm1u,Dm1u] = eig(Lm1u); % zero eigenvalues find connected components of edges through shared triangles
scatter(1:6, angle(Vm1u(:,1)))

%phase = [ pi/3; pi/3; pi/6; pi/6; 0;0];
%(exp(i*phase))'*Lm1u*exp(i*phase)

% substitute the eigenvector into the quadratic form
b = (Vm1u(:,1));
a = b';
L=Lm1u;
a*L*b

for m = 1:6
    for n = 1:6
        m;
        n;
        P(m,n)= abs(A1_up(ceil(m/2),ceil(n/2)))*(b(m)-A1_up(ceil(m/2),ceil(n/2))*T(m,n)*b(n));
        
    end
end

% substitute the normalized  eigenvector
phi = angle(Vm1u(:,1));
b = exp(1i*phi);
a = b';
L=Lm1u;
a*L*b;

for m = 1:6
    for n = 1:6
        m;
        n;
        P(m,n)= abs(A1_up(ceil(m/2),ceil(n/2)))*(b(m)-A1_up(ceil(m/2),ceil(n/2))*T(m,n)*b(n));
        Sum(m,n) = a(m)*Lm1u(m,n)*b(n);
    end
end