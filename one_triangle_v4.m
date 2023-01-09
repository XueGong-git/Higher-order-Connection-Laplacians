%% Directed version {1,2,3}, {1,2}, {2,3}, {1,3}
g = 1/3;
%%% 3 nodes, 3 edges, 1 triangle 
% create boundary operators
%%% B1
B1 = [0 -1 -1; -1 0 1; 1 1 0];

%%% B2
B2 = [1; -1; 1];

%% Undirected Version

% Hodge Laplacian
L0u = B1 * B1'; % 0-up Laplacian
L1d = B1' * B1; % 1-up Laplacian, n_edge*n_edge matrix
L1u = B2 * B2'; % 1-down Laplacian, n_edge*n_edge matrix
L2d = B2' * B2; % 2-down Laplacian, n_triangle*n_triangle matrix

% The undirected infered adjacency matrix
A0_up = diag(diag(L0u)) - L0u; % this is the same as the standard adjacency matrix for standard graph
D1_up = diag(diag(L1u));
A1_up = D1_up - L1u; % infered adjacency between edges

%% First version in Ginestra's notes Section B2
%%% Directed adjacency matrix
A0 = [0 1 1; 0 0 0; 0 1 0]; % node-node
A1 = [0 1 0; 0 0 1; 1 0 0]; % edge-edge

Omega0 = 0.5*[0 1 1; 1 0 1; 1 1 0]; % total flow between nodes
Omega1 = 0.5*[0 1 1; 1 0 1; 1 1 0]; % total flow between edges

%%% Phi - assume direction is the same as the alignment
Phi0 = 2*Omega0;
Phi1 = 2*Omega1;

% plot edges
[V_1,D_1, Phi_1] = maghodge(Omega1,A1, g); % edge
Psi = angle(V_1);
scatter(cos(Psi(:,1)),sin(Psi(:,1)))
set(gca,'XLim',[-1.1 1.1],'YLim',[-1.1 1.1])
edge_name = [23 13 12]'; edge_name = num2str(edge_name); label = cellstr(edge_name);
dx = 0.1; dy = 0.1;% displacement so the text does not overlay the data points
text(cos(Psi(:,1))+dx, sin(Psi(:,1))+dy, label);
set(gca,'fontsize',20);
plt.LineWidth = 2;
ax = gca;
exportgraphics(ax,'plots/triangle4_v1.jpg','Resolution',300) 

%% Second version Ginestra's notes Section III
delta = 2*pi*g;
I2 = ones(2);

sigx = [0 1; 1 0];
sigy = [0 -1i; 1i 0];
sigz = [1 0; 0 -1];
O = zeros(2,2);

Im = exp(-delta*i)*I2;
Ip = exp(delta*i)*I2;
sigy_m = exp(-delta*i)*sigy;
sigy_p = exp(delta*i)*sigy;
sigz_m = exp(-delta*i)*sigz;
sigz_p = exp(delta*i)*sigz;

T =[O sigz_p sigz_m; sigz_m O sigy_m; sigz_p sigy_p O];
Lm1u = kron(D1_up, I2) - T .* kron(A1_up, I2); % matrix product is not Hermitian

% 1-up Laplacian
[Vm1u,Dm1u] = eig(Lm1u); % zero eigenvalues find connected components of edges through shared triangles
scatter(Vm1u(:,2),Vm1u(:,3))
Lm1u == Lm1u'
Psi = angle(Vm1u);
psi1 = Psi(:, 1)';
theta = psi1(1:2:end);
phi = psi1(2:2:end);
scatter(cos(theta),sin(theta));
scatter(cos(phi),sin(phi));
[x,y,z] = sph2cart(theta,phi,1);
scatter3(x,y,z)
set(gca,'XLim',[-1.1 1.1],'YLim',[-1.1 1.1],'ZLim',[-1.1 1.1])
edge_name = [23 13 12]'; edge_name = num2str(edge_name); label = cellstr(edge_name);
dx = 0.1; dy = 0.1; dz = 0.1;% displacement so the text does not overlay the data points
text(x+dx, y+dy, z+dz, label);
set(gca,'fontsize',20);
plt.LineWidth = 2;
ax = gca;
exportgraphics(ax,'plots/triangle4_v2.jpg','Resolution',300) 

%% Third version -- only consider direction of triangle, treat edges as undirected
%solve the eigenvalues
phase = meigenmaps(A1,g);
scatter(cos(phase),sin(phase))
set(gca,'XLim',[-1.1 1.1],'YLim',[-1.1 1.1])
edge_name = [23 13 12]'; edge_name = num2str(edge_name); label = cellstr(edge_name);
dx = 0.1; dy = 0.1;% displacement so the text does not overlay the data points
text(cos(phase)+dx, sin(phase)+dy, label);
set(gca,'fontsize',20);
plt.LineWidth = 2;
ax = gca;
exportgraphics(ax,'plots/triangle4_v3.jpg','Resolution',300) 
