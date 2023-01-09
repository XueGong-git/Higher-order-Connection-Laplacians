% Simpliical complices with two triangles sharing the same edge, directions of edges are against direction of
% triangles

%%% 6 nodes, 8 edges, 3 triangles
n_node = 4; n_edge =5; n_triangle = 2;
node_name = [1 2 3 4]'; node_name = num2str(node_name); node_label = cellstr(node_name);
edge_name = [12 31 23 42 34]'; edge_name = num2str(edge_name); edge_label = cellstr(edge_name);
trg_name = [132 243]'; trg_name = num2str(trg_name); trg_label = cellstr(trg_name);

% create boundary operators
%%% B1
B1 = [-1 1 0 0 0; 1 0 -1 1 0; 0 -1 1 0 -1; 0 0 0 -1 1];

%%% B2
B2 = [-1 0; -1 0; -1 -1; 0 -1; 0 -1];

%% Undirected Version

% Hodge Laplacian
L0u = B1 * B1'; % 0-up Laplacian
L1d = B1' * B1; % 1-up Laplacian, n_edge*n_edge matrix
L1u = B2 * B2'; % 1-up Laplacian, n_edge*n_edge matrix
L2d = B2' * B2; % 2-down Laplacian, n_triangle*n_triangle matrix


%%% Permutation matrix
P = [0 0 0 1; 0 0 1 0; 0 1 0 0; 1 0 0 0];

% The undirected infered adjacency matrix
A0_up = diag(diag(L0u)) - L0u; % this is the same as the standard adjacency matrix for standard graph
D1_up = diag(diag(L1u));
A1_up = D1_up - L1u; % infered adjacency between edges

% solve the eigenvalues
% 0-up Laplacian
[V0u,D0u] = eig(L0u); % same as normal Laplacian, node 1,2,3 on the left, node 4 in the middle, node 5,6 on the right
scatter(V0u(:,2),V0u(:,3))
x = V0u(:,2); y = V0u(:,3); z = V0u(:,4); s = 100.* ones(size(x));
scatter3(x, y, z, s, 'filled')
set(gca,'XLim',[-1.1 1.1],'YLim',[-1.1 1.1],'ZLim',[-1.1 1.1])
dx = 0.1+0.1*rand(size(x)); dy = 0.1+0.1*rand(size(y)); dz = 0.1+0.1*rand(size(z));% displacement so the text does not overlay the data points
text(x+dx, y+dy, z+dz, node_label, 'FontSize', 20);
xlabel('v1');
ylabel('v2');
zlabel('v3');
set(gca,'fontsize',20);
plt.LineWidth = 2;
ax = gca;
exportgraphics(ax,'plots/two_triangles_v2_L0u.eps','Resolution',300) 

% test symmytry
P * L0u* P' == L0u
L0u*(P*V0u) - (P*V0u)*D0u

% 1-down Laplacian
[V1d,D1d] = eig(L1d); % find clusters of eddges by common nodes 
x = V1d(:,1); y = V1d(:,2); z = V1d(:,3); s = 100.* ones(size(x));
scatter3(x, y, z, s, 'filled')
set(gca,'XLim',[-1.1 1.1],'YLim',[-1.1 1.1],'ZLim',[-1.1 1.1])
name = [12 31 23 24 43]'; name = num2str(name); label = cellstr(name);
dx = 0.05+0.2*rand(size(x)); dy = 0.05+0.2*rand(size(y)); dz = 0.05+0.2*rand(size(z));% displacement so the text does not overlay the data points
text(x+dx, y+dy, z+dz, edge_label, 'FontSize', 20);
xlabel('v1');
ylabel('v2');
zlabel('v3');
set(gca,'fontsize',20);
plt.LineWidth = 2;
ax = gca;
exportgraphics(ax,'plots/two_triangles_v2_L1d.eps','Resolution',300) 

% 1-up Laplacian
[V1u,D1u] = eig(L1u); % zero eigenvalues find connected components of edges through shared triangles
x = V1u(:,1); y = V1u(:,2); z = V1u(:,3); s = 100.* ones(size(x));
scatter3(x, y, z, s, 'filled')
set(gca,'XLim',[-1.1 1.1],'YLim',[-1.1 1.1],'ZLim',[-1.1 1.1])
name = [12 31 23 24 43]'; name = num2str(name); label = cellstr(name);
dx = 0.05+0.2*rand(size(x)); dy = 0.05+0.2*rand(size(y)); dz = 0.05+0.2*rand(size(z));% displacement so the text does not overlay the data points
text(x+dx, y+dy, z+dz, edge_label, 'FontSize', 20);
xlabel('v1');
ylabel('v2');
zlabel('v3');
set(gca,'fontsize',20);
plt.LineWidth = 2;
ax = gca;
exportgraphics(ax,'plots/two_triangles_v2_L1u.eps','Resolution',300) 

% 1 Laplacian
[V1,D1] = eig(L1u+L1d); % zero eigenvalues find connected components of edges through shared triangles
x = V1(:,1); y = V1(:,2); z = V1(:,3); s = 100.* ones(size(x));
scatter3(x, y, z, s, 'filled')
set(gca,'XLim',[-1.1 1.1],'YLim',[-1.1 1.1],'ZLim',[-1.1 1.1])
name = [12 31 23 24 43]'; name = num2str(name); label = cellstr(name);
dx = 0.1; dy = 0.1; dz = 0.1;% displacement so the text does not overlay the data points
text(x+dx, y+dy, z+dz, edge_label, 'FontSize', 20);
xlabel('v1');
ylabel('v2');
zlabel('v3');
set(gca,'fontsize',20);
plt.LineWidth = 2;
ax = gca;
exportgraphics(ax,'plots/two_triangles_v2_L1.eps','Resolution',300) 

% 2-down Laplacian
[V2d,D2d] = eig(L2d); % find clusters of triangles by common edges 

%% Directed version {1,2,3}, {1,2}, {2,3}, {3,1}
g = 1/3;

%% First version in Ginestra's notes Section B2

% weight of flow
%%% B1
weight_e = [1 1 1 1 1];
Ph1 = repelem(weight_e,[n_node],[1]);
%%% B2
weight_t = [1 1];
Ph2 = repelem(weight_t,[n_edge],[1]);

%%% Directed adjacency matrix
B1M = 0.5*sqrt(2).*B1.*exp(-1i*pi*g*(Ph1.*B1));
B2M = 0.5*sqrt(2).*B2.*exp(-1i*pi*g*(Ph2.*B2));

L0Mu = B1M*B1M';
L1Md = B1M' * B1M;
L1Mu = B2M * B2M';

% 0-up Laplacian
[V,D] = eig(L0Mu); % zero eigenvalues find connected components of edges through shared triangles
Psi = angle(V);s = 100.* ones(size(Psi(:,1)));
scatter(cos(Psi(:,1)),sin(Psi(:,1)), s, 'filled')
set(gca,'XLim',[-1.1 1.1],'YLim',[-1.1 1.1])
name = [1 2 3 4 5 6]'; name = num2str(name); label = cellstr(name);
dx = 0.1+0.2*rand(n_node,1); dy = 0.1+0.2*rand(n_node,1);% displacement so the text does not overlay the data points
xlabel('cos(\theta)');
ylabel('sin(\theta)');
text(cos(Psi(:,1))+dx, sin(Psi(:,1))+dy, node_label, 'FontSize', 20);
set(gca,'fontsize',20);
plt.LineWidth = 2;
ax = gca;
exportgraphics(ax,'plots/two_triangles_v2_L0Mu.eps','Resolution',300) 

% 1-down Magnetic Laplacian
[V, D] = eig(L1Md); % zero eigenvalues find connected components of edges through shared nodes
Psi = angle(V); s = 100.* ones(size(Psi(:,1)));
scatter(cos(Psi(:,1)),sin(Psi(:,1)), s, 'filled')
set(gca,'XLim',[-1.1 1.1],'YLim',[-1.1 1.1])
name = [12 31 23 24 43 45 64 65]'; name = num2str(name); label = cellstr(name);
xlabel('cos(\theta)');
ylabel('sin(\theta)');
dx = -0.05+0.3*rand(n_edge,1); dy = -0.05+0.3*rand(n_edge,1);% displacement so the text does not overlay the data points
text(cos(Psi(:,1))+dx, sin(Psi(:,1))+dy, edge_label, 'FontSize', 20);
set(gca,'fontsize',20);
plt.LineWidth = 2;
ax = gca;
exportgraphics(ax,'plots/two_triangles_v2_L1Md.eps','Resolution',300) 

% 1-up Magnetic Laplacian
[V,D] = eig(L1Mu); % zero eigenvalues find connected components of edges through shared triangles
Psi = angle(V); s = 100.* ones(size(Psi(:,1)));
scatter(cos(Psi(:,1)),sin(Psi(:,1)), s, 'filled')
set(gca,'XLim',[-1.1 1.1],'YLim',[-1.1 1.1])
xlabel('cos(\theta)');
ylabel('sin(\theta)');
dx = -0.05+0.5*rand(n_edge,1); dy = -0.05+0.5*rand(n_edge,1);% displacement so the text does not overlay the data points
text(cos(Psi(:,1))+dx, sin(Psi(:,1))+dy, edge_label, 'FontSize', 20);
set(gca,'fontsize',20);
plt.LineWidth = 2;
ax = gca;
exportgraphics(ax,'plots/two_triangles_v2_L1Mu.eps','Resolution',300) 

% plot top 2 eigenvectors
psi1 = Psi(:, 1)';
azi = Psi(:, 1)'; %azimuth for the spherical coordinates
ele = Psi(:, 2)'; % elevation for the sperical coordinates
[x,y,z] = sph2cart(azi,ele,1);
scatter3(x,y,z, s, 'filled')
xlabel('x');
ylabel('y');
zlabel('z');
set(gca,'XLim',[-1.1 1.1],'YLim',[-1.1 1.1],'ZLim',[-1.1 1.1])
name = [12 31 23 42 34 45 64 65]'; name = num2str(name); label = cellstr(name);
dx = -0.05+0.3*rand(n_edge,1); dy = -0.05+0.3*rand(n_edge,1); dz = -0.05+0.3*rand(n_edge,1);% displacement so the text does not overlay the data points
text(x+dx', y+dy', z+dz', edge_label, 'FontSize', 20);
set(gca,'fontsize',20);
plt.LineWidth = 2;
ax = gca;
exportgraphics(ax,'plots/two_triangles_v2_L1Mu.eps','Resolution',300) 

% 1 Laplacian
[V1M,D1M] = eig(L1Mu + L1Md); % zero eigenvalues find connected components of edges through shared triangles
Psi = angle(V1M);
scatter(cos(Psi(:,1)),sin(Psi(:,1)), s, 'filled')
set(gca,'XLim',[-1.1 1.1],'YLim',[-1.1 1.1])
name = [12 31 23 42 34 45 64 65]'; name = num2str(name); label = cellstr(name);
xlabel('cos(\theta)');
ylabel('sin(\theta)');
dx = 0.1; dy = 0.1;% displacement so the text does not overlay the data points
text(cos(Psi(:,1))+dx, sin(Psi(:,1))+dy, edge_label, 'FontSize', 20);
set(gca,'fontsize',20);
plt.LineWidth = 2;
ax = gca;
exportgraphics(ax,'plots/two_triangles_v2_L1M.eps','Resolution',300) 

