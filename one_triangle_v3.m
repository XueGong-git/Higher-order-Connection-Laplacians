%% Directed version {1,2,3}, {1,2}, {2,3}, {1,3}
g = 1/3;
%%% 3 nodes, 3 edges, 1 triangle 
% create boundary operators
n_node = 3; n_edge =3; n_triangle = 1;
node_name = [1 2 3]'; node_name = num2str(node_name); node_label = cellstr(node_name);
edge_name = [23 31 12]'; edge_name = num2str(edge_name); edge_label = cellstr(edge_name);
trg_name = [123]'; trg_name = num2str(trg_name); trg_label = cellstr(trg_name);


%%% B0
B1 = [0 1 -1; -1 0 1; 1 -1 0];

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
weight_e = [1 1 1]; % weight of edges
Ph1 = repelem(weight_e,[n_node],[1]);
%%% B2
weight_t = [1]; % weight of triangles
Ph2 = repelem(weight_t,[n_edge],[1]);

%%% Directed boundary operators
B1M = 0.5*sqrt(2).*B1.*exp(-1i*pi*g*(Ph1.*B1));
B2M = 0.5*sqrt(2).*B2.*exp(-1i*pi*g*(Ph2.*B2));

L0Mu = B1M*B1M';
L1Md = B1M' * B1M;
L1Mu = B2M * B2M';

% 0-up Laplacian
[V,D] = eig(L0Mu); % zero eigenvalues find connected components of edges through shared triangles
Psi = angle(V);s = 100.* ones(size(Psi(:,1)));
idx = kmeans([cos(Psi(:,1)), sin(Psi(:,1))], 3);
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
exportgraphics(ax,'plots/triangle_v3_L0Mu.eps','Resolution',300) 

% 1-down Magnetic Laplacian
[V, D] = eig(L1Md); % zero eigenvalues find connected components of edges through shared nodes
Psi = angle(V); s = 100.* ones(size(Psi(:,1)));
idx = kmeans([cos(Psi(:,1:2)), sin(Psi(:,1:2))], 3);
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
exportgraphics(ax,'plots/triangle_v3_L1Md.eps','Resolution',300) 

% 1-up Magnetic Laplacian
[V,D] = eig(L1Mu); % zero eigenvalues find connected components of edges through shared triangles
Psi = angle(V); s = 100.* ones(size(Psi(:,1)));
idx = kmeans([cos(Psi(:,1:3)), sin(Psi(:,1:3))], 3);
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
exportgraphics(ax,'plots/triangle_v3_L1Mu.eps','Resolution',300) 

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
exportgraphics(ax,'plots/triangle_v3_L1Mu.eps','Resolution',300) 

% 1 Laplacian
[V1M,D1M] = eig(L1Mu + L1Md); % zero eigenvalues find connected components of edges through shared triangles
Psi = angle(V1M);
idx = kmeans([cos(Psi(:,1)), sin(Psi(:,1:2))], 3);
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
exportgraphics(ax,'plots/triangle_v3_L1M.eps','Resolution',300)  

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

T =[O sigz_p Ip; sigz_m O sigy_p; Im sigy_m O];
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
exportgraphics(ax,'plots/triangle_v3.jpg','Resolution',300) 

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
exportgraphics(ax,'plots/triangle_v3.jpg','Resolution',300) 
