%%% 3 nodes, 3 edges, 1 triangle 
% create boundary operators
%%% B1
B1 = [0 1 -1; -1 0 1; 1 -1 0];

%%% B2
B2 = [1; 1; 1];

%% Undirected Version

% Hodge Laplacian
L0u = B1 * B1'; % 0-up Laplacian
L1d = B1' * B1; % 1-up Laplacian, n_edge*n_edge matrix
L1u = B2 * B2'; % 1-up Laplacian, n_edge*n_edge matrix
L2d = B2' * B2; % 2-down Laplacian, n_triangle*n_triangle matrix

% The undirected infered adjacency matrix
A0_up = diag(diag(L0u)) - L0u; % this is the same as the standard adjacency matrix for standard graph
D1_up = diag(diag(L1u));
A1_up = D1_up - L1u; % infered adjacency between edges

% solve the eigenvalues
% 0-up Laplacian
[V0u,D0u] = eig(L0u); % same as normal Laplacian, node 1,2,3 on the left, node 4 in the middle, node 5,6 on the right
scatter(V0u(:,2),V0u(:,3))
x = V0u(:,1); y = V0u(:,2); z = V0u(:,3); s = 100.* ones(size(x));
scatter3(x, y, z, s, 'filled')
set(gca,'XLim',[-1.1 1.1],'YLim',[-1.1 1.1],'ZLim',[-1.1 1.1])
name = [1 2 3]'; name = num2str(name); label = cellstr(name);
dx = 0.1; dy = 0.1; dz = 0.1;% displacement so the text does not overlay the data points
text(x+dx, y+dy, z+dz, label, 'FontSize', 20);
xlabel('v1');
ylabel('v2');
zlabel('v3');
set(gca,'fontsize',20);
plt.LineWidth = 2;
ax = gca;
exportgraphics(ax,'plots/triangle1_L0u.eps','Resolution',300) 

% 1-down Laplacian
[V1d,D1d] = eig(L1d); % find clusters of eddges by common nodes 
x = V1d(:,1); y = V1d(:,2); z = V1d(:,3); s = 100.* ones(size(x));
scatter3(x, y, z, s, 'filled')
set(gca,'XLim',[-1.1 1.1],'YLim',[-1.1 1.1],'ZLim',[-1.1 1.1])
label = ["[2,3]", "[3,1]" "[1,2]"];
dx = 0.1; dy = 0.1; dz = 0.1;% displacement so the text does not overlay the data points
text(x+dx, y+dy, z+dz, label, 'FontSize', 20);
xlabel('v1');
ylabel('v2');
zlabel('v3');
set(gca,'fontsize',20);
plt.LineWidth = 2;
ax = gca;
exportgraphics(ax,'plots/triangle1_L1d.eps','Resolution',300) 

% 1-up Laplacian
[V1u,D1u] = eig(L1u); % zero eigenvalues find connected components of edges through shared triangles
x = V1u(:,1); y = V1u(:,2); z = V1u(:,3); s = 100.* ones(size(x));
scatter3(x, y, z, s, 'filled')
set(gca,'XLim',[-1.1 1.1],'YLim',[-1.1 1.1],'ZLim',[-1.1 1.1])
label = ["[2,3]", "[3,1]" "[1,2]"];
dx = 0.1; dy = 0.1; dz = 0.1;% displacement so the text does not overlay the data points
text(x+dx, y+dy, z+dz, label, 'FontSize', 20);
xlabel('v1');
ylabel('v2');
zlabel('v3');
set(gca,'fontsize',20);
plt.LineWidth = 2;
ax = gca;
exportgraphics(ax,'plots/triangle1_L1u.eps','Resolution',300) 

% 1 Laplacian
[V1,D1] = eig(L1u+L1d); % zero eigenvalues find connected components of edges through shared triangles
x = V1(:,1); y = V1(:,2); z = V1(:,3); s = 100.* ones(size(x));
scatter3(x, y, z, s, 'filled')
set(gca,'XLim',[-1.1 1.1],'YLim',[-1.1 1.1],'ZLim',[-1.1 1.1])
label = ["[2,3]", "[3,1]" "[1,2]"];
dx = 0.1; dy = 0.1; dz = 0.1;% displacement so the text does not overlay the data points
text(x+dx, y+dy, z+dz, label, 'FontSize', 20);
xlabel('v1');
ylabel('v2');
zlabel('v3');
set(gca,'fontsize',20);
plt.LineWidth = 2;
ax = gca;
exportgraphics(ax,'plots/triangle1_L1.eps','Resolution',300) 



% 2-down Laplacian
[V2d,D2d] = eig(L2d); % find clusters of triangles by common edges 

%% Directed version {1,2,3}, {1,2}, {2,3}, {3,1}
g = 1/3;

%% First version in Ginestra's notes Section B2
%%% Directed adjacency matrix
A0 = [0 1 0; 0 0 1; 1 0 0]; % node-node
A1 = [0 1 0; 0 0 1; 1 0 0]; % edge-edge

Omega0 = 0.5*[0 1 1; 1 0 1; 1 1 0]; % node-node
Omega1 = 0.5*[0 1 1; 1 0 1; 1 1 0]; % edge-edge

%%% Phi - assume direction is the same as the alignment
Phi0 = 2*Omega0;
Phi1 = 2*Omega1;

%solve the eigenvalues
%[V_0,D_0, Phi_0] = maghodge(Omega0,g); % this should be the same as the Magnetic Laplacian
%phase = meigenmaps(A1,g);

% plot edges
[V_1, Lg] = maghodge(Omega1, A1, g); % edge
Psi = angle(V_1);
scatter(cos(Psi(:,1)),sin(Psi(:,1)), s, 'filled')
set(gca,'XLim',[-1.1 1.1],'YLim',[-1.1 1.1])
label = ["[2,3]", "[3,1]" "[1,2]"];
dx = 0.1; dy = 0.1;% displacement so the text does not overlay the data points
text(cos(Psi(:,1))+dx, sin(Psi(:,1))+dy, label, 'FontSize', 20);
set(gca,'fontsize',20);
plt.LineWidth = 2;
ax = gca;
exportgraphics(ax,'plots/triangle1_v1.eps','Resolution',300) 

% compare with the same definition but calculated with boundary operator
% manual calculation
%B1M = 0.5*sqrt(2)* [0, exp(-1i*pi*g), -exp(1i*pi*g); -exp(1i*pi*g), 0, exp(-1i*pi*g); exp(-1i*pi*g), -exp(1i*pi*g), 0];
%B2M = 0.5*sqrt(2)* [exp(-1i*pi*g); exp(-1i*pi*g); exp(-1i*pi*g)];

B1M = 0.5*sqrt(2).*B1.*exp(-1i*pi*g*B1);
B2M = 0.5*sqrt(2).*B2.*exp(-1i*pi*g*B2);


L0Mu = B1M*B1M';
angle(-L0Mu)/(pi*g)
L1Md = B1M' * B1M;
angle(-L1Md)/(pi*g)
L1Mu = B2M * B2M';
angle(-L1Mu)/(pi*g)

% 0-up Laplacian
[V,D] = eig(L0Mu); % zero eigenvalues find connected components of edges through shared triangles
Psi = angle(V);
scatter(cos(Psi(:,1)),sin(Psi(:,1)), s, 'filled')
set(gca,'XLim',[-1.1 1.1],'YLim',[-1.1 1.1])
name = [1 2 3]'; name = num2str(name); label = cellstr(name);
dx = 0.1; dy = 0.1;% displacement so the text does not overlay the data points
xlabel('cos(\theta)');
ylabel('sin(\theta)');
text(cos(Psi(:,1))+dx, sin(Psi(:,1))+dy, label, 'FontSize', 20);
set(gca,'fontsize',20);
plt.LineWidth = 2;
ax = gca;
exportgraphics(ax,'plots/triangle1_L0Mu.eps','Resolution',300) 

% 1-down Magnetic Laplacian
[V, D] = eig(L1Md); % zero eigenvalues find connected components of edges through shared triangles
Psi = angle(V);
scatter(cos(Psi(:,1)),sin(Psi(:,1)), s, 'filled')
set(gca,'XLim',[-1.1 1.1],'YLim',[-1.1 1.1])
label = ["[2,3]", "[3,1]" "[1,2]"];
xlabel('cos(\theta)');
ylabel('sin(\theta)');
dx = 0.1; dy = 0.1;% displacement so the text does not overlay the data points
text(cos(Psi(:,1))+dx, sin(Psi(:,1))+dy, label, 'FontSize', 20);
set(gca,'fontsize',20);
plt.LineWidth = 2;
ax = gca;
exportgraphics(ax,'plots/triangle1_L1Md.eps','Resolution',300) 

% 1-up Magnetic Laplacian
[V,D] = eig(L1Mu); % zero eigenvalues find connected components of edges through shared triangles
Psi = angle(V);
psi1 = Psi(:, 1)';
azi = Psi(:, 1)'; %azimuth for the spherical coordinates
ele = Psi(:, 2)'; % elevation for the sperical coordinates
[x,y,z] = sph2cart(azi,ele,1);
scatter3(x,y,z, s, 'filled')
xlabel('x');
ylabel('y');
zlabel('z');
set(gca,'XLim',[-1.1 1.1],'YLim',[-1.1 1.1],'ZLim',[-1.1 1.1])
label = ["[2,3]", "[3,1]" "[1,2]"];
dx = 0.1; dy = 0.1; dz = 0.1;% displacement so the text does not overlay the data points
text(x+dx, y+dy, z+dz, label, 'FontSize', 20);
set(gca,'fontsize',20);
plt.LineWidth = 2;
ax = gca;
exportgraphics(ax,'plots/triangle1_L1Mu.eps','Resolution',300) 

[V,D] = eig(L1Mu); % zero eigenvalues find connected components of edges through shared triangles
x = V(:,1); y = V(:,2); z = V(:,3); s = 100.* ones(size(x));
scatter3(x, y, z, s, 'filled')
set(gca,'XLim',[-1.1 1.1],'YLim',[-1.1 1.1],'ZLim',[-1.1 1.1])
name = [1 2 3]'; name = num2str(name); label = cellstr(name);
dx = 0.1; dy = 0.1; dz = 0.1;% displacement so the text does not overlay the data points
text(x+dx, y+dy, z+dz, label, 'FontSize', 20);
xlabel('v1');
ylabel('v2');
zlabel('v3');
set(gca,'fontsize',20);
plt.LineWidth = 2;
ax = gca;
exportgraphics(ax,'plots/triangle1_L1Mu.eps','Resolution',300) 



% 1 Laplacian
[V1M,D1M] = eig(L1Mu + L1Md); % zero eigenvalues find connected components of edges through shared triangles
Psi = angle(V1M);
scatter(cos(Psi(:,1)),sin(Psi(:,1)), s, 'filled')
set(gca,'XLim',[-1.1 1.1],'YLim',[-1.1 1.1])
label = ["[2,3]", "[3,1]" "[1,2]"];
xlabel('cos(\theta)');
ylabel('sin(\theta)');
dx = 0.1; dy = 0.1;% displacement so the text does not overlay the data points
text(cos(Psi(:,1))+dx, sin(Psi(:,1))+dy, label, 'FontSize', 20);
set(gca,'fontsize',20);
plt.LineWidth = 2;
ax = gca;
exportgraphics(ax,'plots/triangle1_L1M.eps','Resolution',300) 

%% Second version Ginestra's notes Section III
delta = 2*pi*g;
I2 = ones(2);
Im = exp(-delta*i)*I2;
Ip = exp(delta*i)*I2;
O = zeros(2,2);

T =[O Im Ip; Ip O Ip; Im Im O];
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
label = ["[2,3]", "[3,1]" "[1,2]"];
dx = 0.1; dy = 0.1; dz = 0.1;% displacement so the text does not overlay the data points
text(x+dx, y+dy, z+dz, label);
set(gca,'fontsize',20);
plt.LineWidth = 2;
ax = gca;
exportgraphics(ax,'plots/triangle1_v2.eps','Resolution',300) 

%% Third version -- only consider direction of triangle, treat edges as undirected
%solve the eigenvalues
phase = meigenmaps(A1,g);
scatter(cos(phase),sin(phase))
set(gca,'XLim',[-1.1 1.1],'YLim',[-1.1 1.1])
label = ["[2,3]", "[3,1]" "[1,2]"];
dx = 0.1; dy = 0.1;% displacement so the text does not overlay the data points
text(cos(phase)+dx, sin(phase)+dy, label);
set(gca,'fontsize',20);
plt.LineWidth = 2;
ax = gca;
exportgraphics(ax,'plots/triangle1_v3.eps','Resolution',300) 


