% Define the nodes
nodes = [1, 1, 1;  % Node 1
         1, -1, 1;  % Node 2
         -1, -1, 1;  % Node 3
         -1, 1, 1;  % Node 4
         1, 0, 0;  % Node 5
         0, -1, 0;  % Node 6
         -1, 0, 0;  % Node 7
         0, 1, 0;  % Node 8
         0, 0, 1];  % Node 9

% Define the triangles
triangles = [1, 5, 9;  % Triangle 1
             5, 2, 6;  % Triangle 2
             6, 3, 7;  % Triangle 3
             7, 4, 8;  % Triangle 4
             8, 1, 9];  % Triangle 5

% Create a triangulation
T = triangulation(triangles, nodes);

% Plot the triangulation
figure
trimesh(T)
axis equal