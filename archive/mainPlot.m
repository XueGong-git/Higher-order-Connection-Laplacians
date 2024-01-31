clear

%% Triangle case 1 
% 1->2, 2->3, 3->1, 1->2->3

% Node list, order does not matter
V = [1 2 3 4];
% Edge list, order represent direction of flow
E = [1 2; 1 3; 2 3; 3 4]; 
% Triangle list, order represent direction of flow
T = [1 2 3]; 

%
[B1, B2] = plotMagneticLaplacian(V, E, T, 96)
L1u = B2*B2'
L1d = B1'*B1
L1 = L1u+L1d 

%% Triangle case 2 1->2->3
% 1->2, 2->3, 3->1 and 3->2->1
% Node list
V = [1 2 3];
% Edge list
E = [1 2; 2 3; 3 1]; 
% Triangle list
T = [3 2 1];
%[B1, B2] = plotMagneticLaplacian(V, E, T, 96);

%% Triangle case 3 1->2->3
% 1->2, 2->3, 1->3 and 1->2->3
% Node list
V = [1 2 3];
% Edge list
E = [1 2; 2 3; 1 3]; 
% Triangle list
T = [1 2 3];

%[B1, B2] = plotMagneticLaplacian(V, E, T, 96)

%% Torus case
% Node list
V = [1 2 3 4 5 6 7 8 9];

% Edge list
E = [2 1; 3 2; 1 3; 
    8 7; 9 8; 7 9;
    5 4; 6 5; 4 6;
    1 4; 4 7; 7 1;
    2 5; 5 8; 8 2;
    3 6; 6 9; 9 3;
    4 2; 5 3; 6 1;
    7 5; 8 6; 9 4;
    1 8; 2 9; 3 7;];
        
% Triangle list
T = [1 4 2; 4 2 5; 2 5 3; 3 6 5; 6 3 1; 1 4 6; 
    5 4 7; 7 5 8; 6 5 8; 8 6 9; 4 6 9; 9 4 7;
    7 1 8; 8 2 1; 2 9 8; 3 2 9; 9 3 7; 3 7 1];

%[B1, B2] = plotMagneticLaplacian(V, E, T, 96);
