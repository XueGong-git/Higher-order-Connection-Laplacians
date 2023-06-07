% Triangle 1->2, 2->3, 3->1 and 1->2->3
clear
%%% B1
% edges {1, 2}, {1, 3}, {2, 3}
B1 = [-1   -1   0 ; 
       1  0   -1 ; 
       0   1   1 ];

%%% B2
B2 = [ 1 ;
       -1 ;
       1 ];

%% Undirected Version

% Hodge Laplacian
L0u = B1 * B1'; % 0-up Laplacian
L1d = B1' * B1; % 1-up Laplacian, n_edge*n_edge matrix
L1u = B2 * B2'; % 1-up Laplacian, n_edge*n_edge matrix
L2d = B2' * B2; % 2-down Laplacian, n_triangle*n_triangle matrix

% The undirected infered adjacency matrix
A0u = diag(diag(L0u)) - L0u; % this is the same as the standard adjacency matrix for standard graph
D1u = diag(diag(L1u));
A1u = D1u - L1u; % infered adjacency between edges
D1d = diag(diag(L1d));
A1d = D1d - L1d; % infered adjacency between edges
i1 = 1*i;

% Triangle 1->2, 2->3, 3->1 and 1->2->3

delta=pi/4;
n = 1
Lup=[1, -exp(i1*delta), exp(-i1*delta);
    -exp(-i1*delta), 1, -exp(i1*delta); 
    exp(i1*delta), -exp(-i1*delta),1];

Ldown=[2,exp(-i1*delta),-exp(i1*delta); 
       exp(i1*delta),2,exp(-i1*delta); 
       -exp(-i1*delta),exp(i1*delta),2];

anti_commutator(n)=sum(sum((Lup*Ldown-Ldown*Lup).^2)); % squared frobenius norm of anti-commutator
lambda_up(n,:)=sort(eig((Lup)));
lambda_down(n,:)=sort(eig((Ldown)));
lambda_L(n,:)=sort(eig((Ldown+Lup)));
hodge_norm1(n) = norm(Lup*Ldown,"fro") % frobenius norm of Lup*Ldown
hodge_norm2(n) = norm(Ldown*Lup,"fro") % calculate  
[Ul, Dl] = eig(Ldown+Lup);
[Uu, Du] = eig(Lup);
[Ud, Dd] = eig(Ldown);


    