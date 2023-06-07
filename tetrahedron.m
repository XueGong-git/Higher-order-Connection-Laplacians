% Triangle 4->1, 4->2, 4->3, 1->2, 2->3, 3->1 and 1->2->3
clear
%%% B1
% edges {1, 2}, {1, 3}, {2, 3}, {1, 4}, {2, 4}, {3, 4}
B1 = [-1 -1 0 -1 0 0; 
       1 0 -1 0 -1 0; 
       0 1 1 0 0 -1;
       0 0 0 1 1 1];

%%% B2 {1 2 3}, {1 3 4}, {1 2 4}, {2 3 4}
B2 = [ 1 0 1 0;
       -1 1 0 0;
       1 0 0 1;
       0 -1 -1 0;
       0 0 1 -1;
       0 1 0 1];

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
a = 1; b = 1; c = 1; d = 1;

% Triangle 1->2, 2->3, 3->1 and 1->2->3
x_grid = (2*pi/100) * linspace(1,100)
Tu = zeros(12, 12); Td = zeros(12, 12);
I = [1 0; 0 1]; sx = [0 1; 1 0]; sy = [0 -1i; 1i 0]; sz = [1 0;0 -1]; One = [1 1; 1 1];
for n=1:100
    delta=2*pi*n/100;
    Tu(1:6, 1:6) = [I a*exp(i1*delta)*I a*exp(-i1*delta)*I;
                     a*exp(-i1*delta)*I I a*exp(i1*delta)*I;
                     a*exp(i1*delta)*I a*exp(-i1*delta)*I I]
    
    Tu(1:6, 7:12) = [ b*exp(-i1*delta)*sx d*exp(i1*delta)*sz I;
                       d*exp(i1*delta)*sz I b*exp(-i1*delta)*sx;
                       I b*exp(-i1*delta)*sx d*exp(i1*delta)*sz];
        
    
    Tu(7:12, 1:6) = Tu(1:6, 7:12)';
        
    Tu(7:12, 7:12) = [I c*exp(-i1*delta)*sy c*exp(i1*delta)*sy;
                       c*exp(i1*delta)*sy I c*exp(-i1*delta)*sy;
                       c*exp(-i1*delta)*sy c*exp(i1*delta)*sy, I];
                  
    Lup = kron(D1u, I) - Tu.* kron(A1u, One);
    
    Td(1:6, 1:6) = [I a*exp(-i1*delta)*I a*exp(i1*delta)*I;
                     a*exp(i1*delta)*I I a*exp(-i1*delta)*I;
                     a*exp(-i1*delta)*I a*exp(i1*delta)*I I];
    
    Td(1:6, 7:12) = [ a*exp(i1*delta)*I d*sz I;
                       d*sz I a*exp(-i1*delta)*I;
                       I a*exp(-i1*delta)*I c*sy];
        
    
    Td(7:12, 1:6) = Td(1:6, 7:12)';
        
    Td(7:12, 7:12) = [I c*sy c*sy;
                         c*sy I c*sy;
                         c*sy c*sy, I];
  
    Ldown = kron(D1d, I) - Td.* kron(A1d, One);
                   
    anti_commutator(n)=sum(sum((Lup*Ldown-Ldown*Lup).^2)); % squared frobenius norm of anti-commutator
    lambda_up(n,:)=sort(eig((Lup)));
    lambda_down(n,:)=sort(eig((Ldown)));
    lambda_L(n,:)=sort(eig((Ldown+Lup)));
    hodge_norm1(n) = norm(Lup*Ldown,"fro") % frobenius norm of Lup*Ldown
    hodge_norm2(n) = norm(Ldown*Lup,"fro") % calculate  

end

figure
subplot(2,2,1);
plot(x_grid, anti_commutator)
title('Anticommutator');
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
% Add a label to the x-axis
xlabel('\delta')

subplot(2,2,2);
plot(x_grid, lambda_up(:,1))
hold on 
for l = 2:12
    plot(x_grid, lambda_up(:,l))
end
hold off
title('Spectrum of L_{[1]}^{M, up}');
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
% Add a label to the x-axis
xlabel('\delta')

subplot(2,2,3);
plot(x_grid, lambda_down(:,1))
hold on 
for l = 2:12
    plot(x_grid, lambda_down(:,l))
end
hold off
title('Spectrum of L_{[1]}^{M, down}');
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
% Add a label to the x-axis
xlabel('\delta')

subplot(2,2,4);
plot(x_grid, lambda_L(:,1))
hold on 
for l = 2:12
    plot(x_grid, lambda_L(:,l))
end
hold off
title('Spectrum of L_{[1]}^{M}');
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
% Add a label to the x-axis
xlabel('\delta')
%set(gca, 'FontSize', 14);
saveas(gcf, 'plots\case1.eps', 'epsc');


%figure
%plot(x_grid, hodge_norm1)
%title('Frobenius Norm of Lup*Ldown');

%figure
%plot(x_grid, hodge_norm2)
%title('Frobenius Norm of Ldown*Lup');
