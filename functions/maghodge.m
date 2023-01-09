function [V, Lg] = maghodge(Omega,A0,g)
%G plot the nodes against the eigenvalue
%eig is the top eigenvalue
%Phi are the phases corresponding to the top eigenvector
%D are the eigenvalues


    Ws = Omega;
    A = 2*(A0-A0'); % directed adjacency matrix
    deg = sum(Ws,2); Deg = diag(deg); % Degree matrix
    Tg = exp(1) .^(2*pi*1i*g*A.');  % Transporter
    Lg = Deg - Ws.*Tg; % Magnetic Laplacian.
    [V,D] = eigs(Lg,size(Ws,1),'smallestabs'); % All eigenvectors ranging from smallest eigenvalue to largest eigenvalue
    D = diag(D);
    Phi = angle(V); %phases corresponding to all eigenvectors

end
