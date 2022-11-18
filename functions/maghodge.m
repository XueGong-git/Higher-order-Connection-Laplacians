function [V,D, Phi] = maghodge(Omega,g)
%G plot the nodes against the eigenvalue
%eig is the top eigenvalue
%Phi are the phases corresponding to the top eigenvector
%D are the eigenvalues


    Ws = Omega;
    A = 2*(Ws-2*tril(Ws));
    deg = sum(Ws,2); Deg = diag(deg); % Degree matrix
    Tg = exp(1) .^(2*pi*1i*g*A.');  % Transporter
    Lg = Deg - Ws.*Tg; % Magnetic Laplacian.
    [V,D] = eigs(Lg,size(Ws,1),'smallestabs'); % All eigenvectors ranging from smallest eigenvalue to largest eigenvalue
    D = diag(D);
    Phi = angle(V); %phases corresponding to all eigenvectors
    %G.Nodes.phase0 = mod(Phi(:,1), 2*pi);%phases corresponding to the top eigenvector
    %G.Nodes.phase1 = mod(Phi(:,2), 2*pi);
    %if (g==0)
    %  G.Nodes.phase0 = V(:,2);%phases corresponding to the second smallest eigenvector
    %end
    %eig = round(D(1,1),4); %top eigenvalue
end
