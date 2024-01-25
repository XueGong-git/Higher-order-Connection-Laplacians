function [L1up, L1down] = constructMagneticLaplacian(V, E, T, delta, L1u, L1d)
T_up = zeros(2*size(E,1), 2*size(E,1));
Theta_up = zeros(size(E,1), size(E,1));
Theta_down = zeros(size(E,1), size(E,1));
P_up = zeros(2*size(E,1), 2*size(E,1));
sx = [0 1; 1 0]; sy = [0 -1i; 1i 0]; sz = [1 0; 0 -1];

% The undirected infered adjacency matrix
D1_up = diag(diag(L1u));
A1_up = diag(diag(L1u)) - L1u; % infered adjacency between edges

D1_down = diag(diag(L1d));
A1_down = diag(diag(L1d)) - L1d; % this is the same as the standard adjacency matrix for standard graph


for n = 1:size(T,1) % triangle
    % find edges
    edgeInd = find(ismember(E(:,2), T(n, :)) == 1 & ismember(E(:,1), T(n, :)) == 1);
    for x = 1:3 % edge
        e1 = edgeInd(x);
        for y = 1:3
            e2 = edgeInd(y);
            if e1 ~= e2
                v1 = find(ismember(T(n, :), E(e1,1)) == 1);
                v2 = find(ismember(T(n, :), E(e1,2)) == 1);
                if E(e1,2) == E(e2,1) % tail of e1 == head of e2
                    % find index of common nodes in triangle
                    
                    if  mod(v2-v1,3) == 1 % case 1
                        Theta_up(e1, e2) = - delta; % ijk; jki; kij
                        P_up(2*e1-1: 2*e1, 2*e2-1: 2*e2) = eye(2);
                        T_up(2*e1-1: 2*e1, 2*e2-1: 2*e2) = exp(-1i*delta)*eye(2);
                    else % case 4
                        Theta_up(e1, e2) = delta; 
                        P_up(2*e1-1: 2*e1, 2*e2-1: 2*e2) = sx;
                        T_up(2*e1-1: 2*e1, 2*e2-1: 2*e2) = exp(1i*delta)*sx;

                    end
                    
                elseif E(e1,1) == E(e2,2) % head of e1 == tail of e2
                    if  mod(v2-v1,3) == 1 % case 2: kji; jik; ikj triangle along the edge
                        Theta_up(e1, e2) =  delta; 
                        P_up(2*e1-1: 2*e1, 2*e2-1: 2*e2) = eye(2);
                        T_up(2*e1-1: 2*e1, 2*e2-1: 2*e2) = exp(1i*delta)*eye(2);

                    else % case 3 : ijk; jki; kij   triangle against edge
                        Theta_up(e1, e2) = -delta; 
                        P_up(2*e1-1: 2*e1, 2*e2-1: 2*e2) = sx;
                        T_up(2*e1-1: 2*e1, 2*e2-1: 2*e2) = exp(-1i*delta)*sx;

                    end
                    
                elseif E(e1,1) == E(e2,1) % head of e1 == head of e2
                    if  mod(v2-v1,3) == 1 % case 5 ijk; jki; kij
                        Theta_up(e1, e2) = - delta; % ijk; jki; kij
                        P_up(2*e1-1: 2*e1, 2*e2-1: 2*e2) = sy;
                        T_up(2*e1-1: 2*e1, 2*e2-1: 2*e2) = exp(-1i*delta)*sy;

                    else % case 6
                        Theta_up(e1, e2) = delta; 
                        P_up(2*e1-1: 2*e1, 2*e2-1: 2*e2) = sy;
                        T_up(2*e1-1: 2*e1, 2*e2-1: 2*e2) = exp(1i*delta)*sy;

                    end
                elseif E(e1,2) == E(e2,2) % tail of e1 == tail of e2
                    if  mod(v2-v1,3) == 1 % case 7 ijk; jki; kij
                        Theta_up(e1, e2) = - delta; % ijk; jki; kij
                        P_up(2*e1-1: 2*e1, 2*e2-1: 2*e2) = sz;
                        T_up(2*e1-1: 2*e1, 2*e2-1: 2*e2) = exp(-1i*delta)*sz;

                    else % case 8
                        Theta_up(e1, e2) = delta; 
                        P_up(2*e1-1: 2*e1, 2*e2-1: 2*e2) = sz;
                        T_up(2*e1-1: 2*e1, 2*e2-1: 2*e2) = exp(1i*delta)*sz;

                    end
                end
            end              
        end
    end
end

L1up = kron(D1_up, eye(2)) - T_up.* kron(A1_up, ones(2));



%% Down Laplacian

for n = 1:size(V,2) % triangle
    % find edges
    [edgeInd, ~] = find(ismember(E(), V(n)) == 1 );
    for x = 1:size(edgeInd,1) % edge
        e1 = edgeInd(x);
        for y = 1:size(edgeInd,1)
            e2 = edgeInd(y);
            if e1 ~= e2
                v1 = find(ismember(E(e1,:), V(n)) == 1);
                v2 = find(ismember(E(e2,:), V(n)) == 1);
                if v1 == 1 & v2 == 1 % j-i, j-k
                        Theta_down(e1, e2) = 0; 
                        P_down(2*e1-1: 2*e1, 2*e2-1: 2*e2) = sy;
                        T_down(2*e1-1: 2*e1, 2*e2-1: 2*e2) = sy;
                elseif v1 == 1 & v2 == 2  % case 2
                        Theta_down(e1, e2) = -delta; 
                        P_down(2*e1-1: 2*e1, 2*e2-1: 2*e2) = eye(2);
                        T_down(2*e1-1: 2*e1, 2*e2-1: 2*e2) = exp(-1i*delta)*eye(2);
                elseif v1 == 2 & v2 == 1  % case 1
                        Theta_down(e1, e2) = delta; 
                        P_down(2*e1-1: 2*e1, 2*e2-1: 2*e2) = eye(2);
                        T_down(2*e1-1: 2*e1, 2*e2-1: 2*e2) = exp(1i*delta)*eye(2);
                elseif v1 == 2 & v2 == 2  % case 4
                        Theta_down(e1, e2) = 0; 
                        P_down(2*e1-1: 2*e1, 2*e2-1: 2*e2) = sz;
                        T_down(2*e1-1: 2*e1, 2*e2-1: 2*e2) = sz;

                end
            end
                    
               

        end
    end
end             

L1down = kron(D1_down, eye(2)) - T_down.* kron(A1_down, ones(2));
end
