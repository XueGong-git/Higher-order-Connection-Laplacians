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
%x0 = [1;0;0]
x0 = randn(6,1);

% The undirected infered adjacency matrix
A0u = diag(diag(L0u)) - L0u; % this is the same as the standard adjacency matrix for standard graph
D1u = diag(diag(L1u));
A1u = D1u - L1u; % infered adjacency between edges
D1d = diag(diag(L1d));
A1d = D1d - L1d; % infered adjacency between edges
i1 = 1*i;

% Triangle 1->2, 2->3, 3->1 and 1->2->3
x_grid = (2*pi/100) * linspace(1,100);


delta=3*pi/4;

Lup(1,:)=[1, 0,                  0, i1*exp(-i1*delta),   0, exp(i1*delta)];
Lup(2,:)=[0, 1,                  -i1*exp(-i1*delta), 0,    exp(i1*delta), 0];

Lup(3,:)=[0, i1*exp(i1*delta),  1,0,                     -exp(-i1*delta), 0];
Lup(4,:)=[-i1*exp(i1*delta), 0,   0,1,                     0, exp(-i1*delta)];

Lup(5,:)=[0, exp(-i1*delta),     -exp(i1*delta), 0,       1,0];
Lup(6,:)=[exp(-i1*delta), 0,     0, exp(i1*delta),        0,1];


Ldown(1,:)=[2,0,               0,-i1,   -exp(i1*delta), 0];
Ldown(2,:)=[0,2,               i1,0,    0, -exp(i1*delta)];

Ldown(3,:)=[0,-i1,             2,0,     1,0];
Ldown(4,:)=[i1,0,              0,2,     0,-1];

Ldown(5,:)=[-exp(-i1*delta),0,   1,0,     2,0];
Ldown(6,:)=[0,-exp(-i1*delta),   0,-1,     0,2];

   
   
lambda_up=max(0, sort(eig((Lup))));
lambda_down=sort(eig((Ldown)));
[Uu, Du] = eig(Lup);
[Ud, Dd] = eig(Ldown);
[UL, DL] = eig(Lup + Ldown);


max_step = 150; steps =  1:max_step; step_size = 1/50;
time_step = steps*step_size;
xtu = zeros(6, max_step); xtd = zeros(6, max_step); xtl = zeros(6, max_step); 
 
    for k = steps
        t = k*step_size;
        
        decay_vector = exp(-t*max(diag(Du), 0)).*(Uu'*x0);
        xtu(:, k) = Uu*decay_vector;
        
        decay_vector = exp(-t*max(diag(Dd), 0)).*(Ud'*x0);
        xtd(:, k) = Ud*decay_vector;
        
        decay_vector = exp(-t*max(diag(DL), 0)).*(UL'*x0);
        xtl(:, k) = UL*decay_vector;
        
    end


%end

figure
subplot(1,3,1);
plot(time_step, angle(xtu(1,:)))
hold on 
plot(time_step, angle(xtu(2,:)))
plot(time_step, angle(xtu(3,:)))
plot(time_step, angle(xtu(4,:)))
plot(time_step, angle(xtu(5,:)))
plot(time_step, angle(xtu(6,:)))
hold off
title('x(t) Up');
xlabel('t')


subplot(1,3,2);
plot(time_step, angle(xtd(1,:)))
hold on 
plot(time_step, angle(xtd(2,:)))
plot(time_step, angle(xtd(3, :)))
plot(time_step, angle(xtd(4,:)))
plot(time_step, angle(xtd(5, :)))
plot(time_step, angle(xtd(6,:)))
hold off
title('x(t) Down');
xlabel('t')

subplot(1,3,3);
plot(time_step, angle(xtl(1,:)))
hold on 
plot(time_step, angle(xtl(2,:)))
plot(time_step, angle(xtl(3, :)))
plot(time_step, angle(xtl(4,:)))
plot(time_step, angle(xtl(5, :)))
plot(time_step, angle(xtl(6,:)))
hold off
title('x(t) Up + Down');
xlabel('t')
