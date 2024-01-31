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
x0 = [1;1;1]/sqrt(3);
%x0 = randn(3,1);

% Set the size of the complex vector
vector_size = 3;

% Generate random real and imaginary parts
real_parts = randn(vector_size, 1);
imaginary_parts = randn(vector_size, 1);

% Combine real and imaginary parts into a complex vector
x0 = complex(real_parts, imaginary_parts);


% The undirected infered adjacency matrix
A0u = diag(diag(L0u)) - L0u; % this is the same as the standard adjacency matrix for standard graph
D1u = diag(diag(L1u));
A1u = D1u - L1u; % infered adjacency between edges
D1d = diag(diag(L1d));
A1d = D1d - L1d; % infered adjacency between edges
i1 = 1i;

% Triangle 1->2, 2->3, 3->1 and 1->2->3
x_grid = (2*pi/100) * linspace(1,100);


delta=pi/3;
Lup=[2, -exp(i1*delta), exp(-i1*delta);
    -exp(-i1*delta), 2, -exp(i1*delta); 
    exp(i1*delta), -exp(-i1*delta),2];

Ldown=[2,exp(-i1*delta),-exp(i1*delta); 
       exp(i1*delta),2,exp(-i1*delta); 
       -exp(-i1*delta),exp(i1*delta),2];
lambda_up=max(0, sort(eig((Lup))));
lambda_down=sort(eig((Ldown)));
[Uu, Du] = eig(Lup);
[Ud, Dd] = eig(Ldown);
[UL, DL] = eig(Lup + Ldown);


max_step = 300; steps =  1:max_step; step_size = 1/50;
time_step = steps*step_size;
xtu = zeros(3, max_step); xtd = zeros(3, max_step); xtl = zeros(3, max_step); 
 
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
hold off
title('Up', 'Interpreter', 'latex');
xlabel('t')
%yticks([-pi -3*pi/4 -pi/2 -pi/4 0])
%yticklabels({'-\pi','-3\pi/4','-\pi/2','-\pi/4', '0'})
ylabel('\psi(t)')
axis square; % makes the plot square



subplot(1,3,2);
plot(time_step, angle(xtd(1,:)))
hold on 
plot(time_step, angle(xtd(2,:)))
plot(time_step, angle(xtd(3, :)))
hold off
title('Down', 'Interpreter', 'latex');
xlabel('t')
%yticks([-pi -3*pi/4 -pi/2 -pi/4 0])
%yticklabels({'-\pi','-3\pi/4','-\pi/2','-\pi/4', '0'})
axis square; % makes the plot square

subplot(1,3,3);
plot(time_step, angle(xtl(1,:)))
hold on 
plot(time_step, angle(xtl(2,:)))
plot(time_step, angle(xtl(3, :)))
hold off
title('Up + Down', 'Interpreter', 'latex');
xlabel('t')
%yticks([-pi -3*pi/4 -pi/2 -pi/4 0])
%yticklabels({'-\pi','-3\pi/4','-\pi/2','-\pi/4', '0'})
axis square; % makes the plot square
legend('$\psi_1$', '$\psi_2$', '$\psi_3$');

%saveas(gcf, 'plots\triangle_1_diffusion.png');
saveas(gcf, 'plots\triangle_1_diffusion_0.5pi.eps', 'epsc');
