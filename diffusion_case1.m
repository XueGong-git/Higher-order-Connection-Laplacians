% Triangle 1->2, 2->3, 3->1 and 1->2->3
clear
close all
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



% Set the size of the complex vector
vector_size = 6;

% Generate random real and imaginary parts
real_parts = randn(vector_size, 1);
imaginary_parts = randn(vector_size, 1);

% Combine real and imaginary parts into a complex vector
%x0 = complex(real_parts, imaginary_parts);

x0 = [exp(i*(-0.5));exp(i*(-0.4));exp(i*(0.1));exp(i*(-0.1));exp(i*(0.5));exp(i*(0.4))];

%x0 = randn(6,1);

% The undirected infered adjacency matrix
A0u = diag(diag(L0u)) - L0u; % this is the same as the standard adjacency matrix for standard graph
D1u = diag(diag(L1u));
A1u = D1u - L1u; % infered adjacency between edges
D1d = diag(diag(L1d));
A1d = D1d - L1d; % infered adjacency between edges
i1 = 1*i;

% Triangle 1->2, 2->3, 3->1 and 1->2->3
x_grid = (2*pi/100) * linspace(1,100);


delta = 2*pi/3;

Lup=[2, -exp(i1*delta), exp(-i1*delta);
    -exp(-i1*delta), 2, -exp(i1*delta); 
    exp(i1*delta), -exp(-i1*delta),2];

Ldown=[2,exp(-i1*delta),-exp(i1*delta); 
       exp(i1*delta),2,exp(-i1*delta); 
       -exp(-i1*delta),exp(i1*delta),2];

Lup = kron(Lup, eye(2)); Ldown = kron(Ldown, eye(2));

lambda_up=max(0, sort(eig((Lup))));
lambda_down=sort(eig((Ldown)));
[Uu, Du] = eig(Lup);
[Ud, Dd] = eig(Ldown);
[UL, DL] = eig(Lup + Ldown);

max_step = 250; steps =  1:max_step; step_size = 1/75; max_time = max_step*step_size;
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
for iter = [3,5]
plot(time_step, angle(xtu(iter,:)))
end
for iter = [2,4,6]
plot(time_step, angle(xtu(iter,:)), 'LineStyle', '--')
end


norm_xtu = abs(xtu);
plot(time_step, atan(norm_xtu([2,4,6], :)./norm_xtu([1,3,5], :)), 'LineStyle', '-.')


hold off
title('Up', 'Interpreter', 'latex');
xlim([0, max_time])
xlabel('t', 'Interpreter', 'latex')
%legend('\theta_1','\psi_1','\phi_1','\theta_2','\psi_2', '\phi_2','\theta_3','\psi_3','\phi_3')
set(gca, 'FontSize', 14);
%ylabel('x(t)')
ax = gca; % Get handle to current axes
lines = ax.Children; % Get all lines in the axes
set(lines, 'LineWidth', 1.5); % Set line width to 2 for all lines

axis square; % makes the plot square

subplot(1,3,2);
plot(time_step, angle(xtd(1,:)))
hold on 
for iter = 3:2:6
plot(time_step, angle(xtd(iter,:)))
end
for iter = 2:2:6
plot(time_step, angle(xtd(iter,:)), 'LineStyle', '--')
end

norm_xtd = abs(xtd);
plot(time_step, atan(norm_xtd([2,4,6], :)./norm_xtd([1,3,5], :)), 'LineStyle', '-.')

hold off
title('Down', 'Interpreter', 'latex');
xlim([0, max_time])
xlabel('t', 'Interpreter', 'latex')
%legend('\theta_1','\psi_1','\phi_1','\theta_2','\psi_2', '\phi_2','\theta_3','\psi_3','\phi_3')
set(gca, 'FontSize', 14);
ax = gca; % Get handle to current axes
lines = ax.Children; % Get all lines in the axes
set(lines, 'LineWidth', 1.5); % Set line width to 2 for all lines

axis square; % makes the plot square


subplot(1,3,3);
plot(time_step, angle(xtl(1,:)))
hold on 
for iter = 3:2:6
plot(time_step, angle(xtl(iter,:)))
end
for iter = 2:2:6
plot(time_step, angle(xtl(iter,:)), 'LineStyle', '--')
end
norm_xtl = abs(xtl);
plot(time_step, atan(norm_xtl([2,4,6], :)./norm_xtl([1,3,5], :)), 'LineStyle', '-.')

hold off
title('Up + Down', 'Interpreter', 'latex');
xlim([0, max_time])

xlabel('t', 'Interpreter', 'latex')
%ylim([-1 1])
legend('$\theta_1$','$\theta_2$','$\theta_3$','$\phi_1$','$\phi_2$','$\phi_3$', '$\psi_1$', '$\psi_2$','$\psi_3$', 'Interpreter', 'latex', 'Location', 'eastoutside','Position',[0.913289468357835 0.347782497525199 0.080966468355549 0.343047616595314]);
set(gca, 'FontSize', 14);

ax = gca; % Get handle to current axes
lines = ax.Children; % Get all lines in the axes
set(lines, 'LineWidth', 1.5); % Set line width to 2 for all lines

axis square; % makes the plot square

saveas(gcf, 'plots\triangle_2_diffusion_0.6pi.eps', 'epsc');


createfigure(time_step, angle(xtu([1, 3, 5, 2, 4, 6, 1, 3, 5],:)), angle(xtu([1, 3, 5, 2, 4, 6, 1, 3, 5],:)), angle(xtu([1, 3, 5, 2, 4, 6, 1, 3, 5],:)))