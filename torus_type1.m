% Triangle 4->1, 4->2, 4->3, 1->2, 2->3, 3->1 and 1->2->3
clear

%%% B1
% edges {1, 2}, {1, 3}, {2, 3}, {1, 4}, {2, 4}, {3, 4}
lup = [1 10 19 4 13; 
    2 13 20 5 16;
    3 16 21 6 10;
    4 11 22 7 14;
    5 14 23 8 17;
    6 17 24 9 11;
    7 12 25 1 15;
    8 15 26 2 18;
    9 18 27 3 12;];

ldown = [3 21 10 1 25 12; 
         1 19 13 2 26 15; 
         2 20 16 3 27 18;
         6 24 11 4 19 10;
         4 22 14 5 20 13;
         5 23 17 6 21 16;
         9 27 12 7 22 11;
         7 25 15 8 23 14;
         8 26 18 9 24 17];

i1 = 1i;
sx = [0 1; 1 0]; sy = [0 -1i; 1i 0]; sz = [1 0; 0 -1];


max_step = 96;
x_grid = (2*pi/max_step) * linspace(1,max_step, max_step);
%L1_up_square = zeros(27, 27); L1_down_square = zeros(27, 27);

for n=1:max_step
%delta = 2*pi/6;    

delta=2*pi*n/max_step;

    
L1_up_square = zeros(27, 27); L1_down_square = zeros(54, 54);


L1_up=[2, -exp(i1*delta), exp(-i1*delta);
        -exp(-i1*delta), 2, -exp(i1*delta); 
        exp(i1*delta), -exp(-i1*delta),2];
    
L1_down=[5*eye(2),exp(-i1*delta)*eye(2),sy,exp(-i1*delta)*eye(2),sy,exp(-i1*delta)*eye(2); 
          exp(i1*delta)*eye(2), 5*eye(2), exp(i1*delta)*eye(2), sz, exp(i1*delta)*eye(2), sz; 
          sy, exp(-i1*delta)*eye(2), 5*eye(2), exp(-i1*delta)*eye(2), sy, exp(-i1*delta)*eye(2);
          exp(i1*delta)*eye(2), sz, exp(i1*delta)*eye(2), 5*eye(2), exp(i1*delta)*eye(2), sz;
          sy, exp(-i1*delta)*eye(2), sy, exp(-i1*delta)*eye(2), 5*eye(2), exp(-i1*delta)*eye(2);
          exp(i1*delta)*eye(2), sz, exp(i1*delta)*eye(2), sz, exp(i1*delta)*eye(2), 5*eye(2)];
       
for row = 1:size(lup,1)
    L1_up_square = L1_square_up_type1(L1_up, lup(row,:), L1_up_square);
end

for row = 1:size(ldown,1)
    L1_down_square = L1_square_down(L1_down, ldown(row,:), L1_down_square);
end

L1_up_square = kron(L1_up_square, eye(2));
%L1_down_square = kron(L1_down_square, eye(2));

Lup =  L1_up_square; Ldown = L1_down_square;
anti_commutator(n) = sum(sum((Lup*Ldown-Ldown*Lup).^2)); % squared frobenius norm of anti-commutator
[V_up, D_up] = eigs(Lup, 54, 'smallestreal');
%V_up(n,:,:) = reshape(V, [1, 27, 27]); 
lambda_up(n,:) = diag(D_up);
[V_down, D_down] = eigs(Ldown, 54, 'smallestreal');
%V_down(n,:,:) = reshape(V, [1, 27, 27]); 
lambda_down(n,:) = diag(D_down);
[V_L, D_L] = eigs(Ldown+Lup, 54, 'smallestreal');
%V_L(n,:,:) = reshape(V, [1, 27, 27]); 
lambda_L(n,:) = diag(D_L);
hodge_norm1(n) = norm(Lup*Ldown,"fro"); % Frobenius norm of Lup*Ldown
hodge_norm2(n) = norm(Ldown*Lup,"fro"); % calculate   



%% Visualize eigenvectors for delta = 2*pi/3

if max_step/n == 3
    
    ind_1 = 1:2:54;
    ind_2 = 2:2:54;

    %find eigenvectors correspondign to non-zero eigenvalues
    k_up = find(abs(diag(D_up))>0.01, 1); k_down = find(abs(diag(D_down))>0.01, 1); k_L = find(abs(diag(D_L))>0.01, 1);
    V_up(37:54,:) = - V_up(37:54,:);
    V_down(37:54,:) = - V_down(37:54,:);
    V_L(37:54,:) = - V_L(37:54,:);
    phi_11 = angle(V_up(ind_1,1)); phi_12 = angle(V_up(ind_2,1)); 
    phi_21 = angle(V_down(ind_1,1)); phi_22 = angle(V_down(ind_2,1)); 
    phi_31 = angle(V_L(ind_1,1)); phi_32 = angle(V_L(ind_2,1)); 

    x = zeros(27, 2); y = zeros(27, 2);

    % calculate coordinate
    for i = 1:9
        y(i,:) = [floor((i-1)/3), floor((i-1)/3)];
    end
    x(1:9, :) = repmat([0 1; 1 2; 2 3;], 3, 1);

    for i = 10:18
        x(i,:) = [floor((i-10)/3), floor((i-10)/3)];
    end
    y(10:18, :) = repmat([0 1; 1 2; 2 3;], 3, 1);

    x(19:27, :) = repmat([0 1; 1 2; 2 3;], 3, 1);
    y(19:27, :) = [1 0; 1 0; 1 0; 2 1; 2 1; 2 1; 3 2; 3 2; 3 2];



    %% Draw colored edges
    cc=colormap(turbo(100));
    mc=-pi;
    Mac=pi;

    figure
    colormap(turbo(100))
    subplot(3,2,1);

    for n=1:27
        c2=(floor((phi_11(n)-mc)*(90)/(Mac-mc)))+1;
        plot([x(n, 1),x(n,2)],[y(n, 1),y(n,2)],'Color',cc(c2,:),'LineWidth',3)
        hold on
    end

    colorbar('Ticks',[0,0.25,0.5,0.75, 1],...
             'TickLabels',{'-\pi','-\pi/2','0','\pi/2','\pi'})
    title('$\theta$ of $L_{[1,up]}^{M}$', 'fontsize', 12, 'Interpreter', 'latex');

    
    subplot(3,2,2);
    for n=1:27
        c2=(floor((phi_12(n)-mc)*(90)/(Mac-mc)))+1;
        plot([x(n, 1),x(n,2)],[y(n, 1),y(n,2)],'Color',cc(c2,:),'LineWidth',3)
        hold on

    end
    title('$\phi$ of $L_{[1,up]}^{M}$', 'fontsize', 12, 'Interpreter', 'latex');
    colorbar('Ticks',[0,0.25,0.5,0.75, 1],...
             'TickLabels',{'-\pi','-\pi/2','0','\pi/2','\pi'})


    subplot(3,2,3);
    for n=1:27
        c2=(floor((phi_21(n)-mc)*(90)/(Mac-mc)))+1;
        plot([x(n, 1),x(n,2)],[y(n, 1),y(n,2)],'Color',cc(c2,:),'LineWidth',3)
        hold on
    end
    title('$\theta$ of $L_{[1,down]}^{M}$', 'fontsize', 12, 'Interpreter', 'latex');
    colorbar('Ticks',[0,0.25,0.5,0.75, 1],...
             'TickLabels',{'-\pi','-\pi/2','0','\pi/2','\pi'})     

    subplot(3,2,4);
    for n=1:27
        c2=(floor((phi_22(n)-mc)*(90)/(Mac-mc)))+1;
        plot([x(n, 1),x(n,2)],[y(n, 1),y(n,2)],'Color',cc(c2,:),'LineWidth',3)
        hold on
    end
    title('$\phi$ of $L_{[1,down]}^{M}$', 'fontsize', 12, 'Interpreter', 'latex');
    colorbar('Ticks',[0,0.25,0.5,0.75, 1],...
             'TickLabels',{'-\pi','-\pi/2','0','\pi/2','\pi'}) 
         
    subplot(3,2,5);
    for n=1:27
        c2=(floor((phi_31(n)-mc)*(90)/(Mac-mc)))+1;
        plot([x(n, 1),x(n,2)],[y(n, 1),y(n,2)],'Color',cc(c2,:),'LineWidth',3)
        hold on
    end
    title('$\theta$ of $L_{[1]}^{M}$', 'fontsize', 12, 'Interpreter', 'latex');
    colorbar('Ticks',[0,0.25,0.5,0.75, 1],...
             'TickLabels',{'-\pi','-\pi/2','0','\pi/2','\pi'})     
         
    subplot(3,2,6);
    for n=1:27
        c2=(floor((phi_32(n)-mc)*(90)/(Mac-mc)))+1;
        plot([x(n, 1),x(n,2)],[y(n, 1),y(n,2)],'Color',cc(c2,:),'LineWidth',3)
        hold on
    end
    title('$\phi$ of $L_{[1]}^{M}$', 'fontsize', 12, 'Interpreter', 'latex');
    colorbar('Ticks',[0,0.25,0.5,0.75, 1],...
             'TickLabels',{'-\pi','-\pi/2','0','\pi/2','\pi'})     

         
    saveas(gcf, 'plots\torus_type1_eigvec_0.6pi.eps', 'epsc');
    %saveas(gcf, 'plots\torus_type1_eigvec.png');

end

end

%{
%% Plot eigenvalues
figure
subplot(2,2,1);
plot(x_grid, lambda_up(:,1))
hold on 
for ind = 2:54
    plot(x_grid, lambda_up(:,ind))
end
%plot(x_grid, 1+2*cos(x_grid), '-.')
hold off
title('Spectrum of $L_{[1]}^{M, up}$', 'Interpreter', 'latex');
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'$0$', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$'})
ylim([-0.5 8.5])
% Add a label to the x-axis
xlabel('\delta')

subplot(2,2,2);
plot(x_grid, lambda_down(:,1))
hold on 
for ind = 2:54
    plot(x_grid, lambda_down(:,ind))
end
hold off
title('Spectrum of $L_{[1]}^{M, down}$', 'Interpreter', 'latex');
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'$0$', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$'})
ylim([0 20])

% Add a label to the x-axis
xlabel('\delta')

subplot(2,2,3);
plot(x_grid, lambda_L(:,1))
hold on 
for ind = 2:54
    plot(x_grid, lambda_L(:,ind))
end
hold off
title('Spectrum of $L_{[1]}^{M}$', 'Interpreter', 'latex');
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'$0$', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$'})
ylim([5 24])

% Add a label to the x-axis
xlabel('\delta')
%set(gca, 'FontSize', 14);

subplot(2,2,4);
plot(x_grid, anti_commutator)
title('Commutator', 'Interpreter', 'latex');
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'$0$', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$'})
% Add a label to the x-axis
xlabel('\delta')
%ylim([-1 1])
%saveas(gcf, 'plots\torus_type1_spectrum.eps', 'epsc');

lines = findall(gcf, 'Type', 'line');
% Loop through each line and set color to blue
for i = 1:numel(lines)
    lines(i).Color =  [0, 0.4470, 0.7410];
end

saveas(gcf, 'plots\torus_type1_spectrum.eps', 'epsc');


figure
plot(x_grid, hodge_norm1)
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
title('Frobenius Norm of Lup*Ldown', 'Interpreter', 'latex');

figure
plot(x_grid, hodge_norm2)
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
title('Frobenius Norm of Ldown*Lup', 'Interpreter', 'latex');

%}
    