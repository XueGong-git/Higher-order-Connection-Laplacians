% Triangle 4->1, 4->2, 4->3, 1->2, 2->3, 3->1 and 1->2->3
clear
%%% B1
% edges {1, 2}, {1, 3}, {2, 3}, {1, 4}, {2, 4}, {3, 4}
l = [1 10 19 4 13; 
    2 13 20 5 16;
    3 16 21 6 10;
    4 11 22 7 14;
    5 14 23 8 17;
    6 17 24 9 11;
    7 12 25 1 15;
    8 15 26 2 18;
    9 18 27 3 12;];
i1 = 1i;


max_step = 96;
x_grid = (2*pi/max_step) * linspace(1,max_step, max_step);
%L1_up_square = zeros(27, 27); L1_down_square = zeros(27, 27);

for n=1:max_step
%delta = 2*pi/6;    

 delta=2*pi*n/max_step;

    
L1_up_square = zeros(27, 27); L1_down_square = zeros(27, 27);


L1_up=[1, -exp(i1*delta), exp(-i1*delta);
        -exp(-i1*delta), 1, -exp(i1*delta); 
        exp(i1*delta), -exp(-i1*delta),1];
    
L1_down=[2,exp(-i1*delta),-exp(i1*delta); 
           exp(i1*delta),2,exp(-i1*delta); 
           -exp(-i1*delta),exp(i1*delta),2];
       
for row = 1:size(l,1)
    [L1_up_square, L1_down_square] = L1_square_b(L1_up, L1_down, l(row,:), L1_up_square, L1_down_square);
    
end

% L1_up_square = kron(L1_up_square, eye(2)); L1_down_square = kron(L1_down_square, eye(2));

Lup =  L1_up_square; Ldown = L1_down_square;
anti_commutator(n) = sum(sum((Lup*Ldown-Ldown*Lup).^2)); % squared frobenius norm of anti-commutator
[V_up, D] = eigs(Lup, 27, 'smallestreal');
%V_up(n,:,:) = reshape(V, [1, 27, 27]); 
lambda_up(n,:) = diag(D);
[V_down, D] = eigs(Ldown, 27, 'smallestreal');
%V_down(n,:,:) = reshape(V, [1, 27, 27]); 
lambda_down(n,:) = diag(D);
[V_L, D] = eigs(Ldown+Lup, 27, 'smallestreal');
%V_L(n,:,:) = reshape(V, [1, 27, 27]); 
lambda_L(n,:) = diag(D);
hodge_norm1(n) = norm(Lup*Ldown,"fro"); % Frobenius norm of Lup*Ldown
hodge_norm2(n) = norm(Ldown*Lup,"fro"); % calculate   

end

%% Plot eigenvalues

figure
subplot(2,2,1);
plot(x_grid, lambda_up(:,1))
hold on 
for ind = 2:27
    plot(x_grid, lambda_up(:,ind))
end
%plot(x_grid, 1+2*cos(x_grid), '-.')
hold off
title('Spectrum of L_{[1]}^{M, up}');
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
% Add a label to the x-axis
xlabel('\delta')

subplot(2,2,2);
plot(x_grid, lambda_down(:,1))
hold on 
for ind = 2:27
    plot(x_grid, lambda_down(:,ind))
end
hold off
title('Spectrum of L_{[1]}^{M, down}');
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
% Add a label to the x-axis
xlabel('\delta')

subplot(2,2,3);
plot(x_grid, lambda_L(:,1))
hold on 
for ind = 2:27
    plot(x_grid, lambda_L(:,ind))
end

hold off
title('Spectrum of L_{[1]}^{M}');
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
% Add a label to the x-axis
xlabel('\delta')
%set(gca, 'FontSize', 14);
saveas(gcf, 'plots\case1.eps', 'epsc');

subplot(2,2,4);
plot(x_grid, anti_commutator)
title('Commutator');
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
% Add a label to the x-axis
xlabel('\delta')
%ylim([-1 1])
saveas(gcf, 'plots\case1.eps', 'eps');


figure
plot(x_grid, hodge_norm1)
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
title('Frobenius Norm of Lup*Ldown');

figure
plot(x_grid, hodge_norm2)
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
title('Frobenius Norm of Ldown*Lup');


%% Visualize eigenvectors
ind_1 = 1:27;
phi_11 = angle(V_up(ind_1,1)); 
phi_12 = angle(V_down(ind_1,1)); 
phi_13 = angle(V_L(ind_1,1)); 

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

cc=colormap(hsv(100));
mc=-pi;
Mac=pi;

figure
colormap(hsv(100))
subplot(1,3,1);

for n=1:27
    c2=(floor((phi_11(n)-mc)*(90)/(Mac-mc)))+1;
    plot([x(n, 1),x(n,2)],[y(n, 1),y(n,2)],'Color',cc(c2,:),'LineWidth',3)
    hold on
    
end

colorbar('Ticks',[0,0.25,0.5,0.75, 1],...
         'TickLabels',{'-\pi','-\pi/2','0','\pi/2','\pi'})
title('\phi_1 of L_{[1,up]}^{M}', 'fontsize', 12);

subplot(1,3,2);
for n=1:27
    c2=(floor((phi_12(n)-mc)*(90)/(Mac-mc)))+1;
    plot([x(n, 1),x(n,2)],[y(n, 1),y(n,2)],'Color',cc(c2,:),'LineWidth',3)
    hold on
    
end
title('\phi_2 of L_{[1,down]}^{M}', 'fontsize', 12);
colorbar('Ticks',[0,0.25,0.5,0.75, 1],...
         'TickLabels',{'-\pi','-\pi/2','0','\pi/2','\pi'})


subplot(1,3,3);
for n=1:27
    c2=(floor((phi_13(n)-mc)*(90)/(Mac-mc)))+1;
    plot([x(n, 1),x(n,2)],[y(n, 1),y(n,2)],'Color',cc(c2,:),'LineWidth',3)
    hold on
    
end
title('\phi_2 of L_{[1]}^{M}', 'fontsize', 12);
colorbar('Ticks',[0,0.25,0.5,0.75, 1],...
         'TickLabels',{'-\pi','-\pi/2','0','\pi/2','\pi'})     
    
    