% visualize eigenvector
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
y(19:27, :) = repmat([0 1; 1 2; 2 3;], 3, 1);


%% Draw colored edges
cc=colormap(hsv(100));
mc=-1;
Mac=1;

figure
subplot(2,2,1);

for n=1:27
    c2=(floor((real(cos(phi(n))-mc)*(90)/(Mac-mc))))+1;
    plot([x(I(n)),x(J(n))],[y(I(n)),y(J(n))],'Color',cc(c2,:),'LineWidth',3)
    hold on
end

