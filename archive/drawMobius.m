% Define parameters
resolution = 50;   % Resolution of the grid

u = linspace(0, 2*pi, resolution);  % Parameter along the length of the strip
v = linspace(-0.25, 0.25, resolution);  % Parameter along the width of the strip
[u, v] = meshgrid(u, v);

% Define the Möbius strip
x = (1 + v.*cos(u./2)).*cos(u);
y = (1 + v.*cos(u./2)).*sin(u);
z = v.*sin(u./2);

% Plot the Möbius strip
figure
mesh(x, y, z);
axis equal
saveas(gcf, 'plots\mobius.png');
