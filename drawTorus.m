% Torus parameters
R = 3;   % Major radius
r = 1;   % Minor radius
resolution = 50;   % Resolution of the grid

% Create the grid
theta = linspace(0, 2*pi, resolution);
phi = linspace(0, 2*pi, resolution);
[theta, phi] = meshgrid(theta, phi);

% Calculate the torus coordinates
x = (R + r*cos(phi)) .* cos(theta);
y = (R + r*cos(phi)) .* sin(theta);
z = r*sin(phi);

% Plot the torus
figure;
surf(x, y, z);
axis equal;
%xlabel('X');
%ylabel('Y');
%zlabel('Z');
%title('Torus');

saveas(gcf, 'plots\torus.png');


