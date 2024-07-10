map = viridis;

% Assuming absorbed_energy is already defined as a 419x282x328 matrix
% Interpolating the absorbed_energy matrix by a factor of 2 in each dimension

% Original size
[x_dim, y_dim, z_dim] = size(absorbed_energy);

% Define the original grid
x = linspace(-x_dim / 2, x_dim / 2, x_dim);
y = linspace(-y_dim / 2, y_dim / 2, y_dim);
z = linspace(-z_dim / 2, z_dim / 2, z_dim);
[X, Y, Z] = meshgrid(y, x, z);  % Note the order for compatibility with meshgrid

% Define the new interpolated grid
xq = linspace(-x_dim / 2, x_dim / 2, 2 * x_dim - 1);
yq = linspace(-y_dim / 2, y_dim / 2, 2 * y_dim - 1);
zq = linspace(-z_dim / 2, z_dim / 2, 2 * z_dim - 1);
[Xq, Yq, Zq] = meshgrid(yq, xq, zq);  % Note the order for compatibility with meshgrid

% Interpolating the absorbed_energy matrix
absorbed_energy_interp = interp3(X, Y, Z, absorbed_energy, Xq, Yq, Zq, 'linear');

% Clipping extreme values (set threshold as needed)
max_value_threshold = prctile(absorbed_energy_interp(:), 99); % 99th percentile
absorbed_energy_interp(absorbed_energy_interp > max_value_threshold) = max_value_threshold;

% Creating a slice plot
figure;
slice(Xq, Yq, Zq, absorbed_energy_interp, 0, 0,50);
shading interp;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Interpolated Absorbed Energy Slice Plot');
colorbar;

% Set color axis limits
clim([min(absorbed_energy_interp(:)), max_value_threshold]);

% Choose a distinct colormap
colormap(map);  % You can choose other colormaps like 'parula', 'hot', etc.






% % % % % % % % 
% % % % % % % % % Define axes
% % % % % % % % x = linspace(-size(absorbed_energy, 1) / 2, size(absorbed_energy, 1) / 2, size(absorbed_energy, 1));
% % % % % % % % y = linspace(-size(absorbed_energy, 2) / 2, size(absorbed_energy, 2) / 2, size(absorbed_energy, 2));
% % % % % % % % f = dataset.frequencies(1:44); % Actual frequency values
% % % % % % % % 
% % % % % % % % % Create finer grids for smoother interpolation
% % % % % % % % [xq, yq, fq] = meshgrid(linspace(min(y), max(y), 200), linspace(min(x), max(x), 200), linspace(min(f), max(f), 200));
% % % % % % % % 
% % % % % % % % % Interpolate the data
% % % % % % % % data_interp = interp3(y, x, f, dataset.spec_all(:, :, 1:44), xq, yq, fq, 'linear');
% % % % % % % % 
% % % % % % % % figure; % Create a new figure for each plot
% % % % % % % % h = slice(xq, yq, fq, data_interp, 0, 0, []);
% % % % % % % % set(h, 'EdgeColor', 'none'); % Remove edge colors
% % % % % % % % colormap(map);
% % % % % % % % colorbar;
% % % % % % % % axis tight
% % % % % % % % grid on;
% % % % % % % % grid minor;
% % % % % % % % xlabel('Y Axis');
% % % % % % % % ylabel('X Axis');
% % % % % % % % zlabel('Frequency');
% % % % % % % % 
