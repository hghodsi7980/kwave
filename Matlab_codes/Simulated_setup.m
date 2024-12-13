close all
clear

addpath(genpath('C:\GitHub\kwave\Software_packages'))
mex   -DUSE_OMP C:\GitHub\kwave\Software_packages\ValoMC-master\ValoMC-master\cpp\2d\MC2Dmex.cpp COMPFLAGS='\$COMPFLAGS /openmp'


%% clot Parameters
matrix_size = 200;       % 250x250x250 matrix
voxel_size = 60e-6;      % Each voxel is 50 um = 50e-6 meters
cylinder_diameter = 9e-3; % Diameter in meters
cylinder_height = 2e-3;   % Height in meters
rotation_angle = 135;    % Rotation around the x-axis in degrees


%% Create cylinderical clot mask
cylinder_radius_voxels = (cylinder_diameter / 2) / voxel_size;
cylinder_height_voxels = cylinder_height / voxel_size;
rotation_matrix = [pi/4 0 0];
clot_mask = makeCylinder(matrix_size,matrix_size,matrix_size,matrix_size/2,matrix_size/2,matrix_size/2,cylinder_radius_voxels,cylinder_height_voxels,rotation_matrix);
voxelPlot(double(clot_mask));




%% generate simulation mesh

x_arr = linspace(round(-matrix_size/2),round(matrix_size/2),round(matrix_size))*voxel_size*1000; %mm
y_arr = linspace(round(-matrix_size/2),round(matrix_size/2),round(matrix_size))*voxel_size*1000;
z_arr = linspace(round(-matrix_size/2),round(matrix_size/2),round(matrix_size))*voxel_size*1000;
vmcmesh = createGridMesh(x_arr, y_arr, z_arr); % function provided by ValoMC
vmcmedium = createMedium(vmcmesh);
[X,Y,Z] = meshgrid(y_arr,x_arr,z_arr); % Matlab function
scattering_coefficients = 0*ones(size(clot_mask));
absorption_coefficients = 0*ones(size(clot_mask));
scattering_anisotropies = 0.9678*ones(size(clot_mask));
refractive_indexes = 1.0*ones(size(clot_mask));
% Homogenius clot case with whole blood 45% RBC 5% platelet Fibrin. 50%
% Plasema
absorption_coefficients(clot_mask == 1) = 10.366;
scattering_coefficients(clot_mask == 1) = 72.78/2;
absorption_coefficients(clot_mask == 2) = 5; %% dummmy value
scattering_coefficients(clot_mask == 2) = 72.78; %% dummmy value
absorption_coefficients(clot_mask == 3) = 5; %% dummmy value
scattering_coefficients(clot_mask == 3) = 72.78; %% dummmy value
refractive_indexes(clot_mask == 1) = 1.4;
refractive_indexes(clot_mask == 2) = 1.4;
refractive_indexes(clot_mask == 3) = 1.4;
vmcmedium.scattering_anisotropy =  repmat(scattering_anisotropies(:),6,1);
vmcmedium.absorption_coefficient = repmat(absorption_coefficients(:),6,1);
vmcmedium.scattering_coefficient = repmat(scattering_coefficients(:),6,1);
vmcmedium.refractive_index = repmat(refractive_indexes(:),6,1);
vmcboundary = createBoundary(vmcmesh, vmcmedium);   % create a boundary for the mesh
clearvars('-except','clot_mask','absorption_coefficients','vmcboundary','vmcmesh','vmcmedium','matrix_size','voxel_size');
%% Create a light sourceaddpath 'C:\GitHub\kwave\Software_packages\k-wave-toolbox-version-1.4\k-Wave'
lightsource = findBoundaries(vmcmesh, 'direction',[0 0 0], [0 -10 10], 5);
vmcboundary.lightsource(lightsource) = {'cosinic'};
solution = ValoMC(vmcmesh, vmcmedium, vmcboundary);
x_arr = linspace(round(-matrix_size/2),round(matrix_size/2),round(matrix_size))*voxel_size*1000; %mm
y_arr = linspace(round(-matrix_size/2),round(matrix_size/2),round(matrix_size))*voxel_size*1000;
z_arr = linspace(round(-matrix_size/2),round(matrix_size/2),round(matrix_size))*voxel_size*1000;
[X,Y,Z] = meshgrid(y_arr,x_arr,z_arr); % Matlab function
TR = triangulation(double(vmcmesh.H),vmcmesh.r); % create a matlab
locations = [X(:) Y(:) Z(:)];              % form a 2D matrix from all
indices = pointLocation(TR,locations);     % query the indices of the
indices(isnan(indices)) = 1;               % set the grid points that
grid_fluence = reshape(solution.element_fluence(indices),size(X));
absorbed_energy = absorption_coefficients .* grid_fluence*1e3; % [J/m3]
clearvars ('vmcboundary','vmcmedium','vmcmesh','absorption_coefficients');
% Determine the slice indices for the central planes
slice_x = round(size(grid_fluence, 1) / 2); % Index for yz-plane
slice_y = round(size(grid_fluence, 2) / 2); % Index for xz-plane
slice_z = round(size(grid_fluence, 3) / 2); % Index for xy-plane

% Create figure and subplots
figure;

% Plot the fluence rate in the XY-plane
subplot(1, 3, 1);
imagesc(x_arr, y_arr, squeeze(grid_fluence(:, :, slice_z))');
axis equal tight;
colorbar;
xlabel('X (mm)');
ylabel('Y (mm)');
title('Fluence Rate in XY-plane');
set(gca, 'YDir', 'normal'); % Adjust to correct orientation

% Plot the fluence rate in the XZ-plane
subplot(1, 3, 2);
imagesc(x_arr, z_arr, squeeze(grid_fluence(:, slice_y, :))');
axis equal tight;
colorbar;
xlabel('X (mm)');
ylabel('Z (mm)');
title('Fluence Rate in XZ-plane');
set(gca, 'YDir', 'normal'); % Adjust to correct orientation

% Plot the fluence rate in the YZ-plane
subplot(1, 3, 3);
imagesc(y_arr, z_arr, squeeze(grid_fluence(slice_x, :, :))');
axis equal tight;
colorbar;
xlabel('Y (mm)');
ylabel('Z (mm)');
title('Fluence Rate in YZ-plane');
set(gca, 'YDir', 'normal'); % Adjust to correct orientation

% Adjust the figure layout
sgtitle('Fluence Rate in Different Planes');


%% Microscale simulations 
% % % xarr_q = linspace(min(X(:)),max(X(:)),3*size(X,1));
% % % Yarr_q = linspace(min(Y(:)),max(Y(:)),3*size(Y,1));
% % % Zarr_q = linspace(min(Z(:)),max(Z(:)),3*size(Z,1));
% % % [Xq,Yq,Zq]=meshgrid(xarr_q,Yarr_q,Zarr_q );
% % % absorbed_energy = interp3(X,Y,Z,absorbed_energy,Xq,Yq,Zq,"linear");
% % % clot_mask = interp3(X,Y,Z,double(clot_mask),Xq,Yq,Zq,"linear"); 
% % % Nx = size(absorbed_energy,1);
% % % Ny = size(absorbed_energy,2);
% % % Nz = size(absorbed_energy,3);
% % % dx = voxel_size/3;
% % % dy = voxel_size/3;
% % % dz = voxel_size/3;
% % % kgrid = kWaveGrid(Nx, dx, Ny, dy,Nz, dz);
% % % gruneisen_parameter = 0.16;
% % % source.p0 = gruneisen_parameter .* absorbed_energy;  % [Pa]
% % % medium.sound_speed = 1540*ones(Nx,Ny,Nz);    % [m/s]
% % % medium.sound_speed(clot_mask == 1) = 2380;   % [m/s]
% % % % % % medium.sound_speed(clot_mask == 2) = 2380;   % [m/s]
% % % % % % medium.sound_speed(clot_mask == 3) = 2380;   % [m/s]
% % % medium.density = 1025*ones(Nx, Ny,Nz);        % [kg/m^3]
% % % medium.density(clot_mask == 1) = 2000;
% % % % % % medium.density(clot_mask == 2) = 2000;
% % % % % % medium.density(clot_mask == 3) = 2000;
% % % T = 1e-5;
% % % dt = 5e-10;
% % % T_array = 0:dt:T;
% % % kgrid.t_array = T_array;
% % % sensor.mask = zeros(Nx, Ny, Nz);
% % % N = round(T / dt) + 1; % Number of time points
% % % detector_size = 200e-6; %  pixelsize of verasonic device or size of the hydrophone head
% % % Num_sensors_y = Ny-2;
% % % Num_sensors_x = detector_size/dx;
% % % sensor_angle = 45;
% % % center_sensor = [0.9*Nx,Ny/2,0.9*Nz];   % position coordinates
% % % detector_size = 200e-6; %  pixelsize of verasonic device or size of the hydrophone head
% % % Num_detectors = Num_sensors_y*dy/detector_size;
% % % % Define rotation matrix for the plane
% % % % Define rotation matrix for 45-degree rotation around the y-axis
% % % sensor_angle_rad = deg2rad(sensor_angle); % Convert 45 degrees to radians
% % % R_y = [cos(sensor_angle_rad), 0, sin(sensor_angle_rad);  % Rotation matrix around y-axis
% % %        0,                  1, 0;
% % %       -sin(sensor_angle_rad), 0, cos(sensor_angle_rad)];
% % % 
% % % % Define the grid for the sensor array (long dimension along y-axis)
% % % [y_sensor, x_sensor] = meshgrid(linspace(-Num_sensors_y/2, Num_sensors_y/2, Num_sensors_y), ...
% % %                                 linspace(-Num_sensors_x/2, Num_sensors_x/2, Num_sensors_x));
% % % z_sensor = zeros(size(x_sensor)); % Plane is initially in the YX plane
% % % sensor_grid = [x_sensor(:), y_sensor(:), z_sensor(:)]'; % Flatten and combine
% % % 
% % % % Rotate the sensor grid around the y-axis
% % % rotated_sensor_grid = R_y * sensor_grid;
% % % 
% % % % Translate the rotated sensor grid to the center position
% % % rotated_sensor_grid(1, :) = rotated_sensor_grid(1, :) + center_sensor(1); % X-coordinates
% % % rotated_sensor_grid(2, :) = rotated_sensor_grid(2, :) + center_sensor(2); % Y-coordinates
% % % rotated_sensor_grid(3, :) = rotated_sensor_grid(3, :) + center_sensor(3); % Z-coordinates
% % % 
% % % % Map the rotated plane to the 3D mask
% % % rotated_x = round(rotated_sensor_grid(1, :));
% % % rotated_y = round(rotated_sensor_grid(2, :));
% % % rotated_z = round(rotated_sensor_grid(3, :));
% % % 
% % % % Ensure indices are within bounds
% % % valid_idx = rotated_x > 0 & rotated_x <= Nx & ...
% % %             rotated_y > 0 & rotated_y <= Ny & ...
% % %             rotated_z > 0 & rotated_z <= Nz;
% % % 
% % % % Assign mask values
% % % sensor.mask(sub2ind(size(sensor.mask), ...
% % %             rotated_x(valid_idx), ...
% % %             rotated_y(valid_idx), ...
% % %             rotated_z(valid_idx))) = 1;

% % % sensor_data = kspaceFirstOrder3DG(kgrid, medium, source, sensor, 'PMLInside', false);

% % % Assume sensor_data is of size (Num_sensors_y * Num_sensors_x) x data_length
% % data_length = size(sensor_data, 2); % Get the number of time samples (data length)
% % 
% % % Reshape sensor_data into a 3D matrix
% % sensor_data_3D = reshape(sensor_data, [Num_sensors_y, Num_sensors_x, data_length]);
% % 
% % % Parameters
% % detector_size_experimental = round(200e-6 / dx); % Detector width in pixels
% % sigma = detector_size_experimental / 2;         % Standard deviation for Gaussian kernel
% % 
% % % Generate a 2D Gaussian kernel
% % kernel_size = 2 * ceil(3 * sigma) + 1; % Kernel size based on 3-sigma rule
% % [x_kernel, y_kernel] = meshgrid(linspace(-3*sigma, 3*sigma, kernel_size));
% % gaussian_kernel = exp(-(x_kernel.^2 + y_kernel.^2) / (2 * sigma^2));
% % gaussian_kernel = gaussian_kernel / sum(gaussian_kernel(:)); % Normalize
% % 
% % % Preallocate for the downsampled 1D array
% % [data_length, ~] = size(sensor_data_3D);
% % downsampled_data = zeros(data_length, 1);
% % 
% % % Apply Gaussian integration and downsample
% % for t = 1:data_length
% %     % Extract the 2D spatial slice at time t
% %     sensor_slice = sensor_data_3D(:, :, t);
% % 
% %     % Convolve with the Gaussian kernel
% %     smoothed_slice = conv2(sensor_slice, gaussian_kernel, 'same');
% % 
% %     % Downsample the smoothed data to the experimental detector size
% %     step_size = detector_size_experimental; % Spacing in pixels
% %     smoothed_subsampled = smoothed_slice(1:step_size:end, 1:step_size:end);
% % 
% %     % Integrate over all downsampled pixels for this time slice
% %     downsampled_data(t) = sum(smoothed_subsampled(:));
% % end
% % % Define parameters for FFT
% % T = 1e-5;              % Total simulation time
% % dt = 5e-10;            % Time step
% % Fs = 1 / dt;           % Sampling frequency
% % N = round(T / dt) + 1; % Number of time points
% % T_array = 0:dt:T;      % Time array
% % f = Fs * (0:(N/2)) / N; % Frequency range
% % 
% % % Original frequency response in MHz and dB
% % frequency_mhz = [2.076, 2.16, 2.295, 2.532, 2.734, 2.97, 3.241, 3.376, 3.612, 4.253, ...
% %                  4.861, 5.367, 5.873, 6.447, 7.021, 7.595, 8.034, 8.641, 9.046, 9.553, ...
% %                  9.823, 10.295, 10.599, 10.937, 11.342, 11.679, 11.882, 12.152, ...
% %                  12.287, 12.489, 12.692, 12.861, 13.165, 13.3, 13.257, 13.907, ...
% %                  13.975, 14.38, 14.549, 14.717, 14.92, 15.19, 15.325, 15.426, 15.544];
% % magnitude_db = [-39.392, -36.96, -34.286, -29.422, -24.681, -19.818, -15.684, ...
% %                 -11.307, -7.538, -4.498, -4.742, -4.498, -3.769, -2.553, -1.824, ...
% %                 -1.337, -0.729, 0.122, -0.122, -0.973, -1.824, -3.04, -4.498, ...
% %                 -6.444, -8.632, -10.456, -12.523, -14.954, -17.143, -19.331, ...
% %                 -21.763, -24.195, -26.626, -28.936, -31.246, -32.948, -33.678, ...
% %                 -34.407, -34.772, -34.529, -34.407, -35.258, -36.474, -38.176, -39.392];
% % 
% % % Convert frequency to Hz and normalize frequency response
% % frequency_hz = frequency_mhz * 1e6;
% % magnitude_linear = 10.^(magnitude_db / 20); % Convert dB to linear scale
% % magnitude_normalized = magnitude_linear / max(magnitude_linear); % Normalize to max = 1
% % 
% % % Interpolate magnitude values to match the frequency range
% % interpolated_magnitude = interp1(frequency_hz, magnitude_normalized, f, 'linear', 'extrap');
% % interpolated_magnitude(interpolated_magnitude < 0) = 0;
% % % Parameters for large pixels
% % large_pixel_size_x = round(detector_size_experimental); % Large pixel size in X (in small pixels)
% % large_pixel_size_y = round(detector_size_experimental); % Large pixel size in Y (in small pixels)
% % 
% % % Determine number of large pixels
% % num_large_pixels_x = floor(size(sensor_data_3D, 2) / large_pixel_size_x);
% % num_large_pixels_y = floor(size(sensor_data_3D, 1) / large_pixel_size_y);
% % 
% % % Preallocate for large pixel time-domain and spectral-domain data
% % large_pixel_time_data = zeros(num_large_pixels_y, num_large_pixels_x, N);
% % large_pixel_spectral_data = zeros(num_large_pixels_y, num_large_pixels_x, floor(N/2+1));
% % 
% % % 2D Gaussian kernel for integration
% % sigma = large_pixel_size_x / 2; % Standard deviation proportional to large pixel size
% % kernel_size = 2 * ceil(3 * sigma) + 1; % Kernel size (3-sigma rule)
% % [x_kernel, y_kernel] = meshgrid(linspace(-3*sigma, 3*sigma, kernel_size));
% % gaussian_kernel = exp(-(x_kernel.^2 + y_kernel.^2) / (2 * sigma^2));
% % gaussian_kernel = gaussian_kernel / sum(gaussian_kernel(:)); % Normalize
% % 
% % % Loop over large pixels
% % for large_y = 1:num_large_pixels_y
% %     for large_x = 1:num_large_pixels_x
% %         % Extract small pixels within the current large pixel
% %         y_range = (large_y-1)*large_pixel_size_y+1 : large_y*large_pixel_size_y;
% %         x_range = (large_x-1)*large_pixel_size_x+1 : large_x*large_pixel_size_x;
% %         small_pixel_data = sensor_data_3D(y_range, x_range, :);
% % 
% %         % Integrate over small pixels using the Gaussian kernel
% %         integrated_signal = zeros(N, 1);
% %         for t = 1:N
% %             % Extract spatial slice for time t
% %             spatial_slice = squeeze(small_pixel_data(:, :, t));
% %             % Convolve with the Gaussian kernel
% %             smoothed_slice = conv2(spatial_slice, gaussian_kernel, 'same');
% %             % Sum the smoothed values
% %             integrated_signal(t) = sum(smoothed_slice(:));
% %         end
% % 
% %         % Time-domain data for the current large pixel
% %         large_pixel_time_data(large_y, large_x, :) = integrated_signal;
% % 
% %         % Perform FFT on the integrated signal
% %         integrated_fft = fft(integrated_signal);
% % 
% %         % Apply frequency response filter
% %         filtered_fft = integrated_fft;
% %         filtered_fft(1:N/2+1) = filtered_fft(1:N/2+1) .* interpolated_magnitude';
% %         filtered_fft(N/2+2:end) = conj(filtered_fft(N/2:-1:2)); % Mirror symmetry
% % 
% %         % Store spectral-domain data (magnitude of filtered FFT)
% %         large_pixel_spectral_data(large_y, large_x, :) = abs(filtered_fft(1:N/2+1));
% % 
% %         % Inverse FFT to get filtered time-domain signal
% %         filtered_signal = ifft(filtered_fft, 'symmetric');
% % 
% %         % Update the time-domain data with the filtered signal
% %         large_pixel_time_data(large_y, large_x, :) = filtered_signal;
% %     end
% % end
% % 
% % % Reshape time-domain data
% % large_pixel_time_data_2D = reshape(large_pixel_time_data, size(large_pixel_time_data, 1), []);
% % 
% % % Reshape spectral-domain data
% % large_pixel_spectral_data_2D = reshape(large_pixel_spectral_data, size(large_pixel_spectral_data, 1), []);
% % 
% % 
% % % Parameters for kspaceLineRecon
% % dy = detector_size_experimental; % Step size in y-direction (spacing between large pixels)
% % dt = 5e-10;                      % Time step
% % c = 1540;                        % Speed of sound in medium (m/s)
% % 
% % % Ensure the data is in the correct order ('ty' by default in kspaceLineRecon)
% % % If your data is not in the default 'ty' order, transpose it:
% % large_pixel_time_data_2D_corrected = large_pixel_time_data_2D.'; % Transpose if needed
% % 
% % % Perform reconstruction using only valid optional inputs
% % reconstructed_image = kspaceLineRecon(large_pixel_time_data_2D_corrected, dy, dt, c, ...
% %                                       'Interp', '*nearest', 'Plot', true, 'PosCond', true);
% % 
% % % Visualize reconstructed image
% % figure;
% % imagesc(reconstructed_image);
% % colormap('hot');
% % colorbar;
% % title('Reconstructed Image');
% % xlabel('Sensor y (pixels)');
% % ylabel('Depth z (pixels)');
% % 
% % % Combine spectral data from all large pixels
% % % Flatten the 3D spectral data to a 2D matrix: (num_large_pixels_y * num_large_pixels_x) x (num_frequency_bins)
% % num_pixels_y = size(large_pixel_spectral_data, 1);
% % num_pixels_x = size(large_pixel_spectral_data, 2);
% % num_frequency_bins = size(large_pixel_spectral_data, 3);
% % 
% % spectral_data_flat = reshape(large_pixel_spectral_data, [], num_frequency_bins);
% % 
% % % Perform PCA on the spectral data
% % [coeff, score, latent] = pca(spectral_data_flat);
% % 
% % % Define frequency array
% % frequency_values = linspace(f(1), f(end), num_frequency_bins); % Frequency corresponding to spectral bins
% % 
% % % Plot the variance explained by each principal component
% % figure;
% % explained_variance = 100 * latent / sum(latent); % Percent variance explained
% % bar(explained_variance);
% % xlabel('Principal Component');
% % ylabel('Variance Explained (%)');
% % title('PCA Variance Explained');
% % grid on;
% % 
% % % Plot the spectral profile for the first principal component (coefficients)
% % figure;
% % plot(frequency_values, coeff(:, 1)); % Coefficients (weights) vs frequency
% % xlabel('Frequency (Hz)');
% % ylabel('PCA Weight');
% % title('PCA Weights for the First Principal Component');
% % grid on;
% % 
% % % Plot the variation of the first principal component (scores) across pixels
% % pixel_indices = 1:size(score, 1); % Pixel indices
% % figure;
% % plot(pixel_indices, score(:, 1)); % Scores vs pixel index
% % xlabel('Pixel Index');
% % ylabel('Amplitude');
% % title('First Principal Component (Variation Across Pixels)');
% % grid on;
