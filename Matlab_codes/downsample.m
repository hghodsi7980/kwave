% ==========================================================
% Assign Parameters and Downsample ClotMatrix
% ==========================================================

clear; clc; close all;

% Load the pre-existing ClotMatrix
load('CroppedClotMatrix_1.mat', 'ClotMatrix');

% Define parameter values for each type
param_empty = 0;  % Value for empty space (0)
param_rbc = 1;   % Value for RBCs (1)
param_fibrin = 2; % Value for Fibrin (2)
param_platelet = 3; % Value for Platelets (3)

% Create a parameter matrix based on ClotMatrix
ParameterMatrix = zeros(size(ClotMatrix));
ParameterMatrix(ClotMatrix == 0) = param_empty;
ParameterMatrix(ClotMatrix == 1) = param_rbc;
ParameterMatrix(ClotMatrix == 2) = param_fibrin;
ParameterMatrix(ClotMatrix == 3) = param_platelet;

% ==========================================================
% Downsample the Matrix
% ==========================================================
downsample_factor = 1; % Factor forv downsampling
[rows, cols, slices] = size(ParameterMatrix);

% Ensure dimensions are divisible by the downsample factor
rows_ds = floor(rows / downsample_factor);
cols_ds = floor(cols / downsample_factor);
slices_ds = floor(slices / downsample_factor);

% Initialize downsampled matrix
DownsampledMatrix = zeros(rows_ds, cols_ds, slices_ds);

% Downsample by averaging within blocks
for i = 1:rows_ds
    for j = 1:cols_ds
        for k = 1:slices_ds
            block = ParameterMatrix((i-1)*downsample_factor+1:i*downsample_factor, ...
                                     (j-1)*downsample_factor+1:j*downsample_factor, ...
                                     (k-1)*downsample_factor+1:k*downsample_factor);
            DownsampledMatrix(i, j, k) = mean(block(:));
        end
    end
end

% ==========================================================
% Visualize the Downsampled Matrix using Slice
% ==========================================================
[x, y, z] = ind2sub(size(DownsampledMatrix), find(DownsampledMatrix ~= 0));
values = DownsampledMatrix(DownsampledMatrix ~= 0);
% Plot the non-zero elements as points in 3D space
figure;
scatter3(x, y, z, 10, values, 'filled'); % 50 is the size of points, values for color
colormap(jet); % Color map for visualization
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Non-Zero Elements of the Matrix in 3D');
grid on;
% ==========================================================
% Save the Downsampled Matrix
% ==========================================================
save('DownsampledClotMatrix.mat', 'DownsampledMatrix');
