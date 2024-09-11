% Load the input and output data
close all;
clear;
load out.mat;
load in.mat;

% inputs = input_data';   % Transpose to fit the required format [features x samples]
inputs = input_data(:, [1:10, 13])';  % Select inputs 1 to 10 and 13, excluding 11 and 12, then transpose
targets = output_data'; % Transpose to fit the required format [outputs x samples]

% --- Scramble the dataset ---
num_samples = size(inputs, 2);  % Number of samples
shuffle_indices = randperm(num_samples);  % Create a random permutation of indices

% Apply the shuffle to both inputs and targets
inputs = inputs(:, shuffle_indices);
targets = targets(:, shuffle_indices);

% --- Standardization (Z-score normalization) of the inputs ---
% Calculate mean and standard deviation for each input feature
input_mean = mean(inputs, 2);  % Mean for each feature (row)
input_std = std(inputs, 0, 2);  % Standard deviation for each feature (row)

% Standardize inputs
inputs = (inputs - input_mean) ./ input_std;

% --- Normalization of the targets ---
% Normalize the target values (for each output variable)
max_val1 = max(targets(1, :));
min_val1 = min(targets(1, :));
max_val2 = max(targets(2, :));
min_val2 = min(targets(2, :));

% If you have a third output, normalize it as well
if size(targets, 1) > 2
    max_val3 = max(targets(3, :));
    min_val3 = min(targets(3, :));
end

% Normalize targets to the range [0, 1]
targets(1, :) = (targets(1, :) - min_val1) / (max_val1 - min_val1);
targets(2, :) = (targets(2, :) - min_val2) / (max_val2 - min_val2);

if size(targets, 1) > 2
    targets(3, :) = (targets(3, :) - min_val3) / (max_val3 - min_val3);
end

% --- Define ranges for hidden layer sizes ---
layer1_sizes = 5:5:30;   % Sizes for the first hidden layer
layer2_sizes = 5:5:30;   % Sizes for the second hidden layer
layer3_sizes = 5:5:30;   % Sizes for the third hidden layer (if applicable)

% Set up the data division for training, validation, and testing
trainRatio = 0.7;
valRatio = 0.15;
testRatio = 0.15;

% Function to calculate R²
calculateR2 = @(trueVals, predVals) 1 - sum((trueVals - predVals).^2) / sum((trueVals - mean(trueVals)).^2);

% --- Train separate networks for each output ---

% --- Optimize for Output 1 (Porosity) ---
fprintf('\nOptimizing Network for Output 1 (Porosity)\n');
bestR2_1 = -Inf;
for size1 = layer1_sizes
    for size2 = layer2_sizes
        for size3 = layer3_sizes
            % Create the hidden layer configuration
            hiddenLayerSize = [size1 size2 size3];

            % Create and configure the network for Output 1
            net = feedforwardnet(hiddenLayerSize);
            net.performFcn = 'mse';  % Mean Squared Error

            % Disable the training window and command-line output for speed
            net.trainParam.showWindow = false;
            net.trainParam.showCommandLine = false;

            % Set up data division
            net.divideParam.trainRatio = trainRatio;
            net.divideParam.valRatio = valRatio;
            net.divideParam.testRatio = testRatio;

            % Train the network for Output 1
            [net, tr] = train(net, inputs, targets(1, :));  % Train on Output 1 only

            % Test the network (predicted outputs for Output 1)
            outputs_1 = net(inputs);

            % Calculate R² for Output 1 (Normalized)
            r2_1 = calculateR2(targets(1, tr.testInd), outputs_1(tr.testInd));

            % Check if this is the best configuration based on R²
            if r2_1 > bestR2_1
                bestR2_1 = r2_1;
                bestConfig_1 = hiddenLayerSize;
                best_net_1 = net; % Save the best network
                best_tr_1 = tr;   % Save the training record
                best_outputs_1 = outputs_1;  % Save best outputs
            end

            % Display current configuration and R² value
            fprintf('Tested Configuration for Output 1: [%d, %d, %d], R² = %.4f\n', size1, size2, size3, r2_1);
        end
    end
end

% --- Optimize for Output 2 (Composition) ---
fprintf('\nOptimizing Network for Output 2 (Composition)\n');
bestR2_2 = -Inf;
for size1 = layer1_sizes
    for size2 = layer2_sizes
        for size3 = layer3_sizes
            % Create the hidden layer configuration
            hiddenLayerSize = [size1, size2, size3];

            % Create and configure the network for Output 2
            net = feedforwardnet(hiddenLayerSize);
            net.performFcn = 'mse';  % Mean Squared Error

            % Disable the training window and command-line output for speed
            net.trainParam.showWindow = false;
            net.trainParam.showCommandLine = false;

            % Set up data division
            net.divideParam.trainRatio = trainRatio;
            net.divideParam.valRatio = valRatio;
            net.divideParam.testRatio = testRatio;

            % Train the network for Output 2
            [net, tr] = train(net, inputs, targets(2, :));  % Train on Output 2 only

            % Test the network (predicted outputs for Output 2)
            outputs_2 = net(inputs);

            % Calculate R² for Output 2 (Normalized)
            r2_2 = calculateR2(targets(2, tr.testInd), outputs_2(tr.testInd));

            % Check if this is the best configuration based on R²
            if r2_2 > bestR2_2
                bestR2_2 = r2_2;
                bestConfig_2 = hiddenLayerSize;
                best_net_2 = net; % Save the best network
                best_tr_2 = tr;   % Save the training record
                best_outputs_2 = outputs_2;  % Save best outputs
            end

            % Display current configuration and R² value
            fprintf('Tested Configuration for Output 2: [%d, %d, %d], R² = %.4f\n', size1, size2, size3, r2_2);
        end
    end
end

% --- Optimize for Output 3 (if applicable) ---
if size(targets, 1) > 2
    fprintf('\nOptimizing Network for Output 3\n');
    bestR2_3 = -Inf;
    for size1 = layer1_sizes
        for size2 = layer2_sizes
            for size3 = layer3_sizes
                % Create the hidden layer configuration
                hiddenLayerSize = [size1, size2, size3];

                % Create and configure the network for Output 3
                net = feedforwardnet(hiddenLayerSize);
                net.performFcn = 'mse';  % Mean Squared Error

                % Disable the training window and command-line output for speed
                net.trainParam.showWindow = false;
                net.trainParam.showCommandLine = false;

                % Set up data division
                net.divideParam.trainRatio = trainRatio;
                net.divideParam.valRatio = valRatio;
                net.divideParam.testRatio = testRatio;

                % Train the network for Output 3
                [net, tr] = train(net, inputs, targets(3, :));  % Train on Output 3 only

                % Test the network (predicted outputs for Output 3)
                outputs_3 = net(inputs);

                % Calculate R² for Output 3 (Normalized)
                r2_3 = calculateR2(targets(3, tr.testInd), outputs_3(tr.testInd));

                % Check if this is the best configuration based on R²
                if r2_3 > bestR2_3
                    bestR2_3 = r2_3;
                    bestConfig_3 = hiddenLayerSize;
                    best_net_3 = net; % Save the best network
                    best_tr_3 = tr;   % Save the training record
                    best_outputs_3 = outputs_3;  % Save best outputs
                end

                % Display current configuration and R² value
                fprintf('Tested Configuration for Output 3: [%d, %d, %d], R² = %.4f\n', size1, size2, size3, r2_3);
            end
        end
    end
end

% --- Display the best configuration and R² values for each output ---
fprintf('\nBest Configuration for Output 1: [%d, %d, %d], Best R²: %.4f\n', bestConfig_1(1), bestConfig_1(2), bestConfig_1(3), bestR2_1);
fprintf('Best Configuration for Output 2: [%d, %d, %d], Best R²: %.4f\n', bestConfig_2(1), bestConfig_2(2), bestConfig_2(3), bestR2_2);
if size(targets, 1) > 2
    fprintf('Best Configuration for Output 3: [%d, %d, %d], Best R²: %.4f\n', bestConfig_3(1), bestConfig_3(2), bestConfig_3(3), bestR2_3);
end

% --- Denormalize predictions for reporting if needed ---
denorm_outputs_1 = best_outputs_1 * (max_val1 - min_val1) + min_val1;
denorm_outputs_2 = best_outputs_2 * (max_val2 - min_val2) + min_val2;

if size(targets, 1) > 2
    denorm_outputs_3 = best_outputs_3 * (max_val3 - min_val3) + min_val3;
end

% --- Plot regression and error histogram for the best networks (Normalized) ---

% Regression Plot for Output 1 (Test Data Only)
custom_plot_regression(targets(1, best_tr_1.testInd), best_outputs_1(best_tr_1.testInd), 'Regression for Output 1 (Test Data Only)');
figure;
ploterrhist(gsubtract(targets(1, best_tr_1.testInd), best_outputs_1(best_tr_1.testInd)));
title('Error Histogram for Output 1 (Test Data Only)');

% Regression Plot for Output 2 (Test Data Only)
custom_plot_regression(targets(2, best_tr_2.testInd), best_outputs_2(best_tr_2.testInd), 'Regression for Output 2 (Test Data Only)');
figure;
ploterrhist(gsubtract(targets(2, best_tr_2.testInd), best_outputs_2(best_tr_2.testInd)));
title('Error Histogram for Output 2 (Test Data Only)');

% Regression Plot for Output 3 (Test Data Only)
if size(targets, 1) > 2
    custom_plot_regression(targets(3, best_tr_3.testInd), best_outputs_3(best_tr_3.testInd), 'Regression for Output 3 (Test Data Only)');
    figure;
    ploterrhist(gsubtract(targets(3, best_tr_3.testInd), best_outputs_3(best_tr_3.testInd)));
    title('Error Histogram for Output 3 (Test Data Only)');
end



% Custom function to plot regression without the fitted line
function custom_plot_regression(true_vals, predicted_vals, title_str)
    % Create the figure
    figure;
    
    % Plot the data points
    scatter(true_vals, predicted_vals, 'filled');
    hold on;
    
    % Plot the y = x line
    max_val = max([true_vals(:); predicted_vals(:)]);
    min_val = min([true_vals(:); predicted_vals(:)]);
    plot([min_val, max_val], [min_val, max_val], 'r--', 'LineWidth', 2); % Plot x = y line
    
    % Set labels and title
    xlabel('True Values');
    ylabel('Predicted Values');
    title(title_str);
    
    % Axis limits for better comparison
    axis([min_val max_val min_val max_val]);
    grid on;
    hold off;
end
