% Sensitivity Analysis: Remove inputs systematically and observe performance
num_inputs = size(inputs, 1);  % Number of input features
sensitivity_R2_results = zeros(num_inputs, 1);  % Store R² results for each input removal
sensitivity_MSE_results = zeros(num_inputs, 1);  % Store MSE results for each input removal

% Baseline: Train with all inputs to get baseline performance
fprintf('\nBaseline Performance with All Inputs (Output 3)\n');
[baseline_R2, baseline_MSE] = train_with_optimized_config(inputs, targets(3, :), [24, 23, 7]);

% Sensitivity Analysis: Remove one input at a time
for i = 1:num_inputs
    fprintf('\nTesting without Input %d\n', i);
    
    % Remove the ith input
    reduced_inputs = inputs;
    reduced_inputs(i, :) = [];  % Remove the ith row (input feature)

    % Train the network with the reduced inputs using the optimized configuration
    [bestR2, bestMSE] = train_with_optimized_config(reduced_inputs, targets(3, :), [24, 23, 7]);

    % Store the R² and MSE results for this input removal
    sensitivity_R2_results(i) = bestR2;
    sensitivity_MSE_results(i) = bestMSE;
    
    % Display results for this configuration
    fprintf('Removed Input %d, R² = %.4f, MSE = %.4f\n', i, bestR2, bestMSE);
end




% --- Plotting MSE and R² Results ---
input_labels = 1:num_inputs;  % Labels for each input feature

% Plot R² values
figure;
hold on;
plot(input_labels, sensitivity_R2_results, '-o', 'DisplayName', 'Sensitivity Analysis (R²)');
yline(baseline_R2, 'r--', 'Baseline R²', 'LineWidth', 2, 'DisplayName', 'Baseline R²');
xlabel('Input Feature Removed');
ylabel('R²');
title('R² Values for Sensitivity Analysis (Output 3)');
legend('show');
grid on;
hold off;

% Plot MSE values
figure;
hold on;
plot(input_labels, sensitivity_MSE_results, '-o', 'DisplayName', 'Sensitivity Analysis (MSE)');
yline(baseline_MSE, 'r--', 'Baseline MSE', 'LineWidth', 2, 'DisplayName', 'Baseline MSE');
xlabel('Input Feature Removed');
ylabel('MSE');
title('MSE Values for Sensitivity Analysis (Output 3)');
legend('show');
grid on;
hold off;

% Function to calculate R²
function r2 = calculateR2(trueVals, predVals)
    r2 = 1 - sum((trueVals - predVals).^2) / sum((trueVals - mean(trueVals)).^2);
end


% --- Function to Train with Optimized Configuration ---
function [bestR2, bestMSE] = train_with_optimized_config(inputs, target, hiddenLayerSize)
    % Set up data division for training, validation, and testing
    trainRatio = 0.7;
    valRatio = 0.15;
    testRatio = 0.15;

    % Create and configure the network using the optimized hidden layer sizes
    net = feedforwardnet(hiddenLayerSize);
    net.performFcn = 'mse';  % Mean Squared Error

    % Disable the training window and command-line output for speed
    net.trainParam.showWindow = false;
    net.trainParam.showCommandLine = false;

    % Set up data division
    net.divideParam.trainRatio = trainRatio;
    net.divideParam.valRatio = valRatio;
    net.divideParam.testRatio = testRatio;

    % Train the network
    [net, tr] = train(net, inputs, target);

    % Test the network
    outputs = net(inputs);

    % Calculate R² and MSE (on test set)
    trueVals = target(tr.testInd);
    predVals = outputs(tr.testInd);
    bestR2 = calculateR2(trueVals, predVals);
    bestMSE = mean((trueVals - predVals).^2);
end
