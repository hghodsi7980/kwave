% Generate some sample data (replace this with your actual data)
data = 1-Porosity;

% Define the number of bins and bin width for the histogram
numBins = 100;
binWidth = (max(data)-min(data))/numBins;

% Plot the histogram
histogram(data, numBins);

% Compute mean and standard deviation of the data
mu = mean(data);
sigma = std(data);

% Generate the normal distribution curve
x = linspace(min(data), max(data), 1000); % Generate x values
y = 1* exp(-(x - mu).^2 / (2 * sigma^2)); % Compute corresponding y values

% Normalize the curve to match the histogram
binCounts = histcounts(data, numBins);
y = y * max(binCounts); 

% Plot the normal curve
hold on;
plot(x, y, 'r', 'LineWidth', 2);
hold off;

% Add legend and labels
legend('Histogram', 'Fitted Normal Curve');
xlabel('Value');
ylabel('Frequency');
title('Porosity ratio distribution');

% Display mean and standard deviation
disp(['Mean: ', num2str(mu)]);
disp(['Standard Deviation: ', num2str(sigma)]);
