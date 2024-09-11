% Clear workspace and close all figures
clear;

% Define folder paths for data and results (update these with your actual paths)
folder_data = 'C:\Users\hghodsi\OneDrive - Delft University of Technology\dataset\Database\data2500';
folder_result = 'C:\Users\hghodsi\OneDrive - Delft University of Technology\dataset\Database\result2500';

% Add the paths to MATLAB's search path
addpath(genpath(folder_result));
addpath(genpath(folder_data));

% Get the list of all .mat files in the data and result folders
fileList1 = dir(fullfile(folder_data, '*.mat'));
fileList2 = dir(fullfile(folder_result, '*.mat'));

% Extract numeric indices from filenames for matching data and results
index1 = zeros(1, length(fileList1));
for i = 1:length(fileList1)
    index1(i) = str2double(cell2mat(regexp(fileList1(i).name, '\d+', 'match')));
end

index2 = zeros(1, length(fileList2));
for i = 1:length(fileList2)
    index2(i) = str2double(cell2mat(regexp(fileList2(i).name, '\d+', 'match')));
end

% Determine the dataset size based on the number of result files
dataset_size = length(fileList2);

% Initialize arrays to store inputs and outputs for the neural network
input_data = [];
output_data = [];

% Loop through each dataset
for i = 1:dataset_size
    fprintf('Processing dataset %d/%d\n', i, dataset_size);

    % Load the volume data and mesh_space from the corresponding data file
    filename1 = fullfile(folder_data, fileList1(index1 == index2(i)).name);
    load(filename1);  % This loads 'mesh_space' and 'volume'

    % Load the spectrum data from the corresponding result file
    filename2 = fullfile(folder_result, fileList2(i).name);
    load(filename2);  % This loads 'dataset' which contains 'spec_all' and 'frequencies'

    % Calculate volume_clot and normalize the spectrum
    volume_clot(i) = calculate_volume(mesh_space);
    spectrum = dataset.spec_all;
    spectrum = spectrum / sum(spectrum(:));

    % Downsample the spectrum along the first two dimensions
    q = 30;  % Downsampling factor to make the pixel size 30
    [m, n, p] = size(spectrum);
    m_padded = ceil(m / q) * q;
    n_padded = ceil(n / q) * q;
    spectrum_padded = padarray(spectrum, [m_padded - m, n_padded - n, 0], 0, 'post');
    spectrum_reshaped = reshape(spectrum_padded, q, m_padded / q, q, n_padded / q, p);
    spectrum_downsampled = squeeze(mean(mean(spectrum_reshaped, 1), 3));

    % Initialize spectrum_final for this iteration
    spectrum_final = zeros(size(spectrum_downsampled, 1), size(spectrum_downsampled, 2), 10 * length(dataset.frequencies));

    % Smoothing and interpolation on the spectrum
    for x = 1:size(spectrum_downsampled, 1)
        for y = 1:size(spectrum_downsampled, 2)
            a = dataset.frequencies;
            b = reshape(spectrum_downsampled(x, y, :), 1, []);

            % Fit a smoothing spline to the spectrum
            [xData, yData] = prepareCurveData(a, b);
            ft = fittype('smoothingspline');
            opts = fitoptions('Method', 'SmoothingSpline');
            opts.SmoothingParam = 1.17036428844444e-16;
            [fitresult, gof] = fit(xData, yData, ft, opts);

            % Store the smoothed spectrum in the final matrix
            spectrum_final(x, y, :) = fitresult(linspace(0, max(dataset.frequencies), 10 * length(dataset.frequencies)));
        end
    end
    dataset.frequencies = linspace(0, max(dataset.frequencies), 10 * length(dataset.frequencies));
    % Reshape spectrum_final to 2D [detectors, frequencies]
    [num_x, num_y, num_frequencies] = size(spectrum_final);
    reshaped_spectrum = reshape(spectrum_final, num_x * num_y, num_frequencies);  % Each row corresponds to a detector

    % Perform PCA on the reshaped spectrum data for this dataset
    [coeff, score, latent] = pca(reshaped_spectrum);

    % The first principal component (score(:,1)) represents the dominant mode of variation
    representative_spectrum = score(:, 1);

    % Make sure representative_spectrum and frequencies are the same length
    representative_spectrum = resample(representative_spectrum, length(dataset.frequencies), length(representative_spectrum));

    % Detrend the representative spectrum for peak detection
    detrended_spectrum = detrend(representative_spectrum);

    % Find the maximum value of the detrended spectrum
    max_value = max(detrended_spectrum);

    % Set the prominence as a percentage of the maximum value (e.g., 5%)
    prominence_percentage = 0.5; % 50% of the max
    min_prominence = prominence_percentage * max_value;

    % Use findpeaks with the detrended signal to find peak locations
    [~, locs] = findpeaks(detrended_spectrum, dataset.frequencies, ...
        'NPeaks', 5, 'SortStr', 'descend', 'MinPeakProminence', min_prominence);

    % Get the original amplitudes from the non-detrended spectrum
    peaks = representative_spectrum(locs);



    % Store the first 5 frequencies and corresponding amplitudes into input_data(1:10) for this dataset
    input_data(i, 1:5) = locs(1:5);
    input_data(i, 6:10) = peaks(1:5);

    %% Compute the average spectrum across all detectors
    avg_spectrum = mean(reshaped_spectrum, 1);

    %% Fit a Gaussian curve to the average spectrum
    [xData, yData] = prepareCurveData(dataset.frequencies, avg_spectrum);
    gaussFit = fit(xData, yData, 'gauss1');

    % Extract the effective bandwidth (FWHM) and slope (derivative at center)
    effective_bandwidth = abs(gaussFit.c1 * 2 * sqrt(2 * log(2)));  % FWHM
    % Calculate the slope at the center of the Gaussian (center at gaussFit.b1)
    center_frequency = gaussFit.b1;  % Center of the Gaussian
    center_idx = find(abs(dataset.frequencies - center_frequency) == min(abs(dataset.frequencies - center_frequency)), 1);
    slope_at_center = diff(gaussFit(dataset.frequencies(center_idx:center_idx+1))) / diff(dataset.frequencies(center_idx:center_idx+1));
    input_data(i, 11) = slope_at_center;  % Store the slope at the center frequency

    % Store effective bandwidth and slope in input_data
    input_data(i, 12) = effective_bandwidth;

    %store the volume as the last input
    input_data(i,13) = volume_clot(i);


    %% Calculate porosity and RBC/Fibrin+Platelet ratio
    num_RBC(i) = sum(mesh_space(:) == 1);
    num_fibrin(i) = sum(mesh_space(:) == 2);
    num_platelet(i) = sum(mesh_space(:) == 3);

    porosity = 1 - ((num_RBC(i) + num_fibrin(i) + num_platelet(i)) / volume_clot(i));
    RBC_ratio = num_RBC(i) / volume_clot(i);
    Platelet_ratio = (num_fibrin(i) + num_platelet(i)) / volume_clot(i);

    % Store the calculated outputs in output_data
    output_data(i, 1) = porosity; %#ok<*SAGROW>
    output_data(i, 2) = RBC_ratio;
    output_data(i, 3) = Platelet_ratio;
end

% Save the input and output data
save('in.mat', 'input_data');
save('out.mat', 'output_data');
