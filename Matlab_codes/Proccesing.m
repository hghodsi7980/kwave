feature_size = 18;
num_vib_modes = 5;
output_size = 2;
% Define the folder containing the .mat files
folder = 'C:\Git\kwave\Simulation_data\'; % Update this with your folder path

% Get a list of all .mat files in the folder
fileList = dir(fullfile(folder, '*.mat'));

% Loop through each file

for i = 1:length(fileList)
    keepvars = {'feature_struct', 'fileList' , 'folder' , 'num_vib_modes' , 'feature_size', 'i', 'output_size' , 'output_struct' , 'outputs'};
    clearvars('-except', keepvars{:});
    fprintf('Progress percent: %f\n',i*100/length(fileList));
    % Get the filename
    filename = fullfile(folder, fileList(i).name);
    feature_struct(i).filename = filename;
    % Display name of the file
    fprintf('Loaded file: %s\n', filename);
    % Load the .mat file
    load(filename);
    Nx = size(output_mask_clot,1);
    Ny = size(output_mask_clot,2);
    frequencies = frequencies(1:150);
    Spectrum = fft_r(Nx/2,1:150);
    Spectrum = Spectrum/sum(Spectrum);
    [peaks, indexes] = findpeaks(Spectrum);
    freqs = frequencies(indexes);
    feature_struct(i).feature = char(zeros(feature_size,20)); %#ok<*SAGROW>
    feature_struct(i).value = zeros(1,feature_size);
    for j = 1:num_vib_modes
        feature_struct(i).feature(j,1:length(sprintf('mode_freq_%d',j))) = sprintf('mode_freq_%d',j);
        feature_struct(i).value(j) = freqs(j);
        feature_struct(i).feature(j+num_vib_modes,1:length(sprintf('mode_value_%d',j))) = sprintf('mode_value_%d',j);
        feature_struct(i).value(j+num_vib_modes) = peaks(j);
        temp_array_spec = Spectrum(indexes(j):indexes(j+1));
        temp_array_freq = frequencies(indexes(j):indexes(j+1));
        Min_value = min(temp_array_spec);
        Freq_min = temp_array_freq(temp_array_spec == Min_value);
        FWHM = 2*(Freq_min-freqs(j));
        feature_struct(i).feature(j+2*num_vib_modes,1:length(sprintf('Q_factor_%d',j))) = sprintf('Q_factor_%d',j);
        feature_struct(i).value(j+2*num_vib_modes) = freqs(j)/FWHM;
    end
    f=fit(freqs',peaks','poly1');
    feature_struct(i).feature(3*num_vib_modes+1,1:length('Slope')) = 'Slope';
    feature_struct(i).value(3*num_vib_modes+1) = f.p1;
    temp_matrix = output_mask_clot;
    temp_matrix(temp_matrix ~=0 ) = 1;
    feature_struct(i).feature(3*num_vib_modes+2,1:length('Size')) = 'Size';
    feature_struct(i).value(3*num_vib_modes+2) = clot_radious;
    feature_struct(i).feature(3*num_vib_modes+3,1:length('Freq_shift')) = 'Freq_shift';
    frequency_shift = sum(frequencies.*Spectrum)/sum(frequencies);
    feature_struct(i).value(3*num_vib_modes+3) =  frequency_shift;
    % Outputs 
    output_struct(i).parameter(1,1:length('Porosity')) = 'Porosity';
    output_struct(i).parameter(2,1:length('Fibrin_ratio')) = 'Fibrin_ratio';
    porosity = 100*(1-(sum(sum(temp_matrix))/(pi*(clot_radious/dx)^2)));
    if porosity > 100
        porosity = 100;
    end
    if porosity < 0
        porosity = 0;
    end
    fibrin = sum(sum(output_mask_clot(output_mask_clot == 2)));
    rbc = sum(sum(output_mask_clot(output_mask_clot == 1)));
    fibrin_ratio = 100*fibrin/(fibrin+rbc);
    output_struct(i).value(1) = porosity;
    disp(fibrin_ratio);
    disp(porosity);
    output_struct(i).value(2) = fibrin_ratio;
end
