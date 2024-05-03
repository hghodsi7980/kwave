clear 
close all
% Define the folder containing the .mat files
folder = 'C:\Users\hghodsi\OneDrive - Delft University of Technology\matfiles'; % Update this with your folder path

% Get a list of all .mat files in the folder
fileList = dir(fullfile(folder, '*.mat'));

% Loop through each file
num_data_samples = 0;
dataset_generated = [];
counter = 1;
for i = 1:length(fileList)
    keepvars = {'fileList' , 'folder' , 'dataset_generated', 'num_data_samples', 'i', 'counter'};
    clearvars('-except', keepvars{:});
    fprintf('Progress percent: %f\n',i*100/length(fileList));
    % Get the filename
    filename = fullfile(folder, fileList(i).name);
    % Display name of the file
    fprintf('Loaded file: %s\n', filename);
    % Load the .mat file
    load(filename);
    file_size = size(dataset2,2);
    num_data_samples = num_data_samples+file_size;
    for j = 1:file_size
        dataset_generated(counter).output = dataset2(j).number_of_spheres; %#ok<*SAGROW>
        dataset_generated(counter).input1 = dataset2(j).coordinate_x(1:10);
        dataset_generated(counter).input2 = dataset2(j).coordinate_y(1:10);
        dataset_generated(counter).input3 = dataset2(j).frequencies(dataset2(j).coordinate_f(1:10));
        frequencies = dataset2(j).frequencies;
        fft_1 = dataset2(j).spec_all(:,:,1:length(frequencies));
        fft_1= fft_1/max(fft_1(:));
        Nx = size(fft_1,1);
        Ny = size(fft_1,2);
        test_fft_1 = ones(1,length(frequencies));
        for m = 1:Nx           
            for n = 1:Ny
                test_fft_1 = test_fft_1.*reshape(fft_1(m,n,:),[1 length(frequencies)]);
                test_fft_1 = test_fft_1/max(test_fft_1);
            end
        end
        indx = test_fft_1 == 1;
        temp_matrix = reshape(fft_1(:,:,indx),[Nx Ny]);
        [~, linearIndex] = max(temp_matrix(:));
        [y, x] = ind2sub(size(temp_matrix), linearIndex);
        Spectrum = reshape(fft_1(x,y,:)/max(max(fft_1(x,y,:))),[1 length(frequencies)]);
        Spectrum = Spectrum/sum(Spectrum);
        [peaks, indexes] = findpeaks(Spectrum);
        freqs = frequencies(indexes);
        dataset_generated(counter).input4 = freqs(1:3);
        f=fit(freqs',peaks','poly1');
        dataset_generated(counter).input5 = f.p1;
        dataset_generated(counter).input6 = sum(frequencies.*Spectrum)/sum(frequencies);
        counter = counter+1;
    end
end