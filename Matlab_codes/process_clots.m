clear
close all
% % Define the folder contmesh_spaceining the .mmesh_spacet files
folder = 'C:\GitHub\kwave\Matlab_codes\newdatafinal\New folder'; % Update this with your folder pmesh_spaceth
% 
% % Get mesh_space list of mesh_spacell .mmesh_spacet files in the folder
% fileList1 = dir(fullfile(folder, 'dataset_*.mat'));
fileList2 = dir(fullfile(folder, '*.mat'));
volume_all = zeros(size(fileList2,1),1);
RBC_ratio = zeros(size(fileList2,1),1);
Fibrin_ratio = zeros(size(fileList2,1),1);
Platelet_ratio = zeros(size(fileList2,1),1);
Porosity = zeros(size(fileList2,1),1);
% input_matrix = zeros(size(fileList1,1),14);
output_matrix = zeros(size(fileList2,1),3);
for i = 1:length(fileList2)
    tic
    fprintf('Progress percent: %f\n',i*100/length(fileList2));
%     filename = fullfile(folder, fileList1(i).name);
%     load(filename);
    filename = fullfile(folder, fileList2(i).name);
    load(filename);
    volume_all(i) = volume;
%     frequencies = dataset.frequencies;
%     input_matrix(i,1) = volume;
%     input_matrix(i,2:6) = frequencies(dataset.coordinate_f(1:5));
%     for j = 1:5
%         peak_value = dataset.spec_all(dataset.coordinate_x(j),dataset.coordinate_y(j),dataset.coordinate_f(j));
%         peak_value = peak_value/max(dataset.spec_all(:));
%         input_matrix(i,6+j:11) = peak_value;
%     end
%     dataset.spec_all = dataset.spec_all/max(dataset.spec_all(:));
%     spec_average_1 = mean(dataset.spec_all,1);
%     spec_average_2 = mean(spec_average_1,2);
%     spec_average = reshape(spec_average_2,1,[]);
%     f=fit(frequencies',spec_average','poly1');
%     input_matrix(i,12) = f.p1;   %slope
%     [~, index] = max(spec_average);
%     input_matrix(i,13) = frequencies(index);   % peak frequency shift
%     a = f(0);
%     BW = (a/2-f.p2)/f.p1;
%     input_matrix(i,14) = BW;
    [m, n, p] = size(mesh_space);
    mesh_space_vector = mesh_space(:);
    % Count occurrences of 1's, 2's, and 3's
    count_ones = sum(mesh_space_vector == 1);
    count_twos = sum(mesh_space_vector == 2);
    count_threes = sum(mesh_space_vector == 3);
    total_nonzero = numel(mesh_space_vector) - sum(mesh_space_vector == 0);
    RBC_ratio(i) = count_ones / total_nonzero;
    Fibrin_ratio(i) = count_twos / total_nonzero;
    Platelet_ratio(i) = count_threes / total_nonzero;
    Porosity(i) = 1-(total_nonzero/volume);
    output_matrix(i,1) = Porosity(i);
    output_matrix(i,2) = RBC_ratio(i)/Fibrin_ratio(i);
    output_matrix(i,3) = Platelet_ratio(i)/Fibrin_ratio(i);
    toc
end
% 
% 
% % % Load your dataset
% % % Assuming 'inputs' is a 14x1000 matrix and 'outputs' is a 3x1000 matrix
% inputs = input_matrix(:,2:14)'./input_matrix(:,1)';
% inputs = normalize(inputs, 2); % Normalize each feature
% outputs = output_matrix(:,2)';
% 
% % Normalize outputs to range [0, 1]
% outputMin = min(outputs, [], 2); % Minimum values for each output
% outputMax = max(outputs, [], 2); % Maximum values for each output
% outputs = (outputs - outputMin) ./ (outputMax - outputMin);
% 
% % 
% % % Define the network with one hidden layer (10 neurons)
% % hiddenLayerSize = 10;
% % net = feedforwardnet(hiddenLayerSize);
% % 
% % % Split data into training and testing sets
% % [trainInd, valInd, testInd] = dividerand(1000, 0.7, 0.15, 0.15);
% % 
% % % Train the network
% % net.divideFcn = 'divideind';
% % net.divideParam.trainInd = trainInd;
% % net.divideParam.valInd = valInd;
% % net.divideParam.testInd = testInd;
% % net = train(net, inputs, outputs);
% % 
% % % Evaluate the network
% % outputsPred = net(inputs);
% % performance = perform(net, outputs, outputsPred);
% % 
% % % Display the performance
% % disp(['Performance: ', num2str(performance)]);
% 
% 
% 
% % Define the network with two hidden layers (20 neurons and 10 neurons)
% hiddenLayerSizes = [20 50 50 20];
% net = feedforwardnet(hiddenLayerSizes);
% 
% % Split data into training and testing sets
% [trainInd, valInd, testInd] = dividerand(1000, 0.7, 0.15, 0.15);
% % 
% % Train the network
% net.divideFcn = 'divideind';
% net.divideParam.trainInd = trainInd;
% net.divideParam.valInd = valInd;
% net.divideParam.testInd = testInd;
% net = train(net, inputs, outputs);
% 
% % Evaluate the network
% outputsPred = net(inputs);
% performance = perform(net, outputs, outputsPred);
% 
% % Display the performance
% disp(['Performance: ', num2str(performance)]);



% % % Convert outputs to a vector of class labels
% % [~, outputs] = max(outputs, [], 1); % This converts one-hot encoding to class labels
% % outputs = outputs'; % Transpose to make it a column vector
% % 
% % % Train an SVM
% % SVMModel = fitcecoc(inputs', outputs);
% % 
% % % Evaluate the SVM
% % predictedLabels = predict(SVMModel, inputs');
% % accuracy = sum(predictedLabels == outputs) / numel(outputs);
% % disp(['SVM Accuracy: ', num2str(accuracy)]);

