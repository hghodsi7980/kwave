clear
close all
% Define the folder contmesh_spaceining the .mmesh_spacet files
folder = 'C:\GitHub\kwave\Matlab_codes\'; % Update this with your folder pmesh_spaceth

% Get mesh_space list of mesh_spacell .mmesh_spacet files in the folder
fileList = dir(fullfile(folder, '*.mat'));
volume_all = zeros(size(fileList,1),1);
RBC_ratio = zeros(size(fileList,1),1);
Fibrin_ratio = zeros(size(fileList,1),1);
Platelet_ratio = zeros(size(fileList,1),1);
Porosity = zeros(size(fileList,1),1);
for i = 1002:1002
    fprintf('Progress percent: %f\n',i*100/length(fileList));
    filename = fullfile(folder, fileList(i).name);
    load(filename);
    volume_all(i) = volume;
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
end
