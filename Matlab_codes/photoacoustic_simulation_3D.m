
%% First simulation
close all
clear
mex   -DUSE_OMP C:\GitHub\kwave\software_packages\ValoMC-master\ValoMC-master\cpp\3d\MC3Dmex.cpp COMPFLAGS='\$COMPFLAGS /openmp'
load mesh_space_1.mat
mesh_space = mesh_space(round(3*end/8):round(5*end/8),round(3*end/8):round(5*end/8),round(3*end/8):round(5*end/8));
x_arr = round(-size(mesh_space,1)/2):round(size(mesh_space,1)/2);
y_arr = round(-size(mesh_space,2)/2):round(size(mesh_space,2)/2);
z_arr = round(-size(mesh_space,3)/2):round(size(mesh_space,3)/2);

vmcmesh = createGridMesh(x_arr, y_arr, z_arr); % function provided by ValoMC
nvoxels_total = length(x_arr)*length(y_arr)*length(z_arr);
voxels_in_a_yx_slice = length(y_arr)*length(x_arr);
vmcmedium = createMedium(vmcmesh);
[X,Y,Z] = meshgrid(x_arr,y_arr,z_arr); % Matlab function
F = 0;
vmcmedium.scattering_coefficient = repmat(F(:),6,1); % repeat six times
vmcmedium.absorption_coefficient = repmat(F(:),6,1); % repeat six times
vmcmedium.scattering_anisotropy = 0.9;        
vmcmedium.refractive_index = 1;
vmcboundary = createBoundary(vmcmesh, vmcmedium);   % create a boundary for the mesh
% Create a light source
lightsource = findBoundaries(vmcmesh, 'direction', [0 0 0], [0 0 10], 1);
vmcboundary.lightsource(lightsource) = {'cosinic'};
solution = ValoMC(vmcmesh, vmcmedium, vmcboundary);
TR = triangulation(double(vmcmesh.H),vmcmesh.r); % create a matlab
locations = [X(:) Y(:) Z(:)];              % form a 2D matrix from all
indices = pointLocation(TR,locations);     % query the indices of the
indices(isnan(indices)) = 1;               % set the grid points that
grid_fluence = reshape(solution.element_fluence(indices),size(X));
slice(X, Y, Z, grid_fluence, 0, 0, 0);
xlabel('x [mm]');
ylabel('y [mm]');
zlabel('z [mm]');

view(125,25);
snapnow;


% % % close all
% % % clear
% % % 
% % % mex   -DUSE_OMP C:\GitHub\kwave\software_packages\ValoMC-master\ValoMC-master\cpp\3d\MC3Dmex.cpp COMPFLAGS='\$COMPFLAGS /openmp'
% % % 
% % % load mesh_space_1.mat
% % % array_size = 50000;
% % % T = 3e-6;
% % % dt = T/array_size;
% % % T_array = linspace(0,T,array_size);
% % % mesh_space = mesh_space(1:round(end/4),1:round(end/4),1:round(end/4));
% % % Nx = size(mesh_space,1);           % number of grid points in the x (row) direction
% % % Ny = size(mesh_space,2);           % number of grid points in the y (column) direction
% % % Nz = size(mesh_space,3);           % number of grid points in the z (height) direction
% % % dx = 1e-6;           % grid point spacing in the x direction [m]
% % % dy = 1e-6;           % grid point spacing in the y direction [m]
% % % dz = 1e-6;           % grid point spacing in the z direction [m]
% % % kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);
% % % RBC_indices = find(mesh_space == 1);
% % % Fibrin_indices = find(mesh_space ==2);
% % % Platelet_indices = find(mesh_space ==3);
% % % medium.sound_speed = 1540*ones(Nx, Ny,Nz);    % [m/s]
% % % medium.sound_speed(RBC_indices) = 2380;   % [m/s]
% % % medium.sound_speed(Fibrin_indices) = 2380;   % [m/s]
% % % medium.sound_speed(Platelet_indices) = 2380;   % [m/s]
% % % medium.density = 1025*ones(Nx, Ny,Nz);        % [kg/m^3]
% % % medium.density(RBC_indices) = 2000;
% % % medium.density(Fibrin_indices) = 2000;
% % % medium.density(Platelet_indices) = 2000;
% % % vmcmesh = createGridMesh(kgrid.x_vec*1e3, kgrid.y_vec*1e3, kgrid.z_vec*1e3); % [m to mm]
% % % vmcmedium = createMedium(vmcmesh);
% % % scattering_coefficients = 0.0001*ones(Nx, Ny, Nz);
% % % absorption_coefficients = 0.0001*ones(Nx, Ny, Nz);
% % % scattering_anisotropies = 0.9678*ones(Nx, Ny, Nz);
% % % refractive_indexes = ones(Nx, Ny, Nz);
% % % absorption_coefficients(RBC_indices) = 1e-4;
% % % scattering_coefficients(RBC_indices) = 1e-4;
% % % absorption_coefficients(Fibrin_indices) = 5; %% dummmy value
% % % scattering_coefficients(Fibrin_indices) = 72.78; %% dummmy value
% % % absorption_coefficients(Platelet_indices) = 5; %% dummmy value 
% % % scattering_coefficients(Platelet_indices) = 72.78; %% dummmy value
% % % refractive_indexes(Fibrin_indices) = 1.4;
% % % refractive_indexes(Platelet_indices) = 1.4;
% % % refractive_indexes(RBC_indices) = 1.4;
% % % vmcmedium.scattering_anisotropy =  repmat(scattering_anisotropies(:),6,1);
% % % vmcmedium.absorption_coefficient = repmat(absorption_coefficients(:),6,1);
% % % vmcmedium.scattering_coefficient = repmat(scattering_coefficients(:),6,1);
% % % vmcmedium.refractive_index = repmat(refractive_indexes(:),6,1);
% % % vmcboundary = createBoundary(vmcmesh, vmcmedium);   % create a boundary for the mesh
% % % % Create a light source
% % % lightsource = findBoundaries(vmcmesh, 'direction', [0 0 0], [0 0 10], 10);
% % % vmcboundary.lightsource(lightsource) = {'cosinic'};
% % % solution = ValoMC(vmcmesh, vmcmedium, vmcboundary);
% % % TR = triangulation(double(vmcmesh.H),vmcmesh.r); 
% % % [X,Y,Z] = meshgrid(1:Nx,1:Ny,1:Nz); % Matlab function
% % % locations = [X(:) Y(:) Z(:)];           
% % % indices = pointLocation(TR,locations);     
% % % indices(isnan(indices)) = 1;             
% % % grid_fluence = reshape(solution.element_fluence(indices),size(X));
% % % slice(X, Y, Z, grid_fluence, 0, 0, 0);
% % % xlabel('x [mm]');
% % % ylabel('y [mm]');
% % % zlabel('z [mm]');
% % % 
% % % view(125,25);
% % % snapnow;