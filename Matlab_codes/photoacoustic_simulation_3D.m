
%% First simulation
close all
clear
mex   -DUSE_OMP C:\GitHub\kwave\software_packages\ValoMC-master\ValoMC-master\cpp\3d\MC3Dmex.cpp COMPFLAGS='\$COMPFLAGS /openmp'
load mesh_space_1.mat
mesh_space = mesh_space(round(3*end/8):round(5*end/8),round(3*end/8):round(5*end/8),round(3*end/8):round(5*end/8));
x_arr = round(-size(mesh_space,1)/2):round(size(mesh_space,1)/2);
y_arr = round(-size(mesh_space,2)/2):round(size(mesh_space,2)/2);
z_arr = round(-size(mesh_space,3)/2):round(size(mesh_space,3)/2);
x_arr = x_arr/1000;
y_arr = y_arr/1000;
z_arr = z_arr/1000;
vmcmesh = createGridMesh(x_arr, y_arr, z_arr); % function provided by ValoMC
nvoxels_total = length(x_arr)*length(y_arr)*length(z_arr);
voxels_in_a_yx_slice = length(y_arr)*length(x_arr);
vmcmedium = createMedium(vmcmesh);
[X,Y,Z] = meshgrid(x_arr,y_arr,z_arr); % Matlab function
scattering_coefficients = 0.0001*ones(Nx, Ny, Nz);
absorption_coefficients = 0.0001*ones(Nx, Ny, Nz);
scattering_anisotropies = 0.9678*ones(Nx, Ny, Nz);
refractive_indexes = ones(Nx, Ny, Nz);
absorption_coefficients(RBC_indices) = 1e-4;
scattering_coefficients(RBC_indices) = 1e-4;
absorption_coefficients(Fibrin_indices) = 5; %% dummmy value
scattering_coefficients(Fibrin_indices) = 72.78; %% dummmy value
absorption_coefficients(Platelet_indices) = 5; %% dummmy value 
scattering_coefficients(Platelet_indices) = 72.78; %% dummmy value
refractive_indexes(Fibrin_indices) = 1.4;
refractive_indexes(Platelet_indices) = 1.4;
refractive_indexes(RBC_indices) = 1.4;
vmcmedium.scattering_anisotropy =  repmat(scattering_anisotropies(:),6,1);
vmcmedium.absorption_coefficient = repmat(absorption_coefficients(:),6,1);
vmcmedium.scattering_coefficient = repmat(scattering_coefficients(:),6,1);
vmcmedium.refractive_index = repmat(refractive_indexes(:),6,1);
vmcboundary = createBoundary(vmcmesh, vmcmedium);   % create a boundary for the mesh
% Create a light source
lightsource = findBoundaries(vmcmesh, 'direction', [0 0 0], [0 0 1], max(x_arr,y_arr,z_arr));
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

