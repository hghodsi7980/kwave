
%% First simulation
close all
clear
mex   -DUSE_OMP C:\GitHub\kwave\software_packages\ValoMC-master\ValoMC-master\cpp\3d\MC3Dmex.cpp COMPFLAGS='\$COMPFLAGS /openmp'
% load mesh_space_1.mat
for counter =1:50
    mesh_space = zeros(250,250,250);
    num_particles = 200;
    for i = 1:num_particles
        mesh_space = mesh_space+makeBall(250,250,250,randi(245)+5,randi(245)+5,randi(245)+5,5);
    end
    mesh_space(mesh_space > 1) = 1;
    x_arr = linspace(round(-size(mesh_space,1)/2),round(size(mesh_space,1)/2),size(mesh_space,1));
    y_arr = linspace(round(-size(mesh_space,2)/2),round(size(mesh_space,2)/2),size(mesh_space,2));
    z_arr = linspace(round(-size(mesh_space,3)/2),round(size(mesh_space,3)/2),size(mesh_space,3));
    x_arr = x_arr/1000;
    y_arr = y_arr/1000;
    z_arr = z_arr/1000;
    vmcmesh = createGridMesh(x_arr, y_arr, z_arr); % function provided by ValoMC
    vmcmedium = createMedium(vmcmesh);
    [X,Y,Z] = meshgrid(y_arr,x_arr,z_arr); % Matlab function
    scattering_coefficients = 0.000001*ones(size(mesh_space));
    absorption_coefficients = 0.000001*ones(size(mesh_space));
    scattering_anisotropies = 0.9678*ones(size(mesh_space));
    refractive_indexes = 1.0*ones(size(mesh_space));
    absorption_coefficients(mesh_space == 1) = 22.48;
    scattering_coefficients(mesh_space == 1) = 72.78;
    absorption_coefficients(mesh_space == 2) = 5; %% dummmy value
    scattering_coefficients(mesh_space == 2) = 72.78; %% dummmy value
    absorption_coefficients(mesh_space == 3) = 5; %% dummmy value
    scattering_coefficients(mesh_space == 3) = 72.78; %% dummmy value
    refractive_indexes(mesh_space == 1) = 1.4;
    refractive_indexes(mesh_space == 2) = 1.4;
    refractive_indexes(mesh_space == 3) = 1.4;
    vmcmedium.scattering_anisotropy =  repmat(scattering_anisotropies(:),6,1);
    vmcmedium.absorption_coefficient = repmat(absorption_coefficients(:),6,1);
    vmcmedium.scattering_coefficient = repmat(scattering_coefficients(:),6,1);
    vmcmedium.refractive_index = repmat(refractive_indexes(:),6,1);
    vmcboundary = createBoundary(vmcmesh, vmcmedium);   % create a boundary for the mesh
    % Create a light source
    lightsource = findBoundaries(vmcmesh, 'direction', [0 0 0], [0 0 1],max([size(mesh_space,1)/2,size(mesh_space,2)/2,size(mesh_space,3)/2])/1000);
    vmcboundary.lightsource(lightsource) = {'cosinic'};
    solution = ValoMC(vmcmesh, vmcmedium, vmcboundary);
    TR = triangulation(double(vmcmesh.H),vmcmesh.r); % create a matlab
    locations = [X(:) Y(:) Z(:)];              % form a 2D matrix from all
    indices = pointLocation(TR,locations);     % query the indices of the
    indices(isnan(indices)) = 1;               % set the grid points that
    grid_fluence = reshape(solution.element_fluence(indices),size(X));
    absorbed_energy = absorption_coefficients .* grid_fluence*1e3; % [J/m3]
    clearvars('-except', 'absorbed_energy', 'mesh_space','counter');
    Nx = size(absorbed_energy,1);
    Ny = size(absorbed_energy,2);
    Nz = size(absorbed_energy,3);
    dx = 1e-6;
    dy = 1e-6;
    dz = 1e-6;
    kgrid = kWaveGrid(Nx, dx, Ny, dy,Nz, dz);
    gruneisen_parameter = 0.16;
    source.p0 = gruneisen_parameter .* absorbed_energy;  % [Pa]
    medium.sound_speed = 1540*ones(Nx,Ny,Nz);    % [m/s]
    medium.sound_speed(mesh_space == 1) = 2380;   % [m/s]
    medium.sound_speed(mesh_space == 2) = 2380;   % [m/s]
    medium.sound_speed(mesh_space == 3) = 2380;   % [m/s]
    medium.density = 1025*ones(Nx, Ny,Nz);        % [kg/m^3]
    medium.density(mesh_space == 1) = 2000;
    medium.density(mesh_space == 2) = 2000;
    medium.density(mesh_space == 3) = 2000;
    array_size = 25000;
    T = 2e-6;
    dt = T/array_size;
    T_array = linspace(0,T,array_size);
    kgrid.t_array = T_array;
    sensor.mask = zeros(Nx, Ny, Nz);
    sensor.mask(:,:,end) = 1;
    sensor_data = kspaceFirstOrder3DG(kgrid, medium, source, sensor, 'PMLInside', false);
    Fs = 1/dt;
    N = length(T_array);
    sensor_data_3D = reshape(sensor_data,[Nx Ny N]);
    fft_result = zeros(Nx,Ny,round(N/2));
    frequencies = (0:N-1) * (Fs / N);
    frequencies = frequencies(1:round(N/2));
    for i = 1:Nx
        for j = 1:Ny
            fft_temp = abs(fft(reshape(sensor_data_3D(i,j,:),[1 N])));
            fft_result(i,j,:) = fft_temp(1:round(N/2));
        end
    end
    filename = fprintf('simulation_%d.mat',counter);
    save(filename,'absorbed_energy','frequencies','fft_result','-v7.3');
end