
%% First simulation
% close all
% run 10 times and average:

% addpath(genpath('/home/hghodsi/Software_packages'))
% addpath(genpath('/home/hghodsi/Matlab_codes'))

simulation_size = 0;
array_size = 100000;
T = 5e-6;
dt = T/array_size;
T_array = linspace(0,T,array_size);
% sensor_data_avr = zeros(average_size,size(T_array,2));
% mex   -DUSE_OMP C:\Git\kwave\software_packages\ValoMC-master\ValoMC-master\cpp\2d\MC2Dmex.cpp COMPFLAGS='\$COMPFLAGS /openmp'
%mex   -DUSE_OMP /home/hghodsi/Software_packages/ValoMC-master/ValoMC-master/cpp/2d/MC2Dmex.cpp COMPFLAGS='\$COMPFLAGS -fopenmp' CXXFLAGS='\$CXXFLAGS -fopenmp' LDFLAGS='\$LDFLAGS -fopenmp'
for index = 1:1
    h = waitbar(0, 'Progress'); % Create a progress bar window
    waitbar(index / 1000, h, sprintf('Overall Progress: %d%%', round((index / 1000) * 100)));
    pause(0.1); % Optional: Add a delay to slow down the progress for demonstration purposes
    % Create the k-Wave grid
    Nx = 2000;           % number of grid points in the x (row) direction
    Ny = 2000;           % number of grid points in the y (column) direction
    dx = 1e-6;        % grid point spacing in the x direction [m]
    dy = 1e-6;        % grid point spacing in the y direction [m]
    kgrid = kWaveGrid(Nx, dx, Ny, dy);
    vessel_mask = ones(Nx,Ny);   % all of the ROI for the test
    clot_radious = 1e-6*250;            % clot radius between 100 and 500
    clot_mask = makeDisc(Nx,Ny,Nx/2,Ny/2,round(clot_radious/dx));
    RBC_radious = round(5e-6/dx);
    rng('shuffle');
    RBC_count = 5000;
    rng('shuffle');
    Fibrin_count = 0;
    average_fibrin_length = round(100e-6/dx);
    output_mask_clot = fill_with_RBC(Nx,Ny,clot_mask,RBC_count,RBC_radious);
    output_mask_clot = output_mask_clot+fill_with_fibrin(Nx, Ny, clot_mask, Fibrin_count, average_fibrin_length);
    output_mask_clot(output_mask_clot == 3) = 1;  % make the overlaps the same as RBC 
    disc_indices = find(output_mask_clot==1);
    line_indices = find(output_mask_clot==2);
    medium.sound_speed = 1540*ones(Nx, Ny);    % [m/s]
    medium.sound_speed(disc_indices) = 2380;   % [m/s]
    medium.sound_speed(line_indices) = 1540;   % [m/s]

    medium.density = 1025*ones(Nx, Ny);        % [kg/m^3]
    medium.density(disc_indices) = 2000;
    medium.density(line_indices) = 2000;

    %% Create a ValoMC mesh
    vmcmesh = createGridMesh(kgrid.y_vec*1e3, kgrid.x_vec*1e3); % [m to mm]

    %% Define optical coefficients
    vmcmedium.scattering_coefficient = 0.0001*ones(Nx, Ny);
    vmcmedium.absorption_coefficient = 0.0001*ones(Nx, Ny);
    % Define the acoustic properties
    vmcmedium.absorption_coefficient(disc_indices) = 22.48;
    vmcmedium.scattering_coefficient(disc_indices) = 72.78;
    vmcmedium.absorption_coefficient(line_indices) = 5; %% dummmy value
    vmcmedium.scattering_coefficient(line_indices) = 72;  %% dummy value
    vmcmedium.scattering_anisotropy = 0.9678;           % scattering anisotropy parameter [unitless]
    vmcmedium.refractive_index = 1.0*ones(Nx, Ny);
    vmcmedium.refractive_index(line_indices) = 1.5;
    vmcmedium.refractive_index(disc_indices) = 1.4;
    % Define the Gruneisen parameter describing photoacoustic efficiency
    vmcmedium.gruneisen_parameter = 0.16*ones(Nx, Ny);      % [unitless]

    %% Create a light source

    % Set a light source with a width of 10u and cosinic directional profile
    % in -x direction
    line_start = [0 ceil(-3/5*Nx)];
    line_end = [0 0];
    line_width = 1;
    boundary_with_lightsource = findBoundaries(vmcmesh, 'direction', ...
        line_start, ...
        line_end,  ...
        line_width);
    vmcboundary.lightsource(boundary_with_lightsource) = {'gaussian'};
    vmcboundary.lightsource_gaussian_sigma = 0.1;



    %% Run the Monte Carlo simulation
    solution = ValoMC(vmcmesh, vmcmedium, vmcboundary);

    %% Compute the initial pressure from the photon fluence
    %
    % Compute the absorbed optical energy density.
    % multiply by
    % 1e6   to convert [J/mm^2] to [J/m^2]
    % 1e-3  to set the total energy in the pulse to 1 mJ
    %
    vmcmedium.absorbed_energy = vmcmedium.absorption_coefficient .* solution.grid_fluence*1e3; % [J/m3]

    % Compute the initial pressure distribution
    source.p0 = vmcmedium.gruneisen_parameter .* vmcmedium.absorbed_energy;  % [Pa]

    %% Define the k-Wave sensor mask

    % define a 2D binary sensor mask in the shape of a line
    width = Ny; % [grid points]
    sensor.mask = zeros(Nx, Ny);
    sensor.mask(ceil(Ny/2 - width/2 + 1):ceil(Ny/2 + width/2),1) = 1;

    %% Move the perfectly matched layer (PML) outside of the computation domain and run the acoustic simulation
    % The PML is a layer that absorbs waves for simulating free regions and
    % is normally contained within the computation  region of k-Wave.
    % For a more straightforward mapping between the
    % computation region of k-Wave and ValoMC, the PML is moved outside
    % of the computation region.
    kgrid.dt = dt;
    kgrid.t_array = T_array;
    sensor_data = kspaceFirstOrder2DG(kgrid, medium, source, sensor, 'PMLInside', false);
    Fs = 1/dt;  % Sampling rate in Hz
    N = length(mean(sensor_data,1));
    center_freq = 5e6;
    bandwidth = 1e7;
    sensor_data_gaussian = gaussianFilter(sensor_data, Fs, center_freq, bandwidth);





    %% process the data
    plot(kgrid.t_array,sensor_data(500,:));
    fft_r = zeros(Nx,N/2);
    for i = 1:Nx
        disp(i)
        fft_result = fft(sensor_data_gaussian(i,:));
        fft_r(i,:) = abs(fft_result(1:N/2));
    end
    frequencies = (0:N-1) * (Fs / N);
    frequencies = frequencies(1:N/2);
    spectrum = fft_r(Nx/2,:);

    %% save the output
    file_name = sprintf("New_data_a%d.mat",index);
    save(file_name,'frequencies','spectrum','output_mask_clot', "dx" , "dy" , "dt");

    p_xy = kspaceLineRecon(sensor_data.', dy, kgrid.dt, medium.sound_speed(1,1), 'Plot', true);

close(h);
end
