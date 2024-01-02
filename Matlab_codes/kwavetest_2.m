
%% First simulation
% close all
% run 10 times and average:

addpath(genpath('/home/hghodsi/Software_packages'))
addpath(genpath('/home/hghodsi/Matlab_codes'))


average_size = 1;
dt =5e-14;   %2MHz
T = 1e-6;
T_array = 0:dt:T;
% sensor_data_avr = zeros(average_size,size(T_array,2));
%mex   -DUSE_OMP C:\Git\kwave\software_packages\ValoMC-master\ValoMC-master\cpp\2d\MC2Dmex.cpp COMPFLAGS='\$COMPFLAGS /openmp'
mex   -DUSE_OMP /home/hghodsi/Software_packages/ValoMC-master/ValoMC-master/cpp/2d/MC2Dmex.cpp COMPFLAGS='\$COMPFLAGS -fopenmp' CXXFLAGS='\$CXXFLAGS -fopenmp' LDFLAGS='\$LDFLAGS -fopenmp'
for index = 1:average_size
    disp(index)
    % Create the k-Wave grid
    Nx = 500;           % number of grid points in the x (row) direction
    Ny = 500;           % number of grid points in the y (column) direction
    dx = 2.5e-6;        % grid point spacing in the x direction [m]
    dy = 2.5e-6;        % grid point spacing in the y direction [m]
    kgrid = kWaveGrid(Nx, dx, Ny, dy);
    vessel_mask = ones(Nx,Ny);   % all of the ROI for the test
    clot_radious = round(500e-6/dx);            % clot radius = 250u
    clot_mask = makeDisc(Nx,Ny,Nx/2,Ny/2,clot_radious);
    RBC_radious = round(5e-6/dx);
    % RBC_ratio_Blood = 0;
    RBC_ratio_clot = 0.8;
    % RBC_count_vessel = round((RBC_ratio_Blood*Nx*Ny)/(pi*RBC_radious^2));
%     RBC_count_clot = round(RBC_ratio_clot*(pi*clot_radious^2)/(pi*RBC_radious^2));
    % Fibrin_count = 0;
    % average_fibrin_length = round(100e-6/dx);
    % output_mask_vessel = fill_with_RBC(Nx,Ny,vessel_mask,RBC_count_vessel,RBC_radious);
    output_mask_clot = fill_with_RBC(Nx,Ny,clot_mask,10000,RBC_radious);   % 80%
    % output_mask_clot = output_mask_clot+fill_with_fibrin(Nx, Ny, clot_mask, Fibrin_count, average_fibrin_length);
    % output_mask_clot(output_mask_clot == 3) = 2;
    disc_indices = find(output_mask_clot==1);
    line_indices = find(output_mask_clot==2);
    %disp(sum(sum(output_mask_clot))/sum(sum(clot_mask)));
    %figure
    %xx = linspace(-Nx*dx/2,+Nx*dx/2,Nx);
    %yy = linspace(-Ny*dy/2,+Ny*dy/2,Ny);
    %imagesc(xx,yy,output_mask_clot);

    xlabel width(meter)
    ylabel depth(meter)

    % Define the acoustic properties

    medium.sound_speed = 1500*ones(Nx, Ny);    % [m/s]
    medium.sound_speed(disc_indices) = 2380;   % [m/s]
    medium.sound_speed(line_indices) = 2380;   % [m/s]

    medium.density = 1025*ones(Nx, Ny);        % [kg/m^3]
    medium.density(disc_indices) = 1400;
    medium.density(line_indices) = 1400;

    %% Create a ValoMC mesh
    vmcmesh = createGridMesh(kgrid.y_vec*1e3, kgrid.x_vec*1e3); % [m to mm]

    %% Define optical coefficients
    vmcmedium.scattering_coefficient = 0.0001*ones(Nx, Ny);
    vmcmedium.absorption_coefficient = 0.0001*ones(Nx, Ny);
    % Define the acoustic properties
    vmcmedium.absorption_coefficient(disc_indices) = 0.2248;
    vmcmedium.scattering_coefficient(disc_indices) = 0.7278;
    vmcmedium.absorption_coefficient(line_indices) = 22;
    vmcmedium.scattering_coefficient(line_indices) = 0.1;
    vmcmedium.scattering_anisotropy = 0.9678;           % scattering anisotropy parameter [unitless]
    vmcmedium.refractive_index = 1.0*ones(Nx, Ny);
    vmcmedium.refractive_index(line_indices) = 1.602;
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
    kgrid.t_array = 0:dt:T;
    sensor_data = kspaceFirstOrder2DC(kgrid, medium, source, sensor, 'PMLInside', false);

end

%% process the data

Fs = T/dt;  % Sampling rate in Hz
N = length(mean(sensor_data,1));
fft_r = zeros(1,N);
for i = 1:Nx
    disp(i)
    fft_r(i,:) = exp(-abs(i-250)*1i)*fft(sensor_data(i,:));
end
fft_result = mean(fft_r);
frequencies = (0:N-1) * (Fs / N);
frequencies = frequencies(1:N/2);
fft_result = fft_result(1:N/2);
fft_result = abs(fft_result);



%% save the output

save k1 frequencies fft_result




