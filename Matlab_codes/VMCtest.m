clear 
close all
figure

% Create a finer mesh

x_arr = -2:0.1:2;
y_arr = -2:0.1:2;
z_arr = -2:0.1:2;

vmcmesh = createGridMesh(x_arr, y_arr, z_arr); % function provided by ValoMC
vmcmedium = createMedium(vmcmesh);

[X,Y,Z] = meshgrid(x_arr,y_arr,z_arr); % Matlab function
F = 1.3+cos(X*3).*cos(Y*3).*cos(Z*3)*0.2+0.2;

vmcmedium.scattering_coefficient = 1.0;
vmcmedium.absorption_coefficient = repmat(F(:),6,1); % repeat six times

vmcmedium.scattering_anisotropy = 0.9;
vmcmedium.refractive_index = 1;

vmcboundary = createBoundary(vmcmesh, vmcmedium);   % create a boundary for the mesh

% Create a light source
lightsource = findBoundaries(vmcmesh, 'direction', [0 0 0], [0 20 20], 1);
vmcboundary.lightsource(lightsource) = {'cosinic'};


solution = ValoMC(vmcmesh, vmcmedium, vmcboundary);

TR = triangulation(double(vmcmesh.H),vmcmesh.r); % create a matlab
                                           % triangulation object
                                           % from the mesh

locations = [X(:) Y(:) Z(:)];              % form a 2D matrix from all
                                           % the grid points

indices = pointLocation(TR,locations);     % query the indices of the
                                           % tetrahedrons at grid
                                           % points

indices(isnan(indices)) = 1;               % set the grid points that
                                           % do not belong to the mesh
                                           % to point at the first
                                           % element

% get the values on a grid
grid_fluence = reshape(solution.element_fluence(indices),size(X));


slice(X, Y, Z, grid_fluence, 2, 2, 2);
xlabel('x [mm]');
ylabel('y [mm]');
zlabel('z [mm]');

view(125,25);
snapnow;




