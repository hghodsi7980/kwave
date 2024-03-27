function cylinder = makeCylinder(Nx, Ny, Nz, cx, cy, cz, radius, height, orientation)
% MAKECYLINDER Create a binary map of a filled cylinder within a 3D grid.

% Define literals
MAGNITUDE = 2;


% Force integer values
Nx = round(Nx);
Ny = round(Ny);
Nz = round(Nz);
cx = round(cx);
cy = round(cy);
cz = round(cz);

% Check for zero values
if cx == 0
    cx = floor(Nx / 2) + 1;
end
if cy == 0
    cy = floor(Ny / 2) + 1;
end
if cz == 0
    cz = floor(Nz / 2) + 1;
end

% Create empty matrix

cylinder = zeros(Nx, Ny, Nz);

% Rotate grid points
[XX, YY, ZZ] = meshgrid(1:Nx, 1:Ny, 1:Nz);
grid_points = [XX(:), YY(:), ZZ(:)];
rot_matrix = makeRotationMatrix(orientation);
rotated_points = (rot_matrix * grid_points')';

% Create cylinder
for i = 1:size(rotated_points, 1)
    x = rotated_points(i, 1) - cx;
    y = rotated_points(i, 2) - cy;
    z = rotated_points(i, 3) - cz;
    if x^2 + y^2 <= radius^2 && abs(z) <= height / 2
        cylinder(rotated_points(i, 1), rotated_points(i, 2), rotated_points(i, 3)) = MAGNITUDE;
    end
end

end

function R = makeRotationMatrix(orientation)
% MAKEROTATIONMATRIX Generate a 3x3 rotation matrix from Euler angles (orientation).

phi = orientation(1); % Rotation about x-axis
theta = orientation(2); % Rotation about y-axis
psi = orientation(3); % Rotation about z-axis

Rx = [1, 0, 0; 0, cos(phi), -sin(phi); 0, sin(phi), cos(phi)];
Ry = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)];
Rz = [cos(psi), -sin(psi), 0; sin(psi), cos(psi), 0; 0, 0, 1];

R = Rz * Ry * Rx; % Combine rotations
end
