function cylinder = makeCylinder(Nx, Ny, Nz, cx, cy, cz, radius, height, orientation)
% MAKECYLINDER Create a binary map of a filled, rotated cylinder within a 3D grid.

% Define literals
MAGNITUDE = 1;

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

% Generate grid points relative to the center
[X, Y, Z] = meshgrid(1:Nx, 1:Ny, 1:Nz);
X = X - cx;
Y = Y - cy;
Z = Z - cz;

% Create cylinder mask in unrotated space
cylinder_mask = (X.^2 + Y.^2 <= radius^2) & ...
                (abs(Z) <= height / 2);

% Rotate the coordinates
rot_matrix = makeRotationMatrix(orientation);
rotated_coords = rot_matrix * [X(:), Y(:), Z(:)]'; % Apply rotation
rotated_X = reshape(rotated_coords(1, :), size(X));
rotated_Y = reshape(rotated_coords(2, :), size(Y));
rotated_Z = reshape(rotated_coords(3, :), size(Z));

% Map rotated coordinates back to the original grid
[Xq, Yq, Zq] = meshgrid(1:Nx, 1:Ny, 1:Nz);

% Interpolate the rotated cylinder onto the grid
rotated_cylinder = interp3(X, Y, Z, double(cylinder_mask), ...
                           rotated_X, rotated_Y, rotated_Z, 'linear', 0);

% Threshold to create binary mask
cylinder = rotated_cylinder >= 0.5;
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
