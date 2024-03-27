% Generate random coordinates
N = 100; % Number of points
coordinates = rand(N, 3); % Generate random 3D coordinates

% Perform Delaunay triangulation
DT = delaunayTriangulation(coordinates);

% Refine the triangulation to ensure each point is connected to 3 or 4 neighbors
p = DT.Points;
t = DT.ConnectivityList;
e = DT.edges;
[p1,e1,t1] = refinemesh("lshapeg",p,e,t);

% Plot the triangulation
tetramesh(TR);
axis equal;
title('Tetrahedral Mesh');

% To access the triangulation connectivity, you can use TR.ConnectivityList
