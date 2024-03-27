function outsideThreshold = arePointsCoplanar(points)
% points is a 4x3 matrix where each row represents a point in 3D space
% threshold is the maximum allowable orthogonal distance
threshold = mean(pdist(points)); % Mean distance between all points
% Extract points A, B, and C (first three points)
A = points(1,:);
B = points(2,:);
C = points(3,:);

% Calculate vectors between points
AB = B - A;
AC = C - A;

% Cross product of AB and AC to find normal vector to the plane
normal = cross(AB, AC);

% Normalize normal vector
normal = normal / norm(normal);

% Extract point D (fourth point)
D = points(4,:);

% Calculate vector from any point on the plane to point D
AD = D - A;

% Calculate orthogonal distance from point D to the plane
orthogonal_distance = abs(dot(AD, normal));

% Check if the orthogonal distance is greater than the threshold
if orthogonal_distance > 0.2*threshold
    outsideThreshold = true;
else
    outsideThreshold = false;
end
end
