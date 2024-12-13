for jj = 1
% Load the figure
hFig = openfig('figure1.fig');

% Get the current axes of the figure
ax = gca;

% Remove any data tips, annotations, or data markers

delete(findall(hFig, 'Type', 'datatip'));
delete(findall(hFig, 'Type', 'text'));
delete(findall(hFig, 'Type', 'annotation'));

% Remove the axes and grid
axis off;
grid off;

% Set the background color to black
set(hFig, 'Color', 'k');

% Set up video writer
outputVideo = VideoWriter(sprintf('rotation_with_zoom_animation_blackkk_bg %d.avi',jj));
outputVideo.FrameRate = 30; % Set the frame rate (adjust as needed)
open(outputVideo);

% Number of frames and rotation increment
nFrames = 360; % Adjust the number of frames for smoother rotation
rotationStep = 360 / nFrames;

% Set initial CameraViewAngle for zoom effect
normalViewAngle = ax.CameraViewAngle; % Start from the default view angle (100% zoom)
zoomInAngle = normalViewAngle / 20; % Zoom in to half the angle (adjust this value as needed)

% Loop through and rotate the view with zoom
for k = 1:nFrames
    disp(k*100/nFrames)
    delete(findall(hFig, 'Type', 'datatip'));
    delete(findall(hFig, 'Type', 'text'));
    delete(findall(hFig, 'Type', 'annotation'));
    % Rotate the view
    view(ax, [k *   rotationStep, 30]); % [Azimuth, Elevation], adjust elevation if needed
    
    % Apply zoom effect: first half zoom in, second half zoom out
    if k <= nFrames / 2
        % Zoom in during the first half
        currentCameraAngle = normalViewAngle - ((normalViewAngle - zoomInAngle) * (k / (nFrames / 2)));
    else
        % Zoom out during the second half
        currentCameraAngle = zoomInAngle + ((normalViewAngle - zoomInAngle) * ((k - nFrames / 2) / (nFrames / 2)));
    end
    
    % Set the camera view angle
    ax.CameraViewAngle = currentCameraAngle;
    
    % Capture the current frame
    frame = getframe(hFig);
    
    % Write the frame to the video
    writeVideo(outputVideo, frame);
end

% Close the video writer
close(outputVideo);

% Close the figure
close(hFig);
end