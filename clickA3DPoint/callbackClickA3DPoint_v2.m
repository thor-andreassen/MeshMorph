function callbackClickA3DPoint(src, eventData, pointCloud)
% CALLBACKCLICK3DPOINT mouse click callback function for CLICKA3DPOINT
%
%   The transformation between the viewing frame and the point cloud frame
%   is calculated using the camera viewing direction and the 'up' vector.
%   Then, the point cloud is transformed into the viewing frame. Finally,
%   the z coordinate in this frame is ignored and the x and y coordinates
%   of all the points are compared with the mouse click location and the 
%   closest point is selected.
%
%   Babak Taati - May 4, 2005
%   revised Oct 31, 2007
%   revised Jun 3, 2008
%   revised May 19, 2009
global points_index
global points
global points_graph
global point_graph_old

point = get(gca, 'CurrentPoint'); % mouse click position
camPos = get(gca, 'CameraPosition'); % camera position
camTgt = get(gca, 'CameraTarget'); % where the camera is pointing to

camDir = camPos - camTgt; % camera direction
camUpVect = get(gca, 'CameraUpVector'); % camera 'up' vector

% build an orthonormal frame based on the viewing direction and the 
% up vector (the "view frame")
zAxis = camDir/norm(camDir);    
upAxis = camUpVect/norm(camUpVect); 
xAxis = cross(upAxis, zAxis);
yAxis = cross(zAxis, xAxis);

rot = [xAxis; yAxis; zAxis]; % view rotation 

% the point cloud represented in the view frame
rotatedPointCloud = rot * pointCloud; 

% the clicked point represented in the view frame
rotatedPointFront = rot * point' ;

% find the nearest neighbour to the clicked point 
pointCloudIndex = dsearchn(rotatedPointCloud(1:2,:)', ... 
    rotatedPointFront(1:2));

point_graph_old = findobj(gca,'Tag','pt'); % try to find the old point
selectedPoint = pointCloud(:, pointCloudIndex);


points_index=[points_index,pointCloudIndex];
points=pointCloud(:, points_index);

if ~isempty(points) && ~isempty(points_graph)
    disp('update graph')
    points_graph.XData=points(1,:);
    points_graph.YData=points(2,:);
    points_graph.ZData=points(3,:);
elseif ~isempty(points)
    points_graph=plot3(points(1,:),points(2,:),points(3,:), 'bx', 'MarkerSize', 30);
end

% if isempty(point_graph_old) % if it's the first click (i.e. no previous point to delete)
%     
%     % highlight the selected point
%     
%     point_graph_old = plot3(selectedPoint(1,:), selectedPoint(2,:), ...
%         selectedPoint(3,:), 'bx', 'MarkerSize', 30); 
%     set(point_graph_old,'Tag','pt'); % set its Tag property for later use   
% 
% else % if it is not the first click
% 
%     delete(point_graph_old); % delete the previously selected point
%     
%     % highlight the newly selected point
%     point_graph_old = plot3(selectedPoint(1,:), selectedPoint(2,:), ...
%         selectedPoint(3,:), 'r.', 'MarkerSize', 20);  
%     set(point_graph_old,'Tag','pt');  % set its Tag property for later use
% 
% end

fprintf('you clicked on point number %d\r\n', pointCloudIndex);
fprintf('the point is located at:\r\n %f, %f, %f \r\n',selectedPoint(1),selectedPoint(2),selectedPoint(3));
