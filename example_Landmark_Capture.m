%% clearing
clear
close all
clc

%% load data
load('input_data.mat','source','target')

%% click point
global points;
points=[];
global points_index;
points_index=[];
global points_graph;
global point_graph_old


%% temp line graphs
points_graph=plot3([],[],[], 'bx', 'MarkerSize', 30);


h = clickA3DPoint_v2(source.nodes(:,2:end)');


%% get points

test_str=input('Delete Point (1), Finish (0): \r\n');
while test_str~=0
   
    switch test_str
        case 0
            break;
        case 1
            delete(point_graph_old)
            points=points(:,1:(end-1));
            points_index=points_index(1:(end-1));
            if ~isempty(points)
                points_graph.XData=points(1,:);
                points_graph.YData=points(2,:);
                points_graph.ZData=points(3,:);
            else
                points_graph.XData=[];
                points_graph.YData=[];
                points_graph.ZData=[];
            end
    end
    
   test_str=input('Delete Point (1), Finish (0): \r\n');
end