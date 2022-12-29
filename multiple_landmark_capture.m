%% clear
clear
close all
clc


%% load data
load('input_data.mat','source','target')



%% choose points source
fig1=figure();
source_nodes=source.nodes(:,2:end);
[~,source_landmark_points]=choose3DPoints(source_nodes);

%% choose points target
fig2=figure();
target_nodes=target.nodes(:,2:end);
[~,target_landmark_points]=choose3DPoints(target_nodes);
%% morph points
beta=-.15;
dist_mat=pdist2(source_landmark_points',target_landmark_points');
c=(diag(dist_mat).^2)';
K=(dist_mat.^2+c).^beta;

P=target_landmark_points-source_landmark_points;

w=K\P';



%% new points
source_target_dist=pdist2(source_nodes,target_landmark_points');
K=(source_target_dist.^2+c).^beta;

source_nodes=source_nodes+K*w;

%% load points
scatter3(target_nodes(:,1),target_nodes(:,2),target_nodes(:,3),'ro');
hold on

scatter3(source_nodes(:,1),source_nodes(:,2),source_nodes(:,3),'bx');

