%% clearing
clear
close all
clc



%% load target mesh
total_time=tic;

stl_path=['C:\Users\Thor.Andreassen\Desktop\Thor Personal Folder\Research\Iterative Alignment Check\MeshMorph\S193761_Morph_bones\Tibia\'];
results_path=['C:\Users\Thor.Andreassen\Desktop\Thor Personal Folder\Research\Iterative Alignment Check\MeshMorph\ICP_Morph_Comparison\'];

target_path=[stl_path,'Target Geom\'];
source_path=[stl_path,'Source Geom\'];
site_path=[stl_path,'Site Geom\'];
landmark_path=[stl_path,'Landmarks\'];
convert_to_mm=0;

%% load target geometries
files=dir([target_path,'*.stl']);
[target.faces,target.nodes]=stlRead2([target_path,files(1).name]);


%% load source goemetries
files=dir([source_path,'*.stl']);
[source.faces,source.nodes]=stlRead2([source_path,files(1).name]);
% load('femur_test')


%% perform initial rigid alignment
Options.Registration='Rigid';

source.nodes_orig=source.nodes;
[source.nodes,M1]=ICP_finite(target.nodes, source.nodes, Options);


% Options.TolX=.0001;
% Options.TolP=.0001;
Options.Registration='Affine';
[source.nodes,M2]=ICP_finite(target.nodes, source.nodes, Options);
source.nodes_affine=source.nodes;

Affine_TransMat=M2*M1;

%% reduce target mesh
target_pc=pointCloud(target.nodes);
source_pc=pointCloud(source.nodes);
[tform,move_pc]=pcregistercpd(source_pc,target_pc,'MaxIterations',2);
source.nodes=move_pc.Location;

time_total=toc(total_time)
%% similarity metrics
inputs.faces=target.faces;
inputs.nodes=target.nodes;
pts=source.nodes;

[distances,project_pts,outside]=fastPoint2TriMesh(inputs,pts,1);
surf_distances=abs(distances);

haus_distance=getHausdorffDistance(source.nodes,target.nodes);
figure()
cdfplot(surf_distances)
hold on
cdfplot(haus_distance)
legend({'Surface Project Distance','Hausdorff Distance'});

save([results_path,'CPD_Tibia.mat'],'surf_distances','haus_distance')