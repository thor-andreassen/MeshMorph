%% clearing
clear
close all
clc



%% load target mesh
total_time=tic;

stl_path=['C:\Users\Thor.Andreassen\Desktop\Thor Personal Folder\Research\Iterative Alignment Check\MeshMorph\ICP_Morph_Comparison\Scapula STL Morphing\'];
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


%% reduce target mesh
[target.faces_reduce,target.nodes_reduce]=reducepatch(target.faces,target.nodes,.05);
[source.faces_reduce,source.nodes_reduce]=reducepatch(source.faces,source.nodes,.05);

%% perform initial rigid alignment
Options.Registration='Rigid';

source.nodes_orig=source.nodes;
[source.nodes_reduce,M1]=ICP_finite(target.nodes_reduce, source.nodes_reduce, Options);


Options.Registration='Affine';
[source.nodes_reduce,M2]=ICP_finite(target.nodes_reduce, source.nodes_reduce, Options);

Affine_TransMat=M2*M1;


source.nodes = transformPts(Affine_TransMat,source.nodes);
source.nodes_affine=source.nodes;

%% reduce target mesh
target_pc=pointCloud(target.nodes);
source_pc=pointCloud(source.nodes);
[tform,move_pc]=pcregistercpd(source_pc,target_pc,'MaxIterations',1,'SmoothingWeight',.1);
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


edge_angles=getAllEdgeAngles(source.faces,source.nodes);
geom_temp.faces=source.faces;
geom_temp.vertices=source.nodes;
aspects=zeros(size(geom_temp.faces,1),1);
skewness=zeros(size(geom_temp.faces,1),1);
for count_face=1:size(geom_temp.faces,1)
    nodel=geom_temp.faces(count_face,:);
    face_nodes=geom_temp.vertices(nodel,:);
    [skewness(count_face),aspects(count_face)]=getMeshQuality2(face_nodes,1);
    
end

node_dist_travel=vecnorm(source.nodes-source.nodes_affine,2,2);


save([results_path,'CPD_Scapula.mat'],'surf_distances','haus_distance',...
    'skewness','aspects','edge_angles','node_dist_travel')



%% final motion figure
figure()
patch('Faces',source.faces,'Vertices',source.nodes,'FaceVertexCData',node_dist_travel,'FaceColor','interp','EdgeAlpha',.3);
c=jet(1000);
colormap(c(125:875,:));
colorbar
caxis([0,25]);
axis off
view ([0,1,0])
axis equal