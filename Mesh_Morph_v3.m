%% clearing
clear
close all
clc



%% load target mesh
total_time=tic;

stl_path=['C:\Users\Thor.Andreassen\Desktop\Thor Personal Folder\Research\Laxity FE Model\S192803_FE_Model\Iterative Alignment Check\MeshMorph\S193761_Morph_bones\Tibia\'];
target_path=[stl_path,'Target Geom\'];
source_path=[stl_path,'Source Geom\'];
site_path=[stl_path,'Site Geom\'];

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
[source.nodes,M]=ICP_finite(target.nodes, source.nodes, Options);
source.nodes_affine=source.nodes;

% Options.TolX=.0001;
% Options.TolP=.0001;
Options.Registration='Affine';
[source.nodes,M]=ICP_finite(target.nodes, source.nodes, Options);


%% reduce target mesh
[target.faces_reduce,target.nodes_reduce]=reducepatch(target.faces,target.nodes,.4);
[source.faces_reduce,source.nodes_reduce]=reducepatch(source.faces,source.nodes,.4);



%% show scaled nodes
figure();
plot3(source.nodes_orig(:,1),source.nodes_orig(:,2),source.nodes_orig(:,3),'go');
hold on
plot3(source.nodes(:,1),source.nodes(:,2),source.nodes(:,3),'bo');
plot3(target.nodes(:,1),target.nodes(:,2),target.nodes(:,3),'ro');

%% morph large f
params.max_iterations=10;
params.want_plot=1;
params.beta=-.99;
params.scale=.95;
params.dist_threshold=25;
params.dist_threshold_scale=.99;
params.scale_scale=.9;
params.knots_scale=1.1;
params.beta_scale=.99;
params.smooth=10;
params.normal_scale=10;
params.smooth_decay=1;
    
[source.nodes_deform]= pointCloudMorph_v3(target.nodes_reduce,source.nodes_reduce,params,target.faces_reduce,source.faces_reduce);


%% smooth mesh
figure();
smooth_mesh.vertices=source.nodes_deform;
smooth_mesh.faces=source.faces_reduce;
FV2=smoothpatch(smooth_mesh,0,30);
patch('Faces',FV2.faces,'Vertices',FV2.vertices,'FaceColor','r','EdgeAlpha',.3);

source.nodes_deform=FV2.vertices;

%% morph small f
params.max_iterations=4;
params.want_plot=1;
params.beta=-.9;
params.scale=.5;
params.dist_threshold=15;
params.dist_threshold_scale=.99;
params.scale_scale=.9;
params.knots_scale=1.1;
params.beta_scale=.99;
params.smooth=1;
    
[source.nodes_deform]= pointCloudMorph_v3(target.nodes_reduce,source.nodes_deform,params,target.faces_reduce,source.faces_reduce);

%% smooth mesh
figure();
smooth_mesh.vertices=source.nodes_deform;
smooth_mesh.faces=source.faces_reduce;
FV2=smoothpatch(smooth_mesh,0,30);
patch('Faces',FV2.faces,'Vertices',FV2.vertices,'FaceColor','r','EdgeAlpha',.3);

source.nodes_deform=FV2.vertices;








%% deform original mesh
node_deform=source.nodes_deform-source.nodes_reduce;
model=newgrnn(source.nodes_reduce',node_deform',10);

new_deform=sim(model,source.nodes');
source.nodes=source.nodes+new_deform';





%% morph small f
params.max_iterations=20;
params.want_plot=1;
params.beta=-.9;
params.scale=.5;
params.dist_threshold=15;
params.dist_threshold_scale=.99;
params.scale_scale=.9;
params.knots_scale=1.1;
params.beta_scale=.99;
params.smooth=10;
params.smooth_decay=.95;
    
[source.nodes]= pointCloudMorph_v3(target.nodes,source.nodes,params,target.faces,source.faces);

%% smooth mesh
figure();
smooth_mesh.vertices=source.nodes;
smooth_mesh.faces=source.faces;


[smooth_mesh.vertices]=improveTriMeshQuality(smooth_mesh.faces,smooth_mesh.vertices,2,2,.01);
patch('Faces',smooth_mesh.faces,'Vertices',smooth_mesh.vertices,'FaceColor','r','EdgeAlpha',.3);

source.nodes=smooth_mesh.vertices;


%% smooth
figure();
smooth_mesh.vertices=source.nodes;
smooth_mesh.faces=source.faces;
FV2=smoothpatch(smooth_mesh,0,1);
patch('Faces',FV2.faces,'Vertices',FV2.vertices,'FaceColor','r','EdgeAlpha',.3);

source.nodes=FV2.vertices;







% % % %% morph small f
% % % params.max_iterations=5;
% % % params.want_plot=1;
% % % params.beta=-.9;
% % % params.scale=.5;
% % % params.dist_threshold=15;
% % % params.dist_threshold_scale=.99;
% % % params.scale_scale=.9;
% % % params.knots_scale=1.1;
% % % params.beta_scale=.99;
% % % params.smooth=5;
% % % params.normal_scale=1;
% % %     
% % % [source.nodes]= pointCloudMorph_v3(target.nodes,source.nodes,params,target.faces,source.faces);

% % % %% smooth mesh
% % % figure();
% % % smooth_mesh.vertices=source.nodes;
% % % smooth_mesh.faces=source.faces;
% % % 
% % % 
% % % [smooth_mesh.vertices]=improveTriMeshQuality(smooth_mesh.faces,smooth_mesh.vertices,2,1,.01);
% % % patch('Faces',smooth_mesh.faces,'Vertices',smooth_mesh.vertices,'FaceColor','r','EdgeAlpha',.3);
% % % 
% % % source.nodes=smooth_mesh.vertices;


%% plot final meshes
figure()
target_geom_orig=patch('Faces',target.faces,'Vertices',target.nodes,'FaceColor','r','EdgeAlpha',.2,'FaceAlpha',.4);
hold on
source_geom_orig=patch('Faces',source.faces,'Vertices',source.nodes_orig,'FaceColor','g','EdgeAlpha',.2);
source_geom_fin=patch('Faces',source.faces,'Vertices',source.nodes,'FaceColor','b','EdgeAlpha',.2,'FaceAlpha',.4);
axis equal



%% plot net change
figure();
source_geom_orig=patch('Faces',source.faces,'Vertices',source.nodes,'FaceColor','b','EdgeAlpha',.2);
hold on
source_geom_fin=patch('Faces',source.faces,'Vertices',source.nodes_affine,'FaceColor','g','EdgeAlpha',.2);
segments=createLineSegments(source.nodes_affine,source.nodes);
plot3(segments(:,1),segments(:,2),segments(:,3),'k','LineWidth',5);

%% plot final geometreis
figure()
target_geom_orig=patch('Faces',target.faces,'Vertices',target.nodes,'FaceColor','r','EdgeAlpha',.5,'FaceAlpha',.9);
hold on
source_geom_fin=patch('Faces',source.faces,'Vertices',source.nodes,'FaceColor','b','EdgeAlpha',.5,'FaceAlpha',.9);
axis equal


%% determine net motion

source.nodes_change=source.nodes-source.nodes_affine;
toc(total_time)
%% animate motion
figure('units','normalized','outerposition',[0 0 1 1]);
num_frames=100;

col=vecnorm(source.nodes_change,2,2);
source_geom_morph=patch('Faces',source.faces,'Vertices',source.nodes_affine,'EdgeAlpha',.6,'FaceVertexCData',col,'FaceColor','interp');

colorbar
colormap jet

view([1,1,1]);
pause(1);
for count_frame=1:num_frames
    new_nodes=source.nodes_affine+(count_frame/num_frames)*source.nodes_change;
    source_geom_morph.Vertices=new_nodes;
    view([1,1,1]);
    axis square
    pause(.01)
    
    
end






%% save final mesh
% stlWrite2([stl_path,target_geom_path,'_morph.stl'],source.faces,source.nodes_deform);