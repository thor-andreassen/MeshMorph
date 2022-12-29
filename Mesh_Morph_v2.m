%% clearing
clear
close all
clc



%% load target mesh
stl_path=['C:\Users\Thor.Andreassen\Desktop\Thor Personal Folder\Research\Laxity FE Model\S192803_FE_Model\Iterative Alignment Check\MeshMorph\mesh_geometries\'];
target_geom_path='VHF_Right_Bone_Femur_smooth.stl';
source_geom_path='VHM_Right_Bone_Femur_smooth_align.stl';
[target.faces,target.nodes]=stlRead2([stl_path,target_geom_path]);
[source.faces,source.nodes]=stlRead2([stl_path,source_geom_path]);
% load('femur_test')


%% perform initial rigid alignment
Options.Registration='Size';

source.nodes_orig=source.nodes;
[source.nodes,M]=ICP_finite(target.nodes, source.nodes, Options);

%% show scaled nodes
figure();
plot3(source.nodes_orig(:,1),source.nodes_orig(:,2),source.nodes_orig(:,3),'go');
hold on
plot3(source.nodes(:,1),source.nodes(:,2),source.nodes(:,3),'bo');
plot3(target.nodes(:,1),target.nodes(:,2),target.nodes(:,3),'ro');

%% morph large f
params.max_iterations=20;
params.want_plot=1;
params.beta=-.9;
params.scale=.5;
params.dist_threshold=50;
params.dist_threshold_scale=.9;
params.scale_scale=.9;
params.knots_scale=1.1;
params.beta_scale=.9;
    
    
[source.nodes_deform]= pointCloudMorph_v2(target.nodes,source.nodes,params);


%% smooth mesh
smooth_mesh.vertices=source.nodes_deform;
smooth_mesh.faces=source.faces;
FV2=smoothpatch(smooth_mesh,0,2);
patch('Faces',FV2.faces,'Vertices',FV2.vertices,'FaceColor','r','EdgeAlpha',.3);

source.nodes_deform=FV2.vertices;

%% morph small f
params.max_iterations=10;
params.dist_threshold=5;
params.beta=-.01;
params.beta_scale=1;
[source.nodes_deform]= pointCloudMorph_v2(target.nodes,source.nodes_deform,params);




%% plot final meshes
figure()
target_geom_orig=patch('Faces',target.faces,'Vertices',target.nodes,'FaceColor','r','EdgeAlpha',.2);
hold on
source_geom_orig=patch('Faces',source.faces,'Vertices',source.nodes_deform,'FaceColor','b','EdgeAlpha',.2);
source_geom_fin=patch('Faces',source.faces,'Vertices',source.nodes,'FaceColor','g','EdgeAlpha',.2);
axis equal



%% plot net change
source_geom_orig=patch('Faces',source.faces,'Vertices',source.nodes_deform,'FaceColor','b','EdgeAlpha',.2);
hold on
source_geom_fin=patch('Faces',source.faces,'Vertices',source.nodes_orig,'FaceColor','g','EdgeAlpha',.2);
segments=createLineSegments(source.nodes_orig,source.nodes_deform);
plot3(segments(:,1),segments(:,2),segments(:,3),'k','LineWidth',5);


%% determine net motion

source.nodes_change=source.nodes_deform-source.nodes_orig;

%% animate motion
figure('units','normalized','outerposition',[0 0 1 1]);
num_frames=10;

col=vecnorm(source.nodes_change,2,2);
source_geom_morph=patch('Faces',source.faces,'Vertices',source.nodes_orig,'EdgeAlpha',.6,'FaceVertexCData',col,'FaceColor','interp');

colorbar
colormap jet

view([0,1,0]);
camroll(90);
pause(1);
for count_frame=1:num_frames
    new_nodes=source.nodes_orig+(count_frame/num_frames)*source.nodes_change;
    source_geom_morph.Vertices=new_nodes;
    view([0,1,0]);
    pause(.1)
    
    
end






%% save final mesh
stlWrite2([stl_path,'VHM_Right_Bone_Femur_smooth_morph.stl'],source.faces,source.nodes_deform);