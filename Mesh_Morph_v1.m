%% clearing
clear
close all
clc



%% load target mesh
stl_path=['C:\Users\Thor.Andreassen\Desktop\Thor Personal Folder\Research\Laxity FE Model\S192803_FE_Model\Iterative Alignment Check\MeshMorph\mesh_geometries\'];
% target_geom_path='VHF_Right_Bone_Femur_smooth.stl';
% source_geom_path='VHM_Right_Bone_Femur_smooth_2.stl';
% [target.faces,target.nodes]=stlRead2([stl_path,target_geom_path]);
% [source.faces,source.nodes]=stlRead2([stl_path,source_geom_path]);
% load('femur_test')



% target_geom_path='S193761_Cart_Femur_align.inp';
% source_geom_path='S192803_Cart_Femur_align.inp';
% [target.nodes,target.elems,target.elems_renum]=READ_MESH_NUMS_AJC([stl_path,target_geom_path]);
% [source.nodes,source.elems,source.elems_renum]=READ_MESH_NUMS_AJC([stl_path,source_geom_path]);
% save('input_data.mat','source','target')
load('input_data.mat','source','target')

%% get mesh surfaces
[target_face_outer_surf,target_face_list,target_nodes_surf,target_nodes_outer_surf_coords,target_inner_nodes,target_nodes_inner_coords]=get3DElementOuterSurface(target.elems_renum(:,2:end),target.nodes(:,2:end));
[source_face_outer_surf,source_face_list,source_nodes_surf,source_nodes_outer_surf_coords,source_inner_nodes,source_nodes_inner_coords]=get3DElementOuterSurface(source.elems_renum(:,2:end),source.nodes(:,2:end));

%% renumber outer faces
[target_faces_renumber,target_nodes_renumber,target_node_correspondance_list]=renumberFacesandNodesSubset(target_face_outer_surf,target.nodes(:,2:end));
target_tri_elems=splitQuadsToTries(target_faces_renumber);
% patch('Faces',target_tri_elems,'Vertices',target_nodes_renumber,'FaceColor','r','EdgeAlpha',.3);

% hold on

[source_faces_renumber,source_nodes_renumber,source_node_correspondance_list]=renumberFacesandNodesSubset(source_face_outer_surf,source.nodes(:,2:end));
source_tri_elems=splitQuadsToTries(source_faces_renumber);
% patch('Faces',source_tri_elems,'Vertices',source_nodes_renumber,'FaceColor','b','EdgeAlpha',.3)

source_faces_plot=getHexorTetFaces(source.elems_renum(:,2:end));
%% reduce source
target.nodes_orig=target.nodes;
target.nodes=target_nodes_renumber;
target.faces=target_tri_elems;

source.nodes_orig=source.nodes;
source.nodes=source_nodes_renumber;
source.faces=source_tri_elems;
source.nodes_reduce=source_nodes_renumber;
source.faces_reduce=source_tri_elems;


%% morph large f
close all
params.f=100000;
params.new_knots_per_iter=20;
params.f_decay=.999;
params.max_iterations=20;
params.want_plot=1;
params.rand_mult=.1;
params.rand_decay=.5;
params.include_rand_knots=0;
params.knot_reset_iter=50;
params.target_source_switch_iter=2;
params.start_target=0;
params.initial_knots=3;
params.d_min=5;
params.beta=-.1;
params.scale=.5;
params.dist_threshold=30;

[source.nodes_deform]= pointCloudMorph_v2(target.nodes,source.nodes,params);


%% smooth mesh
smooth_mesh.vertices=source.nodes_deform;
smooth_mesh.faces=source.faces;
FV2=smoothpatch(smooth_mesh,0,2);
patch('Faces',FV2.faces,'Vertices',FV2.vertices,'FaceColor','r','EdgeAlpha',.3);

source.nodes_deform=FV2.vertices;

%% morph small f
% params.f=10000;
% params.new_knots_per_iter=20;
% params.f_decay=.999;
% params.target_source_switch_iter=2;
% params.start_target=0;
% params.max_iterations=10;
% params.include_rand_knots=0;
% params.rand_mult=0.1;
% params.initial_knots=4;
% params.d_min=1;
% params.dist_threshold=15;
% 
% params.beta=-.9;
% params.scale=.5;
% [source.nodes_deform]= pointCloudMorph_v2(target.nodes,source.nodes_deform,params);




%% plot final meshes
figure()
target_geom_orig=patch('Faces',target.faces,'Vertices',target.nodes,'FaceColor','r','EdgeAlpha',.2);
hold on
source_geom_orig=patch('Faces',source.faces,'Vertices',source.nodes_deform,'FaceColor','b','EdgeAlpha',.2);
source_geom_fin=patch('Faces',source.faces,'Vertices',source.nodes,'FaceColor','g','EdgeAlpha',.2);
axis equal


%% get net change
source.change=source.nodes_deform-source.nodes;
deform_net=newgrnn(source.nodes',source.change');

%% get new points
source.nodes_deform_total=[sim(deform_net,source.nodes_orig(:,2:end)')]'+source.nodes_orig(:,2:end);


%% plot points

% scatter3(source.nodes_deform(:,1),source.nodes_deform(:,2),source.nodes_deform(:,3));
figure()
target_geom_orig=patch('Faces',target.faces,'Vertices',target.nodes,'FaceColor','r','EdgeAlpha',.2);

hold on
source_geom_deform=patch('Faces',source_faces_plot,'Vertices',source.nodes_deform_total,'FaceColor','m');
axis equal

%% deform inner nodes
b_min=min(target.nodes);
b_max=max(target.nodes);

[target_tet_node,target_tet_elem,target_tet_face]=surf2mesh(target.nodes,target.faces,b_min,b_max,.5,0.05);


%% morph full mesh nodes
params.f=5;
params.new_knots_per_iter=20;
params.f_decay=1;
params.target_source_switch_iter=2;
params.start_target=0;
params.max_iterations=5;
params.include_rand_knots=0;
params.rand_mult=0.05;
params.initial_knots=3;
params.knot_reset_iter=5;
params.d_min=0.1;
params.beta=-.01;
params.dist_threshold=5;
[source.nodes_deform_total]= pointCloudMorph_v2(target_tet_node,source.nodes_deform_total,params);

%% plot points

% scatter3(source.nodes_deform(:,1),source.nodes_deform(:,2),source.nodes_deform(:,3));
figure()
target_geom_orig=patch('Faces',target.faces,'Vertices',target.nodes,'FaceColor','r','EdgeAlpha',.2,'FaceAlpha',.3);

hold on
source_geom_deform=patch('Faces',source_faces_plot,'Vertices',source.nodes_deform_total,'FaceColor','m','EdgeAlpha',.2,'FaceAlpha',.3);
axis equal



%% coherent point drift method
% [registered]=nonrigidICPv1(target.nodes,source_nodes,target.faces,source.faces_reduce,50,1);

%% save final mesh
% stlWrite2([stl_path,'VHM_Right_Bone_Femur_smooth_3.stl'],source.faces_reduce,registered);