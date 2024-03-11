%% Example code for 2D mesh Morphing and Prediction of sites
% Written by Thor Andreassen, PhD
% University of Denver
% Created 12/1/22
% Last Edited 3/10/24


% This is the main example script to create a morphing of a 2D mesh
% geometry to another 2D mesh geometry. The morphing function is contained
% in the pointCloudMorph_v4 function.



%%
%
%
%%____________________________CODE______________________________
%% clearing
clear
close all
clc

%% OVERALL PATHS
total_time=tic;

default_path=[pwd,'\Example 2D Structure\'];
base_path=uigetdir(default_path);

%% check valid filepath
if base_path==0
    error('Choose a Valid path location');
elseif base_path(end)~='\'
    base_path=[base_path,'\'];
end

%% Create Storage Paths
results_path=[base_path,'Results\'];
target_path=[base_path,'Target Geom\'];
source_path=[base_path,'Source Geom\'];
site_path=[base_path,'Site Geom\'];
landmark_path=[base_path,'Landmarks\'];



use_known_align=0;



%% load target geometries
files=dir([target_path,'*.stl']);
[target.faces,target.nodes]=stlRead2([target_path,files(1).name]);


%% load source goemetries
files=dir([source_path,'*.stl']);
[source.faces,source.nodes]=stlRead2([source_path,files(1).name]);
% load('femur_test')

%% Step 1 - Create Reduced Geometries
[target.faces_reduce,target.nodes_reduce]=reducepatch(target.faces,target.nodes,.15);
[source.faces_reduce,source.nodes_reduce]=reducepatch(source.faces,source.nodes,.15);

%% Step 2 - Rigid Alignment
target_pts=[125.231,97.5814,482.127;...
    59.1168,44.6367,112.439;...
    112.968,65.8699,113.966];
source_pts=[423.979,491.863,898.497;...
    494.608,490.195,534.043;...
    436.441,463.113,531.81];
if use_known_align==1
    M0=alignKnownPts(target_pts,source_pts);
    M0=rotateTransMat(M0,3);
else
    M0=eye(4);
end
source.nodes_reduce = transformPts(M0,source.nodes_reduce);
Options.Registration='Rigid';
Options.TolX=.001;

source.nodes_orig=source.nodes;
[source.nodes_reduce,M1]=ICP_finite(target.nodes_reduce, source.nodes_reduce, Options);

%% Step 3 - Affine Alignment
Options.Registration='Affine';
[source.nodes_reduce,M2]=ICP_finite(target.nodes_reduce, source.nodes_reduce, Options);

Initial_Alignment_TransMat=M0;
Rigid_TransMat=M1*M0;
Affine_TransMat=M2*M1*M0;

source.nodes_orig_rigid_align=transformPts(M1*M0,source.nodes);
source.nodes = transformPts(Affine_TransMat,source.nodes);
source.nodes_affine=source.nodes;



%% Plotting - Show the Alignment after Rigid and Affine Transformations
figure();
scatter3(source.nodes_orig(:,1),source.nodes_orig(:,2),source.nodes_orig(:,3),1,'g');
hold on
scatter3(source.nodes(:,1),source.nodes(:,2),source.nodes(:,3),3,'b');
scatter3(target.nodes(:,1),target.nodes(:,2),target.nodes(:,3),3,'r');
axis equal
axis off

%% Step 4 - High Smoothing Morphing of Reduced Geometries
params.max_iterations=10; %normally ~10
params.want_plot=1;
params.scale=.95;
params.smooth=10; % normally 10
params.normal_scale=5;
params.normal_scale_decay=.999;
params.use_parallel=1;
params.smooth_decay=1;

[source.nodes_deform]= pointCloudMorph_v4(target.nodes_reduce,source.nodes_reduce,params,target.faces_reduce,source.faces_reduce);
% [source.nodes_deform]= pointCloudMorph_v4(target.nodes_reduce,source.nodes_reduce,params);


%% smooth mesh
% figure();
% smooth_mesh.vertices=source.nodes_deform;
% smooth_mesh.faces=source.faces_reduce;
% FV2=smoothpatch(smooth_mesh,0,30);
% patch('Faces',FV2.faces,'Vertices',FV2.vertices,'FaceColor','r','EdgeAlpha',.3);
%
% source.nodes_deform=FV2.vertices;

%% Step 5 - Low Smoothing Morphing of Reduced Geometries
params.max_iterations=4; %normally ~4
params.want_plot=1;
params.scale=.5;
params.smooth=1; % normally 1

[source.nodes_deform]= pointCloudMorph_v4(target.nodes_reduce,source.nodes_deform,params,target.faces_reduce,source.faces_reduce);

%% smooth mesh
% figure();
% smooth_mesh.vertices=source.nodes_deform;
% smooth_mesh.faces=source.faces_reduce;
% FV2=smoothpatch(smooth_mesh,0,30);
% patch('Faces',FV2.faces,'Vertices',FV2.vertices,'FaceColor','r','EdgeAlpha',.3);
%
% source.nodes_deform=FV2.vertices;

%% Step 6 - Create initial position of dense mesh from reduce solution
node_deform=source.nodes_deform-source.nodes_reduce;
model_orig=newgrnn(source.nodes_reduce',node_deform',10);

new_deform=sim(model_orig,source.nodes');
source.nodes=source.nodes+new_deform';

%% Step 7 - High Smoothing Morphing of Final Dense Geometries
params.max_iterations=10; % normally 10-20
params.want_plot=1;
params.scale=.75;
params.smooth=10; % normally 10
params.smooth_decay=.95;

[source.nodes]= pointCloudMorph_v4(target.nodes,source.nodes,params,target.faces,source.faces);
% [source.nodes]= pointCloudMorph_v4(target.nodes,source.nodes,params);
time_total=toc(total_time)

%% smooth mesh
% figure();
% smooth_mesh.vertices=source.nodes;
% smooth_mesh.faces=source.faces;
%
%
% [smooth_mesh.vertices]=improveTriMeshQuality(smooth_mesh.faces,smooth_mesh.vertices,2,2,.01);
% patch('Faces',smooth_mesh.faces,'Vertices',smooth_mesh.vertices,'FaceColor','r','EdgeAlpha',.3);
%
% source.nodes=smooth_mesh.vertices;

%% Plotting - plot final meshes
final_overlap_fig=figure()
target_geom_orig=patch('Faces',target.faces,'Vertices',target.nodes,'FaceColor','r','EdgeAlpha',.2,'FaceAlpha',.4);
hold on
source_geom_orig=patch('Faces',source.faces,'Vertices',source.nodes_orig,'FaceColor','g','EdgeAlpha',.2);
source_geom_fin=patch('Faces',source.faces,'Vertices',source.nodes,'FaceColor','b','EdgeAlpha',.2,'FaceAlpha',.4);
axis equal
mkdir([results_path,'Images\']);
try
    saveas(final_overlap_fig,[results_path,'Images\','Overlap_Final_Figure.png']);
    saveas(final_overlap_fig,[results_path,'Images\','Overlap_Final_Figure.fig']);
end

%% Plotting - plot net change
figure();
source_geom_orig=patch('Faces',source.faces,'Vertices',source.nodes,'FaceColor','b','EdgeAlpha',.2);
hold on
source_geom_fin=patch('Faces',source.faces,'Vertices',source.nodes_affine,'FaceColor','g','EdgeAlpha',.2);
% segments=createLineSegments(source.nodes_affine,source.nodes);
% plot3(segments(:,1),segments(:,2),segments(:,3),'k','LineWidth',5);

%% Plotting - Plot final geometries
figure()
target_geom_orig=patch('Faces',target.faces,'Vertices',target.nodes,'FaceColor','r','EdgeAlpha',.5,'FaceAlpha',.9);
hold on
source_geom_fin=patch('Faces',source.faces,'Vertices',source.nodes,'FaceColor','b','EdgeAlpha',.5,'FaceAlpha',.9);
axis equal


%% Step 8 - determine net motion for final morphing solution

source.nodes_change=source.nodes-source.nodes_affine;
source.nodes_change_affine=source.nodes_affine-source.nodes_orig_rigid_align;
source.nodes_change_total=source.nodes-source.nodes_orig_rigid_align;
%% Step 9 - Train Final deformation model for morphing solution
model_final=newgrnn(source.nodes_affine',source.nodes_change',1);


%% Plotting - Animate overall morphing motion
v=VideoWriter([results_path,'morph_animation.avi']);
open(v);
fig_anim=figure('units','normalized','outerposition',[0 0 1 1]);
num_frames=100;

col=vecnorm(source.nodes_change_total,2,2);
original_geom=patch('Faces',source.faces,'Vertices',source.nodes_orig_rigid_align,'EdgeAlpha',.15,'FaceColor','k','FaceAlpha',.2);
hold on
source_geom_morph=patch('Faces',source.faces,'Vertices',source.nodes_orig_rigid_align,'EdgeAlpha',.15,'FaceVertexCData',zeros(size(col)),'FaceColor','interp');

colorbar
colormap jet
caxis([0,max(col)]);
view([1,0,0]);

axis equal
axis off
pause(1);

for count_frame=1:num_frames
    new_nodes=source.nodes_orig_rigid_align+(count_frame/num_frames)*source.nodes_change_total;
    source_geom_morph.Vertices=new_nodes;
    source_geom_morph.FaceVertexCData=col*(count_frame/num_frames);
    view([1,1,0]);
    
    frame_val=getframe(fig_anim);
    writeVideo(v,frame_val);
    pause(.01)


end
close(v);

%% load landmarks


files=dir([landmark_path,'*.csv']);
% landmark.orig=csvread([landmark_path,files(1).name]);
temp_node=readtable([landmark_path,files(1).name]);
landmark.orig=table2array(temp_node(:,2:4));

landmark.deform=applyMorphToNodes(landmark.orig,Affine_TransMat,model_final);

new_table=temp_node;
new_table{:,2:4}=landmark.deform;
new_table=renamevars(new_table,1:width(new_table),{'landmark','x','y','z'});
mkdir([results_path,'Landmarks\']);
writetable(new_table,[results_path,'Landmarks\',files(1).name])

%% Plotting/Saving - Site Geometries and Morphing Landmark/Site Figure

morph_fig=figure()
subplot(1,2,2)
target_geom_orig=patch('Faces',target.faces,'Vertices',target.nodes,'FaceColor',[0.3,0.3,0.3],'EdgeAlpha',0,'FaceAlpha',0.3);
hold on
plot3(landmark.deform(:,1),landmark.deform(:,2),landmark.deform(:,3),'kx','MarkerSize',20);



subplot(1,2,1)
source_geom_fin=patch('Faces',source.faces,'Vertices',source.nodes_orig,'FaceColor',[0.3,0.3,0.3],'EdgeAlpha',0,'FaceAlpha',0.3);
hold on
plot3(landmark.orig(:,1),landmark.orig(:,2),landmark.orig(:,3),'kx','MarkerSize',20);



files=dir([site_path,'*.stl']);
colors=jet(length(files));

if ~isempty(files)
    mkdir([results_path,'Site Geom\']);
end
try
    for count_site=1:length(files)
        [site.faces,site.nodes]=stlRead2([site_path,files(count_site).name]);
        subplot(1,2,1)
        hold on
        p_site_orig{count_site}=patch('Faces',site.faces,'Vertices',site.nodes,'FaceColor',colors(count_site,:),'EdgeAlpha',.3);


        site.nodes_deform=applyMorphToNodes(site.nodes,Affine_TransMat,model_final);
        new_geom.faces=site.faces;
        new_geom.vertices=site.nodes_deform;
        subplot(1,2,2)
        hold on
        p_site_new{count_site}=patch('Faces',site.faces,'Vertices',site.nodes_deform,'FaceColor',colors(count_site,:),'EdgeAlpha',.3);

        [~,old_site_name,~]= fileparts([site_path,files(count_site).name]);
        new_site_name=[old_site_name,'_Morph.stl'];


        try
            stlWrite2([[results_path,'Site Geom\'],new_site_name],site.faces,site.nodes_deform);
        catch
            stlwrite([[results_path,'Site Geom\'],new_site_name],new_geom);
        end
    end
end


try
    saveas(morph_fig,[results_path,'Images\','Morph_Figure.png']);
    saveas(morph_fig,[results_path,'Images\','Morph_Figure.fig']);
end
%% Saving - Save Morphed Geometry and Morphing Solution Parameters
files=dir([target_path,'*.stl']);
target_filename=files(1).name;


[~,old_target_name,~]= fileparts([target_path,target_filename]);
new_morphed_name=[old_target_name,'_Morph.stl'];
mkdir([results_path,'Morphed Geometry\']);
try
    stlWrite2([results_path,'Morphed Geometry\',new_morphed_name],source.faces,source.nodes);
catch
    new_geom.faces=source.faces;
    new_geom.vertices=source.nodes;
    stlwrite([results_path,'Morphed Geometry\',new_morphed_name],new_geom);
end


save([results_path,'Morphing_Parameters.mat'],'Affine_TransMat','source','target',...
    'model_final','landmark','M2','M1','M0','Initial_Alignment_TransMat',...
    'Rigid_TransMat');

%% Plotting/Saving - Morphing Acccuracy Metrics

inputs.faces=target.faces;
inputs.nodes=target.nodes;
pts=source.nodes;

try
    [distances,project_pts,outside]=fastPoint2TriMesh(inputs,pts,1);
    metrics.surf_distances=abs(distances);
end

try
    metrics.haus_distance=getHausdorffDistance(source.nodes,target.nodes);
end
try
    figure()
    cdfplot(metrics.surf_distances)
    hold on
    cdfplot(metrics.haus_distance)
    legend({'Surface Project Distance','Hausdorff Distance'});
end

try
    metrics.edge_angles=getAllEdgeAngles(source.faces,source.nodes);
end

try
    geom_temp.faces=source.faces;
    geom_temp.vertices=source.nodes;
    metrics.aspects=zeros(size(geom_temp.faces,1),1);
    metrics.skewness=zeros(size(geom_temp.faces,1),1);
    for count_face=1:size(geom_temp.faces,1)
        nodel=geom_temp.faces(count_face,:);
        face_nodes=geom_temp.vertices(nodel,:);
        [metrics.skewness(count_face),metrics.aspects(count_face)]=getMeshQuality2(face_nodes,1);

    end
end

try
    metrics.node_dist_travel_rig_to_aff=vecnorm(source.nodes_change_affine,2,2);
    metrics.node_dist_travel_GRNN=vecnorm(source.nodes_change,2,2);
    metrics.node_dist_travel_total=vecnorm(source.nodes_change_total,2,2);
end
save([results_path,'Morph_Similarity.mat'],'metrics','time_total');

%% Plotting - Net Accuracy of Morphing

try
    net_surf_figure=figure();
    patch('Faces',source.faces,'Vertices',source.nodes,'EdgeAlpha',.6,'FaceVertexCData',metrics.surf_distances,'FaceColor','interp','EdgeAlpha',.3);
    c=jet(1000);
    colormap(c(125:875,:));
    colorbar
    caxis([0,max(metrics.surf_distances)]);
    axis off
    view ([1,-1,1])
    axis equal
    saveas(net_surf_figure,[results_path,'Images\','Surf_Distance_Figure.png']);
    saveas(net_surf_figure,[results_path,'Images\','Surf_Distance_Figure.fig']);
end

%% Plotting - final bone position
figure()
bone_color=[0.992156863212585,0.917647063732147,0.796078443527222];
patch('Faces',source.faces,'Vertices',source.nodes,'FaceColor',bone_color);
axis off
view ([1,-1,1])
axis equal




%% Plotting - final motion figure
try
    morphed_color_fig=figure();
    subplot(1,3,1);
    patch('Faces',source.faces,'Vertices',source.nodes_affine,'FaceVertexCData',metrics.node_dist_travel_rig_to_aff,'FaceColor','interp','EdgeAlpha',.3);
    c=jet(1000);
    colormap(c(125:875,:));
    colorbar
    caxis([0,max(metrics.node_dist_travel_rig_to_aff)]);
    axis off
    view ([0,1,0])
    axis equal
    title('Bone after Affine Transformation')

    subplot(1,3,3);
    patch('Faces',source.faces,'Vertices',source.nodes,'FaceVertexCData',metrics.node_dist_travel_GRNN,'FaceColor','interp','EdgeAlpha',.3);
    c=jet(1000);
    colormap(c(125:875,:));
    colorbar
    caxis([0,max(metrics.node_dist_travel_GRNN)]);
    axis off
    view ([0,1,0])
    axis equal
    title('GRNN Morph')

    subplot(1,3,3);
    patch('Faces',source.faces,'Vertices',source.nodes,'FaceVertexCData',metrics.node_dist_travel_total,'FaceColor','interp','EdgeAlpha',.3);
    c=jet(1000);
    colormap(c(125:875,:));
    colorbar
    caxis([0,max(metrics.node_dist_travel_total)]);
    axis off
    view ([0,1,0])
    axis equal
    title('Complete Morphed Source')
    saveas(morphed_color_fig,[results_path,'Images\','Morphed_Contour_Figure.png']);
    saveas(morphed_color_fig,[results_path,'Images\','Morphed_Contour_Figure.fig']);
end

%% final motion figure
try
    final_bone_fig=figure();
    subplot(1,3,1)
    patch('Faces',target.faces,'Vertices',target.nodes,'FaceColor',bone_color,'EdgeAlpha',.3);
    axis off
    view ([0,1,0])
    axis equal
    title('Target Mesh');

    subplot(1,3,2)
    patch('Faces',source.faces,'Vertices',source.nodes_affine,'FaceColor',bone_color,'EdgeAlpha',.3);
    axis off
    view ([0,1,0])
    axis equal
    title ('Source Mesh after Affine')

    subplot(1,3,3)
    patch('Faces',source.faces,'Vertices',source.nodes,'FaceColor',bone_color,'EdgeAlpha',.3);
    axis off
    view ([0,1,0])
    axis equal
    title('Morphed Source Mesh')
    saveas(final_bone_fig,[results_path,'Images\','Final_Geometries_Figure.png']);
    saveas(final_bone_fig,[results_path,'Images\','Final_Geometries_Figure.fig']);
end