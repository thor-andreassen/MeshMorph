%% Example code for 2D mesh Morphing and Prediction of sites
% Written by Thor Andreassen
% University of Denver
% 12/1/22

% Edited by Thor E. Andreassen, PhD
% 3/5/24





%%
%
%
%%____________________________CODE______________________________
%% clearing
clear
close all
clc

%% load target mesh
total_time=tic;

stl_path=['C:\Users\Thor.Andreassen\Desktop\Thor Personal Folder\Research\In Vivo Modeling\U01_S05\Morphing\TibFib\'];
results_path=[stl_path,'Results\'];

target_path=[stl_path,'Target Geom\'];
source_path=[stl_path,'Source Geom\'];
site_path=[stl_path,'Site Geom\'];
landmark_path=[stl_path,'Landmarks\'];
use_known_align=0;

%% load target geometries
files=dir([target_path,'*.stl']);
[target.faces,target.nodes]=stlRead2([target_path,files(1).name]);


%% load source goemetries
files=dir([source_path,'*.stl']);
[source.faces,source.nodes]=stlRead2([source_path,files(1).name]);
% load('femur_test')

%% reduce target mesh
[target.faces_reduce,target.nodes_reduce]=reducepatch(target.faces,target.nodes,.25);
[source.faces_reduce,source.nodes_reduce]=reducepatch(source.faces,source.nodes,.25);

% target.faces_reduce=target.faces;
% target.nodes_reduce=target.nodes;
% source.faces_reduce=source.faces;
% source.nodes_reduce=source.nodes;

%% perform initial rigid alignment
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

source.nodes_orig=source.nodes;
[source.nodes_reduce,M1]=ICP_finite(target.nodes_reduce, source.nodes_reduce, Options);


Options.Registration='Affine';
[source.nodes_reduce,M2]=ICP_finite(target.nodes_reduce, source.nodes_reduce, Options);

Affine_TransMat=M2*M1*M0;


source.nodes = transformPts(Affine_TransMat,source.nodes);
source.nodes_affine=source.nodes;





% reduce target mesh
% [target.faces_reduce,target.nodes_reduce]=reducepatch(target.faces,target.nodes,.05);
% [source.faces_reduce,source.nodes_reduce]=reducepatch(source.faces,source.nodes,.05);


%% show scaled nodes
figure();
scatter3(source.nodes_orig(:,1),source.nodes_orig(:,2),source.nodes_orig(:,3),1,'g');
hold on
scatter3(source.nodes(:,1),source.nodes(:,2),source.nodes(:,3),3,'b');
scatter3(target.nodes(:,1),target.nodes(:,2),target.nodes(:,3),3,'r');
axis equal
axis off

%% morph large f
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

%% morph small f
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

%% deform original mesh
node_deform=source.nodes_deform-source.nodes_reduce;
model_orig=newgrnn(source.nodes_reduce',node_deform',10);

new_deform=sim(model_orig,source.nodes');
source.nodes=source.nodes+new_deform';

%% plot deformations
% figure()
% temp_deform=new_deform';
% vector_mesh=figure();
% p6=plot3(source.nodes(:,1),source.nodes(:,2),source.nodes(:,3),'bo','MarkerSize',1);
% hold on
% q2=quiver3(source.nodes(:,1),source.nodes(:,2),source.nodes(:,3),...
%     temp_deform(:,1),temp_deform(:,2),temp_deform(:,3),2.5,'k');
% axis off
% axis equal
% view([0,1,0])
% disp('test')


% temp_deform=vecnorm(new_deform',2,2);
% patch_mesh=figure();
% p10=patch('Faces',source.faces,'Vertices',source.nodes,'FaceVertexCData',temp_deform,'FaceColor','interp','EdgeAlpha',.3);
% axis off
% axis equal
% view([0,1,0])
% c=jet(1000);
% colormap(c(125:875,:));
% colorbar
% caxis([0,max(temp_deform)]);
% axis equal
% disp('test')




%% morph small f
params.max_iterations=20; % normally 10-20
params.want_plot=1;
params.scale=.5;
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


%% smooth
% figure();
% smooth_mesh.vertices=source.nodes;
% smooth_mesh.faces=source.faces;
% FV2=smoothpatch(smooth_mesh,0,1);
% patch('Faces',FV2.faces,'Vertices',FV2.vertices,'FaceColor','r','EdgeAlpha',.3);
% 
% source.nodes=FV2.vertices;







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
% segments=createLineSegments(source.nodes_affine,source.nodes);
% plot3(segments(:,1),segments(:,2),segments(:,3),'k','LineWidth',5);

%% plot final geometreis
figure()
target_geom_orig=patch('Faces',target.faces,'Vertices',target.nodes,'FaceColor','r','EdgeAlpha',.5,'FaceAlpha',.9);
hold on
source_geom_fin=patch('Faces',source.faces,'Vertices',source.nodes,'FaceColor','b','EdgeAlpha',.5,'FaceAlpha',.9);
axis equal


%% determine net motion

source.nodes_change=source.nodes-source.nodes_affine;


%% final deformation model
model_final=newgrnn(source.nodes_affine',source.nodes_change',1);


%% animate motion
v=VideoWriter([results_path,'morph_animation.avi']);
open(v);
fig_anim=figure('units','normalized','outerposition',[0 0 1 1]);
num_frames=100;

col=vecnorm(source.nodes_change,2,2);
source_geom_morph=patch('Faces',source.faces,'Vertices',source.nodes_affine,'EdgeAlpha',.6,'FaceVertexCData',col,'FaceColor','interp');

colorbar
colormap jet

view([1,0,0]);

axis equal
axis off
pause(1);

for count_frame=1:num_frames
    new_nodes=source.nodes_affine+(count_frame/num_frames)*source.nodes_change;
    source_geom_morph.Vertices=new_nodes;
    view([1,1,0]);
    frame_val=getframe(fig_anim);
    writeVideo(v,frame_val);
    pause(.01)
    
    
end
close(v);

%% load landmarks


files=dir([landmark_path,'*.xlsx']);
% landmark.orig=csvread([landmark_path,files(1).name]);
temp_node=readtable([landmark_path,files(1).name]);
landmark.orig=table2array(temp_node(:,1:3));
    landmark.deform=applyMorphToNodes(landmark.orig,Affine_TransMat,model_final);
new_table=temp_node;
new_table{:,1:3}=landmark.deform;
% new_table=renamevars(new_table,1:width(new_table),{'landmark','x','y','z'});
writetable(new_table,[results_path,files(1).name])
% csvwrite([results_path,files(1).name],landmark.deform);

%% save morphing data
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
for count_site=1:length(files)
    [site.faces,site.nodes]=stlRead2([site_path,files(count_site).name]);
    subplot(1,2,1)
    hold on
    p_site_orig{count_site}=patch('Faces',site.faces,'Vertices',site.nodes,'FaceColor',colors(count_site,:),'EdgeAlpha',.3);

    
    site.nodes_deform=applyMorphToNodes(site.nodes,Affine_TransMat,model_final);
%     temp_nodes=[site.nodes,ones(size(site.nodes,1),1)];
%     affine_nodes=[Affine_TransMat*temp_nodes']';
%     site.nodes_deform=affine_nodes(:,1:3);
%     
%     new_deform=sim(model_final,site.nodes_deform');
%     site.nodes_deform=site.nodes_deform+new_deform';
    new_geom.faces=site.faces;
    new_geom.vertices=site.nodes_deform;
    subplot(1,2,2)
    hold on
    p_site_new{count_site}=patch('Faces',site.faces,'Vertices',site.nodes_deform,'FaceColor',colors(count_site,:),'EdgeAlpha',.3);
    try
        stlWrite2([results_path,files(count_site).name,'_Morph.stl'],site.faces,site.nodes_deform);
    catch
        stlwrite([results_path,files(count_site).name,'_Morph.stl'],new_geom);
    end
end

saveas(morph_fig,[results_path,'Morph_Figure.png']);
saveas(morph_fig,[results_path,'Morph_Figure.fig']);
%% save morphing meshes
files=dir([target_path,'*.stl']);
target_filename=files(1).name;
try
    stlWrite2([results_path,target_filename,'_Morph.stl'],source.faces,source.nodes);
catch
    new_geom.faces=source.faces;
    new_geom.vertices=source.nodes;
    stlwrite([results_path,target_filename,'_Morph.stl'],new_geom);
end


save([results_path,'Morphing_Parameters.mat'],'Affine_TransMat','source','target',...
    'model_final','landmark');


%% save final mesh
% stlWrite2([stl_path,target_geom_path,'_morph.stl'],source.faces,source.nodes_deform);

%% computer similarity metrics

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


save([results_path,'Morph_Similarity.mat'],'surf_distances','haus_distance',...
    'skewness','aspects','edge_angles','node_dist_travel')

%% plot net deformation


figure()
patch('Faces',source.faces,'Vertices',source.nodes,'EdgeAlpha',.6,'FaceVertexCData',surf_distances,'FaceColor','interp');
c=jet(1000);
colormap(c(125:875,:));
colorbar
caxis([0,10]);
axis off
view ([1,-1,1])
axis equal

figure()
bone_color=[0.992156863212585,0.917647063732147,0.796078443527222];
patch('Faces',source.faces,'Vertices',source.nodes,'FaceColor',bone_color);
axis off
view ([1,-1,1])
axis equal


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



%% final motion figure
figure()
patch('Faces',target.faces,'Vertices',target.nodes,'FaceColor',bone_color,'EdgeAlpha',.3);
axis off
view ([0,1,0])
axis equal


figure()
patch('Faces',source.faces,'Vertices',source.nodes_affine,'FaceColor',bone_color,'EdgeAlpha',.3);
axis off
view ([0,1,0])
axis equal

figure()
patch('Faces',source.faces,'Vertices',source.nodes,'FaceColor',bone_color,'EdgeAlpha',.3);
axis off
view ([0,1,0])
axis equal