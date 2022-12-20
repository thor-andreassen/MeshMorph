%% clearing
clear
close all
clc

%% constant parameters
% f=10000;
% p=2;
% f=10;
% p=2;


% f=1000;
% p=10;

f=10;
p=2;
rand_mult=1;

%% Initialize Figures
fig_mesh=figure();
fig_error=figure();

%% load target mesh
stl_path=['C:\Users\Thor.Andreassen\Desktop\Thor Personal Folder\Research\Laxity FE Model\S192803_FE_Model\Iterative Alignment Check\MeshMorph\mesh_geometries\'];
% target_geom_path='VHF_Right_Bone_Femur_smooth.stl';
% source_geom_path='VHM_Right_Bone_Femur_smooth_2.stl';
% % [target.faces,target.nodes]=stlRead2([stl_path,target_geom_path]);
% % [source.faces,source.nodes]=stlRead2([stl_path,source_geom_path]);


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
%% plot original bones
% target_geom_orig=patch('Faces',target.faces,'Vertices',target.nodes,'FaceColor','r','EdgeAlpha',.2);
% hold on
% source_geom_orig=patch('Faces',source.faces_reduce,'Vertices',source.nodes_reduce,'FaceColor','b','EdgeAlpha',.2);


%% mesh morphin initialization
target_nodes=target.nodes;
source_nodes_0=source.nodes_reduce;
source_nodes=source_nodes_0;
bb_min=min(source.nodes_reduce);
bb_max=max(source.nodes_reduce);
bb_dif=abs(bb_max-bb_min)*.3;
bb_min=bb_min-bb_dif;
bb_max=bb_max+bb_dif;

x_bb=linspace(bb_min(1),bb_max(1),3);
y_bb=linspace(bb_min(2),bb_max(2),3);
z_bb=linspace(bb_min(3),bb_max(3),3);

[X,Y,Z]=meshgrid(x_bb,y_bb,z_bb);
X=reshape(X,[],1);
Y=reshape(Y,[],1);
Z=reshape(Z,[],1);

bb_dif=bb_max-bb_min;
K0=[X,Y,Z];

K=K0;
% hold on
% scatter3(K(:,1),K(:,2),K(:,3),'go')

phi_func=@(r,sigma) exp(-(r.^2)/(2*sigma.^2));


%% plot original



%% main loop for mesh morphing
tic
counter=1;
error_hist=[];
switcher=0;
while counter<=50
    
    % get sigmas
    [~,D] = knnsearch(K,K,'K',4);
    sigmas=mean(D(:,2:end),2)*f;
    
    
    % determine nearest points
    [Idx_T] = knnsearch(source_nodes,target_nodes,'K',1);
    U=target_nodes;
    d=source_nodes(Idx_T,:);
    
    rbf_mat_test=rbfFuncMat(phi_func,d,K,sigmas);
    
    % fit rbf A coefficients
    A=rbf_mat_test\U;
    
    % determine new node locations
    rbf_mat_full=rbfFuncMat(phi_func,source_nodes,K,sigmas);
    new_source_nodes=rbf_mat_full*A;
    
    figure(fig_mesh)
    clf(fig_mesh)
    scatter3(source_nodes_0(:,1),source_nodes_0(:,2),source_nodes_0(:,3),'go');
    hold on
    scatter3(source_nodes(:,1),source_nodes(:,2),source_nodes(:,3),'bo');
    scatter3(target_nodes(:,1),target_nodes(:,2),target_nodes(:,3),'ro');
    scatter3(K(:,1),K(:,2),K(:,3),'kx')
    axis equal
    pause(.01);
    
    
    %determine new knot lcoations
%     if mod(counter,2)==0
if switcher==0
        [Idx_n,D_n] = knnsearch(source_nodes,target_nodes,'K',1);
        temp_dist=[Idx_n,D_n];
        temp_dist=sortrows(temp_dist,2,'descend');
        new_indices=unique(temp_dist(:,1),'stable');
        new_knot_indices=new_indices(1:p);
        K=[K;source_nodes(new_knot_indices,:)+rand(p,3)*rand_mult];
        
    else
        [Idx_n,D_n] = knnsearch(target_nodes,source_nodes,'K',1);
        temp_dist=[Idx_n,D_n];
        temp_dist=sortrows(temp_dist,2,'descend');
        new_indices=unique(temp_dist(:,1),'stable');
        new_knot_indices=new_indices(1:p);
        K=[K;target_nodes(new_knot_indices,:)+rand(p,3)*rand_mult];
    end
    
%     rand_nodes_x=rand(p,1)*bb_dif(1)+bb_min(1);
%     rand_nodes_y=rand(p,1)*bb_dif(2)+bb_min(2);
%     rand_nodes_z=rand(p,1)*bb_dif(3)+bb_min(3);
%     rand_nodes=[rand_nodes_x,rand_nodes_y,rand_nodes_z];
%     K=[K;rand_nodes];
    
    error=max(D_n)
    error_hist=[error_hist,error];
    
    figure(fig_error);
    plot(error_hist);
    
    source_nodes=new_source_nodes;
    
    counter=counter+1;
    
    rand_mult=rand_mult*1;
%     f=.95*f;
    f=.995*f;
    if mod(counter,10)==0
%         current_rand_ind=randperm(size(K,1));
        K=K0;
%         K=K(current_rand_ind(1:20),:);
    end
    
%     if mod(counter,2)==0
        if mod(counter,20)==0
        if switcher==1
            switcher=0;
        else
            switcher=1;
        end
    end
    
end
toc

%% plot final meshes
figure()
target_geom_orig=patch('Faces',target.faces,'Vertices',target.nodes,'FaceColor','r','EdgeAlpha',.2);
hold on
source_geom_orig=patch('Faces',source.faces_reduce,'Vertices',source_nodes,'FaceColor','b','EdgeAlpha',.2);
source_geom_fin=patch('Faces',source.faces_reduce,'Vertices',source.nodes_reduce,'FaceColor','g','EdgeAlpha',.2);
axis equal


%% get net change
source.change=source_nodes-source.nodes;
deform_net=newrbe(source.nodes',source.change',100);

%% get new points
source.nodes_deform=[sim(deform_net,source.nodes_orig(:,2:end)')]'+source.nodes_orig(:,2:end);


%% plot points

% scatter3(source.nodes_deform(:,1),source.nodes_deform(:,2),source.nodes_deform(:,3));
figure()
target_geom_orig=patch('Faces',target.faces,'Vertices',target.nodes,'FaceColor','r','EdgeAlpha',.2);

hold on
source_geom_deform=patch('Faces',source_faces_plot,'Vertices',source.nodes_deform,'FaceColor','m');
axis equal
%% coherent point drift method
% [registered]=nonrigidICPv1(target.nodes,source_nodes,target.faces,source.faces_reduce,50,1);

%% save final mesh
% stlWrite2([stl_path,'VHM_Right_Bone_Femur_smooth_3.stl'],source.faces_reduce,registered);