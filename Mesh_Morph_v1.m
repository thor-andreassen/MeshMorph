%% clearing
clear
close all
clc

%% constant parameters
% f=10000;
f=10;
p=10;

%% load geom
% [faces,vertices]=stlRead2('Right_Cartilage_FemurDistal.stl');

% % %% create surf mesh
% % b_min=min(vertices);
% % b_max=max(vertices);
% %
% % [node,elem,face]=surf2mesh(vertices,faces,b_min,b_max,.25,0.1);
% %
% % %% plot
% %
% % % hm=plotmesh(node,face,elem)
% % scatter3(node(:,1),node(:,2),node(:,3))
% %


fig_mesh=figure();
fig_error=figure();

%% load target mesh
stl_path=['D:\Files\DU\Research\Mesh morphin code\iso2mesh\test geoms\'];
target_geom_path='VHF_Right_Bone_Femur_smooth.stl';
source_geom_path='VHM_Right_Bone_Femur_smooth_2.stl';
[target.faces,target.nodes]=stlRead2([stl_path,target_geom_path]);
[source.faces,source.nodes]=stlRead2([stl_path,source_geom_path]);


%% reduce source
[source.faces_reduce,source.nodes_reduce]=reducepatch(source.faces,source.nodes,.95);
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
    if mod(counter,2)==0
        [Idx_n,D_n] = knnsearch(source_nodes,target_nodes,'K',1);
        temp_dist=[Idx_n,D_n];
        temp_dist=sortrows(temp_dist,2,'descend');
        new_indices=unique(temp_dist(:,1),'stable');
        new_knot_indices=new_indices(1:p);
        K=[K;source_nodes(new_knot_indices,:)+rand(p,3)*2];
        
    else
        [Idx_n,D_n] = knnsearch(target_nodes,source_nodes,'K',1);
        temp_dist=[Idx_n,D_n];
        temp_dist=sortrows(temp_dist,2,'descend');
        new_indices=unique(temp_dist(:,1),'stable');
        new_knot_indices=new_indices(1:p);
        K=[K;target_nodes(new_knot_indices,:)+rand(p,3)*2];
    end
    
    rand_nodes_x=rand(p,1)*bb_dif(1)+bb_min(1);
    rand_nodes_y=rand(p,1)*bb_dif(2)+bb_min(2);
    rand_nodes_z=rand(p,1)*bb_dif(3)+bb_min(3);
    rand_nodes=[rand_nodes_x,rand_nodes_y,rand_nodes_z];
    K=[K;rand_nodes];
    
    error=max(D_n)
    error_hist=[error_hist,error];
    
    figure(fig_error);
    plot(error_hist);
    
    source_nodes=new_source_nodes;
    
    counter=counter+1;
    
    f=.9*f;
    if mod(counter,50)==0
%         current_rand_ind=randperm(size(K,1));
        K=K0;
%         K=K(current_rand_ind(1:27),:);
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
