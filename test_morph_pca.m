%% clear
clear
close all
clc



%% load geom
stl_path=['C:\Users\Thor.Andreassen\Desktop\Thor Personal Folder\Research\Laxity FE Model\S192803_FE_Model\Iterative Alignment Check\MeshMorph\mesh_geometries\example morphing\morphed\'];



%% get vertices
files=dir([stl_path,'*.stl']);


for count_file=1:numel(files)
    current_name=files(count_file).name;
    if count_file==1
        [source.faces,source.nodes]=stlRead2([stl_path,current_name]);
        mat=zeros(numel(source.nodes),numel(files));
        mat(:,count_file)=reshape(source.nodes,[],1);
    else
        [source.faces,source.nodes]=stlRead2([stl_path,current_name]);
        mat(:,count_file)=reshape(source.nodes,[],1);
    end
end

%% create pca
mat=mat';

[coeff,score,latent,tsquared,explained,mu] = pca(mat,'NumComponents',5,'Economy',true);


%% mean shape
mean_nodes=reshape(mu,[],3);

%% plot geoms
figure();
plot3(mean_nodes(:,1),mean_nodes(:,2),mean_nodes(:,3),'bo');
axis equal


%% animate

PC_index=3;



figure('units','normalized','outerposition',[0 0 1 1]);
PC_scores_mean=mean(score);
PC_scores_std=std(score);

std_range=-3:.1:3;
counter=1;

PC_score=zeros(size(PC_scores_std));

PC_score(PC_index)=1;
PC_scores_new=PC_score.*PC_scores_std+PC_scores_mean;
col=vecnorm(reshape(coeff*PC_scores_new',[],3),2,2);
colorbar
colormap jet
for current_std=std_range
    PC_score(PC_index)=current_std;
    PC_scores_new=PC_score.*PC_scores_std+PC_scores_mean;
    color_data=PC_score.*PC_scores_std;
    new_nodes_tall=mu'+coeff*PC_scores_new';
    new_nodes=reshape(new_nodes_tall,[],3);
    if counter==1
        p1=patch('Faces',source.faces,'Vertices',new_nodes,'FaceVertexCData',col,'FaceColor','interp','EdgeAlpha',.05);
        axis([-200,200,-200,200,-100,700]);
    axis equal
    else
        p1.Vertices=new_nodes;
    end
    
    view([1,1,1]);

    pause(.1)
    counter=counter+1;
end