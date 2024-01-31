%% clearing
clear
close all
clc

%% load data
stl_path=['C:\Users\Thor.Andreassen\Desktop\Thor Personal Folder\Research\Iterative Alignment Check\MeshMorph\S193761_Morph_bones\Tibia\Target Geom\'];

%% load source goemetries
files=dir([stl_path,'*.stl']);
[source.faces,source.nodes]=stlRead2([stl_path,files(1).name]);
% load('femur_test')


%%
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


ea_vals=quantile(edge_angles,[0.5,0.01,0.99])
ar_vals=quantile(aspects,[0.5,0.01,0.99])