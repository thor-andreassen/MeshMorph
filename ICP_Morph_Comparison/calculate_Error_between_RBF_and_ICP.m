%% clearing
clear
close all
clc

%% load RBF values
result_path='C:\Users\Thor.Andreassen\Desktop\Thor Personal Folder\Research\Iterative Alignment Check\MeshMorph\ICP_Morph_Comparison\';


load([result_path,'RBF_Tibia.mat']);
RBF_surf_distance=surf_distances;
RBF_haus_distance=haus_distance;

%% load RBF values
load([result_path,'ICP_Tibia.mat']);
ICP_surf_distance=surf_distances;
ICP_haus_distance=haus_distance;

%% load RBF values
load([result_path,'CPD_Tibia.mat']);
CPD_surf_distance=surf_distances;
CPD_haus_distance=haus_distance;


%% calculate error
RBF_surf_val=quantile(RBF_surf_distance,[.5,.01,.99])';
RBF_haus_val=quantile(RBF_haus_distance,[.5,.01,.99])';

ICP_surf_val=quantile(ICP_surf_distance,[.5,.01,.99])';
ICP_haus_val=quantile(ICP_haus_distance,[.5,.01,.99])';

CPD_surf_val=quantile(CPD_surf_distance,[.5,.01,.99])';
CPD_haus_val=quantile(CPD_haus_distance,[.5,.01,.99])';

vals=table(CPD_haus_val,RBF_haus_val,ICP_haus_val,CPD_surf_val,RBF_surf_val,ICP_surf_val)