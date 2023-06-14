%% clearing
clear
close all
clc

%% load RBF values
result_path='C:\Users\Thor.Andreassen\Desktop\Thor Personal Folder\Research\Iterative Alignment Check\MeshMorph\ICP_Morph_Comparison\';


load([result_path,'RBF_Femur.mat']);
RBF_surf_distance=surf_distances;
RBF_haus_distance=haus_distance;
RBF_aspect=aspects;
RBF_edge_angle=edge_angles;
RBF_node_dist=node_dist_travel;

%% load RBF values
load([result_path,'ICP_Femur.mat']);
ICP_surf_distance=surf_distances;
ICP_haus_distance=haus_distance;
ICP_aspect=aspects;
ICP_edge_angle=edge_angles;
ICP_node_dist=node_dist_travel;

%% load RBF values
load([result_path,'CPD_Femur.mat']);
CPD_surf_distance=surf_distances;
CPD_haus_distance=haus_distance;
CPD_aspect=aspects;
CPD_edge_angle=edge_angles;
CPD_node_dist=node_dist_travel;


%% calculate error
CPD_surf_val=quantile(CPD_surf_distance,[.5,.01,.99])';
CPD_haus_val=quantile(CPD_haus_distance,[.5,.01,.99])';
CPD_aspect_val=quantile(CPD_aspect,[.5,.01,.99])';
CPD_edge_val=quantile(CPD_edge_angle,[.5,.01,.99])';
CPD_node_val=quantile(CPD_node_dist,[.5,.01,.99])';

ICP_surf_val=quantile(ICP_surf_distance,[.5,.01,.99])';
ICP_haus_val=quantile(ICP_haus_distance,[.5,.01,.99])';
ICP_aspect_val=quantile(ICP_aspect,[.5,.01,.99])';
ICP_edge_val=quantile(ICP_edge_angle,[.5,.01,.99])';
ICP_node_val=quantile(ICP_node_dist,[.5,.01,.99])';

RBF_surf_val=quantile(RBF_surf_distance,[.5,.01,.99])';
RBF_haus_val=quantile(RBF_haus_distance,[.5,.01,.99])';
RBF_aspect_val=quantile(RBF_aspect,[.5,.01,.99])';
RBF_edge_val=quantile(RBF_edge_angle,[.5,.01,.99])';
RBF_node_val=quantile(RBF_node_dist,[.5,.01,.99])';

vals=table(CPD_surf_val,CPD_haus_val,CPD_aspect_val,CPD_edge_val,CPD_node_val,...
    ICP_surf_val,ICP_haus_val,ICP_aspect_val,ICP_edge_val,ICP_node_val,...
    RBF_surf_val,RBF_haus_val,RBF_aspect_val,RBF_edge_val,RBF_node_val)