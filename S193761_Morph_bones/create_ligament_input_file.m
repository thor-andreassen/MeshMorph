%% clearing
clear
close all
clc

%% load data
path='C:\Users\Thor.Andreassen\Desktop\Thor Personal Folder\Research\Iterative Alignment Check\MeshMorph\S193761_Morph_bones\';
filename='predicted_ligament_sites.xlsx';
path_name=[path,filename];
vals=readtable(path_name);

%% create input file
fid=fopen('LIGAMENT_SITES.inp','w+');
nodes=vals{:,2:end};

fprintf(fid,'**\r\n');
fprintf(fid,'**\r\n');
fprintf(fid,'**\r\n');
fprintf(fid,'*NODE,NSET=LIGAMENT_NODES\r\n');
counter=1;
for count_lig=1:size(nodes,1)
    fprintf(fid,'%d, %6f, %6f, %6f\r\n',counter,nodes(count_lig,1:3));
    counter=counter+1;
    fprintf(fid,'%d, %6f, %6f, %6f\r\n',counter,nodes(count_lig,4:end));
    counter=counter+1;
end
fprintf(fid,'**\r\n');
fprintf(fid,'**\r\n');
counter=1;
for count_lig=1:size(nodes,1)
    fprintf(fid,'*ELEMENT,TYPE=T2D2,ELSET=%s\r\n',char(vals{count_lig,1}));
    fprintf(fid,'%d, %d, %d\r\n',count_lig,counter,counter+1);
    fprintf(fid,'**\r\n');
    fprintf(fid,'**\r\n');
    counter=counter+2;
end
fprintf(fid,'**\r\n');
fprintf(fid,'**\r\n');



fclose(fid)
