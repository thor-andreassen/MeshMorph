function saveAbaqusInputFile(filename,nodes,elements,elem_category,dim,nset_name,elset_name)

% elem_category=['R','C'];
% dim=[2,3];
if dim==2
    dim_type='2D';
else
    dim_type='3D';
end

num_nodes=size(elements,2)-1;
if testCharPresentInChar(elem_category,'C') && dim==2
    elem_TYPE=['CPS',num2str(num_nodes)];
elseif testCharPresentInChar(elem_category,'C')
    elem_TYPE=[elem_category,dim_type,num2str(num_nodes),'R'];
else
    elem_TYPE=[elem_category,dim_type,num2str(num_nodes)];
end


% elem_TYPE='C3D8R';
% nset_name='FEMUR_NODES';
% elset_name='FEMUR_ELEMS';
fid=fopen(filename,'w+');
fprintf(fid,'**\r\n');
fprintf(fid,'**\r\n');
fprintf(fid,'*NODE, NSET=%s\r\n',nset_name);
for count_node=1:size(nodes,1)
    fprintf(fid,'%d, %12f, %12f, %12f\r\n',...
        nodes(count_node,1),...
        nodes(count_node,2),...
        nodes(count_node,3),...
        nodes(count_node,4));
end
fprintf(fid,'**\r\n');
fprintf(fid,'**\r\n');
fprintf(fid,'*ELEMENT, ELSET=%s, TYPE=%s\r\n',...
    elset_name,elem_TYPE);
for count_row=1:size(elements,1)
    for count_col=1:size(elements,2)
        if count_col<size(elements,2)
            fprintf(fid,'%d,   ',elements(count_row,count_col));
        else
            fprintf(fid,'%d\r\n',elements(count_row,count_col));
        end
    end
end
fclose(fid);
% patch('Faces',elements(:,2:4)-100,'Vertices',nodes(:,2:end),...
%     'FaceColor','r','EdgeAlpha',.3)