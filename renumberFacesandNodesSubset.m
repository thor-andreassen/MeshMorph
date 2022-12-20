function [faces_renumber,nodes_renumber,node_correspondance_list]=renumberFacesandNodesSubset(face_subset,nodes)
    node_correspondance_list=zeros(size(nodes,1),2);
    elem_length=size(face_subset,2);
    faces_renumber=face_subset;
    counter_new_node=1;
    nodes_renumber=nodes;
    for count_elem=1:size(face_subset,1)
        for count_node=1:elem_length
            original_node_number=face_subset(count_elem,count_node);
            if node_correspondance_list(original_node_number,1)==0
                nodes_renumber(counter_new_node,:)=nodes(original_node_number,:);
                node_correspondance_list(original_node_number,2)=counter_new_node;
                node_correspondance_list(original_node_number,1)=original_node_number;
                counter_new_node=counter_new_node+1;
            end
            
            faces_renumber(count_elem,count_node)=node_correspondance_list(original_node_number,2);
                
        end
        
    end
    
    % eliminate unused nodes
    for count_nodes=size(nodes,1):-1:1
        if node_correspondance_list(count_nodes,1)==0
            node_correspondance_list(count_nodes,:)=[];
        end
    end
    
    nodes_renumber=nodes_renumber(1:(counter_new_node-1),:);
    
end