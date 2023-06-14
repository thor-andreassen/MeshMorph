function [source_nodes_fit]= pointCloudMorph_v3(target_nodes,source_nodes,params,target_faces,source_faces)
    max_iterations=params.max_iterations;
    want_plot=params.want_plot;
    beta=params.beta;
    scale=params.scale;
    dist_threshold=params.dist_threshold;
    dist_threshold_scale=params.dist_threshold_scale;
    scale_scale=params.scale_scale;
    knots_scale=params.knots_scale;
    beta_scale=params.beta_scale;
    smooth=params.smooth;
    normal_scale=params.normal_scale;
    smooth_decay=params.smooth_decay;
    
    %% mesh morphin initialization
    source_nodes_0=source_nodes;

    %% main loop for mesh morphing
    counter=1;
    max_knots=500;
    switcher=1;
%     target_scale
vertex_normal_target=target_nodes;
for count_node_target=1:size(target_nodes)
    vertex_normal_target(count_node_target,:)=findVertexNormalFromMesh(target_faces,target_nodes,count_node_target);
end
vertex_normal_target=vertex_normal_target*normal_scale;


    while counter<=max_iterations

        vertex_normal_source=source_nodes;
        for count_node_source=1:size(source_nodes)
            vertex_normal_source(count_node_source,:)=findVertexNormalFromMesh(source_faces,source_nodes,count_node_source);
        end
        vertex_normal_source=vertex_normal_source*normal_scale;

        counter
        [ind_target] = knnsearch([source_nodes,vertex_normal_source],[target_nodes,vertex_normal_target],'K',1);

%         [ind_target] = knnsearch(source_nodes,target_nodes,'K',1);
        point_nodes_target=target_nodes;
        vec_source_to_target1=target_nodes-source_nodes(ind_target,:);
        
        [ind_source] = knnsearch([target_nodes,vertex_normal_target],[source_nodes,vertex_normal_source],'K',1);

%         [ind_source] = knnsearch(target_nodes,source_nodes,'K',1);
        point_nodes_source=source_nodes;
        vec_source_to_target2=target_nodes(ind_source,:)-source_nodes;
        
        all_points=[point_nodes_target;point_nodes_source];
        all_vectors=[vec_source_to_target1;vec_source_to_target2];
        
        
        model=newgrnn(all_points',all_vectors',smooth);
        
        deform_vector=sim(model,source_nodes');
        source_nodes=source_nodes+scale*deform_vector';
        
        smooth=smooth*smooth_decay;
        
        if want_plot==1
            if counter==1
                fig_mesh=figure();
                p1=plot3(source_nodes_0(:,1),source_nodes_0(:,2),source_nodes_0(:,3),'go','MarkerSize',1);
                hold on
                p2=plot3(source_nodes(:,1),source_nodes(:,2),source_nodes(:,3),'bo','MarkerSize',1);
                p3=plot3(target_nodes(:,1),target_nodes(:,2),target_nodes(:,3),'ro','MarkerSize',1);
            else
                p1.XData=source_nodes_0(:,1);
                p1.YData=source_nodes_0(:,2);
                p1.ZData=source_nodes_0(:,3);
                p2.XData=source_nodes(:,1);
                p2.YData=source_nodes(:,2);
                p2.ZData=source_nodes(:,3);
                p3.XData=target_nodes(:,1);
                p3.YData=target_nodes(:,2);
                p3.ZData=target_nodes(:,3);
                pause(.01);                                                                                                                          
            end

        end
        counter=counter+1;
        
        
    end
    source_nodes_fit=source_nodes;

end