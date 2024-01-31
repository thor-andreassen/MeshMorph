function [source_nodes_fit]= pointCloudMorph_v4(target_nodes,source_nodes,params,target_faces,source_faces)
    if nargin <=3
        use_normal=0;
    else
        use_normal=1;
    end
    max_iterations=params.max_iterations;
    want_plot=params.want_plot;
    scale=params.scale;
    smooth=params.smooth;
    normal_scale=params.normal_scale;
    smooth_decay=params.smooth_decay;
    normal_scale_decay=params.normal_scale_decay;
    use_parallel=params.use_parallel;
    
    %% temp figure structure
%        target_nodes=target_nodes(target_nodes(:,3)>16,:);
%        target_nodes(:,3)=target_nodes(:,3)-3.75;
%        source_nodes=source_nodes(source_nodes(:,3)>12,:);
%        source_nodes=source_nodes(source_nodes(:,1)>-1.5,:);
%        source_nodes=source_nodes(source_nodes(:,2)>-.25 &source_nodes(:,2)<.25,:);
% 
%        target_nodes=target_nodes(target_nodes(:,2)>-.25 &target_nodes(:,2)<.25,:);
%        source_nodes=source_nodes(source_nodes(:,2)>0,:);
%    target_nodes(:,2)=0;
%    source_nodes(:,2)=0;
    
    %% mesh morphin initialization

    source_nodes_0=source_nodes;

    %% main loop for mesh morphing
    counter=1;
    
if use_normal==1
    vertex_normal_target=target_nodes;
    for count_node_target=1:size(target_nodes)
        vertex_normal_target(count_node_target,:)=findVertexNormalFromMesh(target_faces,target_nodes,count_node_target);
    end
    vertex_normal_target=vertex_normal_target*normal_scale;
end

    while counter<=max_iterations
        source_nodes_0=source_nodes;
        if use_normal==1
            vertex_normal_source=source_nodes;
            for count_node_source=1:size(source_nodes)
                vertex_normal_source(count_node_source,:)=findVertexNormalFromMesh(source_faces,source_nodes,count_node_source);
            end
            vertex_normal_source=vertex_normal_source*normal_scale;
        end

        counter
        if use_normal==1
            [ind_target] = knnsearch([source_nodes,vertex_normal_source],[target_nodes,vertex_normal_target],'K',1);
        else
            [ind_target] = knnsearch(source_nodes,target_nodes,'K',1);
        end
        
        point_nodes_target=target_nodes;
        vec_source_to_target1=target_nodes-source_nodes(ind_target,:);
        
        if use_normal==1
            [ind_source] = knnsearch([target_nodes,vertex_normal_target],[source_nodes,vertex_normal_source],'K',1);
        else
            [ind_source] = knnsearch(target_nodes,source_nodes,'K',1);
        end
        
        point_nodes_source=source_nodes;
        vec_source_to_target2=target_nodes(ind_source,:)-source_nodes;
        
        all_points=[point_nodes_target;point_nodes_source];
        all_vectors=[vec_source_to_target1;vec_source_to_target2];
        
        
        
        
        % test memory efficiency
        model=newgrnn(all_points',all_vectors',smooth);
        
        
        
        
        
        if use_parallel==1
            deform_vector=sim(model,source_nodes','useParallel','yes');
        else
            deform_vector=sim(model,source_nodes');
        end
        
        source_nodes=source_nodes+scale*deform_vector';
        
        smooth=smooth*smooth_decay;
        normal_scale=normal_scale*normal_scale_decay;
        
        if want_plot==1
            if 1==1
                
                vec_length=0.2;
                
                fig_mesh=figure(212121);
                p1=plot3(source_nodes_0(:,1),source_nodes_0(:,2),source_nodes_0(:,3),'go','MarkerSize',1);
                hold on
                p2=plot3(source_nodes(:,1),source_nodes(:,2),source_nodes(:,3),'bo','MarkerSize',1);
                p3=plot3(target_nodes(:,1),target_nodes(:,2),target_nodes(:,3),'ro','MarkerSize',1);
                axis off
                axis equal
                
                
%                 fig_mesh2=figure('units','normalized','outerposition',[0 0 1 1]);
% % %                 p1=plot3(source_nodes_0(:,1),source_nodes_0(:,2),source_nodes_0(:,3),'bo','MarkerSize',7.5);
% % %                 hold on
% % % %                 p2=plot3(source_nodes(:,1),source_nodes(:,2),source_nodes(:,3),'bo','MarkerSize',5);
% % %                 p3=plot3(target_nodes(:,1),target_nodes(:,2),target_nodes(:,3),'ro','MarkerSize',7.5);
% % %                 q1=quiver3(source_nodes_0(:,1),source_nodes_0(:,2),source_nodes_0(:,3),...
% % %                     vec_source_to_target2(:,1),vec_source_to_target2(:,2),vec_source_to_target2(:,3),vec_length,'k');
% % % %                 q2=quiver3(target_nodes(:,1),target_nodes(:,2),target_nodes(:,3),...
% % % %                     -vec_source_to_target1(:,1),-vec_source_to_target1(:,2),-vec_source_to_target1(:,3),vec_length,'k');
% % % %                 q3=quiver3(target_nodes(:,1),target_nodes(:,2),target_nodes(:,3),...
% % % %                     vec_source_to_target1(:,1),vec_source_to_target1(:,2),vec_source_to_target1(:,3),vec_length,'k');
% % % %                 temp_vec=deform_vector';
% % % %                 q4=quiver3(source_nodes_0(:,1),source_nodes_0(:,2),source_nodes_0(:,3),...
% % % %                     temp_vec(:,1),temp_vec(:,2),temp_vec(:,3),vec_length,'k');


% figure()
% p1=plot(source_nodes_0(:,1),source_nodes_0(:,3),'bo','MarkerSize',10,'MarkerFaceColor','#0000FF','LineWidth', 4);
% % p2=plot3(source_nodes(:,1),source_nodes(:,2),source_nodes(:,3),'bo','MarkerSize',25,'MarkerFaceColor','#0000FF');
% hold on
% p3=plot(target_nodes(:,1),target_nodes(:,3),'ro','MarkerSize',10,'MarkerFaceColor','#FF0000','LineWidth', 4);
% 
% 
%                 q1=quiver(source_nodes_0(:,1),source_nodes_0(:,3),...
%                     vec_source_to_target2(:,1),vec_source_to_target2(:,3),0,'LineWidth',2.5,'Color','b',...
%                     'MaxHeadSize',0.1,'ShowArrowHead','on');
%                 q2=quiver(target_nodes(:,1),target_nodes(:,3),...
%                     -vec_source_to_target1(:,1),-vec_source_to_target1(:,3),0,'LineWidth',2.5,'Color','r',...
%                     'MaxHeadSize',0.1,'ShowArrowHead','on');
%                 q3=quiver(target_nodes(:,1),target_nodes(:,3),...
%                     vec_source_to_target1(:,1),vec_source_to_target1(:,3),0,'LineWidth',2.5,'Color','r',...
%                     'MaxHeadSize',0.1,'ShowArrowHead','on');
%                 temp_vec=deform_vector';
%                 q4=quiver(source_nodes_0(:,1),source_nodes_0(:,3),...
%                     temp_vec(:,1),temp_vec(:,3),0,'LineWidth',2.5,'Color','k',...
%                     'MaxHeadSize',0.1,'ShowArrowHead','on');
% 
%                 
%                 [X,Y,Z]=meshgrid(linspace(-2,2,10),[0],linspace(12,17,10));
%                 X=reshape(X,[],1);
%                 Y=reshape(Y,[],1);
%                 Z=reshape(Z,[],1);
%                 test_pts=[X,Y,Z];
%                 deform_vector_field=sim(model,test_pts','useParallel','yes');
%                 temp_vec=deform_vector_field';
%                 q4=quiver(test_pts(:,1),test_pts(:,3),...
%                     temp_vec(:,1),temp_vec(:,3),'LineWidth',2.5,'Color','k',...
%                     'MaxHeadSize',0.1,'ShowArrowHead','on');
%                 
%                 
% 
% axis off
% 
% axis([-2,2,12,17]);
% % axis equal
% % view([0,1,0]);
% %                 axis equal
% disp('test')
% % camzoom(1.5)
%   disp('test');              
                
% % % %                 deform_mesh=figure();
% % % %                 p4=plot3(source_nodes(:,1),source_nodes(:,2),source_nodes(:,3),'bo','MarkerSize',1);
% % % %                 hold on
% % % %                 p5=plot3(target_nodes(:,1),target_nodes(:,2),target_nodes(:,3),'ro','MarkerSize',1);
% % % %                 q1=quiver3(source_nodes(:,1),source_nodes(:,2),source_nodes(:,3),...
% % % %                     vec_source_to_target2(:,1),vec_source_to_target2(:,2),vec_source_to_target2(:,3),vec_length,'k');
% % % %                 axis off
% % % %                 axis equal
% % % %                 view([0,1,0])
% % % %                 disp('test')
% % % %                 
% % % %                 
% % % %                 deform_mesh2=figure();
% % % %                 p6=plot3(source_nodes(:,1),source_nodes(:,2),source_nodes(:,3),'bo','MarkerSize',1);
% % % %                 hold on
% % % %                 p7=plot3(target_nodes(:,1),target_nodes(:,2),target_nodes(:,3),'ro','MarkerSize',1);
% % % %                 q2=quiver3(target_nodes(:,1),target_nodes(:,2),target_nodes(:,3),...
% % % %                     vec_source_to_target1(:,1),vec_source_to_target1(:,2),vec_source_to_target1(:,3),vec_length,'k');
% % % %                 axis off
% % % %                 axis equal
% % % %                 view([0,1,0])
% % % %                 disp('test')
% % % % 
% % % %                 temp_deform=deform_vector';
% % % %                 vector_mesh=figure();
% % % %                 p6=plot3(source_nodes(:,1),source_nodes(:,2),source_nodes(:,3),'bo','MarkerSize',1);
% % % %                 hold on
% % % %                 q2=quiver3(source_nodes(:,1),source_nodes(:,2),source_nodes(:,3),...
% % % %                     temp_deform(:,1),temp_deform(:,2),temp_deform(:,3),vec_length,'k');
% % % %                 axis off
% % % %                 axis equal
% % % %                 view([0,1,0])
% % % %                 disp('test')
% % % %                 
% % % %                 temp_deform=vecnorm(deform_vector',2,2);
% % % %                 patch_mesh=figure();
% % % %                 p6=patch('Faces',source_faces,'Vertices',source_nodes,'FaceVertexCData',temp_deform,'FaceColor','interp','EdgeAlpha',.3);
% % % %                 axis off
% % % %                 axis equal
% % % %                 view([0,1,0])
% % % %                 c=jet(1000);
% % % %                 colormap(c(125:875,:));
% % % %                 colorbar
% % % %                 caxis([0,max(temp_deform)]);
% % % %                 axis equal
%                 disp('test')
                

                
                
                
            else
                fig_mesh;
                p1.XData=source_nodes_0(:,1);
                p1.YData=source_nodes_0(:,2);
                p1.ZData=source_nodes_0(:,3);
                p2.XData=source_nodes(:,1);
                p2.YData=source_nodes(:,2);
                p2.ZData=source_nodes(:,3);
                p3.XData=target_nodes(:,1);
                p3.YData=target_nodes(:,2);
                p3.ZData=target_nodes(:,3);
                pause(.001);
            end

        end
        counter=counter+1;
        
        
    end
    source_nodes_fit=source_nodes;

end