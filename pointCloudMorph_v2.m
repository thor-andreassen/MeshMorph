function [source_nodes_fit]= pointCloudMorph_v2(target_nodes,source_nodes,params)
    max_iterations=params.max_iterations;
    want_plot=params.want_plot;
    beta=params.beta;
    scale=params.scale;
    dist_threshold=params.dist_threshold;
    dist_threshold_scale=params.dist_threshold_scale;
    scale_scale=params.scale_scale;
    knots_scale=params.knots_scale;
    beta_scale=params.beta_scale;
    
    %% mesh morphin initialization
    source_nodes_0=source_nodes;

    %% main loop for mesh morphing
    counter=1;
    max_knots=500;
    switcher=1;
    while counter<=max_iterations
        
        if switcher==1
            % use target nodes as knot
            [ind,D] = knnsearch(source_nodes,target_nodes,'K',1);
            index_mat=[(1:size(target_nodes,1))',ind,D];
            index_mat=sortrows(index_mat,3,'descend');
            
            count_knot=1;
            K=target_nodes(index_mat(1,1),:);
            count_pot=2;
            index_list=index_mat(1,:);
            while count_knot<max_knots && count_pot<size(target_nodes,1)
                dist_vals=pdist2(K,target_nodes(index_mat(count_pot,1),:));
                if min(dist_vals)>dist_threshold
                    K=[K;target_nodes(index_mat(count_pot,1),:)];
                    index_list=[index_list;index_mat(count_pot,:)];
                    count_knot=count_knot+1;
                end
                count_pot=count_pot+1;
            end
            switcher=1;
            target_points=target_nodes(index_list(:,1),:);
            source_points=source_nodes(index_list(:,2),:);
        else
            [ind,D] = knnsearch(target_nodes,source_nodes,'K',1);
            index_mat=[(1:size(source_nodes,1))',ind,D];
            index_mat=sortrows(index_mat,3,'descend');
            
            count_knot=1;
            K=target_nodes(index_mat(1,2),:);
            count_pot=2;
            index_list=index_mat(1,:);
            while count_knot<max_knots && count_pot<size(source_nodes,1)
                dist_vals=pdist2(K,source_nodes(index_mat(count_pot,1),:));
                if min(dist_vals)>dist_threshold
                    K=[K;source_nodes(index_mat(count_pot,2),:)];
                    index_list=[index_list;index_mat(count_pot,:)];
                    count_knot=count_knot+1;
                end
                count_pot=count_pot+1;
            end
            switcher=1;
            
            target_points=target_nodes(index_list(:,2),:);
            source_points=source_nodes(index_list(:,1),:);
        end
        
        
        knot_points=target_points;
        dist_test=pdist2(source_points,knot_points);
        c=diag(dist_test).^3;
        K_test=(dist_test.^3+(c)').^beta;
        P=(target_points-source_points);
        w=K_test\P;


        source_dist=pdist2(source_nodes,knot_points);
        K_new=(source_dist.^3+(c)').^beta;
        source_nodes=source_nodes+scale*K_new*w;


        scale=scale*scale_scale;
        dist_threshold=dist_threshold*dist_threshold_scale;
        max_knots=round(max_knots*knots_scale)
        num_knots=size(K,1)
        beta=beta*beta_scale;
        
        if want_plot==1
            if counter==1
                fig_mesh=figure();
            else
                fig_mesh;
            end
            clf(fig_mesh);
            subplot(1,2,1);
            plot3(source_nodes_0(:,1),source_nodes_0(:,2),source_nodes_0(:,3),'go','MarkerSize',1);
            hold on
            plot3(source_nodes(:,1),source_nodes(:,2),source_nodes(:,3),'bo','MarkerSize',1);
            plot3(target_nodes(:,1),target_nodes(:,2),target_nodes(:,3),'ro','MarkerSize',1);
            plot3(K(:,1),K(:,2),K(:,3),'kx','LineWidth',5)
            axis equal
            view([0 0 1])
            pause(.01);
        end
        counter=counter+1;
        
        
    end
    source_nodes_fit=source_nodes;

end