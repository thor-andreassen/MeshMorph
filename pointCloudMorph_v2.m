function [source_nodes_fit]= pointCloudMorph_v2(target_nodes,source_nodes,params)
    %% clearing
%     clear
%     close all
%     clc
%     
%     %% load data
%     load('temp_test.mat');
    
    %%data
    f=params.f;
    f_decay=params.f_decay;
    max_iterations=params.max_iterations;
    want_plot=params.want_plot;
    rand_mult=params.rand_mult;
    rand_decay=params.rand_decay;
    include_rand_knots=params.include_rand_knots;
    knot_reset_iter=params.knot_reset_iter;
    target_source_switch_iter=params.target_source_switch_iter;
    start_target=params.start_target;
    p=params.new_knots_per_iter;
    initial_knots=params.initial_knots;
    d_min=params.d_min;
    beta=params.beta;
    scale=params.scale;
    dist_threshold=params.dist_threshold;
    
    
    %% mesh morphin initialization
    source_nodes_0=source_nodes;

    %% main loop for mesh morphing
    counter=1;
    error_hist=[];
    switcher=start_target;
    max_count=500;
    switcher=1;
    

    
    
    while counter<=max_iterations
        
        
        % get sigmas
        if switcher==1
            [ind,D] = knnsearch(source_nodes,target_nodes,'K',1);
            index_mat=[(1:size(target_nodes,1))',ind,D];
            
            
            index_mat=sortrows(index_mat,3,'descend');
            
            count_knot=1;
            
            K=target_nodes(index_mat(1,1),:);
            count_pot=2;
            index_list=index_mat(1,:);
            while count_knot<max_count && count_pot<size(target_nodes,1)
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
%             knot_points=source_points;
        else
            [ind,D] = knnsearch(target_nodes,source_nodes,'K',1);
            index_mat=[(1:size(source_nodes,1))',ind,D];
            
            
            index_mat=sortrows(index_mat,3,'descend');
            
            count_knot=1;
            
            K=target_nodes(index_mat(1,2),:);
            count_pot=2;
            index_list=index_mat(1,:);
            while count_knot<max_count && count_pot<size(source_nodes,1)
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
%             knot_points=target_points;
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


        scale=scale*.9;
        dist_threshold=dist_threshold*.9
        max_count=round(max_count*1.05)
        num_knots=size(K,1)
        beta=beta*.95
        
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