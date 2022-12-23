function [source_nodes_fit]= pointCloudMorph(target_nodes,source_nodes,params)
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
    
    %% mesh morphin initialization
    source_nodes_0=source_nodes;
    bb_min=min(source_nodes_0);
    bb_max=max(source_nodes_0);
    bb_dif=abs(bb_max-bb_min)*.25;
    bb_min=bb_min-bb_dif;
    bb_max=bb_max+bb_dif;

    x_bb=linspace(bb_min(1),bb_max(1),initial_knots);
    y_bb=linspace(bb_min(2),bb_max(2),initial_knots);
    z_bb=linspace(bb_min(3),bb_max(3),initial_knots);

    [X,Y,Z]=meshgrid(x_bb,y_bb,z_bb);
    X=reshape(X,[],1);
    Y=reshape(Y,[],1);
    Z=reshape(Z,[],1);

    bb_dif=bb_max-bb_min;
    K0=[X,Y,Z];

    K=K0;
    phi_func=@(r,sigma) exp(-(r.^2)/(2*sigma.^2));

    %% main loop for mesh morphing
    counter=1;
    error_hist=[];
    switcher=start_target;
    while counter<=max_iterations
        % get sigmas
        [~,D] = knnsearch(K,K,'K',4);
        sigmas=mean(D(:,2:end),2)*f;

        % determine nearest points
        [Idx_T] = knnsearch(source_nodes,target_nodes,'K',1);
        U=target_nodes;
        d=source_nodes(Idx_T,:);
        rbf_mat_test=rbfFuncMat(phi_func,d,K,sigmas);

        % fit rbf A coefficients
        A=rbf_mat_test\U;

        % determine new node locations
        rbf_mat_full=rbfFuncMat(phi_func,source_nodes,K,sigmas);
        new_source_nodes=rbf_mat_full*A;

        if want_plot==1
            if counter==1
                fig_mesh=figure();
            else
                fig_mesh;
            end
            clf(fig_mesh);
            subplot(1,2,1);
            scatter3(source_nodes_0(:,1),source_nodes_0(:,2),source_nodes_0(:,3),'go','LineWidth',.25);
            hold on
            scatter3(source_nodes(:,1),source_nodes(:,2),source_nodes(:,3),'bo','LineWidth',.25);
            scatter3(target_nodes(:,1),target_nodes(:,2),target_nodes(:,3),'ro','LineWidth',.25);
            scatter3(K(:,1),K(:,2),K(:,3),'kx','LineWidth',5)
            axis equal
            view([0 0 1])
            pause(.01);
        end


        %determine new knot lcoations
        %     if mod(counter,2)==0
        if switcher==0
            [Idx_n,D_n] = knnsearch(source_nodes,target_nodes,'K',1);
            temp_dist=[Idx_n,D_n];
            temp_dist=sortrows(temp_dist,2,'descend');
            new_indices=unique(temp_dist(:,1),'stable');
            
            pot_knots=source_nodes(new_indices,:);
            count_p=1;
            count_pot=1;
            while count_p<p && count_pot<size(pot_knots,1)
                temp_knots=[K;pot_knots(count_pot,:)];
                [~,new_potential_knots_D]=...
                    knnsearch(temp_knots,pot_knots(count_pot,:),'K',2);
                if new_potential_knots_D(2) >= d_min
                    K=[K;pot_knots(count_pot,:)];
                    count_p=count_p+1;
                end
                count_pot=count_pot+1;
            end

            
%             new_knot_indices=new_indices(1:p);
%             K=[K;source_nodes(new_knot_indices,:)+rand(p,3)*rand_mult];

        else
            [Idx_n,D_n] = knnsearch(target_nodes,source_nodes,'K',1);
            temp_dist=[Idx_n,D_n];
            temp_dist=sortrows(temp_dist,2,'descend');
            new_indices=unique(temp_dist(:,1),'stable');
            
            pot_knots=target_nodes(new_indices,:);
            count_p=1;
            count_pot=1;
            while count_p<p && count_pot<size(pot_knots,1)
                temp_knots=[K;pot_knots(count_pot,:)];
                [~,new_potential_knots_D]=...
                    knnsearch(temp_knots,pot_knots(count_pot,:),'K',2);
                if new_potential_knots_D(2) >= d_min
                    K=[K;pot_knots(count_pot,:)];
                    count_p=count_p+1;
                end
                count_pot=count_pot+1;
            end
            
            
            
%             new_knot_indices=new_indices(1:p);
%             K=[K;target_nodes(new_knot_indices,:)+rand(p,3)*rand_mult];
        end

        if include_rand_knots==1
        
            rand_nodes_x=rand(p,1)*bb_dif(1)+bb_min(1);
            rand_nodes_y=rand(p,1)*bb_dif(2)+bb_min(2);
            rand_nodes_z=rand(p,1)*bb_dif(3)+bb_min(3);
            rand_nodes=[rand_nodes_x,rand_nodes_y,rand_nodes_z];
            K=[K;rand_nodes];
        end

        
        error=max(D_n);
        if want_plot==1
            subplot(1,2,2);
            error_hist=[error_hist,error];
            plot(error_hist);
            pause(.05);
        end

        
        source_nodes=new_source_nodes;
        counter=counter+1;

        rand_mult=rand_mult*rand_decay;
        f=f_decay*f;
        
        if mod(counter,knot_reset_iter)==0
            K=K0;
        end

        %     if mod(counter,2)==0
        if mod(counter,target_source_switch_iter)==0
            if switcher==1
                switcher=0;
            else
                switcher=1;
            end
        end

    end
    source_nodes_fit=source_nodes;

end