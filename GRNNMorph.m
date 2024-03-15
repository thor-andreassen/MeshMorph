function [source_nodes]= GRNNMorph(target_nodes,source_nodes,params,target_faces,source_faces)
%% main
% Created by Thor E. Andreassen, PhD
% University of Denver
% Last Edited 3/5/2024

% This code morphs one geometry to another. In particular, this code morphs a set
% of nodes of one geometry to those of the other. Additional speedup and mesh
% geometry can be taken into account by using a set of faces from a triangular
% mesh to determine relative mesh similairity based on normal directions.
% The code accepts as inputs meshes of 2D types with triangular faces.
% The function that calls this code, can perform the morphing on 3D
% element types or quadrilateral face types by reducing or converting
% to triangular, and then applying the morph to the original geometry.
% In the end the source geometry nodes will be moved to approximately
% line up with the resulting shape of the target geometry. If plotting is
% included, the initial source mesh at the start of these steps will be
% plotted in green. The target mesh will be plotted in red, and the morphed
% source nodes will be plotted in blue.
%These are structures that have the folloiwng variables:
    % target_nodes: (m x 3) A list of nodes as vertices of the
            % elements with the row corresponding to the assumed node number
            % corresponding to the list of nodes in the elements. The values
            % represent the location in cartesian x, y, z
            % coordinate space. The nodes given are those for the
            % stationary (unchanged) target geometry
    % source_nodes: (m x 3) A list of nodes as vertices of the
            % elements with the row corresponding to the assumed node number
            % corresponding to the list of nodes in the elements. The values
            % represent the location in cartesian x, y, z
            % coordinate space. The nodes given are those for the
            % unstationay (moving) source geometry

    % the params input is a structure of configuration controls to define
    % the method by which the overclosure adjustment occurs.
        % params.use_parallel (0 or 1) - Default = 0
            % A binary value used to determine if the algorithm
            % will use parallel computing for calculation of the
            % GRNN and the overclosure distances, or if only a
            % single core will be used. NOTE this functionality
            % requires the parallel computing toolbox.
            % 1 = use parallel computing
            % 0 = use single thread
        % params.max_iterations (x > 1) - Default = 20
            % The number of iterations to perform for this morphing
            % step. Values represent the number of times the
            % required displacement field will be calculated and
            % applied to the source nodes.
        % params.want_plot (0 or 1) - Default = 1
            % This parameter contorls whether graphs will be shown
            % to the user, showing the current node graphs at each
            % iteration. The green points represent the initial
            % location of the source points. The red points
            % represent the target points. The blue points
            % represent the current location of the source points
            % after each iteration.
        % params.smooth (X > 0) - Default = 10
            % This is the value that controls the amount of
            % smoothing applied to the GRNN to allow for smooth
            % contours of the moprphing. This parameter is the most
            % important one to test to create good results. A greater
            % number applies a larger smoothing and will usually
            % require more iterations to complete (Useful for morphing
            % gross anatomy). A smaller value will apply
            % less smoothing (useful for ensuring morphing of small
            % localized details). T Values between 1
            % and 100 work reasonably well for biomechanical
            % structures of the human body as meshes in mm.
            % Approximate values for the smoothing can be obtained
            % as the approximate size between the features that are
            % currently desired to be morphed.
        % params.scale (0 < X < 1) - Default = 0.5
            % This parameter represensts the scalar quantity
            % applied to the displacment field prior to application
            % to the orignal nodes. This helps smooth the motion
            % and acts as a learning rate for the algorithm.
        % params.normal_scale (0 <= X) - Default = 1
            % This parameter is used to help morph nodes to the correct
            % side of the geometries by including nnormal directions
            % of the points in the presence of surface mesh
            % information given (inclusion of faces above) as part
            % of the nearest neighbor distance calculation for each
            % mesh.
        % params.smooth_decay (0 <= x <= 1) - Default =
        % 0.995
            % This is a parameter that will adjust the smoothing
            % parameter each iteration so as to guarantee that the
            % algorithm is able to pickup smaller details. The value
            % will create a geometric series of diminishing
            % smoothing each iteration applied to the overall
            % smoothing factor.
        % params.normal_scale_decay (0 <= X <= 1) - Default =
        % 0.95
            % This is a parameter that will adjust the scale of the
            % normal values of the vertices, if the faces are
            % provided of the surfaces.So as to guarantee that the
            % algorithm is able to pickup smaller details. The value
            % will create a geometric series of diminishing
            % smoothing each iteration applied to the overall
            % normal scale factor.
        % params.memory_splits (1 <= X) - Default =1
            % this is a parameter that allows the computation to split up
            % the matrix for the GRNN calculation to allow the solution to
            % complete even in the presence of very many nodes, where
            % insufficient RAM is found. If as the code is run an error
            % similar to the following is found:
            %"Requested 100000x100000 (74.5GB) array exceeds maximum array
            % size preference (63.9GB). This might cause MATLAB to become
            % unresponsive."
            % Then this line allows the user to split the computation into
            % X number of splits, with a larger number splitting up the
            % computation that many times, at the cost of that many times
            % increase in computational time. I.E. params.memory_splits=10
            % means 10 times less RAM is required, but the computation will
            % take approximately 10 times as long. 

    % The following loop is used to determine if the face information has
    % been included and therefore the normal directions can be used for the
    % alignment.
    if nargin <=3
        use_normal=0;
    else
        use_normal=1;
    end
    
    % the following lines are used to set the defualt values for the
    % parameters, in the event they were not provided by the user.
    params=setDefaultParamValue(params,'use_parallel',0);
    params=setDefaultParamValue(params,'max_iterations',20);
    params=setDefaultParamValue(params,'want_plot',1);
    params=setDefaultParamValue(params,'smooth',10);
    params=setDefaultParamValue(params,'scale',0.5);
    params=setDefaultParamValue(params,'normal_scale',1);
    params=setDefaultParamValue(params,'smooth_decay',0.995);
    params=setDefaultParamValue(params,'normal_scale_decay',0.95);
    params=setDefaultParamValue(params,'memory_splits',1);
    
    
    
    % the following lines are used to set the vlaues for the internal
    % variables based on the information provided in the parameters and the
    % default values if not chosen by the user.
    use_parallel=params.use_parallel;
    max_iterations=params.max_iterations;
    want_plot=params.want_plot;
    smooth=params.smooth;
    scale=params.scale;
    normal_scale=params.normal_scale;
    smooth_decay=params.smooth_decay;
    normal_scale_decay=params.normal_scale_decay;
    memory_splits=params.memory_splits;
    
    %% mesh morphin initialization
    % the following lines initialize the counter, and set the values for
    % initial source nodes
    source_nodes_0=source_nodes;
    counter=1;
    
    % the following lines calculate the values for the target surface
    % normals, if the faces are provided, since these values do not change
    % for each iteration.
    if use_normal==1
        vertex_normal_target=target_nodes;
        for count_node_target=1:size(target_nodes)
            vertex_normal_target(count_node_target,:)=findVertexNormalFromMesh(target_faces,target_nodes,count_node_target);
        end
        vertex_normal_target=vertex_normal_target*normal_scale;
    end

    %% main loop for mesh morphing
    while counter<=max_iterations
    disp(counter)
    % the following lines calculate the values for the target surface
    % normals, if the faces are provided, since these values do not change
    % for each iteration.
        if use_normal==1
            vertex_normal_source=source_nodes;
            for count_node_source=1:size(source_nodes)
                vertex_normal_source(count_node_source,:)=findVertexNormalFromMesh(source_faces,source_nodes,count_node_source);
            end
            vertex_normal_source=vertex_normal_source*normal_scale;
        end
    
        % the following lines determine the nearest neighbor from the
        % target to the corresponding points in the source. The
        % distances are either calculated with the inclusion of the vertex
        % normals or not. 
        if use_normal==1
            [ind_target] = knnsearch([source_nodes,vertex_normal_source],[target_nodes,vertex_normal_target],'K',1);
        else
            [ind_target] = knnsearch(source_nodes,target_nodes,'K',1);
        end
        point_nodes_target=target_nodes;
        vec_source_to_target1=target_nodes-source_nodes(ind_target,:);


        % the following lines determine the nearest neighbor from the
        % source to the corresponding points in the target. The
        % distances are either calculated with the inclusion of the vertex
        % normals or not.
        if use_normal==1
            [ind_source] = knnsearch([target_nodes,vertex_normal_target],[source_nodes,vertex_normal_source],'K',1);
        else
            [ind_source] = knnsearch(target_nodes,source_nodes,'K',1);
        end
        point_nodes_source=source_nodes;
        vec_source_to_target2=target_nodes(ind_source,:)-source_nodes;


        % the following lines create the concatenated set of points and
        % displacement field vectors for the GRNN training. NOTE: all of
        % the points are the original cartesian coordinates. The vectors
        % are all created such that they are in a uniform coordinate system
        % (as if all directed form the source to the target)
        all_points=[point_nodes_target;point_nodes_source];
        all_vectors=[vec_source_to_target1;vec_source_to_target2];
    
    
        % the following line creates the Generalized Regression Neural
        % Network based on the input of the cartesian points in space, and
        % the resulting displacement vector field.
        model=newgrnn(all_points',all_vectors',smooth);
    
    
    
    
        % the following determines the displacement field vectors to apply
        % to the source nodes for the current iteration. This can be done
        % either in parallel or individually for the original points. The
        % index ranges are used to split the memory up so as to not full up
        % the ram of the computer. For example, the size of the GRNN will
        % be approximately:
        % [# of source_nodes + # of target_nodes,# of source_nodes]

        % as such, this matrix can become too large for the available
        % memory of the computer. To combat this, the program allows for
        % the matrix to be split into several smaller matrices of:
        % [# of source_nodes + # of target_nodes,(# of source_nodes)/splits] 

        % allowing for the computation to finish (even for small amounts of
        % RAM) at the cost of decreased computational speed. 
        index_ranges=createSplitsOfTotal(size(source_nodes,1),memory_splits);
        deform_vector=zeros(size(source_nodes'));
        for count_memory_split=1:memory_splits
            indices_current=index_ranges(count_memory_split,:);
            if use_parallel==1
                deform_vector(:,indices_current(1):indices_current(2))=sim(model,source_nodes(indices_current(1):indices_current(2),:)','useParallel','yes');
            else
                deform_vector(:,indices_current(1):indices_current(2))=sim(model,source_nodes(indices_current(1):indices_current(2),:)');
            end
        end
        
        % the following lines are to implement the momentum term
        if counter==1
            previous_vector=zeros(size(deform_vector,2),size(deform_vector,1));
            momentum=.5;
        end



        % the following line is used to apply the determined displacement
        % vector field to the current source nodes.
        source_nodes=source_nodes+scale*deform_vector'+momentum*previous_vector;
        previous_vector=scale*deform_vector'+momentum*previous_vector;
        momentum=momentum*0.9;
    
        % the following lines update the smoothing and the normal scale to
        % decay the values over time. This ensures that the morphing starts
        % with gross/large geometries, and eventually focuses on
        % smaller/fine components of the meshes with localized differences.
        smooth=smooth*smooth_decay;
        normal_scale=normal_scale*normal_scale_decay;
    
        % the following lines plot the state of the current morphing after
        % the current iteration.
        if want_plot==1
            fig_mesh=figure(212121);
            p1=plot3(source_nodes_0(:,1),source_nodes_0(:,2),source_nodes_0(:,3),'go','MarkerSize',1);
            hold on
            p2=plot3(source_nodes(:,1),source_nodes(:,2),source_nodes(:,3),'bo','MarkerSize',1);
            p3=plot3(target_nodes(:,1),target_nodes(:,2),target_nodes(:,3),'ro','MarkerSize',1);
            axis off
            axis equal
        end

        % the following line updates the counter with the current
        % iteration. The morphing stops when the current iteration exceeds
        % the maximum number allowed.
        counter=counter+1;
    end
end