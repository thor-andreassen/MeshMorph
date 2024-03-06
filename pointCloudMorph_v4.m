function [source_nodes_fit]= pointCloudMorph_v4(target_nodes,source_nodes,params,target_faces,source_faces)
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
    % line up with the resulting shape of the target geometry.
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
                % params.use_parallel_loops (0 or 1) - Default = 0
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
                    %

                % params.relative_gap_weight (0 <= x <= 1) - Default = 0.5
                        % 1 = fixed geom 1 moving geom 2.
                        % 0 = moving geom 1, fixed geom 2.
                        % 0.5 = equal deformation
                % params.element_3d_type ([1,0]) - Value is Required
                    % A binary array of two values containing flags for if
                    % the geometry being input is of 2D or 3D type for geom
                    % 1, and then geom 2, respectively.
                % params.smooth_2D_surface (0 or 1) - Default = 0
                    % A binary value that is only used if the element type
                    % is a 2D triangular face element. When enabled this
                    % allows smoothing of the mesh between iterations to
                    % improve the final quality of the mesh.
                    % 1 = enable smoothing
                    % 0 = no smoothing

                % params.smoothing_reduction (0 <= x <= 1) - Default =
                % 0.995
                    % This is a parameter that will adjust the smoothing
                    % parameter each iteration so as to guarantee that the
                    % algorithm is able to remove all overclosure as they
                    % become smaller and less frequent in number. The value
                    % will create a geometric series of diminishing
                    % smoothing each iteration.

                % params.stop_tolerance (x > 0) - Default = 1E-5
                    % The parameter that controls the point that the
                    % algorithm will stop when it is within is tolerance of
                    % the given desired gap.
                % params.geom1_mesh_reduction_factor (0 < x < 1) - Default
                % = 0.01
                    % This parameter controls the initial amount to reduce
                    % the nubmer of elements and nodes of the first
                    % geometry by to improve the speed of the algorithm. A
                    % value of 1 means no reduction, and a value of 0.01
                    % means to reduce the mesh to approximately 1% of its
                    % initial number of nodes and elements. Note, this has
                    % no effect on the final mesh which will have the same
                    % number of nodes and elements as the original
                    % geometries. 
                % params.geom2_mesh_reduction_factor (0 < x < 1) - Default
                % = 0.01
                    % This parameter controls the initial amount to reduce
                    % the nubmer of elements and nodes of the second
                    % geometry by to improve the speed of the algorithm. A
                    % value of 1 means no reduction, and a value of 0.01
                    % means to reduce the mesh to approximately 1% of its
                    % initial number of nodes and elements. Note, this has
                    % no effect on the final mesh which will have the same
                    % number of nodes and elements as the original
                    % geometries.     


    if nargin <=3
        use_normal=0;
    else
        use_normal=1;
    end

    params=setDefaultParamValue(params,'use_parallel_loops',0);
    params=setDefaultParamValue(params,'max_iterations',20);
    params=setDefaultParamValue(params,'want_plot',1);
    params=setDefaultParamValue(params,'smooth',10);
    params=setDefaultParamValue(params,'scale',0.5);

    use_parallel=params.use_parallel;
    max_iterations=params.max_iterations;
    want_plot=params.want_plot;
    smooth=params.smooth;
    
    scale=params.scale;
    
    normal_scale=params.normal_scale;
    smooth_decay=params.smooth_decay;
    normal_scale_decay=params.normal_scale_decay;
    
    
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
                fig_mesh=figure(212121);
                p1=plot3(source_nodes_0(:,1),source_nodes_0(:,2),source_nodes_0(:,3),'go','MarkerSize',1);
                hold on
                p2=plot3(source_nodes(:,1),source_nodes(:,2),source_nodes(:,3),'bo','MarkerSize',1);
                p3=plot3(target_nodes(:,1),target_nodes(:,2),target_nodes(:,3),'ro','MarkerSize',1);
                axis off
                axis equal
                

        end
        counter=counter+1;
        
        
    end
    source_nodes_fit=source_nodes;

end