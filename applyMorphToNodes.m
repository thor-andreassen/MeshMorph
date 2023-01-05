function deform_nodes=applyMorphToNodes(nodes_orig,Affine_TransMat,model_final,use_parallel)
    if nargin<=3
        use_parallel=0;
    end
    temp_nodes=[nodes_orig,ones(size(nodes_orig,1),1)];
    affine_nodes=[Affine_TransMat*temp_nodes']';
    deform_nodes=affine_nodes(:,1:3);
    if use_parallel==1
        new_deform=sim(model_final,deform_nodes','useParallel','yes');
    else
        new_deform=sim(model_final,deform_nodes');
    end
    deform_nodes=deform_nodes+new_deform';
end