function deform_nodes=applyMorphToNodes(nodes_orig,Affine_TransMat,model_final)
    temp_nodes=[nodes_orig,ones(size(nodes_orig,1),1)];
    affine_nodes=[Affine_TransMat*temp_nodes']';
    deform_nodes=affine_nodes(:,1:3);
    new_deform=sim(model_final,deform_nodes');
    deform_nodes=deform_nodes+new_deform';
end