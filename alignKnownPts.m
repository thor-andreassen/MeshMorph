function source_to_target_TransMat=alignKnownPts(target,source)
%% Calculate Center of Mass of coordinates
    source_com=mean(source,1);
    source_atzero=bsxfun(@minus,source,source_com);

    target_com=mean(target,1);
    target_atzero=bsxfun(@minus,target,target_com);

    %% Calculate Rotation from Laser Scan to Optotrak
    H=source_atzero'*(target_atzero);
    [U,~,V]=svd(H);
    rotation_mat=(V*U');

    if det(rotation_mat)<0
        [U,~,V]=svd(rotation_mat);
        V(:,3)=V(:,3)*-1;
        rotation_mat=V*U';
        disp('det less than zero');
    end

    %% Calculate Translation vector
    trans_vector=-rotation_mat*source_com'+target_com';

    %% Create Transformation Matrix
    source_to_target_TransMat=zeros(4,4);
    source_to_target_TransMat(1,1)=1;
    source_to_target_TransMat(2:4,1)=trans_vector;
    source_to_target_TransMat(2:4,2:4)=rotation_mat;

end