function rbf_mat=rbfFuncMat(rbf_func,d,K,sigmas)
    rbf_mat=zeros(size(d,1),size(K,1));
    for count_d=1:size(d,1)
        for count_K=1:size(K,1)
            dist_val=norm(d(count_d,:)-K(count_K,:));
            rbf_mat(count_d,count_K)=rbf_func(dist_val,sigmas(count_K));
        end
    end
end