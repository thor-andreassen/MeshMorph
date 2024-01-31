function pts = transformPts(TransMat,pts_orig)

    pts_temp=[pts_orig,ones(size(pts_orig,1),1)]';
    pts_temp=TransMat*pts_temp;
    pts=pts_temp(1:3,:)';
end