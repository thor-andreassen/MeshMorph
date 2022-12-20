function tri_elems=splitQuadsToTries(quad_elems)
    
    tri_elems=zeros(size(quad_elems,1)*2,3);
    
    counter=1;
    for count_elem=1:size(quad_elems,1)
        tri_elems(counter,:)=quad_elems(count_elem,[1,2,3]);
        counter=counter+1;
        tri_elems(counter,:)=quad_elems(count_elem,[3,4,1]);
        counter=counter+1;
    end
    
end