function index_ranges=createSplitsOfTotal(total_nums,splits)
    indices=linspace(1,total_nums,splits+1);
    indices=round(indices);
    index_ranges=[indices(1:end-1)',indices(2:end)'-1];
    index_ranges(end,end)=index_ranges(end,end)+1;
end