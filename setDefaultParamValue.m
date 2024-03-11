function params=setDefaultParamValue(params,name,value)
    if isfield(params, name)
        NaN;
    else
        params.(name)=value;
    end
end