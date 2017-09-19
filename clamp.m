% clamp value M to range specified by minval and maxval
function C = clamp(M, minval, maxval) 
    C = min(max(M, minval), maxval);
end