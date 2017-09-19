%mean squared error distance function
function D = MSQdistance(patchA, patchB)
    distance = patchA - patchB;
    distance = distance(~isnan(distance));
    D = sum(distance.^2) / numel(distance);
end