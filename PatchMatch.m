% 
% PatchMatch.m
%
% the Matlab code of PatchMatch algorithm
% PatchMatch returns approximate nearest neighbor field (NNF).

% Usage: [NNF, intermediate] = PatchMatch(imgTarget, imgSource, patchSize, maxIterations)
% 
% Inputs: 
% - imgTarget: An image (usually masked by NaN. NaN is lost domain)
% - imgSource: An image from which patches are extracted, same size as
% targetImg.
% - NNF: nearest neighbor field
% - patchSize: patch size (patchSize x patchSize). Default is 9. 
% - maxIterations: maximum number of iterations during PatchMatch iteration. Default is 5. 
% 
% Outputs:
% - NNF: Nearest Neighbor Field, which contains indices of imgSource for 
%        each corresponding indices of imgTarget
% - intermediate: debugging information.
%

function NNF = PatchMatch(imgTarget, imgSource, offsetW, NNF, patchDistances, patchSize, reverseOrder)

sizeSource = size(imgSource);%[size(imgSource, 1), size(imgSource, 2)];

epsilonDistance = 10;

% define decreasing window sizes for random search steps
radius = sizeSource(1)/4;
alpha = .5;
radii = round(radius * alpha.^(0:(-floor(log(radius) / log(alpha)))));
randomSteps = length(radii);

if mod(patchSize, 2) == 0
    error('patchSize must be odd.');
end

%pD equals to distance from patch center to patch edge (half patch size)
pD = (patchSize - 1) / 2;

numPatches = size(NNF, 1) * size(NNF, 2);

 %progress report
progress = 0;
msg = sprintf('PatchMatch: %3.0f percent completed.', progress);
fprintf(msg);
reverseStr = repmat(sprintf('\b'), 1, length(msg));

yAxis = 1:size(NNF, 1); 
xAxis = 1:size(NNF, 2);

if reverseOrder % even iteration
    %reverse raster scan order
    yAxis = flip(yAxis); 
    xAxis = flip(xAxis);
end

for py = yAxis
  for px = xAxis
    y = (py - 1) * offsetW + 1;
    x = (px - 1) * offsetW + 1;
    
    targetPatch = imgTarget(y : y+2*pD, x : x+2*pD, :);

    currentNNF = NNF(py, px, :);%squeeze(NNF(py, px, :));

    currentDistance = patchDistances(py, px);
    % CHECK IF THIS HELPS: if patchDistance(py, px) is lower than predefined threshold, continue.
    if(currentDistance < epsilonDistance) 
        continue;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%            Propagation           %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if reverseOrder     %propagate from top and left
        dy = +1; dx = +1;
    else                %propagate from bottom and right
        dy = -1; dx = -1;
    end

    neighborDistances(1) = currentDistance; %current
    neighborDistances(2) = patchDistances(clamp(py+dy, 1, size(NNF, 1)), px); %top/bottom
    neighborDistances(3) = patchDistances(py, clamp(px+dx, 1, size(NNF, 2))); %left/right
    [~, neighborPatch] = min(neighborDistances);
    checkNeighbor = false;
    neighborNNF = currentNNF;
    switch neighborPatch
        case 2 %vertical neighbor 
            if inRange(NNF(py+dy, px, 1) - dy, 1 + pD, sizeSource(1) - pD)
                neighborNNF = NNF(py+dy, px, :);
                neighborNNF(1) = neighborNNF(1) - dy;
                checkNeighbor = true;
            end
        case 3 %horizontal neighbor
            if inRange(NNF(py, px+dx, 2) - dx, 1 + pD, sizeSource(2) - pD)
                neighborNNF = NNF(py, px+dx, :);
                neighborNNF(2) = neighborNNF(2) - dx;
                checkNeighbor = true;
            end
    end

    if checkNeighbor      
        neighborCandidatePatch = imgSource(neighborNNF(1)-pD : neighborNNF(1)+pD, neighborNNF(2)-pD : neighborNNF(2)+pD, :);    
        neighborDistance = MSQdistance(targetPatch, neighborCandidatePatch);
        if(neighborDistance < currentDistance)
            %use mean squared error as distance measure
            currentDistance = neighborDistance;
            currentNNF = neighborNNF;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%           Random Search          %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    minRandY = max(1 + pD, currentNNF(1) - radii);
    maxRandY = min(currentNNF(1) + radii, sizeSource(1)-pD);
    minRandX = max(1 + pD, currentNNF(2) - radii);
    maxRandX = min(currentNNF(2) + radii, sizeSource(2)-pD);

    randY = minRandY + floor((maxRandY-minRandY+1) .* rand(1, randomSteps));
    randX = minRandX + floor((maxRandX-minRandX+1) .* rand(1, randomSteps));

    randomDistances = zeros(randomSteps, 1);
    for randomStep = 1:randomSteps
        rY = randY(randomStep);
        rX = randX(randomStep);

        randomPotentialPatch = imgSource(rY-pD : rY+pD, rX-pD : rX+pD, :);
        %use mean squared error as distance measure
        randomDistances(randomStep) = MSQdistance(targetPatch, randomPotentialPatch);        
    end

    [minDistance, bestRandomPatch] = min(randomDistances);
    %if found better patch through randomized search, update distance
    %and NNF to reflect finding
    if(minDistance < currentDistance)
        currentDistance = randomDistances(bestRandomPatch);
        currentNNF = [randY(bestRandomPatch);  randX(bestRandomPatch)];
    end

    NNF(py, px, :) = currentNNF;
    if(currentDistance > patchDistances(py, px))
        fprintf('something is fishy >{|||*>');
    else 
        patchDistances(py, px) = currentDistance;
    end
    
    %keep track of progress
    if reverseOrder
        currentProgress = round(100 * (1 - ((py-1)*size(NNF, 2)+px) / numPatches));
    else
        currentProgress = round(100 * ((py-1)*size(NNF, 2)+px) / numPatches);
    end

    %display progress in console
    if(currentProgress > progress)
        progress = currentProgress;
        msg = sprintf('PatchMatch: %3.0f percent completed.', progress);
        fprintf([reverseStr, msg]);
        %clever trick: \b is the backspace character
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
  end
end

fprintf('\npatch distance between %d and %d', min(patchDistances(:)), max(patchDistances(:)));
fprintf('\n');
end % end of function

