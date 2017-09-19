function [ target, NNF ] = textureSynthesis(source, target, NNF, patchSize, maxIterations, offsetW)
    targetHeight = size(target, 1);
    targetWidth  = size(target, 2);
    sizeTarget = [targetHeight, targetWidth];
    
    pD = floor(patchSize / 2);
    if(nargin < 6)
        offsetW = floor(patchSize / 4);
    end

    if( isempty(target) )
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%                Initialize target                %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %random initialization at given resolution
        target = zeros(targetHeight, targetWidth, size(source, 3));

        for y = 1 : patchSize : targetHeight
            for x = 1 : patchSize : targetWidth     
                randSourceY = randi([1 + pD, size(source, 1) - pD], 1);
                randSourceX = randi([1 + pD, size(source, 2) - pD], 1);
                %pull a random patch from source and assign to target texture
                target(y : y + patchSize - 1, x : x + patchSize - 1, :) = source(randSourceY - pD : randSourceY + pD, randSourceX - pD : randSourceX + pD, :);
            end
        end
    end
    
    if( isempty(NNF) )
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%                   Initialize NNF                %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        NNFsize = floor((sizeTarget - 2*pD) ./ offsetW);% - 2*pD;
        % initialize NNF with random patch coordinates 
        NNF = cat(3, randi([1 + pD, size(source, 1) - pD], NNFsize ), randi([1+  pD, size(source, 2)-pD], NNFsize ));

        % CONSIDER: is it necessary to calculate offset in advance? 
        patchDistances = inf(NNFsize);
        for py = 1:NNFsize(1)
           y = (py - 1) * offsetW + 1 + pD;
           for px = 1:NNFsize(2)
               NNy = NNF(py, px, 1);
               NNx = NNF(py, px, 2);

               x = (px - 1) * offsetW + 1 + pD;

               targetPatch = target(y - pD : y + pD, x - pD : x + pD, :);
               %target(y - pD : y + pD, x - pD : x + pD, :) = targetPatch;%target(y - pD : y + pD, x - pD : x + pD, :) 
               %initialize distance
               currentNearestSourcePatch = source(NNy - pD : NNy + pD, NNx - pD : NNx + pD, :);
               %use mean squared error as distance measure
               patchDistances(py, px) = MSQdistance(targetPatch, currentNearestSourcePatch);
           end
        end
    end

    subplot(1, 4, 4); imshow (uint8(target));
    subplot(1, 4, 2); imagesc(NNF(:, :, 1));
    subplot(1, 4, 3); imagesc(NNF(:, :, 2));

    % start iterations
    for iteration = 1 : maxIterations
        fprintf('Iteration %d | ', iteration);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% STEP 1: find nearest neighbors using PatchMatch %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        t1 = tic;
        scanLineOrder = mod(iteration, 2) == 1;
        NNF = PatchMatch(target, source, offsetW, NNF, patchDistances, patchSize, scanLineOrder);

        t2 = toc(t1); 
        t3 = tic;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% STEP 2: minimize the global texture energy      %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        target_current = target;
        target = zeros(targetHeight, targetWidth, size(source, 3));
        Xf = zeros(targetHeight, targetWidth);
        targetWeight = zeros(targetHeight, targetWidth);
        fs = fspecial('gaussian', patchSize, patchSize);

        %step through all patches
        for py = 1:NNFsize(1)
           y = (py - 1) * offsetW + 1 + pD;
           for px = 1:NNFsize(2)
                x = (px - 1) * offsetW + 1 + pD;

                NNy = NNF(py, px, 1);
                NNx = NNF(py, px, 2);

                % patch distance
                dx = target_current(y - pD : y + pD, x - pD : x + pD) - source(NNy - pD : NNy + pD, NNx - pD : NNx + pD, :);
                dx = dx.^2;

                weight = (1.0e-5 + sum(dx(:)))^(-0.6);%(-1.2);%

                % gaussian blur
                dz = source(NNy - pD : NNy + pD, NNx - pD : NNx + pD, :);
                for color = 1:size(source, 3)
                    dz(:, :, color) = dz(:, :, color) .* fs;
                end
                target(y - pD : y + pD, x - pD : x + pD, :) = target(y - pD : y + pD, x - pD : x + pD, :) + dz.*weight;
                targetWeight(y - pD : y + pD, x - pD : x + pD) = targetWeight(y - pD : y + pD, x - pD : x + pD) + weight*fs;
                Xf(y - pD : y + pD, x - pD : x + pD) = ones(patchSize);
            end
        end		

        % normalize
        for y = 1:targetHeight
            for x = 1:targetWidth
                if Xf(y, x) ~= 0
                    target(y, x, :) =  target(y, x, :) ./ targetWeight(y, x);
                else
                    target(y, x, :) = target_current(y, x, :);
                end                    
            end
        end

        t4 = toc(t3);
        % result
        str1 = num2str(t2);
        str2 = num2str(t4);
        disp(['timing nearest neighbors: ' str1]);
        disp(['timing minimization: ' str2]);

        subplot(1, 4, 2); imagesc(NNF(:, :, 1));
        subplot(1, 4, 3); imagesc(NNF(:, :, 2));
        subplot(1, 4, 4); imshow (uint8(target));
        drawnow;
    end
end