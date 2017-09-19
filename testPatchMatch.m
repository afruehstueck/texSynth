clear all;
close all;
clc;

SaveFolderName = datestr(now,'yymmdd_HHMM');
mkdir('results', SaveFolderName);

%diary(fullfile('results', SaveFolderName, 'log.txt'));

pathSource = 'data/BrickOldMixed.jpg';
pathReconst = 'data/BrickRound.jpg';

imgTarget = imread(pathReconst); 
imgSource = imread(pathSource);

if size(imgTarget, 3) == 3
    disp('convert input to grayscale');
    imgTarget = rgb2gray(imgTarget);
end

if size(imgSource, 3) == 3
    disp('convert source to grayscale');
    imgSource = rgb2gray(imgSource);
end

numFigs = 4;
figIndex = 1;

figure;
screensize = get( groot, 'Screensize' );
set(gcf, 'Position', [200, 500, screensize(3) *0.75, screensize(4) *0.4]);

subplot(1, numFigs, figIndex);
figIndex = figIndex+1;
imshow(imgTarget);
title('Original');
drawnow;

imgTarget = double(imgTarget);
imgSource = double(imgSource);

sizeSource = [size(imgSource, 1), size(imgSource, 2)];
sizeTarget = [size(imgTarget, 1), size(imgTarget, 2)];

patchSize = 15;
offset_w = 3;

maxIterations = 5;
pad imgTarget with pD
pD = (patchSize - 1) / 2;
padTarget = nan(sizeTarget(1) + 2*pD, sizeTarget(2) + 2*pD);
padTarget(1+pD : sizeTarget(1)+pD, 1+pD : sizeTarget(2)+pD) = imgTarget;
imgTarget = padTarget;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%          Initialize NNF          %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NNFsize = floor(sizeTarget / offset_w);
% initialize NNF with random patch coordinates 
NNF = cat(3, randi([1+pD, sizeSource(1)-pD], NNFsize ), randi([1+pD, sizeSource(2)-pD], NNFsize ));

% CONSIDER: is it necessary to calculate offset in advance? 
patchDistances = inf(NNFsize);
for py = 1:NNFsize(1)
   for px = 1:NNFsize(2)
       y = (py - 1) * offset_w + 1;
       x = (px - 1) * offset_w + 1;
       targetPatch = imgTarget(y : y + 2 * pD, x : x + 2 * pD);
       %initialize distance
       currentNearestSourcePatch = imgSource(NNF(py, px, 1) - pD : NNF(py, px, 1) + pD, NNF(py, px, 2) - pD : NNF(py, px, 2) + pD);
       %use mean squared error as distance measure
       patchDistances(py, px) = MSQdistance(targetPatch, currentNearestSourcePatch);
   end
end

intermediate.distances = patchDistances;
intermediate.NNF = NNF;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%             Main loop            %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
for iteration = 1:maxIterations
    fprintf('Iteration %d | ', iteration);
    
    scanLineOrder = mod(iteration, 2) == 1;
     
    NNF = PatchMatch(imgTarget, imgSource, offset_w, NNF, patchDistances, patchSize, scanLineOrder);



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%          Reconstruction          %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    sizeReconst = size(padTarget);
    imgReconst1 = zeros(sizeReconst);
    imgReconst2 = zeros(sizeReconst);
    imgReconst3 = zeros(sizeReconst);
    countPatches = zeros(sizeReconst); %count number of patches that influence any pixel value in imgReconst3
    for py = 1:NNFsize(1)
       for px = 1:NNFsize(2)
    % for y = 1+pD : patchSize : sizeTarget(1)-pD
    %     for x = 1+pD : patchSize : sizeTarget(2)-pD
            y = (py - 1) * offset_w + 1;
            x = (px - 1) * offset_w + 1;
            if(mod(y - pD, patchSize) == 0 && mod(x - pD, patchSize) == 0) 
                imgReconst1(y : y + 2*pD, x : x + 2*pD) = imgSource(NNF(py, px, 1) - pD : NNF(py, px, 1) + pD, NNF(py, px, 2) - pD : NNF(py, px, 2) + pD);
            end
            dl = floor(offset_w/2);
            dr = ceil(offset_w/2);
            imgReconst2(y +pD - dl : y +pD + dr, x+pD - dl : x+pD + dr) = imgSource(NNF(py, px, 1) - dl : NNF(py, px, 1) + dr, NNF(py, px, 2) - dl : NNF(py, px, 2) + dr);
            imgReconst3(y : y + 2*pD, x : x + 2*pD) = imgReconst3(y : y + 2*pD, x : x + 2*pD) + double(imgSource(NNF(py, px, 1) - pD : NNF(py, px, 1) + pD, NNF(py, px, 2) - pD : NNF(py, px, 2) + pD));
            countPatches(y : y + 2*pD, x : x + 2*pD) = countPatches(y : y + 2*pD, x : x + 2*pD) + 1;
        end
    end
    imgReconst3 = imgReconst3 ./ countPatches; %divide by number of influencing patches

    imgReconst1 = uint8(imgReconst1);%
    imgReconst2 = uint8(imgReconst2);%
    imgReconst3 = uint8(imgReconst3);%

    figIndex = 2;   
    subplot(1, numFigs, figIndex);
    figIndex = figIndex+1;
    imshow(imgReconst1);

    subplot(1, numFigs, figIndex);
    figIndex = figIndex+1;
    imshow(imgReconst2);
    
    subplot(1, numFigs, figIndex);
    figIndex = figIndex+1;
    imshow(imgReconst3);
    drawnow;
end
toc
disp('PatchMatch done!');
% 
% if(size(imgReconst1) == size(imgSource)) 
%     PSNR = psnr(double(imgReconst1), double(imgSource), 255);
% 
%     psnrout = sprintf('Copy patches, no blending\nPSNR is %.4f\n',PSNR);
%     fprintf(psnrout);
%     title(psnrout);
% end
% if(size(imgReconst2) == size(imgSource))
%     PSNR = psnr(double(imgReconst2), double(imgSource), 255);
% 
%     psnrout = sprintf('Copy center pixel of each patch\nPSNR is %.4f\n',PSNR);
%     fprintf(psnrout);
%     title(psnrout);
% end
% if(size(imgReconst3) == size(imgSource)) 
%     PSNR = psnr(double(imgReconst3), double(imgSource), 255);
% 
%     psnrout = sprintf('Average patches\nPSNR is %.4f\n',PSNR);
%     fprintf(psnrout);
%     title(psnrout);
% end
%     

%imwrite(patchReconst,fullfile('results', SaveFolderName, 'reconstImg_ini.bmp'),'BMP');
imwrite(imgReconst1,fullfile('results', SaveFolderName, 'reconstructedImg.png'), 'PNG');

%diary off

