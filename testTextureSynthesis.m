clear all;
close all;
clc;

%pathSource = 'data/bricks.jpg';
pathSource = 'data/Texture-01.png';
% pathSource = 'data/water.jpg';
% pathSource = 'data/weave.jpg';
% pathSource = 'data/grapes.jpg';
source = imread(pathSource);

% if size(source, 3) == 3
%     disp('convert source to grayscale');
%     source = rgb2gray(source);
% end

numFigs = 4;
figIndex = 1;

figure;
screensize = get( groot, 'Screensize' );
set(gcf, 'Position', [200, 500, screensize(3) *0.75, screensize(4) *0.4]);

subplot(1, numFigs, figIndex);
figIndex = figIndex+1;
imshow(source);
title('Source texture');
drawnow;

%convert to double
source = double(source);
sizeSource = size(source);

patchSize = 15;
pD = floor(patchSize / 2);
offsetW = floor(patchSize / 4);

maxIterations = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%     Initialize target texture    %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
resX = 200;
resY = 200;

targetHeight = resY / 2;
targetWidth  = resX / 2;
target = double.empty(targetHeight, targetWidth, 0);
NNF = double.empty(0, 0);

[ target, NNF ] = textureSynthesis(source, target, NNF, maxIterations, patchSize);

%%
%% upsample texture and repeat
% initialize X for level 2

maxIterations = 30;

new_target = zeros(resY, resX, size(source, 3));
new_target(1:2:end, 1:2:end, :) = target;
target = 4*imfilter(new_target, fspecial('gaussian'), 'replicate');

targetHeight = resY;
targetWidth  = resX;
sizeTarget = [targetHeight, targetWidth];

NNF = double.empty(0, 0);
[ target, NNF ] = textureSynthesis(source, target, NNF, patchSize, maxIterations);