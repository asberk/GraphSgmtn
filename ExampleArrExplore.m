% This script is an example of how to use the arrExplore function to
% explore feature vectors. Note that the getWeights functions is
% computationally intensive to run on large image sizes NN.


%% Clear workspace
clear all; close all; clc;

%% Define initial parameter values
% must run this before any of the other sections work
NN = 100; % domain length (square)
N = 3; % neighbourhood radius
tau = .2;

[xx, yy] = meshgrid(linspace(-1,1,NN));
inputImg = xx.^2 + yy.^2 < .5;
szu0 = size(inputImg);

%% meat space
% depends on "Define..." section
FtrArr = getFeatureArr(inputImg, N); % Neumann, 2D matrix
W = getWeights(FtrArr, @gaussNorm, 1, 1, tau);
W_arr = reshape(W, [szu0, szu0]);
%arrExplore(W_arr, 1, 1);
%pause;

%% freq space
% depends on "Define..." section
myFFT = @(x) abs(fftshift(fft2(x)));
FreqArr = getFeatureArr(inputImg, N, myFFT);
W2 = getWeights(FreqArr, @gaussNorm, 1, 1, tau);


%% wvlt space
% depends on "Define..." section
wname = 'haar';
wMaxLev = 2; % wmaxlev([7 7], 'haar') = 2
FtrArr = getFeatureArr(inputImg, N); % Neumann, 2D matrix
[asdf, S] = wavedec2(reshape(FtrArr(:, 1), [2*N+1, 2*N+1]), wMaxLev, wname);
WvltArr = zeros(length(asdf), NN.^2);
for j = 1:NN.^2
    WvltArr(:, j) = wavedec2(reshape(FtrArr(:,j), [2*N+1, 2*N+1]), wMaxLev, wname);
    if (mod(j,100)==0) 
        display(j);
    end
end
W3 = getWeights(WvltArr, @gaussNorm, 1, 1, tau);
W3_arr = reshape(W3, [szu0, szu0]);
