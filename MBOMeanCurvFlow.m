%% Implementation of MBO scheme for mean curvature evolution
% It's worth noting that for certain initial conditions, that 
% (heat flux)(time) must be larger than some threshold to "kick" the
% evolution out of certain local minima. 
%% 
clear all; close all; clc;
%% Initial grid, time step
alpha = 1;
N = 255; h = 1./N;
xv = h*(0:N);
yv = xv.';
[x,y] = meshgrid(xv,yv);

dt = .01;
%% Initial surface
A0 = zeros(size(x));
% A0(x< .333 & x > .1 & y < .667 & y > .333) = 1;
% A0(x > .2 & x < .8 & y > .4 & y < .5) = 1;
% A0(x < .9 & x > .667 & y < .667 & y > .333) = 1;

% A(x<.667 & x > .333 & y < .667 & y > .333) = 1;
% A(x < .75 & x > .25 & y < .75 & y > .25) = .5;
% imshow(A);

for j = .125:.25:.875
    for k = .125:.25:.875
        A0((x-j).^2 + (y-k).^2 < (.102)^2) = 1;
    end
end

%% MBO scheme for evolution by mean curvature flow
% A = mbomcf(input,flux*time,Niter,SubSteps,Cut-Off)
A = mbomcf2(A0, alpha*dt, 1000, 'vp', 1);

%% Visualize
subplot(2,2,1);
imshow(A0, []);
subplot(2,2,2);
imshow(A, []);