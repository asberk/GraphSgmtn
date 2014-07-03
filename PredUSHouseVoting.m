% % PredUSHouseVoting - Designed to test Bertozzi and Flenner's example, US House of Representatives Voting Records from 1984.
%     This script file loads a data set used by Bertozzi and Flenner to
%     verify that our implementation of the diffuse interface method for
%     arbitrary graphs yields the same results as that of Bertozzi and
%     Flenner, section 4.1. See the help file for PredImgSgmtn.m for more 
%     detail. 
%
% This script depends on the file(s):
% * gaussNorm.m
% * getGraphLaplacian.m
% * getImage.m
% * getWeights.m
% * imgnmz.m

%% clear the workspace
clear all; close all; clc;

%% Load and format dataset
fid = fopen('data/house-votes-84.data');
tmp = textscan(fid, '%s', 'Delimiter', ',');
tmp = reshape(tmp{1}, [17, 435]);
fclose(fid);
hv84 = strcmp('y', tmp) - strcmp('n', tmp);
hv84 = hv84(2:end, :);
DRlabels = strcmp('democrat', tmp(1,:)) - strcmp('republican',tmp(1,:));
clear tmp fid

%% Initialize parameters, variables
M = 750; % Number of iterations for convergence
tau = .2;
%% Get weights

W = getWeights(hv84, @gaussNorm, 1, 1, tau);

%% Check if W is positive semi-definite
[R, p] = chol(W);

if p ~= 0 % then W is not positive semi definite:
    warning('W is not positive semi-definite.');
end

%% Compute graph Laplacian

[GL, Deg] = getGraphLaplacian(W); % use default settings to return symmetric graph Laplacian

%% Get EigenInformation of Graph Laplacian
[Phi,Lambda] = eig(GL);
Phi = real(Phi);
Lambda = real(diag(Lambda));

%% Convex splitting for Graph Laplacian
c = 1; 
dt = .1; 
epsilon = 2; 

% initialize (all are 1-by-L vectors)
u0 = zeros(size(DRlabels));
u0(1:5) = DRlabels(1:5);

display('Running convex splitting scheme for graph Laplacian...');
tic;
a = u0*Phi; 
b = (u0.^3)*Phi; 
d = zeros(1, 435);
D = 1 + dt*(epsilon*Lambda.' + c);
eta = abs(u0); %fidelity term

for m = 1:M
    a = ((1+dt/epsilon + c*dt)*a - dt/epsilon*b - dt*d)./D;
    u = a*Phi.';
    b = (u.^3)*Phi;
    d = (eta.*(u-u0(:).'))*Phi;

end
toc;


%%

DRpred = 2*double(u >= 0)-1;

display(sum(DRpred == DRlabels)/length(DRlabels));

%% 
figure; plot(u)