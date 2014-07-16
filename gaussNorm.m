function y = gaussNorm(x, tau)
% GAUSSNORM computes the Gaussian of the squared 2-norm of a vector x, with scaling coefficient tau.
%     Notice that (gaussNorm^{-1} - 1) is a metric. 

% fix weird compatibility issue with getWeights.m usage
if iscell(tau)
    tau = cell2mat(tau);
end
y = exp(-sum(x.^2, 1)./tau); % sum over rows / sum each feature matrix
% sum(...,1) means "collapse rows in sum"
end