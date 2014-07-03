function y = gaussNorm(x, tau)
if iscell(tau)
    tau = cell2mat(tau);
end
y = exp(-sum(x.^2, 1)./tau); % sum over rows / sum each feature matrix
end