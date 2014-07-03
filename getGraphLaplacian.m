function [GL, varargout] = getGraphLaplacian(W, varargin)
% GETGRAPHLAPLACIAN computes the Graph Laplacian; returns symmetric graph Laplacian by default. 
% varargin{1} = ReturnSymmetricLaplacian (true by default)
% varargout{1} = Degree Matrix - diagonal matrix whose entries correspond
%                to the degree of each vertex
D = diag(sum(W,1) + sum(W,2).'); % Degree matrix for graph
L = D - W; % Laplacian
if (nargin>2) && varargin{1}==0 % not symmetric
    GL = L;
else
    invsqrtD = D^(-1/2); % takes a short while to invert D
    GL = invsqrtD*L*invsqrtD; % "Symmetric" Laplacian
end


nout = max(nargout,1)-1;
if nout == 1
    varargout(1) = {D};
end
if nout >= 2
    error('Too many output arguments');
end

end