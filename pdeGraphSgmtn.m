function u = pdeGraphSgmtn(u0, Phi, Lambda, Parms, varargin)
% PDEGRAPHSGMTN(u0, Phi, Lambda, Parms) performs a graph cutting algorithm
% using a nonlinear PDE solved by gradient descent with a convex splitting
% scheme. 
% varargin{1} : verbose?

switch nargin-4
    case 0
        verbose = 1; % default
    case 1
        verbose = varargin{1};
    otherwise
        error('Whoops! Wrong number of arguments to pdeGraphSgmtn.');
end

% extract Parms
dt = Parms.dt;
epsilon = Parms.epsilon;
L = Parms.L;
M = Parms.M;

if ~isfield(Parms, 'c') || isempty(Parms.c)
    c = 1; % this option not recommended due to potential for convexity issues. 
else
    c = Parms.c;
end
if ~isfield(Parms, 'eta') || isempty(Parms.eta)
    eta = abs(u0);
else
    eta = Parms.eta;
end

% check condition on convex splitting scheme
if c < 0 || -2/epsilon + c - 1 > 0
    warning('functionals for convex splitting may not be convex! consider re-choosing c and/or epsilon');
end

% initialization step (all are 1-by-L vectors)
a = u0*Phi; 
b = (u0.^3)*Phi; 
d = zeros(1, L);
D = 1 + dt*(epsilon*Lambda.' + c);

if verbose
    tic;
end

for m = 1:M
    a = ((1+dt/epsilon + c*dt)*a - dt/epsilon*b - dt*d)./D;
    u = a*Phi.';
    b = (u.^3)*Phi;
    d = (eta.*(u-u0))*Phi;

end

if verbose
    toc;
end
