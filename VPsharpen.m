function lambda = VPsharpen(A, varargin)
% VPsharpen computes the cut-off value, lambda, for the MBO scheme for volume-preserving evolution by mean curvature.
%           A : input matrix
% varargin{1} : either Fourier transform of A or [preferably] first coefficient, fft2(A)(1,1)
% varargin{2} : a/b = initial guess for threshold
% varargin{3} : maxIter = maximum number of iterations
% f(t) := m({\chi_A > t}) is an increasing function of t. So let's search
% for the value that matches what we want, realizing that t \in [0, 1] (no
% need to look outside of this range)

if (nargin - 1) >= 1
    volA = varargin{1};
    if size(volA,1)==0 || size(volA,2)==0
        error('Whoops! Fourier coefficient can''t have size zero!');
    end    
    volA = volA(1);
else
    volA = sum(A(:));
end
if (nargin - 2) >= 2
    [a,b] = rat(varargin{2}); % varargin{2} needs to be a number, not a matrix, but we won't check for this. 
else
    a = 1; b = 2;
end
if (nargin - 3) >= 3
    maxIter = varargin{3};
else
    maxIter = 500;
end

volCO = sum(A(:) >= a/b);
nIter = 0;

% TO FIX: What about the case where inf |volCO - volA| > 1? Do we always
% want to run to maxIter? That seems expensive...
while (abs(volCO - volA) > 1) && nIter<=maxIter
    nIter = nIter + 1;
    switch sign(volCO - volA)
        case 1 % volCO is greater
            a = 2*a + 1;
            b = 2*b;
            volCO = sum(A(:) >= a/b);
        case -1
            a = 2*a -1;
            b = 2*b;
            volCO = sum(A(:) >= a/b);
        otherwise
            % this means they're equal, but something weird happened for the loop to have continued.
            break;
    end
end
lambda = a/b;
end
    
