function [Phi, Lambda] = nystrom(Wxx, Wxy, varargin)
% NYSTROM approximates eigenvectors, eigenvalues for a matrix W = [Wxx, Wxy; Wyx, Wyy]
% varargin{1} : displaying timing? 

switch nargin-2
    case 0
        verbose = 1;
    case 1
        verbose = varargin{1};
    otherwise
        error('nystrom takes min 2 args and max 3.');
end

if verbose
    display('Computing eigenvectors of sampled subspace, components of orthognal approximation...');
    tic;
end
[Bx, Gamma] = eig(Wxx);
sqrtGamma = Gamma^(1/2);
S = Bx*(sqrtGamma\Bx.');
Q = Wxx + S*(Wxy*Wxy.')*S;
[A Xi] = eig(Q);
if verbose
    toc;
end

%%
if verbose
    display('Approximating leading eigenvectors of graph Laplacian on orthogonal subspace...');
    tic;
end
Phi = [Bx*sqrtGamma; Wxy.'*Bx/sqrtGamma]*Bx.'*(A/(Xi^(1/2)));
Lambda = 1- diag(Xi);
if verbose
    toc;
end
