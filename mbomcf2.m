function Aout = mbomcf2(A, alphadt, varargin)
% MBOMCF2 implements the Merriman-Bence-Osher algorithm for mean curvature flow of a curve C = \partial\Sigma where \Sigma is given as the characteristic function \chi_{int C}: Reals^2 \to Reals of the interior of C.
%     Note that:
%     * This algorithm is for two-dimensional stuff only (uses fft2)
%     Input parameters:
%     alphadt = alpha*dt = (heat flux parameter)*(discrete time step)
%     varargin{1} = T = number of time-steps = 500 (default)
%     varargin{2} = co = cutoff value = 1/2 (default)
%     varargin{3} = visualize = whether to show visualization

if nargin-2 > 3
    error('Whoops! Too many input arguments.');
end
if nargin-2 > 0
    T = varargin{1};
else
    T = 500;
end
if nargin-2 > 1
    co = varargin{2};
else
    co = .5
end
if nargin-2 > 2
    visualize = varargin{3};
else 
    visualize = 0;
end

N = size(A);
Nx = N(1);
Ny = N(2);

kx = (1i*[0:(Nx+1)/2-1 0 -(Nx+1)/2+1:-1]);
ky = (1i*[0:(Ny+1)/2-1 0 -(Ny+1)/2+1:-1].');
k2x = kx.^2;
k2y = ky.^2;
[kxx, kyy] = meshgrid(k2x, k2y);

for tau = 1:T
    % diffusion step
    % implicit spectral method for heat equation with Neumann BC
    A_hat = fft2(A);
    Anew = A_hat./(1 - (kxx + kyy)*alphadt);
    A = ifft2(Anew);

    % cut-off step
    if isnumeric(co)
        A = double(A>=co);
    elseif ischar(co)
        switch co
            case {'volume-preserving', 'volume preserving', 'vp', 'VP'}
                lambda = VPsharpen(A, A_hat, 1/2, 20);
                A = double(A >=lambda);
            otherwise
                error('co must be numeric or ''volume-preserving'' at present.');
        end
    end
    if visualize
        imshow(A); title(['Time ', num2str(tau)]); drawnow;
    end
end
Aout = A;
end

