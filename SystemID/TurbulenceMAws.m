function [ phi ] = TurbulenceMAws(pxSize, m, r0, L0, wsX,wsY, Nt)
%TurbulenceMAws Generate nonstationary wind speed turbulence according to the Moving Average method [1].
%   phi = TurbulenceMAws(pxSize,m,r0,L0,wsX,wsY,Nt), returns a 
%   square turbulence phase screen of size (m^2 x Nt) which is generated
%   based on the Moving Average method in [1] detailed in the section about
%   generation of low resolution turbulence.
%
%   For the parameters:
%
%       pxSize      The individual pixel size [m]
%       m           Resolution of the wavefront in one dimension
%       r0          The Fried parameter [m]
%       L0          Outer scale of spatial coherence [m]
%       dt          The sample time [s]
%       wsX         Either a vector or a scalar containing the piece-wise 
%                   constant wind speeds in x direction [pixels/sample]
%       wsY         Either a vector or a scalar containing the piece-wise 
%                   constant wind speeds in y direction [pixels/sample]
%       Nt          The number of time samples to generate
%
%   This method bases the generation of turbulence on a moving average (MA)
%   filter where a random process e is convoluted with some parameters
%   theta that are calculated based on the inverse FFT of the square root
%   of the FFT of the covariance C_phi.
%
%   The nonstationary turbulence is generated by calculating an oversized
%   phase screen and moving a smaller aperture over it. Wind speed is
%   simulated as piece-wise constant wind speeds varying over the entire
%   simulation according to the distribution in wsX/wsY. You must specify
%   the wind speed in pixels/sample.
%
%   References:
%   [1] Beghi et. al. "Multiscale phase screen synthesis based on local
%       principal component analysis". Opt. Soc. Am. A (2013).
%
% G. Monchen <guido.monchen@gmail.com>
% DCSC, TU Delft (2017)

nX = length(wsX);

if nX > 1
    niterpart = ceil(Nt / nX);

    xArr = ones(niterpart * wsX(1),1) * wsX(1);
    for i = 2:nX
        xArr = [xArr; ones(niterpart * wsX(i),1) * wsX(i)];
    end
    m01 = Nt*wsY + m;
    n01 = length(xArr) + m - 1;
else
    m01 = Nt*wsY + m;
    n01 = Nt*wsX + m - 1;
end


% Neighbourhood
de = 500;
d = 5;
dn = (-de:1:de)*pxSize(1);

C_phi0 = spatCov(dn,r0(1),L0(1));

% The process x0 can be defined as:
% x0(u,v) = SUM(teta_0(u, v) * e_0(u - k_u, v - k_v))
% Using fft2 we can calculate the parameters teta_0
C_phi_FT = fft2(C_phi0);
C_teta_FT = sqrt(C_phi_FT);
C_teta = ifft2(C_teta_FT);

%e0 = randn(dn+m0+1);
e0 = randn(m01+2*d,n01+2*d);

% This describes the turbulence at the lowest resolution
x = filter2(C_teta,e0);

phi = zeros(m^2,Nt);

fg = 0;
j = 0;
q = 1;
for i = 1:Nt
    
    if nX > 1
        if mod(i,niterpart) == 0 && q ~= nX
            fg = fg + j * wsX(q);
            j = 0;
            q = q+1;
        end
    end
    
    phi_temp = x(i*wsY+1+d:i*wsY+1+m+d-1,fg+j*wsX(q)+d+1:fg+j*wsX(q)+m+d);

    phi(:,i) = phi_temp(:);
    j = j + 1;
end

end