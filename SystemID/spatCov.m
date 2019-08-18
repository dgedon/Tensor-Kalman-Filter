function [ C ] = spatCov( d, r0, L0 )
%SPATCOV Create a spatial covariance matrix
%   C = spatCov(d,r0,L0), returns the matrix C which is a spatial
%   covariance matrix based on the von Karman theory.
%
%   For the parameters:
%
%       d   A vector which contains the physical location of the pixels in
%           meters.
%       r0  The Fried parameter
%       L0  The outer scale parameter
%
%   This function has only been implemented for square phase screens, so
%   the vector d represents the location in both the u and v directions.
%
%   References:
%   Beghi et al. "Stochastic realization approach to the efficient
%   simulation of phase screens", Opt. Soc. Am. A, 2008.
%


dn = size(d,2);

% The constant value for calculating the spatial covariance according to
% the von Karman theory
c = (2^(1/6)*gamma(11/6))/(pi^(8/3))*((24/5)*gamma(6/5))^(5/6);

C = zeros(dn,dn);
    
% Calculate the spatial covariance matrix for 
for k = 1:dn
    for l = 1:dn        
        % We want to calculate the spatial covariance for -d0:+d0 of x0,
        % but only from the point (0,0), so:
        % C_phi0(u,v) = E[ x0(0,0) x0(u,v) ]
        dx = d(k);
        dy = d(l);        
        r = sqrt(dx^2 + dy^2);
        r = r+eps;

        % Calculate the value for a specific r according to 
        % Von Karman theory
        C(k,l) = (L0/r0)^(5/3)*c*0.5*(2*pi*r/L0)^(5/6) ...
                      * besselk(5/6,2*pi*r/L0);
    end
end
end

