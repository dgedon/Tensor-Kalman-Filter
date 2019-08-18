function [ G,m ] = SHSensor( m )
%SHSENSOR Create a Wavefront Sensor matrix G
%   G = SHSensor(m,n), returns the matrix G in the model for the
%   Shack-Hartmann output measurement given as:
%
%   s(k) = G (phi(k) - phi_dm(k)) + v(k)
%
%   For the parameters:
%
%       m       Resolution of the wavefront in one dimension
%
%   Note that the number of phase points in one dimension is defined as:
%       p = m - 1
%   Since the phase is evaluated at the intersections of the wavefront.
%   
%   This function has only been implemented for square phase screens and
%   evenly distributed phase points.
%
% G. Monchen <guido.monchen@gmail.com>
% DCSC, TU Delft (2017)

% Generate the G matrix which defines the wavefront sensor
n = m;
G = zeros(2 * (m-1) * (n-1), m*n);

for i=1:m-1
    for j=1:n-1
        % Calculate the horizontal gradient first
        Gtemp = zeros(m,n);
        Gtemp(i,j) = -0.5;
        Gtemp(i,j+1) = 0.5;
        Gtemp(i+1,j) = -0.5;
        Gtemp(i+1,j+1) = 0.5;
        temp = 2*(i+(j-1)*(m-1))-1;
        G(temp,:) = Gtemp(:)';
        
        % For the vertical gradient
        Gtemp = zeros(m,n);
        Gtemp(i,j) = -0.5;
        Gtemp(i,j+1) = -0.5;
        Gtemp(i+1,j) = 0.5;
        Gtemp(i+1,j+1) = 0.5;
        temp = 2*(i+(j-1)*(m-1));
        G(temp,:) = Gtemp(:)';
    end
end

end

