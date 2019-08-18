function dat= getAOsystemMat(iter,cLensletsIn,diameterIn,varargin)
% Run basic AO simulator and get main matrices
% AOreshape.m
% Date:             24.04.2019
% Authors:          Daniel Gedon, 4735226
% Description:      This function runs the basic AO simulator and does a
%                   basic system identifcation such that it can be used for
%                   a Kalman filter of the system
%% variable input

Aid= 0;
if nargin > 3
    Aid= varargin{1};
end

%% basic setup of parameters

%   Sensor
% Number of pixels in one dimensions in the sensor
cLenslets = cLensletsIn;
% SNR = 15;
%   Telescope
% [m] Diameter of the telescope
cDiameter = diameterIn; % 1

%   Sampling frequency
% % [Hz] Frequency at which the loop is executed
% cFs = 500;
%   Turbulence
% % Number of phase points to consider for calculating the slopes (default = 2)
% cLensletPhasepoints = 4;
% Number of random samples to consider for changing R0 and L0
cR0 = 0.5;
cL0 = 25;   10;
% [px/sample] Speed at which the wind moves in horizontal direction
cWindX = 1; [1 0];
% [px/sample] Speed at which the wind moves in vertical direction
cWindY = 0; [0 1];
%  Iterations
cIdent = iter;

% Sensor
[G,cResolution] = SHSensor(cLenslets);
ng = size(G,1);
m2 = size(G,2);
% Pixel size
nPxSize = cDiameter/cResolution;
% % Time length of each sample
% nDt = 1/cFs;
% % Number of phase points in one dimension (number of columns of s(k))
% s = cLenslets;
% % Number of rows of s(k)
% s2 = 2 * s;

%   Deformable mirror
% % Factor the number of actuators with respect to number of pixels
% % (Fried geometry ->  cGeometry = 1)
% cGeometry = 1;
% % Gaussian falloff of the actuators in the DM
% cCoupling = 0.2;
% % Number of actuators
% nActuators = 5;%sqrt(size(G,2));
% % Matrix relating actuators to phase points
% H = deformableMirror(cResolution,nActuators,cDiameter,cCoupling);

%% Identify the matrices

phiTemp= cell(1,length(cWindX));
for i = 1:length(cWindX)
    phiTemp{i} = TurbulenceMAws(nPxSize, cResolution, cR0, cL0, cWindX(i), cWindY, cIdent);
end
phiIdent = zeros(size(phiTemp{1}));
for i = 1:length(cWindX)
    phiIdent = phiIdent+phiTemp{i};
end

% Initiate the variables we need
% epsilon_k = zeros(m2,cIdent);
% sIdent = zeros(ng,cIdent);

% Do a simple open loop simulation so we can get data for identification
% for k=1:cIdent-1
%     epsilon_k(:,k+1) = phiIdent(:,k+1);
%     sIdent(:,k+1) = G * epsilon_k(:,k+1) + sigma*randn(ng,1);
% end

% AR1 estimation
if Aid
    a= 0.99;
    A= a*eye(size(G,2));
    %     Q= (1-a^2)*cov(phiIdent');
else
    A = phiIdent(:,2:end) / phiIdent(:,1:end-1);
end
% identify Q
if 1
    % always use this method. The analytical one does not consider the
    % cross terms which are actually present in case that A=a*I.
    e = phiIdent(:,2:end) - A * phiIdent(:,1:end-1);
    Q = cov(e');
else
    Cphi= cov(phiIdent');
    Qana = Cphi - A*Cphi*A';
    Q= 0.5*(Qana+Qana');
    % Q(Q<1e-3)= 0;
end

%% Output

dat.A= A;
dat.C= G;
dat.Q= Q;