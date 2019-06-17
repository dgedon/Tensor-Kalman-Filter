function [yhat,xhat,timeKF]= myETKF_FullSim(A,C,Q,R,N,t,y,varargin)
%% Function Kalman Filter without input
% myKF_FullSim.m
% Date:             13.11.2018
% Authors:          Daniel Gedon, 4735226
% Description:      The function does a recursive Kalman state estimation.
% Inputs:           A - state evolution matrix
%                   C - measurement matrix
%                   Q - process noise covariance matrix
%                   R - measurement noise covariance matrix
%                   t - time vector
%                   y - measurement vector
% Outputs:          yhat - estimated measurements
%                   xhat - estimated state
%                   timeKF - time for computation

%% Check Input

timeLimit= inf;
if nargin == 8
    timeLimit= varargin{1};
end

% state and output dimension
n= size(A,1);
p= size(C,1);

% time
kMax= max(t);

%% Algorithm

% initialization and allocation
X_= randn(n,N);
Z_= zeros(n,N);
xhat= zeros(n,kMax);
yhat= zeros(p,kMax);

timeKF= 0;
% simulation
try
    for k= 1:length(t)
        
        % set timer
        tic
        % execute kalman filter
        [yhat(:,k),xhat(:,k),X_,Z_]= myETKF(A,C,Q,R,N,y(:,k),Z_,X_);
        
        % timer
        timeKF= timeKF+toc;
        % limit for computation time
        if toc > timeLimit
            error('timeout');
        end
    end
catch ME
    if strcmpi(ME.message,'timeout')
        %fprintf('\tException caught timeout (%4.2e > %i) conventional KF!\n',toc,timeLimit);
    end
    %fprintf('\tVariables assigned NaN\n');
    timeKF= nan;
    yhat= nan(p,kMax);
    xhat= nan(n,kMax);
end
