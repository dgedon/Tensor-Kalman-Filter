function [yhat,xhat,timeKF,Pnorm]= myKF_FullSim(A,C,Q,R,t,y,varargin)
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
if nargin == 7
    timeLimit= varargin{1};
end

% state and output dimension
n= size(A,1);
p= size(C,1);

% time
kMax= max(t);

%% Algorithm

% initialization and allocation
xhat_= zeros(n,1);
xhat= zeros(n,kMax);
yhat= zeros(p,kMax);
P0= eye(size(Q));
P_= P0; %0.1*

timeKF= 0;
% simulation
% timeLimit= 3;
try
    for k= 1:length(t)
        % prevent computation when timeout error occurs anyways
        if (timeLimit < 10^6) && (n>1300) % for timeLimit=inf and for c>40
            error('preventing computation where not necessary')
        end
        
        % set timer
        tic
        % execute kalman filter
        if 0
            [yhat(1:p/2,k),xhatint,P_int]= myKF_MU(C(1:p/2,:),R(1:p/2,1:p/2),y(1:p/2,k),P_,xhat_);
            [yhat(p/2+1:end,k),xhat(:,k),P]= myKF_MU(C(p/2+1:end,:),R(1:p/2,1:p/2),y(p/2+1:end,k),P_int,xhatint);
        else
            [yhat(:,k),xhat(:,k),P]= myKF_MU(C,R,y(:,k),P_,xhat_);
        end
        % TU
        [xhat_,P_]= myKF_TU(A,Q,P,xhat(:,k));
        % timer
        timeKF= timeKF+toc;
        % limit for computation time
        if toc > timeLimit
            error('timeout');
        end
        %Pnorm.P(k)= norm(P);
        %Pnorm.P_(k)= norm(P_);
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
