function [yhat,xhat,timeKF]= massioniKF_FullSim(A,C,Q,R,t,y,varargin)
%% Function Kalman Filter without input
% massioniKF_FullSim.m
% Date:             15.05.2019
% Authors:          Daniel Gedon, 4735226
% Description:      Massionis KF full simulation and initialization
% Inputs:           A - state transition matrix
%                   C - measurement matrix
%                   Q - process noice covariance
%                   R - measurement noise covariance
%                   t - time vector
%                   y - measurement vector for all time
% Outputs:          yhat - estimated measurement vector
%                   xhat - estiamted state vector
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

timeKF= 0;
% simulation
% % timeLimit= 3;
try
    for k= 1:length(t)
        % prevent computation when timeout error occurs anyways
        if (timeLimit < 10^6) && (n>1500) % for timeLimit=inf and for c>40
            error('preventing computation where not necessary')
        end
        
        % set timer
        tic
        %%% Precompute Kalman gain since no recursion necessary
        % assumption A= a*I
        a= A(1);
        % assumption R= eq*I
        sigma= R(1);
        
        % only observable part of C
        [~,~,S]=svd(C);
        Cs= C*S;
        Cs=Cs(:,1:end-2);
        CC1= (Cs'*Cs)^-1;
        
        %%% Kalman gain
        
        % partition Q
        Qvs= S'*Q*S;
        Qvs1= Qvs(1:end-2,1:end-2);
        Qvs12= Qvs(1:end-2,end-1:end);
        
        % covariance first order approximation
        P11= Qvs1+a^2*sigma*CC1;
        P12= Qvs12+sigma*a^2*CC1/(Qvs1)*Qvs12; % originially 1/Qvs1 but det(Qvs1)=~0
        
        K= S*[P11;P12']*Cs'/(Cs*P11*Cs'+sigma*eye(size(Cs,1))); % originally not pinv but ^-1
        % K= K(1:end-2,:);
        
        
        %%%% execute kalman filter
        [yhat(:,k),xhat(:,k),xhat_]= massioniKF(A,C,K,y(:,k),xhat_);
        timeKF= timeKF+toc;
        % limit for computation time
        if toc > timeLimit
            error('timeout');
        end
    end
catch ME
    if strcmpi(ME.message,'timeout')
        %fprintf('\tException caught timeout (%4.2e > %i) Massioni KF!\n',toc,timeLimit);
    end
    %fprintf('\tVariables assigned NaN\n');
    timeKF= nan;
    yhat= nan(p,kMax);
    xhat= nan(n,kMax);
end
