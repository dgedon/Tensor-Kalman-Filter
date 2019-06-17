function [yhatOut,xhatOut,timeTKF]= ...
    myMIMOTensorKF_FullSim(A,C,Q,R,t,y,tol,varargin)
%% Function Tensor Kalman Filter with vector output
% myMatrixTensorKF_FullSim.m
% Date:             24.04.2019
% Authors:          Daniel Gedon, 4735226
% Description:      The function does a Kalman state estimation in tensor
%                   representation.
% Inputs:           A - state evolution matrix in TT format
%                   C - measurement matrix in TT format
%                   Q - process noise covariance tensor in TT format
%                   R - measurement noise covariance tensor 
%                   t - time vector
%                   y - measurement vector in TT format
%                   tol - tensor network conversion and rounding tolerance
% Outputs:          yhatOut - estimated observations (k|k)
%                   xhatOut - estimated state (k|k)
%                   timeTKF - Time for computation
%% Input assignment

% defaul values
PrankMax= Inf;
SrankMax= Inf;

for i=1:2:length(varargin)-1
    if (~isempty(varargin{i+1}))
        switch lower(varargin{i})
            case 'prankmax'
                PrankMax=varargin{i+1};
            case 'srankmax',
                SrankMax=varargin{i+1};
            otherwise
                error('Unrecognized option: %s\n',varargin{i});
        end
    end
end

%% Initialization of Variables

% get initialization values
% sizes
n= size(contract(A),1);
ni= A.n(:,2);
p= size(contract(C),1);
d= size(A.n,1);
% simulation time
tMax= numel(t);

%%%% Output variables
xhatOut= zeros(n,tMax);
yhatOut= zeros(p,tMax);

%%%% initialize Xhat, P and K in correct format
% initialize Xhat:
xhat_.n= [ones(d,1) ni ones(d,1)];
xhat_.core= cell(1,d);
xhat_.core{1}= zeros(1,ni(1));
for i= 2:d
    xhat_.core{1,i}= zeros(1,ni(i));
end

% initialize P: 
sigmaP= 1;
P_.n= [ones(d,1) ni ni ones(d,1)];
P_.core{1}= reshape(sigmaP*eye(ni(1)),P_.n(1,:));
for i= 2:d
    P_.core{1,i}= reshape(eye(ni(i)),P_.n(i,:));
end
P_= roundTN(P_,tol,PrankMax);

%% Time Simulation of Tensor KF 


% % symmetry, positive definiteness
% tempPpos= 1;
% tempPsymm= 1;
% tempP_pos= 1;
% tempP_symm= 1;

timeTKF= 0;
for k= t
    % convert y to TT
    if 0    % use TT_SVD
        ytt= myTTSVD(y(:,k),C.n(:,2));
        ytt= roundTN(ytt,tol);
    else % use TT_ALS
        if k== t(1)
            % first time step: need TT_SVD conversion
            % necessary orthogonalization
            ytt= myTTSVD(y(:,k),C.n(:,2));
            ytt= roundTN(ytt,tol);
            % orthogonalization
            for id= d:-1:2
                % move norm to first core
                ytt= TT_orth_i(ytt,id,-1);
            end
        else
            % otherwise use TT_ALS with previous yTT as init
            % only do one half sweep (sufficient)
            ytt= myVec2TT_ALS(y(:,k),ytt,'sweepBackward',0);
            % orthogonalization (norm currently in last core)
            for id= d:-1:2
                % move norm to first core
                ytt= TT_orth_i(ytt,id,-1);
            end
        end
    end  
    
    % set timer
    tic
    % Execute matrix tensor KF
    %%% Measurement update
    [xhat,yhat,P]= myMIMOTensorKF_MU(C,R,xhat_,ytt,P_,tol,'PrankMax',PrankMax,...
        'SrankMax',SrankMax);
    
    %%% Time update
    [xhat_,P_]= myMIMOTensorKF_TU(A,Q,xhat,P,tol,'Prankmax',PrankMax);

    % get timer
    timeTKF= timeTKF+ toc;
    
%     % check for symmetry and positive-definitness of P (MU), P_ (TU)
%     tempPpos= min(tempPpos,min(eig(contract(P))));
%     temp= contract(P);
%     tempPsymm= min(tempPsymm,norm(temp-temp'));
%     % P_
%     tempP_pos= min(tempP_pos,min(eig(contract(P_))));
%     temp= contract(P_);
%     tempP_symm= min(tempP_symm,norm(temp-temp'));
    
    %---------- Function Output Variables ----------
    % State xhat_out
    xhatOut(:,k)= contract(xhat);
    
    % Output yHat_Out
    yhatOut(:,k)= contract(yhat);
end


% % symmetry positive definitenss
% fprintf('\n\tmin P eig val\t\t %2.3e\n',tempPpos)
% fprintf('\tP symmetry up to\t %2.3e\n',tempPsymm)
% fprintf('\tmin P_ eig val\t\t %2.3e\n',tempP_pos)
% fprintf('\tP_ symmetry up to\t %2.3e\n',tempP_symm)