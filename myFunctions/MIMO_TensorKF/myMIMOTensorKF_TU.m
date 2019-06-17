function [xhat_,P_]= myMIMOTensorKF_TU(A,Q,xhat,P,tol,varargin)
%% Tensor Kalman Filter Time Update
% myMIMOTensorKF_MU.m
% Date:             01.05.2019
% Authors:          Daniel Gedon, 4735226
% Description:      The function does a Kalman state estimation in tensor
%                   representation for the time update in MIMO case.
% Inputs:           A - state transition matrix in TT-format
%                   Q - Process noise covariance in TT-format
%                   xhat - state from MU in TT-format
%                   P - Covariance matrix from MU in TT-format
%                   tol - TT rounding tolerance
% Outputs:          xhat_ - propagated state in TT-format
%                   P_ - Covariance matrix in TT-format
%% Limitation on TN-rank of P and P_

PrankMax= Inf;

for i=1:2:length(varargin)-1
    if (~isempty(varargin{i+1}))
        switch lower(varargin{i})
            case 'prankmax'
                PrankMax=varargin{i+1};
            otherwise
                error('Unrecognized option: %s\n',varargin{i});
        end
    end
end

%% ---------- Time Update ----------
%%%% A priori estimate of the current state %%%%
% xhat(t+1|t) = A * xhat(t|t)
xhat_= contractab(xhat,A,[2,3]);
xhat_= roundTN(xhat_,tol);

%%%% A priori estimate of the state covariance matrix %%%%
% P(t+1|t) = A * P(t|t) * A' + Q
P_= contractab(P,A,[2,3]);
P_= roundTN(P_,tol);
P_= contractab(P_,A,[3,3]);
P_= roundTN(P_,tol);
P_= addTN(P_,Q);
P_= roundTN(P_,tol,PrankMax);