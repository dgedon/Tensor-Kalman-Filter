function [xhat,yhat,P,K]= myMIMOTensorKF_MU(C,R,xhat_,y,P_,tol,varargin)
%% Tensor Kalman Filter Measurement Update
% myMIMOTensorKF_MU.m
% Date:             01.05.2019
% Authors:          Daniel Gedon, 4735226
% Description:      The function does a Kalman state estimation in tensor
%                   representation for the measurement update in MIMO case.
% Inputs:           C - Measurement matrix in TT-format
%                   R - Measurement noise covariance in TT-format
%                   xhat_ - state from TU in TT-format
%                   y - measurements in TT-format
%                   P_ - Covariance matrix from TU in TT-format
%                   tol - TT rounding tolerance
% Outputs:          xhat - estimated state in TT-format
%                   yhat - estimated output in TT-format
%                   P - Covariance matrix in TT-format
%                   K - Kalman gain in TT-format
%% Limitation on TN-rank of P and P_

PrankMax= Inf;
SrankMax= Inf;

for i=1:2:length(varargin)-1
    if (~isempty(varargin{i+1}))
        switch lower(varargin{i})
            case 'prankmax'
                PrankMax=varargin{i+1};
            case 'srankmax'
                SrankMax=varargin{i+1};
            otherwise
                error('Unrecognized option: %s\n',varargin{i});
        end
    end
end

%% get input

d= size(C.n,1);

%% ---------- Measurement Update ----------

%%%% Riccati part %%%%
% S = C * P(t|t-1) * C' + R
S= contractab(P_,C,[2,3]);
S= roundTN(S,tol);
S= contractab(S,C,[3,3]);
S= roundTN(S,tol);
S= addTN(S,R);
S= roundTN(S,tol);

if ~isequal([S.n(:,1) S.n(:,end)], ones(d,2))
    S= roundTN(S,tol,SrankMax);
    fprintf('\tTrunction of S to ttr(S)=%i done\n',SrankMax);
end

%%%% Riccati part inverse %%%%
% if all TN ranks are one then the cores are matrices and inversion
% is simply the inversion of the matrices
if 1
    Sinv.n= S.n;
    for i= 1:d
        Sinv.core{i}= inv(squeeze(S.core{i})); % pinv
        Sinv.core{i}= reshape(Sinv.core{i},Sinv.n(i,:));
    end
else
    Sinv= myTTSVD(pinv(contract(S)),S.n(:,2:3));
end
    
%%%% Kalman filter gain %%%%
% K = P(t|t-1) * C' * Sinv
K= contractab(P_,C,[3,3]);
K= roundTN(K,tol);
K= contractab(Sinv,K,[2,3]);
K= roundTN(K,tol);

%%%% Measurement residual error or innovation error %%%%
% v = y(t) - C * xhat(t|t-1)
yhat_= contractab(xhat_,C,[2,3]);
yhat_= roundTN(yhat_,tol);
yhat_min.n= yhat_.n;    % minus
yhat_min.core= yhat_.core;
yhat_min.core{d}= -yhat_.core{d};
v= addTN(y,yhat_min);

%%%% A posterioi (updated) estimate of the current state %%%%
% xhat(t|t) = xhat(t|t-1) + K * v
xhat= contractab(v,K,[2,3]);
xhat= roundTN(xhat,tol);
xhat= addTN(xhat_,xhat);

%%%% A posteriori (updated) state covariance matrix %%%%
% P(t|t) = P(t|t-1) + K * C * P(t|t-1)
Ktilde= contractab(C,K,[2,3]);
Ktilde= roundTN(Ktilde,tol);
Ktilde= contractab(P_,Ktilde,[2,3]);
Ktilde= roundTN(Ktilde,tol);
% negative of Ktilde for subtraction
Ktilde_min.n= Ktilde.n;
Ktilde_min.core= Ktilde.core;
Ktilde_min.core{d}= -Ktilde_min.core{d};
% addition with P
P= addTN(P_,Ktilde_min);
P= roundTN(P,tol,PrankMax);

%%%% A posteriori output estimate %%%%
% y(t|t) = C * xhat(t|t)
yhat= contractab(xhat,C,[2,3]);
yhat= roundTN(yhat,tol);