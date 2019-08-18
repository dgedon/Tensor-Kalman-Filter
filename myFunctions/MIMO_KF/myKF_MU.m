function [yhat,xhat,P,K]= myKF_MU(C,R,y,P_,xhat_)
%% Function Kalman Filter Measurement Update
% myKF_MU.m
% Date:             22.11.2018
% Authors:          Daniel Gedon, 4735226
% Description:      The function does the Kalman state measurement update.
% Inputs:           C - measurement matrix
%                   R - measurement noise covariance matrix
%                   y - plant output vector
%                   P_ - a priori state covariance matrix at time k
%                   xhat_ - a priori state vector at time k
% Outputs:          yhat - estimated observations
%                   xhat - estimated state vector 
%                   P - state covariance 
%% Algorithm

%---------- Measurement Update ----------
% Kalman filter coefficient ( K(t) = P(t|t-1) * C' * inv(C*P(t|t-1) * C' + R) )
% if det(C*P_*C'+R) < 1e-7
%     stopvar= 1;
% end
S= (C*P_*C' + R);
% [Us,Ss,Vs]= svds(S,11);
K = P_ * C'*pinv(S);%/S;% *pinv(Us*Ss*Vs');%*(Vs*diag(1./diag(Ss))*Us'); %/S;

% Estimated observation ( y(t|t-1) = C*x(t|t-1) )
yhat_ = C * xhat_;

% Measurement residual error or innovation error ( y(t) - y(t|t-1) )
inov = y - yhat_;

% A posterioi (updated) estimate of the current state ( x(t|t) = x(t|t-1) + K(t)*(y(t)-y(t|t-1)) )
xhat = xhat_ + K * inov;

% A posteriori (updated) state covariance matrix ( P(t|t) = (I - K(t)*C) * P(t|t-1) )
P = (eye(size(P_)) - K*C)*P_;
% P = (eye(size(A)) - K*C)*P_*(eye(size(A)) - K*C)' + K*R*K';

% A posteriori output estimate ( y(t|t) = C*x(t|t) )
yhat = C*xhat;

% stopvar= 1;