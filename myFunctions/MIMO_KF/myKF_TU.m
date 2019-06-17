function [xhat_,P_]= myKF_TU(A,Q,P,xhat)
%% Function Kalman Filter Time update
% myKF_TU.m
% Date:             22.11.2018
% Authors:          Daniel Gedon, 4735226
% Description:      The function does a Kalman state estimation.
% Inputs:           A - state evolution matrix
%                   Q - process noise covariance matrix
%                   P - covariance matrix
%                   xhat - estimated state
% Outputs:          xhat_ - estaimted state 
%                   P_ - covariance matrix
%% Algorithm

%---------- Time Update ----------
% A priori estimate of the current state ( x(t+1|t) = A*x(t|t) + G*u(t))
xhat_ = A * xhat;

% A priori estimate of the state covariance matrix ( P(t+1|t) = A*P(t|t)*A' + Q )
P_ = A*P*A' + Q;