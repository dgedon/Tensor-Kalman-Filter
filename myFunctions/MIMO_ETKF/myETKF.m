function [yhat,xhat,X_,Z_]= myETKF(A,C,Q,R,N,y,Z_,X_)
%% Function Stochastic Ensemble Kalman Filter Test
% myEnKF_Det_Trans.m
% Date:             25.09.2018
% Authors:          Daniel Gedon, 4735226
% Description:      The function does a Kalman state estimation.
% Inputs:           A - state evolution matrix
%                   C - measurement matrix
%                   Q - process noise covariance matrix
%                   R - measurement noise covariance matrix
%                   N - number of ensembles
%                   y - measurement vector
% Outputs:          TODO
%% Get inputs

% state size
n= size(A,1);
% p= size(C,1);

%% Alogrithm

% ones
ONE= ones(N,1);

%% ---------- Measurement Update ----------

% Estimated observation ( yhat(k|k-1) = C*xhat(k|k-1) )
xhat_= X_*ONE/N;mean(X_,2);
y_= C * xhat_;

% Measurement residual error or innovation error ( y(k) 1' - Y(k|k-1) )
inov = y - y_;

% Innovation form ( Sinov =  inv(sqrt(R)) * inov )
Sinov= sqrt(R) \ inov;

% CZ form ( Scz = inv(Sqrt(R)) * C * Z(k|k-1) )
Scz = sqrt(R) \ C * Z_;

% eigenvalue decomposition
[Qk,Gamma] = eig( eye(N) + Scz'*Scz );

% A posterioi (updated) estimate of the current state ( xhat(k|k) = xhat(k|k-1) + ... )
xhat = xhat_ + Z_*Scz'* ( Sinov - Scz * Qk/Gamma*Qk' * Scz'*Sinov );

% Ensemble from ( X(k|k) = sqrt(N-1)* Z(k|k-1)*Qk*inv(sqrt(Gamma))*Qk' + xhat*1' )
X = Z_ * Qk*(sqrt(N-1)*inv(sqrt(Gamma)))*Qk' + xhat*ONE';

% A posteriori output estimate ( y(k|k) = C*x(k|k) )
Y = C * X;
yhat= Y*ONE/N; %mean(Y,2);

%% ---------- Time Update  ----------
% A priori state estimate ( X(k+1|k) = A*X(k|k) + G*(V(k)+U(K)) )
V= (mvnrnd(zeros(1,n),Q,N))';
X_= A*X + V;

% get mean of a priori state estimate
% xhat_ = X_*ONE/N;%mean(X_,2);

% A priori state anomalies ( Z(k|k-1) = 1/sqrt(N-1)* X(k|k-1)*(I - 1/N*1*1') )
Z_= 1/sqrt(N-1) * X_ * ( eye(N) - 1/N*ONE*(ONE') );