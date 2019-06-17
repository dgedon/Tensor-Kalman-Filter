function [yhat,xhat,xhat_]= massioniKF(A,C,K,y,xhat_)
%% Function Kalman Filter without input
% massioniKF.m
% Date:             15.05.2019
% Authors:          Daniel Gedon, 4735226
% Description:      Massionis KF according to 2015 IEEE transaction paper
% Inputs:           A - state transition matrix
%                   C - measurement matrix
%                   Q - process noice covariance
%                   R - measurement noise covariance
%                   y - measurement vector for one time step
%                   xhat_ - previously estimated state vector
% Otuputs           yhat - MU estimated measurement vector
%                   xhat - MU estimated state vector
%                   xhat_ - TU estiamted state vector
%% inputs

% assumption A= a*I
a= A(1);
% % assumption R= eq*I
% sigma= R(1);

%% Kalman gain
% 
% % partition Q
% Qvs= S'*Q*S;
% Qvs1= Qvs(1:end-2,1:end-2);
% Qvs12= Qvs(1:end-2,end-1:end);
% 
% % covariance first order approximation
% P11= Qvs1+a^2*sigma*CC1; 
% P12= Qvs12+sigma*a^2*CC1/(Qvs1)*Qvs12; % originially 1/Qvs1 but det(Qvs1)=~0
% 
% K= S*[P11;P12']*Cs'/(Cs*P11*Cs'+sigma*eye(size(Cs,1))); % originally not pinv but ^-1
% % K= K(1:end-2,:);

%% measurement update

xhat= xhat_ + K * (y - C*xhat_);
yhat= C*xhat;

%% Time update

xhat_= a * xhat;