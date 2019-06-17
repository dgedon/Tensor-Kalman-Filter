function [x,y]= myMIMOPlantSimulation(A,C,Q,R,t,x0)
%% Function Plant Evoluation for Extended State
% myMIMOPlantsimulation.m
% Date:             20.11.2018
% Authors:          Daniel Gedon, 4735226
% Description:      The function evaluates the LTI plant with extended
%                   system state to matrix output of size pxl
% Inputs:           A - state evolution matrix
%                   C - measurement matrix
%                   Q - process noise covariance matrix
%                   R - measurement noise covariance matrix
%                   t - time vector
%                   x0 - initial state vector
% Outputs:          x - state vector
%                   y - output vector
%% Check Input

% state and output dimension
n= size(A,1);
p= size(C,1);

% time
tMax= max(t);

%% Noise

% allocation
w= zeros(n,tMax);
v= zeros(p,tMax);
% store noise values
w= (mvnrnd(zeros(n,1),Q,tMax))';%*sqrt(Q)
v= (mvnrnd(zeros(p,1),R,tMax))';%*sqrt(R)

%% plant simulation

% allocation
x= zeros(n,tMax);
y= zeros(p,tMax);
% init
x(:,1)= x0;
% simulation loop
for i= 1:tMax
    x(:,i+1)= A*x(:,i)+w(:,i);
    y(:,i)= C*x(:,i)+v(:,i);
end
x= x(:,1:tMax);