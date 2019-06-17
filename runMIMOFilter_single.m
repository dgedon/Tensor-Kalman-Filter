% Run of a random MIMO Filter for a single system
% runMIMOFilter_single.m
% Date:             08.05.2019
% Authors:          Daniel Gedon, 4735226
% Description:      Runs a random system simulation for ONE specific system
%                   size
%% General Initialization

clear
close all
clc

warning off

% set random number generator
rng(1234)

% include subfunction
addpath(genpath([pwd,'/myFunctions']))
addpath(genpath([pwd,'/TensorFunctions']))

% plotting options
showOutputEvolution= 0;
showStateEvolution= 1;
showFullComparison= 0;

%% System Initialization

d= 8;           % tensor order
n= 2;           % tensor mode
p= n^d;         % output size

% TT ranks
rA= 5;
rC= 1;
PrankMax= 1;

% TT rounding tolerance
tol= 1e-6;

% noise variance values
Qcov= 5e-1^2;
Rcov= 1e-6^2;   %1e-8^2;   % 1e-4^2 not working

% Simulation constraints
tMax=  50;   300;
t= 1:tMax;

% Initial state
x0= randn(n^d,1);
x0= x0-mean(x0);

%% Init system matrices (TT and Matrix)

if 0
    A_tt= genStableTTrank(n,d,rA);
else
    A_tt.n= [ones(d,1) n*ones(d,2) ones(d,1)];
    A_tt.core{1}= reshape(0.99*eye(n),A_tt.n(1,:));
    for i= 2:d
        A_tt.core{i}= reshape(eye(n),A_tt.n(i,:));
    end
end
A= contract(A_tt);

% this does not really work for the correct outputs
temp= 2*rand(n,n,d,rC); %
C_tt= getKroneckerTT(n,d,rC,temp);
C= contract(C_tt);
clear temp

%% Init Noise Matrix

% Noise characteristics
Q= Qcov*eye(n^d);
R= Rcov*eye(p);

% test for numeric feasibility
if  abs(det(C*Q*C'+R)) < 1e-10
    warning('det(...)= %2.3e',det(C*Q*C'+R))
end

%% Init noise TT

Q_tt.n= [ones(d,1) n*ones(d,1) n*ones(d,1) ones(d,1)];
Q_tt.core= cell(1,d);
Q_tt.core{1,1}= Qcov*eye(n);
for i= 2:d
    Q_tt.core{1,i}= eye(n);
end

R_tt.n= [ones(d,1) C_tt.n(:,2) C_tt.n(:,2) ones(d,1)];
R_tt.core= cell(1,d);
R_tt.core{1,1}= Rcov*eye(n);
for i= 2:d
    R_tt.core{1,i}= eye(R_tt.n(i,2:3));
end

%% Plant evolution

fprintf('\n\tState size:\t%i\n',n^d);
fprintf('\tOutput size:\t%i\n\n',p);

fprintf('Plant simulation started...')
[x,y]= myMIMOPlantSimulation(A,C,Q,R,t,x0);
fprintf('\t...done\n\n')

%% Tensor Kalman Filter

try
    % Tensor Kalman filter
    fprintf('Tensor KF started...');
    [yhat_ten,xhat_ten,time_ten]= myMIMOTensorKF_FullSim...
        (A_tt,C_tt,Q_tt,R_tt,t,y,tol,'PrankMax',PrankMax);
    
    fprintf('\t\t...Finished in t=%4.3f [sec/step]\n',time_ten/tMax);
catch ME
    % ended with exception
    fprintf('\t...Ended with exception error!\n');
    
    % set output variables of tensor KF to NaN
    xhat_ten= NaN(n^d,tMax);
    yhat_ten= NaN(p,tMax);
    time_ten= NaN(1);
end
%% Conventional KF

fprintf('Conventional KF started...');
try
    % run filter
    [yhat_mat,xhat_mat,time_mat]= myKF_FullSim(A,C,Q,R,t,y);
    fprintf('\t...Finished in t=%4.3f [sec/step]\n',time_mat/tMax);
catch ME
    % ended with exception
    fprintf('\t...Ended with exception error!\n');
    
    % set output variables of tensor KF to NaN
    xhat_mat= NaN(n^d,tMax);
    yhat_mat= NaN(p,tMax);
    time_mat= NaN(1);
end

%% Massioni KF

fprintf('Massioni KF started...');
try
    % run filter
    [yhat_mass,xhat_mass,time_mass]= massioniKF_FullSim(A,C,Q,R,t,y);
    fprintf('\t\t...Finished in t=%4.3f [sec/step]\n',time_mass/tMax);
catch ME
    % ended with exception
    fprintf('\t...Ended with exception error!\n');
    
    % set output variables of tensor KF to NaN
    xhat_mass= NaN(n^d,tMax);
    yhat_mass= NaN(p,tMax);
    time_mass= NaN(1);
end

%% Simulink Kalman Filter

if 0
    SimInY.time= t;
    SimInY.signals.values= y';
    SimInY.signals.dimensions= p;
    
    fprintf('Simulink KF started...');
    tic;
    sim('KalmanSimComparison')
    fprintf('\t\t...Finished in t=%4.3f [sec/step] (including open Simulink)\n',toc/tMax);
    
    xhat_sim= simXhat.signals.values(2:end,:)';
    yhat_sim= simYhat.signals.values(2:end,:)';
end

%% Norm of state

t0= 15; % zero time for accuracy evaluation (default: 1)

fprintf('\n\nRelative State Norm:\n');

% tensor filter
temp= xhat_ten(:,t0:end)-x(:,t0:end);
temp= norm(temp(:)).^2/norm(vec(x(:,t0:end))).^2;
fprintf('\tTensor KF:\t\t%4.3e\n',temp);
% conventional filter
temp= xhat_mat(:,t0:end)-x(:,t0:end);
temp= norm(temp(:)).^2/norm(vec(x(:,t0:end))).^2;
fprintf('\tConventional KF:\t%4.3e\n',temp);
% Massioni filter
temp= xhat_mass(:,t0:end)-x(:,t0:end);
temp= norm(temp(:)).^2/norm(vec(x(:,t0:end))).^2;
fprintf('\tMassioni KF:\t\t%4.3e\n',temp);

%% Norm of output

fprintf('\n\nRelative Output Norm:\n');

% tensor filter
temp= y(:,t0:end)-yhat_ten(:,t0:end);
temp= norm(temp(:)).^2/norm(vec(y(:,t0:end))).^2;
fprintf('\tTensor KF:\t\t%4.3e\n',temp);
% conventional filter
temp= yhat_mat(:,t0:end)-y(:,t0:end);
temp= norm(temp(:)).^2/norm(vec(y(:,t0:end))).^2;
fprintf('\tConventional KF:\t%4.3e\n',temp);
% Massioni filter
temp= yhat_mass(:,t0:end)-y(:,t0:end);
temp= norm(temp(:)).^2/norm(vec(y(:,t0:end))).^2;
fprintf('\tMassioni KF:\t\t%4.3e\n',temp);

%% Graphical Output

% plot one random output
nrOut= randi(p);

%%%% System output || system state
if showOutputEvolution
    
    figure('position',[100 100 640 480]);
    hold on
    % real system
    plot(t(t0:end),y(nrOut,t0:end),'kx-','LineWidth',2)
    % KF
    plot(t(t0:end),yhat_mat(nrOut,t0:end),'b--','LineWidth',2)
    % TKF
    plot(t(t0:end),yhat_ten(nrOut,t0:end),'r:','LineWidth',2)
    % Mass
    plot(t(t0:end),shat_mass(nrOut,t0:end),'g-.','LineWidth',2);
    % settings
    grid on
    xlabel('time')
    ylabel('system output')
    title(['System output evolution of output ',num2str(nrOut)])
    legend('Noisy system',['KF est.; time=',num2str(time_mat/tMax),'sec/step'],...
        ['Tensor KF est.; time=',num2str(time_ten/tMax),'sec/step'],...
        ['Mass. KF est.; time=',num2str(time_mass/tMax),'sec/step'],'Location','best');
end
%%

if showStateEvolution
    figure('position',[100 100 640 480]);
    hold on
    % real system
    plot(t(t0:end),x(nrOut,t0:end),'kx-','LineWidth',2)
    hold on;
    % KF
    plot(t(t0:end),xhat_mat(nrOut,t0:end),'b--','LineWidth',2)
    % TKF
    plot(t(t0:end),xhat_ten(nrOut,t0:end),'r:','LineWidth',2)
    % Mass
    plot(t(t0:end),xhat_mass(nrOut,t0:end),'g-.','LineWidth',2)
    % simulink
    % plot(t(t0:end),xhat_sim(nrOut,t0:end),'m','Marker','o');
    % settings
    grid on
    xlabel('time')
    ylabel('system state')
    title(['System state evolution of state ',num2str(nrOut)])
    legend('Noisy system',['KF est.; time=',num2str(time_mat/tMax),'sec/step'],...
        ['Tensor KF est.; time=',num2str(time_ten/tMax),'sec/step'],...
        ['Mass. KF est.; time=',num2str(time_mass/tMax),'sec/step'],'Location','best');
end
%%
if showFullComparison
    fig= figure('position',[100 100 640 1.3*480]);
    subplot(3,1,1)
    imagesc(x(:,t0:end)-xhat_mat(:,t0:end)), colorbar, title('x-xhat_{mat}'), xlabel('time'), ylabel('n')
    subplot(3,1,2)
    imagesc(x(:,t0:end)-xhat_ten(:,t0:end)), colorbar, title('x-xhat_{ten}'), xlabel('time'), ylabel('n')
    subplot(3,1,3)
    imagesc(x(:,t0:end)-xhat_mass(:,t0:end)), colorbar, title('x-xhat_{mass}'), xlabel('time'), ylabel('n')
end