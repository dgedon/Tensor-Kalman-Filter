% example_randMIMO_loop.m
% Date:             18.08.2019
% Authors:          Daniel Gedon, 4735226
% Description:      Runs a random simulation for several specific system
%                   sizes in a monte carlo way. The average time and 
%                   relative accuracy is plotted in the end.
%% General Initialization

clear
close all
clc

% set random number generator
rng(12345)

% include subfunction
addpath(genpath([pwd,'/myFunctions']))
addpath(genpath([pwd,'/TensorFunctions']))

% plotting options
showOutputEvolution= 0;
showStateEvolution= 1;
showFullComparison= 1;

%% System Initialization

d= 2:9;                 % tensor order
n= 2*ones(size(d));     % tensor mode

%  TT ranks
rA= 5;          % rank of A
rC= 1;          % rankd of C
PrankMax= 1;    % max rank of P
SrankMax= 1;    % max rank of S

% TT rounding tolerance
tol= 1e-4;

% simulation time
tMax= 150;
t= 1:tMax;
timeLimit= 10; % limit for computation. If exceeded just nan

% Monte Carlo runs
N= 2;

% covariance values
Qcov= 0.1^2;
Rcov= 1e-6^2;

%% Simulation loop

% allocation
normX_ten= zeros(length(n),N);
normX_mat= zeros(length(n),N);
normX_mass= zeros(length(n),N);
normY_ten= zeros(length(n),N);
normY_mat= zeros(length(n),N);
normY_mass= zeros(length(n),N);
time_mat=  zeros(length(n),N);
time_ten=  zeros(length(n),N);
time_mass=  zeros(length(n),N);

% workspace output
fprintf('|  N \t| n, d\t| sim\t| KF mat\t| KF mass\t| KF ten |\n')
fprintf('------------------------------------------------------------------\n')

% Monte Carlo loop
for iN= 1:N    
    
    fprintf('\n');
    
    % system size loop
    for in= 1:length(n)
        % workspace output
        fprintf('| %i\t',iN);
        
        nloop= n(in);
        dloop= d(in);
        
        % workspace output
        fprintf('| %i,%i\t',nloop,dloop);
        
        % inital state
        x0= randn(nloop^dloop,1);
        x0= x0-mean(x0);
        
        %% TT system
        
        % A-TT
        A_tt= genStableTTrank(nloop,dloop,rA);
        
        % C-TT
        temp= 2*rand(nloop,nloop,dloop,rC); %
        C_tt= getKroneckerTT(nloop,dloop,rC,temp);
        
        % Process noise Covariance
        Q_tt.n= [ones(dloop,1) nloop*ones(dloop,2) ones(dloop,1)];
        for i= 1:dloop
            Q_tt.core{i}= reshape(eye(Q_tt.n(i,2:3)),Q_tt.n(i,:));
        end
        Q_tt.core{1}= Qcov*Q_tt.core{1};
        
        % Measurement noise Covariance
        R_tt.n= [ones(dloop,1) nloop*ones(dloop,2) ones(dloop,1)];
        for i= 1:dloop
            R_tt.core{i}= reshape(eye(R_tt.n(i,2:3)),R_tt.n(i,:));
        end
        R_tt.core{1}= Rcov*R_tt.core{1};
        
        %% Matrix system
        
        % A-matrix
        A= contract(A_tt);
        
        % C-matrix
        C= contract(C_tt);
        
        % Process noise covariance
        Q= Qcov*eye(nloop^dloop);
        
        % Measurement noise covariance
        R= Rcov*eye(nloop^dloop);
        
        %% Open loop simulation
        
        [x,y]= myMIMOPlantSimulation(A,C,Q,R,t,x0);
        
        % workspace output
        fprintf('| x\t');
            
        %% Matrix Kalman filter
        
        % Klaman filter run
        [yhat_mat,xhat_mat,time_mat(in,iN)]= myKF_FullSim(A,C,Q,R,t,y,timeLimit);
        
        % workspace output
        if isnan(time_mat(in,iN))
            fprintf('| o\t\t');
        else
            fprintf('| x\t\t');
        end
        
        % norm of state
        temp= x-xhat_mat;
        normX_mat(in,iN)= norm(temp(:)).^2/norm(x(:)).^2;
        
        % norm of output
        temp= y-yhat_mat;
        normY_mat(in,iN) = norm(temp(:)).^2/norm(y(:)).^2;
        
        %% Massioni Kalman filter
        
        % Klaman filter run
        [yhat_mass,xhat_mass,time_mass(in,iN)]= massioniKF_FullSim(A,C,Q,R,t,y,timeLimit);
        
        % workspace output
        if isnan(time_mass(in,iN))
            fprintf('| o\t\t');
        else
            fprintf('| x\t\t');
        end
        
        % norm of state
        temp= x-xhat_mass;
        normX_mass(in,iN)= norm(temp(:)).^2/norm(x(:)).^2;
        
        % norm of output
        temp= y-yhat_mass;
        normY_mass(in,iN) = norm(temp(:)).^2/norm(y(:)).^2;
        
        %% Tensor Kalman filter
        
        % Tensor Kalman filter
        [yhat_ten,xhat_ten,time_ten(in,iN)]= myMIMOTensorKF_FullSim(A_tt,...
            C_tt,Q_tt,R_tt,t,y,tol,'Prankmax',PrankMax,'Srankmax',SrankMax);
        
        % workspace output
        fprintf('| x      |\n');
        
        % norm of state
        temp= x-xhat_ten;
        normX_ten(in,iN) = norm(temp(:)).^2/norm(x(:)).^2;
        
        % norm of output
        temp= y-yhat_ten;
        normY_ten(in,iN) = norm(temp(:)).^2/norm(y(:)).^2;
    end
    
end

% postprocessing: time per time step
time_mat= time_mat/tMax;
time_mass= time_mass/tMax;
time_ten= time_ten/tMax;

%% Graphical output: Computation time

%%% times
fig1= figure('position',[100 100 640 480]);
%subplot(2,1,1)
hold on
% matrix
semilogy(n.^d,mean(time_mat,2),'k-','LineWidth',2);
% tensor
semilogy(n.^d,mean(time_ten,2),'r--','LineWidth',2);
% massioni
semilogy(n.^d,mean(time_mass,2),'g-.','LineWidth',2);
% setttings
grid on
xlabel('states, outputs')
ylabel('time [sec/step]')
title('Computation time per step')
legend('Matrix filter','Tensor filter','Massioni filter','location','best');
set(gca,'XScale','log','YScale','log')

%% Graphical output: Norms (Accuracy)

fig2= figure('position',[100 100 640 1.5*480]);
%%% state norms
subplot(2,1,1)
hold on
% matrix
semilogy(n.^d,mean(normX_mat,2),'k-','LineWidth',2);
% tensor
semilogy(n.^d,mean(normX_ten,2),'r--','LineWidth',2);
% Massioni
semilogy(n.^d,mean(normX_mass,2),'g-.','LineWidth',2);
% setttings
grid on
xlabel('states, outputs')
ylabel('rel. 2-Norm')
title('Relative 2-Norm of est. state')
legend('Matrix filter','Tensor filter','Massioni filter','location','best');
set(gca,'XScale','log','YScale','log')

%%% output norms
subplot(2,1,2)
hold on
% matrix
semilogy(n.^d,mean(normY_mat,2),'k-','LineWidth',2);
% tensor
semilogy(n.^d,mean(normY_ten,2),'r--','LineWidth',2);
% tensor
semilogy(n.^d,mean(normY_mass,2),'g-.','LineWidth',2);
% setttings
grid on
xlabel('states, outputs')
ylabel('rel. 2-Norm')
title('Relative 2-Norm of est. output')
% legend('Matrix filter','Tensor filter','location','best');
set(gca,'XScale','log','YScale','log')

%% saving

if 0
% figures
saveas(fig1,['data/',tempSelID,'/LoopTime_',tempSelID,'_',date])
saveas(fig1,['data/',tempSelID,'/LoopTime_',tempSelID,'_',date],'epsc')

saveas(fig2,['data/',tempSelID,'/LoopNormX_',tempSelID,'_',date])
saveas(fig2,['data/',tempSelID,'/LoopNormX_',tempSelID,'_',date],'epsc')

% data
temp= ['AOfullLoop_',tempSelID,'_',date];
save(temp)
end