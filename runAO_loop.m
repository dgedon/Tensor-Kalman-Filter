% AO filter Loop
% runAO_loop.m
% Date:             13.05.2019
% Authors:          Daniel Gedon, 4735226
% Description:      Runs an AO simulation for several specific system size
%                   in a monte carlo way. The average time and relative
%                   accuracy is plotted in the end.
%% General Initialization
clear
close all
clc

%start= tic;

warning on

% set random number generator
rng(45378)

% include subfunction
addpath(genpath([pwd,'/myFunctions']))
addpath(genpath([pwd,'/TensorFunctions']))

% plotting options
showOutputEvolution= 0;
showStateEvolution= 1;
showFullComparison= 1;

%% System Initialization

% TT rounding tolerance
tol= 1e-4;

% simulation time
t0= 1; % zero time for accuracy evaluation (default:1)
tMax= 150;
t= 1:tMax;
timeLimit= 10; % limit for computation. If exceeded just nan
fprintf('Simulation time: %i [steps]\n',tMax);
fprintf('Simulation limit: %i [sec]\n\n',timeLimit);

% number of lenslets (loop)
numLenslets= 4:2:70; %70
% Monte Carlo runs
N= 1;4;   10;

% AO settings
sigma= 3.33e-3; % make it a bit lower than this!! 1e-4;    45e-9;  1e-5;  1e-8;
iter= 5e4;
samplingLens= 2;

% Select System Identification (1: identified, 0: A=Id)
selectSysID= 0;

% TT-rank truncations
PrankMax= 1;
SrankMax= 1;
rQ= 10;
rA= 10;
% for identified system ttr(Q)=1, ttr(A)=c+1
% for A=1 ttr(Q)=c^2, ttr(A)= 1;

%% Loop

% allocation
normPhihat_mat= zeros(length(numLenslets),N);
normPhihat_mass= zeros(length(numLenslets),N);
normPhihat_ens= zeros(length(numLenslets),N);
normPhihat_ten= zeros(length(numLenslets),N);

normShat_mat= zeros(length(numLenslets),N);
normShat_mass= zeros(length(numLenslets),N);
normShat_ens= zeros(length(numLenslets),N);
normShat_ten= zeros(length(numLenslets),N);

time_mat= zeros(length(numLenslets),N);
time_mass= zeros(length(numLenslets),N);
time_ens= zeros(length(numLenslets),N);
time_ten= zeros(length(numLenslets),N);

% workspace output
fprintf('| N,  d,  c \t| sim\t| KF mat\t| KF mass\t| ETKF \t| KF ten |\n')
fprintf('-------------------------------------------------------\n')

% loop over Monte Carlo runs
for iN= 1:N
    
    fprintf('\n')
    
    % loop over lenselets
    for CNT= 1:length(numLenslets)
        % workspace output
        fprintf('| %i,  ',iN);
        
        cLenslets= numLenslets(CNT);
        diameter= cLenslets/samplingLens;
        
        % workspace output
        fprintf('%i, %i \t',diameter, cLenslets)
        
        % simple system identification
        if selectSysID
            tempSelID= 'Identified';
        else
            tempSelID= 'AIdentity';
        end
        try
            temp= ['00_SystemID/',tempSelID,'/AOSystemID_d',num2str(diameter),...
                '_iter',num2str(iter),'.mat'];
            load(temp);
        catch ME
            error('System identification setting not yet stored!');
        end
        A= data.A;
        G= data.C;
        Q= data.Q;
        n= size(A,1);
        p= size(G,1);
        
        % adjustment of tMax for removal of initial condition
        % tMax= tMax+2*n;
        
        % initial state
        phi0= randn(n,1);
        phi0= phi0-mean(phi0);
        % only because of unobservable piston/waffle mode one obtains a bias by the
        % mean of the initial condition. Hence remove mean from initial condition.
        
        %% Matrix system
        
        %%% A matrix
        A;
        
        %%% Rearrange C such that s=[sx,sy] instead of [sx1,sy1,sx2,sy2,...]
        % separate sx and sy
        C= zeros(size(G));
        for i= 1:p
            if mod(i,2)==1
                C((i+1)/2,:)= G(i,:);
            else
                C(p/2+i/2,:)= G(i,:);
            end
        end
        
        %%% Covariances
        Q;
        R= sigma^2*eye(p);
        
        %% Open loop simulations
        
        [phi,s]= myMIMOPlantSimulation(A,C,Q,R,t,phi0);
        
        % workspace output
        fprintf('| x\t');
        
        if ~selectSysID% only for A=a*I
            [~,~,V1]= svd(C);
            V1= V1(:,1:end-2);
            
            phi= V1'*phi;
        end
        
        %% TT Representation
        
        % number of TT-cores
        d= 2;
        
        % A matrix
        Att= myTTSVD(A,[cLenslets cLenslets;cLenslets cLenslets],'maxrank',rA);
        %Att= roundTN(Att,tol,rA);
        
        % C1 matrix
        % Inits
        G1= C(1:cLenslets-1,1:cLenslets);
        E1= eye(cLenslets-1,cLenslets);
        temp= zeros(cLenslets-1,cLenslets);
        temp(:,2:end)= eye(cLenslets-1);
        E1= E1-temp;
        % TT
        Ctt1.n= [ones(d,1) [cLenslets-1 cLenslets;cLenslets-1 cLenslets] ones(d,1)];
        Ctt1.core{1}= reshape(G1,Ctt1.n(1,:));
        Ctt1.core{2}= reshape(E1,Ctt1.n(2,:));
        
        %%%% C2 matrix
        % Inits
        G2= C(p/2+1:p/2+cLenslets-1,1:cLenslets);
        E2= eye(cLenslets-1,cLenslets);
        temp= zeros(cLenslets-1,cLenslets);
        temp(:,2:end)= eye(cLenslets-1);
        E2= E2+temp;
        % TT
        Ctt2.n= [ones(d,1) [cLenslets-1 cLenslets;cLenslets-1 cLenslets] ones(d,1)];
        Ctt2.core{1}= reshape(G2,Ctt2.n(1,:));
        Ctt2.core{2}= reshape(E2,Ctt2.n(2,:));
        
        % Q matrix
        Qtt= myTTSVD(Q,[cLenslets cLenslets;cLenslets cLenslets],'maxrank',rQ);
        % Qtt= roundTN(Qtt,tol,rQ);
        
        % R matrix
        Rtt.n= [ones(d,1) [cLenslets-1 cLenslets-1;cLenslets-1 cLenslets-1] ones(d,1)];
        for i= 1:d
            Rtt.core{i}= sigma^(2/d)*reshape(eye(Rtt.n(i,2:3)),Rtt.n(i,:));
        end
        
        %% Matrix Kalman filter
        
        %%% Matrix system
        % Kalman filter run
        [shat_mat,phihat_mat,temp]= myKF_FullSim(A,C,Q,R,t,s,timeLimit);
        time_mat(CNT,iN)= temp/tMax;
        
        % workspace output
        if isnan(temp)
            fprintf('| o\t\t');
        else
            fprintf('| x\t\t');
        end
        
        if ~selectSysID% only for A=a*I
            phihat_mat= V1'*phihat_mat;
        end
        
        % norm of state
        temp= phi(:,t0:end)-phihat_mat(:,t0:end);
        normPhihat_mat(CNT,iN)= norm(vec(temp)).^2/norm(vec(phi(:,t0:end))).^2;
        
        % norm of output
        temp= s(:,t0:end)-shat_mat(:,t0:end);
        normShat_mat(CNT,iN) = norm(vec(temp)).^2/norm(vec(s(:,t0:end))).^2;
        
        %% Massioni Kalman filter
        
        %%% Matrix system
        % Kalman filter run
        [shat_mass,phihat_mass,temp]= massioniKF_FullSim(A,C,Q,R,t,s,timeLimit);
        time_mass(CNT,iN)= temp/tMax;
        
        % workspace output
        if isnan(temp)
            fprintf('| o\t\t');
        else
            fprintf('| x\t\t');
        end
        
        if ~selectSysID% only for A=a*I
            phihat_mass= V1'*phihat_mass;
        end
        
        % norm of state
        temp= phi(:,t0:end)-phihat_mass(:,t0:end);
        normPhihat_mass(CNT,iN)= norm(vec(temp)).^2/norm(vec(phi(:,t0:end))).^2;
        
        % norm of output
        temp= s(:,t0:end)-shat_mass(:,t0:end);
        normShat_mass(CNT,iN) = norm(vec(temp)).^2/norm(vec(s(:,t0:end))).^2;
        
        %% Ensemble Transformation Kalman
        
        Ens= ceil(p/2); %n , p
        
        %%% Matrix system
        % ETKF filter run
        [shat_ens,phihat_ens,temp]= myETKF_FullSim(A,C,Q,R,Ens,t,s,timeLimit);
        time_ens(CNT,iN)= temp/tMax;
        
        % workspace output
        if isnan(temp)
            fprintf('| o\t');
        else
            fprintf('| x\t');
        end
        
        if ~selectSysID% only for A=a*I
            phihat_ens= V1'*phihat_ens;
        end
        
        % norm of state
        temp= phi(:,t0:end)-phihat_ens(:,t0:end);
        normPhihat_ens(CNT,iN)= norm(vec(temp)).^2/norm(vec(phi(:,t0:end))).^2;
        
        % norm of output
        temp= s(:,t0:end)-shat_ens(:,t0:end);
        normShat_ens(CNT,iN) = norm(vec(temp)).^2/norm(vec(s(:,t0:end))).^2;
        
        %% Tensor Kalman filter
        
        cdata.G1= G1;
        cdata.E1= E1;
        cdata.G2= G2;
        cdata.E2= E2;
        
        Rtt2.n= ones(2,4);
        for i= 1:d
            Rtt2.core{i}= sigma^(2/d)*reshape(eye(Rtt2.n(i,2:3)),Rtt2.n(i,:));
        end
        
        % Tensor Kalman filter
        [shat_ten,phihat_ten,temp]= myMIMOTensorKF_FullSim_AO(Att,Ctt1,Ctt2,...
            Qtt,Rtt2,t,s,tol,'PrankMax',PrankMax,'SrankMax',SrankMax,...
            'cdata',cdata);
        time_ten(CNT,iN)= temp/tMax;
        % workspace output
        fprintf('| x      |\n');
        
        if ~selectSysID% only for A=a*I
            phihat_ten= V1'*phihat_ten;
        end
        
        % norm of state
        temp= phi(:,t0:end)-phihat_ten(:,t0:end);
        normPhihat_ten(CNT,iN) = norm(vec(temp)).^2/norm(vec(phi(:,t0:end))).^2;
        
        % norm of output
        temp= s(:,t0:end)-shat_ten(:,t0:end);
        normShat_ten(CNT,iN) = norm(vec(temp)).^2/norm(vec(s(:,t0:end))).^2;
    end
end

%% Graphical output

%%% times
fig1= figure('position',[100 100 640 480]);
%subplot(2,1,1)
hold on
% matrix
semilogy(numLenslets/samplingLens,mean(time_mat,2),'k-','LineWidth',2);
% tensor
semilogy(numLenslets/samplingLens,mean(time_ten,2),'r--','LineWidth',2);
% massioni
semilogy(numLenslets/samplingLens,mean(time_mass,2),'g-.','LineWidth',2);
% ETKF
semilogy(numLenslets/samplingLens,mean(time_ens,2),'m:','LineWidth',2);
% setttings
grid on
xlabel('Diameter [m]')
ylabel('time [sec/step]')
title('Computation time per step')
legend('Matrix filter','Tensor filter','Massioni filter','ET filter','location','best');
set(gca,'XScale','lin','YScale','log')
% ylim([5e-5 3e0])

%%% state norms
fig2= figure('position',[100 100 640 480]);
%subplot(2,1,2)
hold on
% matrix
semilogy(numLenslets/samplingLens,mean(normPhihat_mat,2),'k-','LineWidth',2);
% tensor
semilogy(numLenslets/samplingLens,mean(normPhihat_ten,2),'r--','LineWidth',2);
% massioni
semilogy(numLenslets/samplingLens,mean(normPhihat_mass,2),'g-.','LineWidth',2);
% ETKF
semilogy(numLenslets/samplingLens,mean(normPhihat_ens,2),'m:','LineWidth',2);
% setttings
grid on
xlabel('Diameter [m]')
ylabel('rel. 2-Norm')
title('Relative 2-Norm of est. wavefront')
% legend('Matrix filter','Tensor filter','location','best');
set(gca,'XScale','lin','YScale','log')
% ylim([1e-7 5e-1])

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

% toc(start)