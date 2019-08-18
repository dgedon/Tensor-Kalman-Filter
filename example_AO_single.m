% example_AO_single.m
% Date:             18.08.2019
% Authors:          Daniel Gedon, 4735226
% Description:      Runs an AO simulation for ONE specific system size
%% General Initialization

clear
close all
clc

% set random number generator
rng(2308)

% include subfunction
addpath(genpath([pwd,'/myFunctions']))
addpath(genpath([pwd,'/TensorFunctions']))

% plotting options
showOutputEvolution= 0;
showStateEvolution= 1;
showFullComparison= 1;

% SISO batches or full MIMO batch in the measurement update
doMU_SISO= 1;

%% System Initialization

% system size (available between 2 and 35)
diameter= 5;
sigma= 3.33e-3;

samplingLens= 2;
cLenslets= samplingLens*diameter; 

% TT ranks
PrankMax= 6;
SrankMax= inf;
% for identified system ttr(Q)=1, ttr(A)=c+1
% for A=I ttr(Q)=c^2, ttr(A)= 1;
rA= 10;
rQ= 10;

% TT rounding tolerance
tol= 1e-4;

% Select System Identification (1: identified, 0: A=Id)
selectSysID= 0;

% Simulation constraints
tMax= 50;
t= 1:tMax;  
% this is the limit for the simulation time. If exceeded the simulation gets cancelled
timeLimit= inf;
% zero time for accuracy evaluation (default: 1)
t0= 1;
fprintf('Simulation time: %i [steps]\n',tMax);
fprintf('Simulation limit: %i [sec]\n\n',timeLimit);

%% Matrix system

% simple system identification
if selectSysID
    tempSelID= 'Identified';
else
    tempSelID= 'AIdentity';
end
try
    iter= 5e4;
    % try to load system from storage
    temp= ['SystemID/',tempSelID,'/AOSystemID_d',num2str(diameter),...
            '_iter',num2str(iter),'.mat'];
    load(temp);
catch ME
    % if it is not in storage, do system identification
    fprintf('System identification setting not stored! Will be identified!\n');
    addpath(genpath([pwd,'/SystemID']))
    data= getAOsystemMat(iter,cLenslets,diameter,double(not(selectSysID)));
end
% get matrices from loaded data
A= data.A;
G= data.C;
Q= data.Q;
n= size(A,1);
p= size(G,1);
fprintf('\tLenslets:\t%i\n',cLenslets)
fprintf('\tState size:\t%i\n',n)
fprintf('\tOutput size:\t%i\n\n',p)

% Rearrange C such that s=[sx,sy] instead of [sx1,sy1,sx2,sy2,...]
% separate sx and sy
C= zeros(size(G));
for i= 1:p
    if mod(i,2)==1
        C((i+1)/2,:)= G(i,:);
    else
        C(p/2+i/2,:)= G(i,:);
    end
end

% Measurement covariance matrix
R= sigma^2*eye(p);

%% TT representation

d= 2;

%%%% A matrix
if selectSysID
    % TT-decomposition of A-matrix
    Att= myTTSVD(A,[cLenslets cLenslets;cLenslets cLenslets],'maxrank',rA);
    % Att= roundTN(Att,tol,rA);
else
    % A-matrix is diagonal, direct initialization possible
    Att.n= [ones(d,1) cLenslets*ones(d,2) ones(d,1)];
    Att.core{1}= reshape(A(1)*eye(cLenslets),Att.n(1,:));
    Att.core{2}= reshape(eye(cLenslets),Att.n(2,:));
end

%%%% C1 matrix
% Inits, see paper in thesis for explanation
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

%%%% Q matrix
Qtt= myTTSVD(Q,[cLenslets cLenslets;cLenslets cLenslets],'maxrank',rQ);
% Qtt= roundTN(Qtt,tol,rQ);

%%%% R1 matrix
Rtt.n= [ones(d,1) [cLenslets-1 cLenslets-1;cLenslets-1 cLenslets-1] ones(d,1)];
for i= 1:d
    Rtt.core{i}= sigma^(2/d)*reshape(eye(Rtt.n(i,2:3)),Rtt.n(i,:));
end

%% Plant evolution

% Initial state
phi0= randn(n,1);
phi0= phi0-mean(phi0);

fprintf('Plant simulation started...')
[phi,s]= myMIMOPlantSimulation(A,C,Q,R,t,phi0);
fprintf('\t...done\n\n')

if ~selectSysID % only for A=a*I
    [U1,~,V]= svd(C);
    U1= U1(:,1:end-2);
    V1= V(:,1:end-2);
    
    phi= V1'*phi;
end

%% Conventional KF

fprintf('Conventional KF started...');
try
    % run filter
    [shat_mat,phihat_mat,time_mat,K_mat]= myKF_FullSim(A,C,Q,R,t,s,timeLimit); %Pnorm_mat
    
    if ~selectSysID % only for A=a*I
        phihat_mat= V1'*phihat_mat;
    end
    
    fprintf('\t...Finished in t=%4.3f [sec/step]\n',time_mat/tMax);
catch ME
    % ended with exception
    fprintf('\t...Ended with exception error!\n');
    
    % set output variables of tensor KF to NaN
    phihat_mat= NaN(n,tMax);
    shat_mat= NaN(p,tMax);
    time_mat= NaN(1);
end

%% Tensor Kalman Filter

cdata.G1= G1;
cdata.E1= E1;
cdata.G2= G2;
cdata.E2= E2;

Rtt2.n= ones(2,4);
for i= 1:d
    Rtt2.core{i}= sigma^(2/d)*reshape(eye(Rtt2.n(i,2:3)),Rtt2.n(i,:));
end

if doMU_SISO
    Rtt= Rtt2;
end

try
    % Tensor Kalman filter
    fprintf('Tensor KF started...');
    [shat_ten,phihat_ten,time_ten,K_ten]= myMIMOTensorKF_FullSim_AO(Att,Ctt1,... % Pnorm_ten
        Ctt2,Qtt,Rtt,1:tMax,s,tol,'PrankMax',PrankMax,'SrankMax',SrankMax,...
        'cdata',cdata);
    
    if ~selectSysID % only for A=a*I
        phihat_ten= V1'*phihat_ten;
    end
    
    fprintf('\t\t...Finished in t=%4.3f [sec/step]\n',time_ten/tMax);
catch ME
    % ended with exception
    fprintf('\t...Ended with exception error!\n');
    
    % set output variables of tensor KF to NaN
    phihat_ten= NaN(n,tMax);
    shat_ten= NaN(p,tMax);
    time_ten= NaN(1);
end

%% Massioni KF

fprintf('Massioni KF started...');
try
    % run filter
    [shat_mass,phihat_mass,time_mass]= massioniKF_FullSim(A,C,Q,R,t,s,timeLimit);
    
    if ~selectSysID % only for A=a*I
        phihat_mass= V1'*phihat_mass;
    end
    
    fprintf('\t\t...Finished in t=%4.3f [sec/step]\n',time_mass/tMax);
catch ME
    % ended with exception
    fprintf('\t...Ended with exception error!\n');
    
    % set output variables of tensor KF to NaN
    phihat_mass= NaN(n,tMax);
    shat_mass= NaN(p,tMax);
    time_mass= NaN(1);
end
clear temp timeKFi

%% ETKF

% number of ensembles
Ens= ceil(p/2);

fprintf('Ensemble KF started...');
try
    % run filter
    [shat_ens,phihat_ens,time_ens]= myETKF_FullSim(A,C,Q,R,Ens,t,s,timeLimit);
    
    if ~selectSysID % only for A=a*I
        phihat_ens= V1'*phihat_ens;
    end
    
    fprintf('\t\t...Finished in t=%4.3f [sec/step]\n',time_ens/tMax);
catch ME
    % ended with exception
    fprintf('\t...Ended with exception error!\n');
    
    % set output variables of tensor KF to NaN
    phihat_ens= NaN(n,tMax);
    shat_ens= NaN(p,tMax);
    time_ens= NaN(1);
end
clear temp timeKFi

%% Norm of state

fprintf('\n\nRelative State Norm:\n');

% tensor filter
temp= phihat_ten(:,t0:end)-phi(:,t0:end);
temp= norm(vec(temp)).^2/norm(vec(phi(:,t0:end))).^2;
fprintf('\tTensor KF:\t\t%4.3e\n',temp);
% conventional filter
temp= phihat_mat(:,t0:end)-phi(:,t0:end);
temp= norm(vec(temp)).^2/norm(vec(phi(:,t0:end))).^2;
fprintf('\tConventional KF:\t%4.3e\n',temp);
% Massioni filter
temp= phihat_mass(:,t0:end)-phi(:,t0:end);
temp= norm(vec(temp)).^2/norm(vec(phi(:,t0:end))).^2;
fprintf('\tMassioni KF:\t\t%4.3e\n',temp);
% Ensemble transformation filter
temp= phihat_ens(:,t0:end)-phi(:,t0:end);
temp= norm(vec(temp)).^2/norm(vec(phi(:,t0:end))).^2;
fprintf('\tETKF:\t\t\t%4.3e\n',temp);


%% Norm of output

fprintf('\n\nRelative Output Norm:\n');

% tensor filter
temp= s(:,t0:end)-shat_ten(:,t0:end);
temp= norm(vec(temp)).^2/norm(vec(s(:,t0:end))).^2;
fprintf('\tTensor KF:\t\t%4.3e\n',temp);
% conventional filter
temp= shat_mat(:,t0:end)-s(:,t0:end);
temp= norm(vec(temp)).^2/norm(vec(s(:,t0:end))).^2;
fprintf('\tConventional KF:\t%4.3e\n',temp);
% Massioni filter
temp= shat_mass(:,t0:end)-s(:,t0:end);
temp= norm(vec(temp)).^2/norm(vec(s(:,t0:end))).^2;
fprintf('\tMassioni KF:\t\t%4.3e\n',temp);
% Ensemble Transformation filter
temp= shat_ens(:,t0:end)-s(:,t0:end);
temp= norm(vec(temp)).^2/norm(vec(s(:,t0:end))).^2;
fprintf('\tETKF:\t\t\t%4.3e\n',temp);

%% Graphical Output

% plot one random output
nrOut= randi(p);

%%%% System output || system state
if showOutputEvolution
    figure('position',[100 100 640 480]);
    hold on
    % real system
    plot(t,s(nrOut,t0:end),'kx-','LineWidth',2)
    % KF
    plot(t,shat_mat(nrOut,t0:end),'b--','LineWidth',2)
    % TKF
    plot(t,shat_ten(nrOut,t0:end),'r:','LineWidth',2)
    % Mass
    plot(t,shat_mass(nrOut,t0:end),'g-.','LineWidth',2);
    % ETKF
    plot(t,shat_ens(nrOut,t0:end),'m:','LineWidth',2);
    % settings
    grid on
    xlabel('time')
    ylabel('system output')
    title(['System output evolution of output ',num2str(nrOut)])
    legend('Noisy system',['KF est.; time=',num2str(time_mat/tMax),'sec/step'],...
        ['Tensor KF est.; time=',num2str(time_ten/tMax),'sec/step'],...
        ['Mass. KF est.; time=',num2str(time_mass/tMax),'sec/step'],...
        ['ETKF est.; time=',num2str(time_ens/tMax),'sec/step'],'Location','best');
end
%%

nrOut= randi(n);

if showStateEvolution
    figure('position',[100 100 640 480]);
    hold on
    % real system
    plot(t,phi(nrOut,t0:end),'kx-','LineWidth',2)
    hold on;
    % KF
    plot(t,phihat_mat(nrOut,t0:end),'b--','LineWidth',2)
    % TKF
    plot(t,phihat_ten(nrOut,t0:end),'r:','LineWidth',2)
    % Mass
    plot(t,phihat_mass(nrOut,t0:end),'g-.','LineWidth',2)
    % ETKF
    plot(t,phihat_ens(nrOut,t0:end),'m:','LineWidth',2)
    % settings
    grid on
    xlabel('time')
    ylabel('system state')
    title(['System state evolution of state ',num2str(nrOut)])
    legend('Noisy system',['KF est.; time=',num2str(time_mat/tMax),'sec/step'],...
        ['Tensor KF est.; time=',num2str(time_ten/tMax),'sec/step'],...
        ['Mass. KF est.; time=',num2str(time_mass/tMax),'sec/step'],...
        ['ETKF est.; time=',num2str(time_ens/tMax),'sec/step'],'Location','best');
end
%%
if showFullComparison
    fig= figure('position',[100 100 640 1.6*480]);
    subplot(4,1,1)
    imagesc(phi(:,t0:end)-phihat_mat(:,t0:end)), colorbar, title('x-xhat_{mat}'), xlabel('time'), ylabel('n')
    subplot(4,1,2)
    imagesc(phi(:,t0:end)-phihat_ten(:,t0:end)), colorbar, title('x-xhat_{ten}'), xlabel('time'), ylabel('n')
    subplot(4,1,3)
    imagesc(phi(:,t0:end)-phihat_mass(:,t0:end)), colorbar, title('x-xhat_{mass}'), xlabel('time'), ylabel('n')
    subplot(4,1,4)
    imagesc(phi(:,t0:end)-phihat_ens(:,t0:end)), colorbar, title('x-xhat_{ens}'), xlabel('time'), ylabel('n')
end
