function [yhatOut,xhatOut,timeTKF,K]= ... % Pnorm
    myMIMOTensorKF_FullSim_AO(A,C1,C2,Q,R,t,y,tol,varargin)
%% Function Tensor Kalman Filter with vector output
% myMatrixTensorKF_FullSim_AO.m
% Date:             24.04.2019
% Authors:          Daniel Gedon, 4735226
% Description:      The function does a Kalman state estimation in tensor
%                   representation for the AO system
% Inputs:           A - state evolution matrix in TT format
%                   C1 - measurement matrix in TT format for MU
%                   C2 - measurement matrix in TT format for TU
%                   Q - process noise covariance tensor in TT format
%                   R - measurement noise covariance tensor; size C1 or C2
%                   t - time vector
%                   y - measurement vector in TT format
%                   tol - tensor network conversion and rounding tolerance
% Outputs:          yhatOut - estimated observations (k|k)
%                   xhatOut - estimated state (k|k)
%                   timeTKF - Time for computation
%%


doK= 1;
K.norm= 0;
K.K= 0;

% defaul values
PrankMax= Inf;
SrankMax= Inf;
timeLimit= Inf;

for i=1:2:length(varargin)-1
    if (~isempty(varargin{i+1}))
        switch lower(varargin{i})
            case 'prankmax'
                PrankMax=varargin{i+1};
            case 'srankmax'
                SrankMax=varargin{i+1};
            case 'cdata'
                cdata= varargin{i+1};
            case 'timelimit'
                timeLimit= varargin{i+1};
            otherwise
                error('Unrecognized option: %s\n',varargin{i});
        end
    end
end

%% Initialization of Variables

% get initialization values
% sizes
n= size(contract(A),1);
ni= C1.n(:,3);
p= 2*size(contract(C1),1);
%pi= C.n(:,2);
d= size(C1.n,1);
% simulation time
tMax= numel(t);

%%%% Output variables
xhatOut= zeros(n,tMax);
yhatOut= zeros(p,tMax);

%%%% initialize Xhat, P and K in correct format
% initialize Xhat:
xhat_.n= [ones(d,1) ni ones(d,1)];
xhat_.core= cell(1,d);
xhat_.core{1}= zeros(1,ni(1));
for i= 2:d
    xhat_.core{1,i}= zeros(1,ni(i));
end

% initialize P:
sigmaP= 1;
P_.n= [ones(d,1) ni ni ones(d,1)];
P_.core{1}= reshape(sigmaP*eye(ni(1)),P_.n(1,:));
for i= 2:d
    P_.core{1,i}= reshape(eye(ni(i)),P_.n(i,:));
end
P_= roundTN(P_,eps,PrankMax);

%% Time Simulation of Tensor KF

timeTKF= 0;
try
    for k= t
        % set timer
        tic
        if ~isequal(R.n,ones(size(R.n))) % MIMO MU updates
            
            % convert y to TT
            if 0    % use TT_SVD
                ytt= myTTSVD(y(:,k),[C1.n(:,2);2]);
                ytt= roundTN(ytt,tol);
            else % use TT_ALS
                if k== t(1)
                    % first time step: need TT_SVD conversion
                    % necessary orthogonalization
                    ytt= myTTSVD(y(:,k),[C1.n(:,2);2]);
                    ytt= roundTN(ytt,tol);
                    % orthogonalization
                    di= d+1; % since d=2 and actually y has d=3
                    for id= di:-1:2
                        % move norm to first core
                        ytt= TT_orth_i(ytt,id,-1);
                    end
                else
                    % otherwise use TT_ALS with previous yTT as init
                    % only do one half sweep (sufficient)
                    ytt= myVec2TT_ALS(y(:,k),ytt,'sweepBackward',0);
                    % orthogonalization (norm currently in last core)
                    for id= di:-1:2
                        % move norm to first core
                        ytt= TT_orth_i(ytt,id,-1);
                    end
                end
            end
            
            % select first and second part of y
            ytt1= selectLastCoreAO(ytt,1);
            ytt2= selectLastCoreAO(ytt,2);
            
            %%%% Measurement update X with intermediate results
            [xhatint,yhat1,P_int,temp]= myMIMOTensorKF_MU(C1,R,xhat_,ytt1,P_,tol,...
                'PrankMax',PrankMax,'SrankMax',SrankMax);
            
            K_MU= zeros(n,p);
            K_MU(:,1:p/2)= contract(temp);

            
            %%%% Measurement update Y
            [xhat,yhat2,P,temp]= myMIMOTensorKF_MU(C2,R,xhatint,ytt2,P_int,tol,...
                'PrankMax',inf,'SrankMax',SrankMax);
            
            K_MU(:,p/2+1:end)= contract(temp);
            
        else % SISO measurement update
            
            K_MU= zeros(n,p);
            
            %%%% Measurement update X with intermediate results
            tempsize= C1.n(1,2);
            for row= 1:p/2
                i= ceil(row/tempsize);
                j= row-(i-1)*tempsize;
                
                C1new.n= C1.n;
                C1new.n(:,2)= ones(2,1);
                C1new.core{1}= reshape(cdata.G1(j,:),C1new.n(1,:));
                C1new.core{2}= reshape(cdata.E1(i,:),C1new.n(2,:));
                
                ytt1new.n= ones(2,3);
                ytt1new.core{1}= reshape(y(row,k),ytt1new.n(1,:));
                ytt1new.core{2}= reshape(1,ytt1new.n(2,:));
                
                [xhatint,yhat1,P_int,temp]= myMIMOTensorKF_MU_SISO(C1new,R,xhat_,ytt1new,P_,tol,...
                    'PrankMax',inf);%2*PrankMax);
                xhat_= xhatint;
                P_= P_int;

                if k==max(t) && doK== 1
                    K_MU(:,row)= contract(temp);
                end
            end
            
            %%%% Measurement update Y
            tempsize= C2.n(1,2);
            for row= 1:p/2
                i= ceil(row/tempsize);
                j= row-(i-1)*tempsize;
                
                C2new.n= C2.n;
                C2new.n(:,2)= ones(2,1);
                C2new.core{1}= reshape(cdata.G2(j,:),C2new.n(1,:));
                C2new.core{2}= reshape(cdata.E2(i,:),C2new.n(2,:));
                
                ytt2new.n= ones(2,3);
                ytt2new.core{1}= reshape(y(row+p/2,k),ytt2new.n(1,:));
                ytt2new.core{2}= reshape(1,ytt2new.n(2,:));
                
                [xhat,yhat2,P,temp]= myMIMOTensorKF_MU_SISO(C2new,R,xhatint,ytt2new,P_int,tol,... % P_int
                    'PrankMax',inf);%2*PrankMax);
                xhatint= xhat;
                P_int= P;

                if k==max(t) && doK==1
                    K_MU(:,row+p/2)= contract(temp);
                end
            end
        end
        
        %%%% Time update
        [xhat_,P_]= myMIMOTensorKF_TU(A,Q,xhat,P,tol,'PrankMax',PrankMax);
        
        % timer count
        timeTKF= timeTKF+ toc;
        
        % limit for computation time
        if toc > timeLimit
            error('timeout');
        end
        
        if k==max(t) && doK==1
            K.norm= norm(K_MU); % (k)
            K.K= K_MU; % only for last one
        end
        
        %---------- Function Output Variables ----------
        % State xhat_out
        xhatOut(:,k)= contract(xhat);
        
        % Output yHat_Out
        yhatOut(:,k)= 1;    %[contract(yhat1);contract(yhat2)];
    end
catch ME
    if strcmpi(ME.message,'timeout')
        %fprintf('\tException caught timeout (%4.2e > %i) conventional KF!\n',toc,timeLimit);
    end
    %fprintf('\tVariables assigned NaN\n');
    timeTKF= nan;
    yhatOut= nan(p,max(t));
    xhatOut= nan(n,max(t));
end