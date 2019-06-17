function TT= myVec2TT_ALS(vec,TT_0,varargin)
%% Function conversion from vector to TT using ALS
% myVec2TT_ALS.m
% Date:             20.05.2019
% Authors:          Daniel Gedon, 4735226
% Description:      The function converts a vector of any size in a TT
%                   using the ALS scheme. It is initialized with a vector
%                   TT_0 in TT-form. The output vector will have the same
%                   form as TT_0 with the same TT-ranks.
% Inputs:           vec - vector which should be in TT
%                   TT_0 - initialization in TT form
%                   varargin - for only forward or forward and backward
%                   sweep
% Outputs:          TT - vector vec in TT form.
%% inputs

% defaul values
sweepback= 1;
for i=1:2:length(varargin)-1
    if (~isempty(varargin{i+1}))
        switch lower(varargin{i})
            case 'sweepbackward'
                sweepback=varargin{i+1};
            otherwise
                error('Unrecognized option: %s\n',varargin{i});
        end
    end
end

% augmentation
X_TT= TT_0;
y1= vec;

% tensor size
d= size(TT_0.n,1);

%% Forward sweep

% for all core to be optimized
for i=1:d-1
    % contract all cores above i
    if i== 1
        newCore= y1;
    else
        newCore= contractAbove(y1,X_TT,i);
    end
    
    % contract all cores below i
    newCore= contractBelow(newCore,X_TT,i);
    
    % new ith core of X_TT
    X_TT.core{i}= reshape(newCore,X_TT.n(i,:));
    
    % orthogonalization of first core. Move norm to second core
    X_TT= TT_orth_i(X_TT,i,1);
end

% comparison
% varin(count)= norm(y1-contract(X_TT));

%% Backward sweep

% only do backward sweep if desired
if sweepback
    % for all cores to be optimized
    for i= d:-1:2
        % contract all cores above i
        newCore= contractAbove(y1,X_TT,i);
        
        % contract all cores below i
        if i == d
            % do nothing, newCore= newCore;
        else
            newCore= contractBelow(newCore,X_TT,i);
        end
        
        % new it core of X_TT
        X_TT.core{i}= reshape(newCore,X_TT.n(i,:));
        
        % orthogonalization of second core. Move norm to first core
        X_TT= TT_orth_i(X_TT,i,-1);
    end
    
    % comparison
    % var(count)= norm(y1-contract(X_TT));
end

%% output assignment
TT= X_TT;

end

%% Function contract cores above the core to be optimized

function Ci= contractAbove(A_mat,B_TT,id)

% contract all above
coreB= contractCores(B_TT,id,1);
tempSize= size(coreB);
if length(tempSize) == 2
    tempSize= [tempSize 1];
end
coreB= reshape(coreB,[prod(tempSize(1:2)) tempSize(3)]);
coreB= coreB';

% reshape vector for contraction
coreA= reshape(A_mat,[tempSize(2) numel(A_mat)/tempSize(2)]);

% contraction
newCore= coreB*coreA;

% vectorization of newCore
Ci= vec(newCore);

end

%% Function contract cores below the core to be optimized

function Ci= contractBelow(A_mat,B_TT,id)

% contract all below
coreB= contractCores(B_TT,id,-1);
coreB= coreB';

% reshape vector for contraction
tempSize= size(coreB);
coreA= reshape(A_mat,[numel(A_mat)/tempSize(1) tempSize(1)]);

% contraction
newCore= coreA*coreB;

% vectorization
Ci= vec(newCore);

end

%% Function contract all Cores above/below core id

function core= contractCores(A,id,dir)


if dir == 1 % left -> right (top -> bottom)
    % initialization
    ni1= [1 1 1];
    core= reshape(1,ni1);
    for i= 0:id-2
        % modes of i+1
        ni2= A.n(i+1,:);
        % reshape core i
        temp1= reshape(core,[prod(ni1(1:end-1)) ni1(end)]);
        % reshape core i+1
        temp2= reshape(A.core{i+1},[ni2(1) prod(ni2(2:end))]);
        % combine to new core
        core= temp1*temp2;
        % reshape correctly
        ni1= [ni1(1) ni1(2)*ni2(2) ni2(3)];
        core= reshape(core,ni1);
    end
elseif dir == -1 % right -> left (bottom -> top)
    % tensor train length
    d= size(A.n,1);
    % initialization
    %     ni1= A.n(d,:);
    %     core= A.core{d};
    ni1= [1 1 1];
    core= reshape(1,ni1);
    for i= d+1:-1:id+2
        % modes of i-1
        ni2= A.n(i-1,:);
        % reshape core i
        temp1= reshape(core,[ni1(1) prod(ni1(2:end))]);
        % reshape core i-1
        temp2= reshape(A.core{i-1},[prod(ni2(1:end-1)) ni2(end)]);
        % combine to new core
        core= temp2*temp1;
        % reshape correctly
        ni1= [ni2(1) ni2(2)*ni1(2) ni1(3)];
        core= reshape(core,ni1);
    end
else
    error('Variable dir not correctly assigned');
end

end