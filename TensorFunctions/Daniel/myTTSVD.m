function TT= myTTSVD(A,modes,varargin)
%% TT SVD function
% myTTSVD.m
% Date:             13.05.2019
% Authors:          Daniel Gedon, 4735226
% Description:      The function computes the tensor train svd algorithm
% Inputs:           A - matrix to be decomposed
%                   modes - desired mode sizes
%                   varargin - e.g. tolerance
% Outputs:          TT - decomposed A into a tensor train
%% Inputs

% defaul values
tol= eps;
maxrank= inf;

for i=1:2:length(varargin)-1
    if (~isempty(varargin{i+1}))
        switch lower(varargin{i})
            case 'tolerance'
                tol=varargin{i+1};
            case 'maxrank'
                maxrank= varargin{i+1};
            otherwise
                error('Unrecognized option: %s\n',varargin{i});
        end
    end
end

% tensor order and modes
d= size(modes,1);
n= prod(modes,2);

% if matrix input
ismatrix= 0;
if ~isvector(A)
    ismatrix= 1;
end
A= A(:);

%% initialization

TT.n= [ones(d,1) n ones(d,1)];
for i= 1:d
    TT.core{i}= zeros(TT.n(i,:));
end

r= [TT.n(:,1);TT.n(end,end)];

%% algorithm

% truncation parameter
delta= tol*norm(A(:))/sqrt(d-1);

% temporary tensor
C= A;

for i= 1:d-1
    % for matrix case: permutation needed
    if ismatrix
        temp= [modes(i,1) prod(modes(i+1:end,1)) ...
            modes(i,2) prod(modes(i+1:end,2))];
        C= reshape(C,[r(i) temp r(i+1)]);
        C= permute(C,[1 2 4 3 5 6]);
    end
    
    % reshape core i
    temp= r(i)*prod(n(i,:));
    C= reshape(C,[temp,numel(C)/temp]);
    
    % economic SVD
    [U,S,V]=svd(C,'econ');
    
    % choose TT rank
    s= diag(S);
    ri=0;
    while norm(s(end:-1:end-ri)) < delta 
        ri=ri+1;
    end
    ri=length(s)-ri;
    ri= min(ri,maxrank);
    % update ranks
    r(i+1)= ri;
    TT= updateTTranks(TT,r);
    
    % new core i
    TT.core{i}= reshape(U(:,1:ri),TT.n(i,:));
    
    % new temporary tensor 
    C= S(1:ri,1:ri)*V(:,1:ri)';
end

% last core
TT.core{d}= reshape(C,TT.n(d,:));

% reshape to correct modes if matrix input
if size(modes,2) > 1
    % original modes
    TT.n= [TT.n(:,1) modes TT.n(:,end)];
    % reshape cores
    for i= 1:d
        TT.core{i}= reshape(TT.core{i},TT.n(i,:));
    end
end
end


%% Subfunctions

function TT= updateTTranks(TT,r)

TT.n(:,1)= r(1:end-1);
TT.n(:,end)= r(2:end);

end