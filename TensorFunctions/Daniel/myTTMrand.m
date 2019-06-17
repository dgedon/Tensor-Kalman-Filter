function TT= myTTMrand(n,d,r)
%% TT rand function
% myTTrand.m
% Date:             13.05.2019
% Authors:          Daniel Gedon, 4735226
% Description:      The function generates a random TT in matrix form
% Inputs:           n - mode size
%                   d - tensor train order
%                   r - tensor train rank
% Outputs:          TT - random output tensor train
%% Initialization

% ranks
r= r*ones(d+1,1);
r(1)= 1;
r(end)= 1;

% modes (ttm/MPO specific)
n= n*ones(d,2);

TT.n= [ones(d,1) n ones(d,1)];
for i= 1:d
    TT.core{i}= zeros(TT.n(i,:));
end

%% algorithm

for i= 1:d
    cr= randn(r(i)*prod(n(i,:)),r(i+1));
    [cr1,~]= qr(cr,0);
    r(i+1)= size(cr1,2);
    TT= updateTTranks(TT,r);
    TT.core{i}= reshape(cr1,TT.n(i,:));
end


end

%% Subfunctions

function TT= updateTTranks(TT,r)

TT.n(:,1)= r(1:end-1);
TT.n(:,end)= r(2:end);

end