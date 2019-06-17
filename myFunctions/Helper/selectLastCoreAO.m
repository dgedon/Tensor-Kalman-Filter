function TTout= selectLastCoreAO(TTin,sel1,varargin)
%% Selects a specific sub-tensor from the last core (specific for AO)
% selectLastCoreAO.m
% Date:             29.04.2019
% Authors:          Daniel Gedon, 4735226
% Description:      Function selects parts of last TT core or 3-way tensor
%                   contracts with previous core
% Inputs:           TTin - TT to be selected
%                   sel1 - selection of last core
% Outputs:          TTout - selected subtensor
%% Algorithm

b= reshape(TTin.core{2},[prod(TTin.n(2,1:end-1)) TTin.n(2,end)]);
if size(TTin.n,2) == 3
    b= b*TTin.core{3}(:,sel1);
else
    b= b*TTin.core{3}(:,sel1,varargin{1});
end

TTout.n= TTin.n(1:2,:);
TTout.n(2,end)= 1;
core= reshape(b,TTout.n(2,:));

TTout.core{1}= TTin.core{1};
TTout.core{2}= core;

