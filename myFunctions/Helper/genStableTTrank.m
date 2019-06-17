function Out= genStableTTrank(n,d,r)
% Random Stable TT with given TN-Rank Test
% genStableTTrank.m
% Date:             16.11.2018
% Authors:          Daniel Gedon, 4735226
% Description:      The function generates a stable random TT with
%                   specified TN-rank which resembles a matrix of size (n^d
%                   x n^d).
% Inputs:           n - mode size
%                   d - tensor order
%                   r - TN-rank
% Outputs:          Out - Output TT/TTm in Kims format

%% Alogrithm

%%%%%%% For TTm %%%%%%%
% check for rank consistency
if r > n^2
    %warning('r > n^2. Rank is artificially increased to r=%i but rounding decreases to minmal rank n^2=%i.',r,n^2);
end

% random TT-tensor
ttm= myTTMrand(n,d,r);
% increase TN-rank by one until the desired one is achived
% ATTENTION: By rounding this reduces again to a TN-rank of n^2
if max(ttm.n(:,1)) < r
    temp= r-ttm.n(2:end,1);
    tt_temp= myTTMrand(n,d,temp(end));
    ttm= addTN(ttm,tt_temp);
end

% stabilize it by reducing its spectral radius
while ~all(abs(eig(contract(ttm))) <1)
    ttm.core{end}= ttm.core{end}/1.2;
end

% Output
Out= ttm;
