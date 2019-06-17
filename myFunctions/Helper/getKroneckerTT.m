function [A_tt]= getKroneckerTT(n,d,r,A1)
%% get a TT in Kronecker representation
% getKroneckerTT.m
% Date:             20.05.2019
% Authors:          Daniel Gedon, 4735226
% Description:      The function uses each slice of the input order-3
%                   tensor (A1) and generates a TT that is equal to a
%                   Kronecker system with A11 x A12 x A13...
% Inputs:           n - mode size
%                   d - tensor order
%                   r - tensor train rank
%                   A - order-3 tensor containing the slices (optional)
% Outputs:          A_TT - A in TT format equal to Kronecker representation
%% Algorithm

% Kronecker rank
A_tt.n= [[1;r*ones(d-1,1)] n*ones(d,2) [r*ones(d-1,1);1]];
% initialization of cores
for id= 1:d
    A_tt.core{id}= zeros(A_tt.n(id,:));
end

% get numbers
% get random values, multiplier to ensure det()~= 0
if nargin == 3
    A1= abs(randn(n,n,d,r))+1*ones(n,n,d,r);%5*rand(n,n,d,r);
elseif nargin == 4
    % do nothing
end

% sum over all ranks
for ir= 1:r
    %%%%%% get TT rank 1 model A1_tt %%%%%%
    A1_tt.n= [ones(d,1) n*ones(d,2) ones(d,1)];
    % get rank 1 TT model
    for id= 1:d
        % check determinants
        if abs(det(A1(:,:,id,ir))) < 1e-7
            error('Determinant too low: %2.3e',det(A1(:,:,id,ir)))
        end
        % check ranks
        if rank(A1(:,:,id,ir)) ~= n
            error('Submatrix not fulll rank')
        end
        % reshape to correct rank 1 model
        A1_tt.core{id}= reshape(A1(:,:,id,ir),A1_tt.n(id,:));
    end
    
    %%%%%% Combine to TT rank r model A_tt %%%%%%
    % first core
    A_tt.core{1}(1,:,:,ir)= A1_tt.core{1};
    % cores 2 to d-1
    for id= 2:d-1
        A_tt.core{id}(ir,:,:,ir)= A1_tt.core{id};
    end
    % last core
    A_tt.core{d}(ir,:,:,1)= A1_tt.core{d};
end