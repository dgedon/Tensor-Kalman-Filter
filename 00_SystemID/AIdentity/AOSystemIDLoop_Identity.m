% AO System identification
% AOSystemIDLoop.m
% Date:             03.05.2019
% Authors:          Daniel Gedon, 4735226
% Description:      This script runs the system identification for the AO
%                   system and stores the identified matrices for later
%                   use.
clear
close all
clc

% set random number generator
rng(12345)

%%

AIdentity= 1;
samplingLens= 2;
iter= 5e4;

for diameter= 2:42
    cLenslets= samplingLens*diameter; 

    fprintf('Identification for d=%i, c=%i\n',diameter,cLenslets);
    % simple system identification
    data= getAOsystemMat(iter,cLenslets,diameter,AIdentity);
    % store
    temp= ['AOSystemID_d',num2str(diameter),'_iter',num2str(iter),'.mat'];
    save(temp,'data')
end

%% 
% theoretically one would need about two times more iterations than entries
% in the larges matrix to identify. The largest matrix for c>4 is the
% measurement matrix C with pn=2c^4-4c^4+2c^2 entries. Hence with 5e4
% iterations only a reliable system identification up to c=13 is possible.
% However analysis has shown that for c=40 the A and Q matrix still look
% reasonably well. This could accure due to the high sparsity and
% regularity of the matrices. Hence all systems here are identified with
% 5e4 iterations. For more iterations, MATLAB crashes due to the inversion
% of a matrix of size c^2 x iterations which is demanding.