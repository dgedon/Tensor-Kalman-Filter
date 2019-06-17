function a= TT_orth_i(a,i,dir)
%% Function orthogonalization of ith core of TT a in given direction
% TT_orth_i.m
% Date:             20.05.2019
% Authors:          Daniel Gedon, 4735226
% Description:      The function moves the non-orthogonal part from core i
%                   depending on the direction either to core i+1 (for
%                   dir==1) or to core i-1 (for dir==-1).
% Inputs:           a - complete TT 
%                   i - core to be orhogonalized
%                   dir - direction where to move non-orthogonal part
% Outputs:          a - new TT representation
%% Algorithm

if dir == -1 % right to left orthogonalization via QR
    [Q,R]=qr(reshape(a.core{i},[a.n(i,1),prod(a.n(i,2:end))])',0);
    % ith core, make orthogonal
    a.core{i}=reshape(Q(:,1:a.n(i,1))',[a.n(i,1),prod(a.n(i,2:end-1)),size(Q,1)/(prod(a.n(i,2:end-1)))]);
    % i-1th core, move norm there
    temp= reshape(a.core{i-1},[prod(a.n(i-1,1:end-1)),a.n(i-1,end)]);
    a.core{i-1}=reshape(temp*R(1:a.n(i,1),:)',[a.n(i-1,1),prod(a.n(i-1,2:end-1)),a.n(i-1,end)]);
elseif dir == 1 % left to right orthogonalization via QR
    % qr decomposition (left to right)
    sizeResh= [prod(a.n(i,1:end-1)), a.n(i,end)];
    [Q,R]= qr(reshape(a.core{i},sizeResh),0);
    % new first core
    a.core{i}= reshape(Q(:,1:a.n(i,end)),a.n(i,:));
    % new next core
    temp= reshape(a.core{i+1},[a.n(i+1,1) prod(a.n(i+1,2:end))]);
    a.core{i+1}= reshape(R(1:a.n(i,end),:)*temp,a.n(i+1,:));
else
    error('Variable dir not correctly assigned');
end