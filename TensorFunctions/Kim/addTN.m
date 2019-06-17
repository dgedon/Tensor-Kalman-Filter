function c=addTN(a,b)
% c=addTN(a,b)
% ------------
% Adds two Tensor Networks a,b together in c. This increases the TN-ranks,
% so you might want to do a rounding step next.
%
% c         =   Tensor Network, sum of a with b,
%
% a         =   Tensor Network, first summand,
%
% b         =   Tensor Network, second summand.
%
% Reference
% ---------
%
% A Tensor Network Kalman filter with an application in recursive MIMO Volterra system identification
%
% 2016, Kim Batselier, Zhongming Chen, Ngai Wong

[d,n]=size(a.n);
c.core=cell(1,d);
c.n=a.n;

% first do borders
c.n(1,end)=c.n(1,end)+b.n(1,end);
c.core{1}=zeros(1,prod(a.n(1,2:end-1)),a.n(1,end)+b.n(1,end));
c.core{1}(:,:,1:a.n(1,end))=reshape(a.core{1},[1,prod(a.n(1,2:end-1)),a.n(1,end)]);
c.core{1}(:,:,a.n(1,end)+1:end)=reshape(b.core{1},[1,prod(b.n(1,2:end-1)),b.n(1,end)]);
% modification Daniel Gedon 15.11.2018: for correct size between c.core
% and c.n
c.core{1}= reshape(c.core{1},c.n(1,:));

c.n(d,1)=c.n(d,1)+b.n(d,1);
c.core{d}=zeros(a.n(d,1)+b.n(d,1),prod(a.n(d,2:end-1)),1);
c.core{d}(1:a.n(d,1),:,:)=reshape(a.core{d},[a.n(d,1),prod(a.n(d,2:end-1)),1]);
c.core{d}(a.n(d,1)+1:end,:,:)=reshape(b.core{d},[b.n(d,1),prod(b.n(d,2:end-1)),1]);
% modification Daniel Gedon 15.11.2018: for correct size between c.core
% and c.n
c.core{d}= reshape(c.core{d},c.n(d,:));

for i=2:d-1
	c.n(i,end)=c.n(i,end)+b.n(i,end);
	c.n(i,1)=c.n(i,1)+b.n(i,1);
    c.core{i}=zeros(a.n(i,1)+b.n(i,1),prod(a.n(i,2:end-1)),a.n(i,end)+b.n(i,end));
    c.core{i}(1:a.n(i,1),:,1:a.n(i,end))=reshape(a.core{i},[a.n(i,1),prod(a.n(i,2:end-1)),a.n(i,end)]);
    c.core{i}(a.n(i,1)+1:end,:,a.n(i,end)+1:end)=reshape(b.core{i},[b.n(i,1),prod(b.n(i,2:end-1)),b.n(i,end)]);   
    % modification Daniel Gedon 15.11.2018: for correct size between c.core
    % and c.n
    c.core{i}= reshape(c.core{i},c.n(i,:));
end


end
