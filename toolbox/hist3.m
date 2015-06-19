function [nn,x1,x2] = hist3(x,n,m)
%HIST3  simple 3-D Histogram.
%   X must be a 2 by * matrix.
%   N must be a scalar which specifys the number of bins 
%     for the first culumn of X. (optional)
%   M must be a scalar which specifys the number of bins 
%     for the second culumn of X. (optional)
%
%   NN = HIST(X) bins the elements of X into 10 equally spaced containers
%   and returns the number of elements in each container. 
%
%   NN = HIST(X,N), where N is a scalar, uses N bins for both columns.
%
%   NN = HIST(X,N,M), where N and M are scalars, uses N bins for the
%   first column and M for the second.
%
%   [NN,X1,X2] = HIST(...) also returns the position of the bin centers in X.
%
if nargin == 0
    error('Requires one or two or three input arguments.')
end
if nargin == 1
    n = 10;
    m = 10;
end
if nargin == 2
    m = n;
end
if isstr(x) | isstr(n) | isstr(m)
    error('Input arguments must be numeric.')
end
%
minx1 = min(x(1,:));
maxx1 = max(x(1,:));
bw1 = (maxx1 - minx1) ./ n;
x1 = minx1 + bw1*(1:n);
x1(length(x1)) = maxx1;
%
minx2 = min(x(2,:));
maxx2 = max(x(2,:));
bw2 = (maxx2 - minx2) ./ m;
x2 = minx2 + bw2*(1:m);
x2(length(x2)) = maxx2;
%
nn=zeros(length(x1),length(x2));
for k=1:length(x)

  [y,j]=min(abs(ceil(x(1,k)-x1)));
  [z,i]=min(abs(ceil(x(2,k)-x2)));
  nn(i,j)=nn(i,j)+1;
end
if nargout == 0
   [X,Y]=meshgrid(x1,x2);
%   mesh(X,Y,nn);
  bar3(x1-bw1/2,nn,'hist');
end
