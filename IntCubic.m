function [y] = IntCubic(x,p,q)
% given x(n), produces y(n) = x(n*p/q)  
% using cubic interpolation to fill in the missing pixels

% pad zeros for the boundaries
x = x(:)';
s = [0 x 0];

c1 = -2 + sqrt(3);
c2 = -2 - sqrt(3);
y1 = filter(1,[1 -c1],s,[],1);
y2 = filter(-1/c2,[1 -1/c2],y1(end:-1:1,:),[],1);
c = y2(end:-1:1,:)*6;

cub = @(x) double(abs(x)<1).*(2/3 - (abs(x)).^2 + (abs(x)).^3 / 2) + double((abs(x)>=1).*(abs(x)<2)).*(2 - abs(x)).^3/6;

h = cub(-2 + 1/q:1/q:2 - 1/q);
y = upfirdn(c,h,p,q);