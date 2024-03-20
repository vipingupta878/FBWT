function [w] = RAnDwt(x,p,q,r,s,J,F,alfa)

% make sure 'x' is a row vector of even length
x = x(:);x = x.';
L = length(x); 
N = L + mod(L,2); % x should be even length
x = [x zeros(1,N-L)];
clear L;

if (N * ((p/q)^J))*r/(s) < 2,
    error('Too many subbands -- Reduce ''J''');
end
    
X = FB1(x,alfa);


for n = 1:J,
           
    G = F{n,1};
    f = F{n,2};
    
    % positive frequencies
    dd = length(G);
    if n==1
    sub = G(1:length((f{1}:length(X)))).*X(f{1}:end);
     sub2 = iFB(sub,alfa(1:length(sub)));%*sqrt(length(1+(f{1}:length(X)-1)));
    else
    sub = G(1:end).*X(f{1}:f{1}+dd-1);
     sub2 = iFB(sub,alfa(1:length(sub)));%*sqrt(N1);
    end

w{n,1} = (sub2);   
end
H = F{J+1};
%N1 = PQ(J+1,1); 
sub = X(1+(0:f{2}-1)).*H(1+(0:f{2}-1));
sub2 = iFB(sub,alfa(1:length(sub)));%*sqrt(N1);
w{J+1,1} = (sub2);