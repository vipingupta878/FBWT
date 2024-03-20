function [x] = iRAnDwt(w,N,F,alfa)
% make sure 'x' is a row vector of even length
N = N + mod(N,2);  
X = zeros(1,N);
J = size(w,1)-1;


for n = 1:J,
       
    G = conj(F{n,1});
    f = F{n,2};  
    
    
    if n==1
     sub = FB1(w{n,1},alfa(1:size(w{n,1},2)));%/sqrt(length(1+(f{1}:length(X)-1)));    
    X(f{1}:end) = X(f{1}:end) + G(1:length(f{1}:length(X))).*sub; 
     else
    sub = FB1(w{n,1},alfa(1:size(w{n,1},2)));%/sqrt(N1);
%     sub = circshift(sub.',-d);
%    sub = sub.';    
    dd = length(G);
    X(f{1}:f{1}+dd-1) = X(f{1}:f{1}+dd-1) + G.*sub;
   end
    
    
end
sub = FB1(w{J+1,1},alfa(1:size(w{J+1,1},2)));%/sqrt(length(w{J+1,1}));
H = F{J+1,1};  
X(1+(0:f{2}-1)) = X(1+(0:f{2}-1)) + sub.*H;
x = iFB(X,alfa);%*sqrt(N);