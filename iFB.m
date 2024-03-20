function [F1_r] = iFB(a3,alfa)
%a3=[zeros(1,1000) a3 0];
a=length(a3);
g1_r=(alfa(1:end))/a; F1_r=zeros(1,a);
for mm=1:a
    F1_r(mm)=sum(bsxfun(@times,a3,besselj(0,bsxfun(@times,g1_r,mm))));
end
