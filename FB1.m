function [a3] = FB1(f,alfa)
N=length(f);
nb=(1:N);
a3=zeros(1,N);
for m1=1:N
%a3(m1)=(2/(N^2*(besselj(1,alfa(m1))).^2))*sum(nb.*f.*besselj(0,alfa(m1)/N*nb));
a3(m1)=bsxfun(@times,(2/(bsxfun(@times,N^2,(besselj(1,alfa(m1))).^2))),sum(bsxfun(@times,bsxfun(@times,nb,f),besselj(0,alfa(m1)/N*nb))));
end