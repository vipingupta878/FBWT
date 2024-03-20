function [d] = make2D(w,p,q);

N = length(w);
len = length(w{1,1});

d = zeros(N-1,len);
d(1,1:len) = w{1,1};

for k = 1 : N-2,
    c = zeros(1,len);
    cc = w{k+1,1}+1i*w{k+1,2};
    ss = round(len*((p/q)^k));
    ss = min(length(cc),ss);
    in = round(1:len/ss:len);
    c(in) = abs(cc(1:ss));
    for n = 1:ss-1,
        cin = in(n);
        cin2 = in(n+1);
        dif = cin2 - cin;
        for nn = cin+1:cin2-1,
            dif1 = (nn - cin)/dif;
            dif2 = 1-dif1;
            c(nn) = dif1*c(cin) + dif2*c(cin2);
        end
    end
    d(k,1:len) = c;
end
