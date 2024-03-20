function [F] = CreateFilters(N,p,q,r,s,bet,J,varargin)

% make sure 'x' is a row vector of even length
N = N + mod(N,2); % N should be even

if (N * ((p/q)^J))*r/(2*s) < 2,
    error('Too many subbands -- Reduce ''J''');
end

%the sampling factors

PQ = zeros(J,2);
RS = zeros(J,2);
for k = 1:J
    p0 = ceil(N * ((p/q)^k) ); p0 = p0 + mod(p0,2); % make sure p0 is even
    PQ(k,1) = p0;
    RS(k,1) = round(N * ( (p/q)^(k-1) ) * r/s);    
end
PQ(1,2) = N;
PQ(2:end,2) = PQ(1:end-1,1);
RS(1:end,2) = PQ(1:end,2);%/2;
 
PQ = [1 1;PQ];
PQk = 1;

Hk = ones(1,2*N);

for n = 1:J

    PQk = PQk*PQ(n,1)/PQ(n,2);
    
    pp = PQ(n+1,1);qq = PQ(n+1,2);
    rr = RS(n,1);ss = RS(n,2);
        
    epsi = (N/32)*(pp-qq + bet*qq)/(pp+qq);
    
    f{1} = ((1-bet)*N + epsi);
    f{2} = (N*pp/qq);
    f{3} = (N - epsi);
    f{4} = (N + epsi);
    
    for k = 1:4,
        %f{k} = round(f{k}*PQk);
        f{k} = (f{k}*PQk);
    end
    [H,G] = MakeFilters(f);
    
    f{1} = ceil(f{1});  
    f{2} = floor(f{2});
    f{4} = floor(f{4});
    % filter for positive frequencies
    Gk = Hk(1+(f{1}:f{4})).*G(1:end);

    dd = min(rr,length(Gk));
    F{n,1} = Gk(1:dd);
    F{n,2} = f;
    
    %update the lowpass filter
    Hk(1+(f{1}:f{2})) = Hk(1+(f{1}:f{2})).*H(1+(f{1}:f{2})); 
    Hk(1+(f{2}+1:f{4}))=0;
    
end
F{J+1,1} = Hk(1:f{2});

%%%%%

function [H, G] = MakeFilters(f)
% Make frequency responses

% MAKE H0
w = (0:floor(f{2}));
k_pass = (w < f{1});                 % pass-band
k_trans = (w >= f{1});    % transition-band
b = (f{2}-f{1})/pi;
w_scaled = (w - f{1})/b;

H = zeros(size(w));
H(k_pass) = 1;%pass band of low pass filter
H(k_trans) = (1+cos(w_scaled(k_trans))) .* sqrt(2-cos(w_scaled(k_trans)))/2;%transition band of low pass filter
seq = sqrt(1 - H(k_trans).^2);

w = (ceil(f{1}):floor(f{4}));
k_pass = (w <= f{3}) & (w >= f{2});                 % pass-band 
k_trans1 = (w <= f{2});    % transition-band1 
k_trans2 = (w >= f{3});    % transition-band2

b = (f{4}-f{3})/pi;
if b > 0,
    w_scaled = (w - f{3})/b;
else
    w_scaled = 0*w;
end
   
G = zeros(size(w));
G(k_pass) = 1; % pass-band of high pass filter
G(k_trans1) = seq; %transition-band1 of high pass filter
G(k_trans2) = (1+cos(w_scaled(k_trans2))) .* sqrt(2-cos(w_scaled(k_trans2)))/2;  % transition-band2 of high pass filter