clear all;clc;
[x,Fs] = audioread('QASK.m4a');
for n = 1: 442 : 43218
    mess      = x(n:n+881,1);
    for s     = 1 : 882
    w         = 1 - cos(2*pi*s/881);        %window
    mess(s,1) = w.*mess(s,1)/2;
    end
    Rxx_sum   = lagmatrix(xcorr(mess)',-881)';
    Rxx_sum(isnan(Rxx_sum)) = 0;
    a(0+1)    = 1;
    a(1+1)    = -(Rxx_sum(2)/Rxx_sum(1));
%Declaration
    gam       = zeros(1,31);
    e_min     = zeros(1,31);
    g         = zeros(30,63);
    gam(1)    = -a(1+1);
    e_min(1)  = Rxx_sum(1)*(1-a(1+1)^2);
    e_min(2)  = e_min(1).*(1-gam(1)^2);
%g(1,k) with origin = g(1,32)

for b         = -31 : 31
    for q     = 0 : 1
    g(1,b+31+1) = g(1,b+31+1) + Rxx_sum(abs(q-b)+1)*a(1+q);
    end
    end
    gam(2)    = -g(1,31+3)/g(1,31+1);
%finding gamma and e_min    
for m        = 2 : 30
K            = lagmatrix(fliplr(g(m-1,:)),m)';
K(isnan(K))  = 0;
g(m,:)       = g(m-1,:) + gam(m).*K;
gam(m+1)     = -g(m,31+m+2)/g(m,31+1);
e_min(m+1)   = e_min(m).*(1-gam(m)^2);  
end
    e        = zeros(1,882);
%generating error
for h        = 1 : 31
if(h == 1)
hjk           = lagmatrix(mess',1)';
else
hjk           = lagmatrix(er,1)';
end
hjk(isnan(hjk)) = 0;
e             = mess' - gam(h).*hjk; 
er            = hjk - mess'*gam(h);
mess          = e';
end
if(n ==1)
syn = 0;
error = 0;
end
error = [error e];
end


op = zeros(1,43129);
kl = size(error)/2;
sn = zeros(1,43129);
for h         = 1 : 43129
    sn(h) = error(2*h);
end
%sound(100*sn,44100)
error2 = randn(1,43129);
%Synthesize
for(h = 1 : 31)
     hjk = lagmatrix(op,1)';
     hjk(isnan(hjk)) = 0;
     op = error2 + gam(h).*hjk;
     error2 = op;
end
    syn = [syn op];
