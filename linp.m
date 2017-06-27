clear all;clc;
[x,Fs] = audioread('QASK.m4a');%sound(x,44100);
theta = [0 : 0.0001 : 2*pi];z = exp(1i*theta);
for(n = 1: 92 : 44008)
    mess = x(n:n+91,1);
    Rxx_sum = lagmatrix(xcorr(mess)',-91)';
    Rxx_sum(isnan(Rxx_sum)) = 0;
    gam = zeros(1,11);
    e_min = zeros(1,11);
    a = zeros(1,11);
    a(0+1) = 1;
    if(Rxx_sum(1) == 0)
        a(2) = 0;
    else
        a(1+1) = -(Rxx_sum(2)/Rxx_sum(1));
    end
%     [a,g] = lpc(mess,10);
    gam(1) = -a(1+1);
    e_min(1) = Rxx_sum(1)*(1-a(1+1)^2);
    g = zeros(9,25);
    %g(1,:)
    for(b = -10 : 14)
        for(q = 0 : 2)
            g(1,b+10+1) = g(1,b+10+1) + Rxx_sum(abs(q-b)+1)*a(1+q);
        end%g(1,k)
    end
    if(g(1,11) == 0)
        gam(2) = 0;
    else
        gam(2) = -g(1,10+3)/g(1,10+1);
    end
    e_min(2) = e_min(1).*(1-gam(1)^2);
    for(m = 2 : 10)
        K = fliplr(lagmatrix(g(m-1,:),m-(2*(m-2)))');
        K(isnan(K)) = 0;
        g(m,:) = g(m-1,:) + gam(m).*K;
        if(g(m,11) == 0)
            gam(m+1) = 0;
        else
            gam(m+1) = -g(m,10+m+2)/g(m,10+1);
        end
        e_min(m+1) = e_min(m).*(1-gam(m)^2);  
    end
%    
%     A(1,1) = (z+a(2))/z;
%     A(2,1) = (z.*z.*a(2)+z)/z;


A = [1 0 0 0 0 0 0 0 0 0];
Ar = A;
    for h = 1:10
    jh = lagmatrix(Ar,1)';
    jh(isnan(jh)) = 0;
    A = A + gam(h+1).*jh
    Ar = fliplr(A);
    end
 %e = conv2(mess,A);
    e = zeros(1,92);

   % if(h == 1)
%     hjk = lagmatrix(mess',1)';
%     hjk(isnan(hjk)) = 0;
%     e = mess' + gam(1)*hjk;
%     er = hjk - mess'*gam(1);
%     mess = e';
%     h = 2;
%     e = mess' + gam(h)*hjk;
%     er = hjk - mess'*gam(h);

%     end
%     hjk = lagmatrix(er,1)';
%      hjk(isnan(hjk)) = 0;
%      er = hjk - gam(h)*e;
%      e =  e - gam(h)*hjk;    

    [c,g] = lpc(mess,10);
    for h=1:10
    hjk = lagmatrix(mess',h-1)';
    hjk(isnan(hjk)) = 0;
    e = e + c(h)*hjk;
    mess = e';
    end
    
    
    if(n == 1)
    error = 0;
    end
    error = [error e];
end
 sound(error,44100);       
        
        
        
        
