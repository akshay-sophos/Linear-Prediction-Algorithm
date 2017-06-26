[x,Fs] = audioread('QASK.m4a');
sound(x,44100);
eps = 0.000001;
theta = [0 : 0.0001 : 2*pi];
z = exp(1i*theta);
%segment the audio signal
for(n = 1: 92 : 44008)
    temp = x(n:n+91,1);
    temp = [temp ; zeros(92,1)];
    for L = 0:92%start of correlation
        for j = 1 : 92
            Rxx(j) = (temp(j)*temp(j+L));             
        end
        Rxx_sum(L+1) = sum(Rxx);
        Rxx = 0;
    end% of correlation
      
    a = zeros(1,11);
    a(0+1) = 1;
    a(1+1) = -(Rxx_sum(2)/Rxx_sum(1));
    gam(1) = -a(1+1);
    e_min(1) = Rxx_sum(1)*(1-a(1+1)^2);
    g = zeros(9,25);
    %g(1,:)
    for(x = -10 : 14)
        for(q = 0 : 2)
            g(1,x+10+1) = g(1,x+10+1) + Rxx_sum(abs(q-x)+1)*a(1+q);
        end%g(1,k)
    end
    stem([-10:14],g(1,:));
    gam(2) = -g(1,10+2)/g(1,10);
    e_min(2) = e_min(1).*(1-gam(1)^2);
    
    for(m = 2 : 10)
        for(k = -10 : 14)
            K = fliplr(lagmatrix(g(1,:),m)');
            K(isnan(K)) = 0;
            g(m,10+k+1) = g(m-1,10+k+1) + gam(m).*K;
        end
        gam(m+1) = -g(m,10+m+2)/g(m,10+1);
        e_min(m+1) = e_min(m).*(1-gam(m)^2);  
    end
    
end
        
        
        
        
        
