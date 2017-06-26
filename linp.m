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
    for(b = -10 : 14)
        for(q = 0 : 2)
            g(1,b+10+1) = g(1,b+10+1) + Rxx_sum(abs(q-b)+1)*a(1+q);
        end%g(1,k)
    end
    gam(2) = -g(1,10+3)/g(1,10+1);
    e_min(2) = e_min(1).*(1-gam(1)^2);
   
    for(m = 2 : 10)
        K = fliplr(lagmatrix(g(m-1,:),m-(2*(m-2)))');
        K(isnan(K)) = 0;
        g(m,:) = g(m-1,:) + gam(m).*K;
        gam(m+1) = -g(m,10+m+2)/g(m,10+1);
        e_min(m+1) = e_min(m).*(1-gam(m)^2);  
    end
    if(n > 93+92+92*5)
        plot(e_min)
    end
end
        
        
        
        
        
