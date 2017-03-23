function y = manualFFT001(ft)

N = length(ft);



for k=1:N
   fw(k) = 0;
   for n=1:N
      fw(k) = fw(k) + ft(n)*exp(-j*2*pi*(k-1)*(n-1)/N);
   end
end

y = fw;


