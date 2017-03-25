function y = manualFBT001(ft, w)

N = length(ft);


for k=1:N
   fw(k)=0;
   for n=1:length(w)
      fw(k)= fw(k) + ft(w(n))*exp(j*2*pi*(k-1)*(n-1)/N);
   end
end




y = (1/N)*fw;

