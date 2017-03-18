function y = maualFBT001(fw, w)

N = length(fw);
M = length(w);


for n=1:N
   ft(n) = 0;
   for k=1:M
      ft(n) = ft(n) + fw(w(k))*exp(j*2*pi*(w(k)-1)*(n-1)/N);
   end
end

y = (1/())*ft;

