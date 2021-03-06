function y = fullFBT001(fw)

%Ay = Area y

N = length(fw);

for n=1:N
   ft(n) = 0;
   for k=1:N
      ft(n) = ft(n) + fw(k)*exp(j*2*pi*(k-1)*(n-1)/N);
   end
end

y = (1/N)*ft;

