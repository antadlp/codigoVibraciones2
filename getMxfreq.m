function [mx y] = getMxfreq(yf, n)

N = length(yf)

ypow = abs(yf/N);

for i=1:n
   [mx(i) y(i)] = max(ypow);
   ypow(y(i)) = -Inf;
end


