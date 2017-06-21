function [mx y] = getMxfreq(yf, n)

N = length(yf);

ypow = abs(yf/N);

ypow2 = ypow((N/2):N);

for i=1:n
   [mx(i) y(i)] = max(ypow2);
   ypow2(y(i)) = -Inf;
end


