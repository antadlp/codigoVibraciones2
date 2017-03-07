function [x2 y2 Pf2 Pf] = analisisfft003(fs, Z1, X1, Y1, S)

[M N] = size(Z1);

Z1f = fft2(Z1);
Z1sh = fftshift(Z1f);
Z1P = Z1sh.*conj(Z1sh);
x2 = (-M/2:M/2-1)*(fs/M);
y2 = (-N/2:N/2-1)*(fs/N);


numFreq = 10;
dummie = Z1P;
for i=1:numFreq
   
   [aa bb cc] = maxValueMatrix2(dummie);
   Pf(i,:) = [x2(bb) y2(cc)];
   Pf2(i) = aa;
   dummie(bb, cc) = -inf;

end














