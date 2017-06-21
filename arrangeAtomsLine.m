function arr = arrangeAtomsLine(xlinea, ylinea, A, B, dx)
%function distancia = distanciaPuntoLinea001(xlinea, ylinea, P, dx)


for i=1:length(B)

   P = [A(B(i), 1) A(B(i), 2)];
   d(i) = distanciaPuntoLinea001(xlinea, ylinea, P);

end

ymin = min(d);
ymax = max(d);

ypMin = min(A(B, 2));

yArr = ypMin + (ymax - ymin)/2;

for i=1:length(B)

   arr(i, :) = [A(B(i), 1) yArr];

end




