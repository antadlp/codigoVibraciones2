function arr = arrangeAtomsLineY(xlinea, ylinea, A, B, C)
%function distancia = distanciaPuntoLinea001(xlinea, ylinea, P, dx)


for i=1:length(B)

   P = [A(B(i), 1) A(B(i), 2)];
   d(i) = distanciaPuntoLinea002(xlinea, ylinea, P);

end

xmin = min(d);
xmax = max(d);

xpMin = min(A(B, 1));

xArr = xpMin + (xmax - xmin)/2;

for i=1:length(B)

   [ii jj] = ismember(B(i), C);

   if (ii == 1)

       arr(i, :) = [xArr C(jj, 3)];

   elseif (ii == 0)

       arr(i, :) = [xArr A(B(i), 2)];

    end


end


%para alinearlos horizontal
%for i=1:length(B)
%
%   arr(i, :) = [A(B(i), 1) yArr];
%
%end
%




