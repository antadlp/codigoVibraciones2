function lista = crossLine2dAtom(Atoms, xline, yline, radii, ep, orden)

%orden 1 = horizontal 
%orden 2 = vertical 

dx = abs(radii/ep);

xmin = min(xline);
xmax = max(xline);
ymin = min(yline);
ymax = max(yline);


x = xmin:dx:xmax;
m = (yline(length(yline)) - yline(1))/(xline(length(xline)) - xline(1));

y = m*(x - xline(1)) + yline(1);

%y = m*(x - P1(1)) + P1(2); 

A = Atoms;
sz = size(Atoms);
r = radii;

l = 1;
for i=1:sz(1)
   for j=1:length(x)

      xA = A(i, 1);
      yA = A(i, 2);

      rp = (x(j) - xA)^2 + (y(j) - yA)^2;

      if ( rp <= radii^2)

         lista(l) = i;
         l = l + 1;

      end

   end
end

lista = unique(lista);

if (orden == 1)
   
   B = A(lista, 1);

elseif (orden == 2)

   B = A(lista, 2);

end

[C D] = sort(B);

F = lista(D);

lista = F;





      




