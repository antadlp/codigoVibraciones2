function distancia = distanciaPuntoLinea002(xlinea, ylinea, P)


Rz = rotz(90);

M = length(xlinea);
N = length(ylinea);

for i=1:M
      
   Ai = [xlinea(i); ylinea(i); 0];
   Bi = Rz*Ai;
   xR(i) = Bi(1);
   yR(i) = Bi(2);

end

P0 = [P(1); P(2); 0];
Pr = Rz*P0;


x = xR;
y = yR;

m = (yR(M) - yR(1))/(xR(M) - x(1));

b = yR(1) -m*xR(1);

A = -m;
B = 1;
C = -b;

distancia = abs(A*Pr(1) + B*Pr(2) + C)/sqrt(A^2 + B^2);

yP = m*Pr(1) + b;

if (yP > Pr(2))

   distancia = -distancia;

end


