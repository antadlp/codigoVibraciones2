function distancia = distanciaPuntoLinea001(xlinea, ylinea, P, dx)


x = xlinea;
y = ylinea;

m = (y(length(y)) - y(1))/(x(length(x)) - x(1));

b = y(1) -m*x(1);

A = -m;
B = 1;
C = -b;



distancia = abs(A*P(1) + B*P(2) + C)/sqrt(A^2 + B^2);


yP = m*P(1) + b;

if (yP > P(2))

   distancia = -distancia;

end


