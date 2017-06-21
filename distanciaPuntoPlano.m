function y = distanciaPuntoPlano(P, N)


d = N(4);

Zp = (-N(1)*P(1) - N(2)*P(2) -d)/N(3);

y = (abs(N(1)*P(1) + N(2)*P(2) ...
+ N(3)*P(3) + d))/(sqrt(N(1)^2 + N(2)^2 + N(3)^2));

if (Zp > P(3))

   y = -y;

end



