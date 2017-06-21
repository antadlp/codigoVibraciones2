% malla de prueba para transformada
%


x = -4*pi:1:4*pi;
y = x;
[X, Y] = meshgrid(x, y);
Z = (0.3)*sin(X);
figure(1), surface(X, Y, Z), view(3)

x = -2:0.2:2;
y = x;
[X, Y] = meshgrid(x, y);
Z = X.*exp(-X.^2 - Y.^2);
figure(2), surface(X, Y, Z), view(3)



