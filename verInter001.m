[X Y] = meshgrid(-3:3);
V = peaks(X, Y);

figure
surf(X, Y, V)

x = [-3 -1 -0.5 0 1 2];
y = [-4 -1  0 1 2 3]; 
[X1 Y1] = meshgrid(x,y);
V1 = peaks(X1, Y1);

figure
surf(X1, Y1, V1)

[Xq Yq] = meshgrid(-3:0.25:3);
Vq = interp2(X, Y, V, Xq, Yq);
figure
surf(Xq, Yq, Vq)

Vq1 = interp2(X1, Y1, V1, Xq, Yq);
figure
surf(Xq, Yq, Vq1)

x = -2*pi:0.3:2*pi;
y = x;

[X Y] = meshgrid(x, y);

Z = sin(X) + cos(Y);

figure
surf(X, Y, Z)







