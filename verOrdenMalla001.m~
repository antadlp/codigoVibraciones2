filename = 'ordenMalla001.dat';
Zfile = importdata(filename);
zz = Zfile;

x = 1:23;
y = 1:8;

%x columnas
%y renglones

zz2 = reshape(zz, [23, 8]);

[X, Y] = meshgrid(x,y);

X = X';
Y = Y';


figure('Name', 'mesh')
mesh(X, Y, zz2);
axis([-inf inf -inf inf -2 2])

figure('Name', 'surface')
surface(X, Y, zz2, 'EdgeColor', 'none'), view(3)
axis([-inf inf -inf inf -2 2])


fs = 10;
[v, w, vv, ww, Pf] = analisisfft002(zz2, X, Y,  'malla500');


