function gridXY001

% Dado que las posiciones de la malla del grafeno forman un arreglo hexagonal y
% es necesaria una malla rectangular para la interpolacion 2d, este script
% ademas de alinear las posiciones, interpola en los huecos donde deberia ir
% una posicion del atomo, transformando la malla hexagonal del grafeno en una
% malla rectangular como se necesita.

close all
clear all
filename = 'frame1-BN-rota.dat'; % Primer frame del movie a analizar
XYZFile = importdata(filename);

format short e 

MA = XYZFile;
X = XYZFile(:,1);
Y = XYZFile(:,2);
Z = XYZFile(:,3);

fs = 10;

X = floor(X*fs)/fs;
Y = floor(Y*fs)/fs;

xmax = max(X);
xmin= min(X);
x = xmin:1/fs:xmax;

ymax = max(Y);
ymin= min(Y);
y = ymin:1/fs:ymax;

fileIDX = fopen('verX2.dat', 'w');
fileIDY = fopen('verY2.dat', 'w');

Zl = zeros(length(x), length(y));
for i=1:length(x)
   for j=1:length(y)

      Zl(i,j) = NaN;

   end
end

for i=1:length(X)
   
   fprintf(fileIDX, '\nVERIFICANDO %f\n', X(i));

   xCheck = abs(x - X(i));
   [xx I(i)] = min(xCheck);
     
   fprintf(fileIDX, '%f\t%f\n\n', X(i), x(I(i)));
   
   yCheck = abs(y - Y(i));
   [yy J(i)] = min(yCheck);

   fprintf(fileIDY, '\nVERIFICANDO %f\n', Y(i));
   fprintf(fileIDY, '%f\t%f\n\n', Y(i), y(J(i)));

   Zl(I(i), J(i)) = Z(i);
   z2(i,3) = Zl(I(i), J(i));
   z2(i,1) = x(I(i));
   z2(i,2) = y(J(i));

end
fclose(fileIDX);
fclose(fileIDY);



r = 0.3;
ZrC = zeros(length(x), length(y));
for i=1:length(x)
   for j=1:length(y)

      ZrC(i,j) = NaN;

   end
end



for l=1:length(X)
   for i=1:length(x)
      for j=1:length(y)

         ecC = (x(i) - x(I(l)))^2 + (y(j) - y(J(l)))^2;

         if (ecC <= r^2)
        
            Zr(i,j,l) = sqrt(r^2 - ecC);
            ZrC(i,j)= sqrt(r^2 - ecC) + Zl(I(l),J(l)); 

         else

            Zr(i,j,l) = NaN;

         end

      end
   end
end

D = importdata('frame1-mil-rota.dat'); %no cambiar nombre -mil!!
nAtoms = 142;
fileMap0 = fopen('map0-BN.dat','w');
for i=1:nAtoms

   F = abs(X-D(i,1));
   G = abs(Y-D(i,2));

   H = abs(F + G);

   [I J] = min(H);

   fprintf(fileMap0, '%f  %f\n', i, J);
   m0(i,:) = [i, J];

end
fclose(fileMap0)




H11 = [MA(m0(11,2), 1) MA(m0(11,2), 2) MA(m0(11,2), 3)];
H07 = [MA(m0(7,2), 1) MA(m0(7,2), 2) MA(m0(7,2), 3)];
H22 = [MA(m0(22,2), 1) MA(m0(22,2), 2) MA(m0(22,2), 3)];
H04 = [MA(m0(4,2), 1) MA(m0(4,2), 2) MA(m0(4,2), 3)];
H15 = [MA(m0(15,2), 1) MA(m0(15,2), 2) MA(m0(15,2), 3)];
H14 = [MA(m0(14,2), 1) MA(m0(14,2), 2) MA(m0(14,2), 3)];
H01 = [MA(m0(1,2), 1) MA(m0(1,2), 2) MA(m0(1,2), 3)];
H05 = [MA(m0(5,2), 1) MA(m0(5,2), 2) MA(m0(5,2), 3)];
H69 = [MA(m0(69,2), 1) MA(m0(69,2), 2) MA(m0(69,2), 3)];
H47 = [MA(m0(47,2), 1) MA(m0(47,2), 2) MA(m0(47,2), 3)];
H91 = [MA(m0(91,2), 1) MA(m0(91,2), 2) MA(m0(91,2), 3)];
H48 = [MA(m0(48,2), 1) MA(m0(48,2), 2) MA(m0(48,2), 3)];
C88 = [MA(m0(88,2), 1) MA(m0(88,2), 2) MA(m0(88,2), 3)];
C42 = [MA(m0(42,2), 1) MA(m0(42,2), 2) MA(m0(42,2), 3)];
C111 = [MA(m0(111,2), 1) MA(m0(111,2), 2) MA(m0(111,2), 3)];
C43 = [MA(m0(43,2), 1) MA(m0(43,2), 2) MA(m0(43,2), 3)];
C108 = [MA(m0(108,2), 1) MA(m0(108,2), 2) MA(m0(108,2), 3)];
C44 = [MA(m0(44,2), 1) MA(m0(44,2), 2) MA(m0(44,2), 3)];
C124 = [MA(m0(124,2), 1) MA(m0(124,2), 2) MA(m0(124,2), 3)];
C45 = [MA(m0(45,2), 1) MA(m0(45,2), 2) MA(m0(45,2), 3)];
C122 = [MA(m0(122,2), 1) MA(m0(122,2), 2) MA(m0(122,2), 3)];
C46 = [MA(m0(46,2), 1) MA(m0(46,2), 2) MA(m0(46,2), 3)];
C120 = [MA(m0(120,2), 1) MA(m0(120,2), 2) MA(m0(120,2), 3)];
C40 = [MA(m0(40,2), 1) MA(m0(40,2), 2) MA(m0(40,2), 3)];
C117 = [MA(m0(117,2), 1) MA(m0(117,2), 2) MA(m0(117,2), 3)];
C41 = [MA(m0(41,2), 1) MA(m0(41,2), 2) MA(m0(41,2), 3)];
C114 = [MA(m0(114,2), 1) MA(m0(114,2), 2) MA(m0(114,2), 3)];
C34 = [MA(m0(34,2), 1) MA(m0(34,2), 2) MA(m0(34,2), 3)];
C134 = [MA(m0(134,2), 1) MA(m0(134,2), 2) MA(m0(134,2), 3)];
C35 = [MA(m0(35,2), 1) MA(m0(35,2), 2) MA(m0(35,2), 3)];
C31 = [MA(m0(31,2), 1) MA(m0(31,2), 2) MA(m0(31,2), 3)];
C132 = [MA(m0(132,2), 1) MA(m0(132,2), 2) MA(m0(132,2), 3)];
C66 = [MA(m0(66,2), 1) MA(m0(66,2), 2) MA(m0(66,2), 3)];
C129 = [MA(m0(129,2), 1) MA(m0(129,2), 2) MA(m0(129,2), 3)];
C54 = [MA(m0(54,2), 1) MA(m0(54,2), 2) MA(m0(54,2), 3)];
C128 = [MA(m0(128,2), 1) MA(m0(128,2), 2) MA(m0(128,2), 3)];
C126 = [MA(m0(126,2), 1) MA(m0(126,2), 2) MA(m0(126,2), 3)];
C52 = [MA(m0(52,2), 1) MA(m0(52,2), 2) MA(m0(52,2), 3)];
C141 = [MA(m0(141,2), 1) MA(m0(141,2), 2) MA(m0(141,2), 3)];
C125 = [MA(m0(125,2), 1) MA(m0(125,2), 2) MA(m0(125,2), 3)];
C136 = [MA(m0(136,2), 1) MA(m0(136,2), 2) MA(m0(136,2), 3)];
C115 = [MA(m0(115,2), 1) MA(m0(115,2), 2) MA(m0(115,2), 3)];
C140 = [MA(m0(140,2), 1) MA(m0(140,2), 2) MA(m0(140,2), 3)];
C112 = [MA(m0(112,2), 1) MA(m0(112,2), 2) MA(m0(112,2), 3)];
C94 = [MA(m0(94,2), 1) MA(m0(94,2), 2) MA(m0(94,2), 3)];
C137 = [MA(m0(137,2), 1) MA(m0(137,2), 2) MA(m0(137,2), 3)];
C139 = [MA(m0(139,2), 1) MA(m0(139,2), 2) MA(m0(139,2), 3)];
C92 = [MA(m0(92,2), 1) MA(m0(92,2), 2) MA(m0(92,2), 3)];
C138 = [MA(m0(138,2), 1) MA(m0(138,2), 2) MA(m0(138,2), 3)];
C70 = [MA(m0(70,2), 1) MA(m0(70,2), 2) MA(m0(70,2), 3)];
C142 = [MA(m0(142,2), 1) MA(m0(142,2), 2) MA(m0(142,2), 3)];
C135 = [MA(m0(135,2), 1) MA(m0(135,2), 2) MA(m0(135,2), 3)];
C90 = [MA(m0(90,2), 1) MA(m0(90,2), 2) MA(m0(90,2), 3)];
C133 = [MA(m0(133,2), 1) MA(m0(133,2), 2) MA(m0(133,2), 3)];
C68 = [MA(m0(68,2), 1) MA(m0(68,2), 2) MA(m0(68,2), 3)];
C131 = [MA(m0(131,2), 1) MA(m0(131,2), 2) MA(m0(131,2), 3)];
C65 = [MA(m0(65,2), 1) MA(m0(65,2), 2) MA(m0(65,2), 3)];
C30 = [MA(m0(30,2), 1) MA(m0(30,2), 2) MA(m0(30,2), 3)];
C47 = [MA(m0(47,2), 1) MA(m0(47,2), 2) MA(m0(47,2), 3)];
C31 = [MA(m0(31,2), 1) MA(m0(31,2), 2) MA(m0(31,2), 3)];
C48 = [MA(m0(48,2), 1) MA(m0(48,2), 2) MA(m0(48,2), 3)];
C35 = [MA(m0(35,2), 1) MA(m0(35,2), 2) MA(m0(35,2), 3)];


[xY01 lineY01 ] = makeAtomsLine(H11, H07);
[xY02 lineY02 ] = makeAtomsLine(H15, H14);


extraLimitA = abs(H11(1) - H22(1));
extraLimitB = abs(H04(1) - H15(1));

extraLimit = max([extraLimitA extraLimitB]);
extraLimit2 = extraLimit/2;

limitXIzq = xmin - extraLimit2;
limitXDer = xmax + extraLimit2;

limitYArr = ymax + extraLimit2;
limitYAbj = ymin - extraLimit2;

%P = H01;
%Q = H22;
%R = H11;
%PQ = P - Q;
%PR = P - R;
%normal = cross(PQ, PR);

P = H15;
Q = H14;
R = H04;
PQ = P - Q;
PR = P - R;
normal = cross(PQ, PR);

xp = limitXIzq:01:limitXDer;
yp = limitYAbj:01:limitYArr;

[Xp Yp] = meshgrid(xp, yp);
[Xl Yl] = meshgrid(x, y);
[Zp N] = makePlano001(Xp, Yp, H15, normal');
%Zd = distanciaPuntoPlano(P, N)



r = .3;
ep = 20;
offset = 0;

sHor = '/home/toshiba/gridXY-BN/horizontales-BN-';
sVert = '/home/toshiba/gridXY-BN/verticales-BN-';
[xX02 lineaX02] = makeAtomsLine(H69, H47);
listaHor02 = crossLine2dAtom(MA, xX02, lineaX02, r, ep, 1);
arrHor02 = arrangeAtomsLine(xX02, lineaX02, MA, listaHor02, 0.001);
namefile = strcat(sHor, '2');
fileID = fopen(namefile, 'w');
for i=1:length(listaHor02)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrHor02(i,1), arrHor02(i, 2), listaHor02(i));
end
fprintf(fileID, '\n');
fclose(fileID);
horizontalesL(2, 1) = arrHor02(1,2);
for i=1:length(arrHor02) 

   nMA(i + offset, :) = [listaHor02(i) arrHor02(i, 1) arrHor02(i, 2)];

end
offset = length(nMA);
A1 = [arrHor02(1, 1) arrHor02(1, 2)];
A2 = [arrHor02(length(arrHor02), 1) arrHor02(length(arrHor02), 2)];
[xX02t lineaX02t] = makeAtomsLine(A1, A2);
plot(xX02t, lineaX02t)
%plot(nMA(:,2), nMA(:,3), 'ro')
hold on
axis([limitXIzq limitXDer limitYAbj limitYArr]);
grid on




[xX03 lineaX03] = makeAtomsLine(H91, H48);
listaHor03 = crossLine2dAtom(MA, xX03, lineaX03, r, ep, 1);
arrHor03 = arrangeAtomsLine(xX03, lineaX03, MA, listaHor03, 0.001);

namefile = strcat(sHor, '3');
fileID = fopen(namefile, 'w');
for i=1:length(listaHor03)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrHor03(i,1), arrHor03(i, 2), listaHor03(i));
end
fprintf(fileID, '\n');
fclose(fileID);
horizontalesL(3, 1) = arrHor03(1,2);

for i=1:length(arrHor03) 

   nMA(i + offset, :) = [listaHor03(i) arrHor03(i, 1) arrHor03(i, 2)];

end
offset = length(nMA);
%plot(nMA(:,2), nMA(:,3), 'ro')
A1 = [arrHor03(1, 1) arrHor03(1, 2)];
A2 = [arrHor03(length(arrHor03), 1) arrHor03(length(arrHor03), 2)];
[xX03t lineaX03t] = makeAtomsLine(A1, A2);
plot(xX03t, lineaX03t)


[xX04 lineaX04] = makeAtomsLine(C88, C42);
listaHor04 = crossLine2dAtom(MA, xX04, lineaX04, r, ep, 1);
arrHor04 = arrangeAtomsLine(xX04, lineaX04, MA, listaHor04, 0.001);

namefile = strcat(sHor, '4');
fileID = fopen(namefile, 'w');
for i=1:length(listaHor04)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrHor04(i,1), arrHor04(i, 2), listaHor04(i));
end
fprintf(fileID, '\n');
fclose(fileID);
horizontalesL(4, 1) = arrHor04(1,2);

for i=1:length(arrHor04) 

   nMA(i + offset, :) = [listaHor04(i) arrHor04(i, 1) arrHor04(i, 2)];

end
offset = length(nMA);
%plot(nMA(:,2), nMA(:,3), 'ro')
A1 = [arrHor04(1, 1) arrHor04(1, 2)];
A2 = [arrHor04(length(arrHor04), 1) arrHor04(length(arrHor04), 2)];
[xX04t lineaX04t] = makeAtomsLine(A1, A2);
plot(xX04t, lineaX04t)



[xX05 lineaX05] = makeAtomsLine(C111, C43);
listaHor05 = crossLine2dAtom(MA, xX05, lineaX05, r, ep, 1);
arrHor05 = arrangeAtomsLine(xX05, lineaX05, MA, listaHor05, 0.001);

namefile = strcat(sHor, '5');
fileID = fopen(namefile, 'w');
for i=1:length(listaHor05)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrHor05(i,1), arrHor05(i, 2), listaHor05(i));
end
fprintf(fileID, '\n');
fclose(fileID);
horizontalesL(5, 1) = arrHor05(1,2);

for i=1:length(arrHor05) 

   nMA(i + offset, :) = [listaHor05(i) arrHor05(i, 1) arrHor05(i, 2)];

end
offset = length(nMA);
%plot(nMA(:,2), nMA(:,3), 'ro')
A1 = [arrHor05(1, 1) arrHor05(1, 2)];
A2 = [arrHor05(length(arrHor05), 1) arrHor05(length(arrHor05), 2)];
[xX05t lineaX05t] = makeAtomsLine(A1, A2);
plot(xX05t, lineaX05t)




[xX06 lineaX06] = makeAtomsLine(C108, C44);
listaHor06 = crossLine2dAtom(MA, xX06, lineaX06, r, ep, 1);
arrHor06 = arrangeAtomsLine(xX06, lineaX06, MA, listaHor06, 0.001);

namefile = strcat(sHor, '6');
fileID = fopen(namefile, 'w');
for i=1:length(listaHor06)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrHor06(i,1), arrHor06(i, 2), listaHor06(i));
end
fprintf(fileID, '\n');
fclose(fileID);
horizontalesL(6, 1) = arrHor06(1,2);

horizontales( :, :, 6) =  [arrHor06(:,1) arrHor06(:, 2) MA(listaHor06, 3)];
for i=1:length(arrHor06) 

   nMA(i + offset, :) = [listaHor06(i) arrHor06(i, 1) arrHor06(i, 2)];

end
offset = length(nMA);
%plot(nMA(:,2), nMA(:,3), 'ro')
A1 = [arrHor06(1, 1) arrHor06(1, 2)];
A2 = [arrHor06(length(arrHor06), 1) arrHor06(length(arrHor06), 2)];
[xX06t lineaX06t] = makeAtomsLine(A1, A2);
plot(xX06t, lineaX06t)




[xX07 lineaX07] = makeAtomsLine(C124, C45);
listaHor07 = crossLine2dAtom(MA, xX07, lineaX07, r, ep, 1);
arrHor07 = arrangeAtomsLine(xX07, lineaX07, MA, listaHor07, 0.001);

namefile = strcat(sHor, '7');
fileID = fopen(namefile, 'w');
for i=1:length(listaHor07)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrHor07(i,1), arrHor07(i, 2), listaHor07(i));
end
fprintf(fileID, '\n');
fclose(fileID);
horizontalesL(7, 1) = arrHor07(1,2);

for i=1:length(arrHor07) 

   nMA(i + offset, :) = [listaHor07(i) arrHor07(i, 1) arrHor07(i, 2)];

end
offset = length(nMA);
%plot(nMA(:,2), nMA(:,3), 'ro')
A1 = [arrHor07(1, 1) arrHor07(1, 2)];
A2 = [arrHor07(length(arrHor07), 1) arrHor07(length(arrHor07), 2)];
[xX07t lineaX07t] = makeAtomsLine(A1, A2);
plot(xX07t, lineaX07t)




[xX08 lineaX08] = makeAtomsLine(C122, C46);
listaHor08 = crossLine2dAtom(MA, xX08, lineaX08, r, ep, 1);
arrHor08 = arrangeAtomsLine(xX08, lineaX08, MA, listaHor08, 0.001);

namefile = strcat(sHor, '8');
fileID = fopen(namefile, 'w');
for i=1:length(listaHor08)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrHor08(i,1), arrHor08(i, 2), listaHor08(i));
end
fprintf(fileID, '\n');
fclose(fileID);
horizontalesL(8, 1) = arrHor08(1,2);

for i=1:length(arrHor08) 

   nMA(i + offset, :) = [listaHor08(i) arrHor08(i, 1) arrHor08(i, 2)];

end
offset = length(nMA);
%plot(nMA(:,2), nMA(:,3), 'ro')
A1 = [arrHor08(1, 1) arrHor08(1, 2)];
A2 = [arrHor08(length(arrHor08), 1) arrHor08(length(arrHor08), 2)];
[xX08t lineaX08t] = makeAtomsLine(A1, A2);
plot(xX08t, lineaX08t)





[xX09 lineaX09] = makeAtomsLine(C120, C40);
listaHor09 = crossLine2dAtom(MA, xX09, lineaX09, r, ep, 1);
arrHor09 = arrangeAtomsLine(xX09, lineaX09, MA, listaHor09, 0.001);

namefile = strcat(sHor, '9');
fileID = fopen(namefile, 'w');
for i=1:length(listaHor09)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrHor09(i,1), arrHor09(i, 2), listaHor09(i));
end
fprintf(fileID, '\n');
fclose(fileID);
horizontalesL(9, 1) = arrHor09(1,2);

for i=1:length(arrHor09) 

   nMA(i + offset, :) = [listaHor09(i) arrHor09(i, 1) arrHor09(i, 2)];

end
offset = length(nMA);
%plot(nMA(:,2), nMA(:,3), 'ro')
A1 = [arrHor09(1, 1) arrHor09(1, 2)];
A2 = [arrHor09(length(arrHor09), 1) arrHor09(length(arrHor09), 2)];
[xX09t lineaX09t] = makeAtomsLine(A1, A2);
plot(xX09t, lineaX09t)



[xX10 lineaX10] = makeAtomsLine(C117, C41);
listaHor10 = crossLine2dAtom(MA, xX10, lineaX10, r, ep, 1);
arrHor10 = arrangeAtomsLine(xX10, lineaX10, MA, listaHor10, 0.001);

namefile = strcat(sHor, '10');
fileID = fopen(namefile, 'w');
for i=1:length(listaHor10)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrHor10(i,1), arrHor10(i, 2), listaHor10(i));
end
fprintf(fileID, '\n');
fclose(fileID);
horizontalesL(10, 1) = arrHor10(1,2);

for i=1:length(arrHor10) 

   nMA(i + offset, :) = [listaHor10(i) arrHor10(i, 1) arrHor10(i, 2)];

end
offset = length(nMA);
%plot(nMA(:,2), nMA(:,3), 'ro')
A1 = [arrHor10(1, 1) arrHor10(1, 2)];
A2 = [arrHor10(length(arrHor10), 1) arrHor10(length(arrHor10), 2)];
[xX10t lineaX10t] = makeAtomsLine(A1, A2);
plot(xX10t, lineaX10t)




[xX11 lineaX11] = makeAtomsLine(C114, C34);
listaHor11 = crossLine2dAtom(MA, xX11, lineaX11, r, ep, 1);
arrHor11 = arrangeAtomsLine(xX11, lineaX11, MA, listaHor11, 0.001);

namefile = strcat(sHor, '11');
fileID = fopen(namefile, 'w');
for i=1:length(listaHor11)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrHor11(i,1), arrHor11(i, 2), listaHor11(i));
end
fprintf(fileID, '\n');
fclose(fileID);
horizontalesL(11, 1) = arrHor11(1,2);

for i=1:length(arrHor11) 

   nMA(i + offset, :) = [listaHor11(i) arrHor11(i, 1) arrHor11(i, 2)];

end
offset = length(nMA);
%plot(nMA(:,2), nMA(:,3), 'ro')
A1 = [arrHor11(1, 1) arrHor11(1, 2)];
A2 = [arrHor11(length(arrHor11), 1) arrHor11(length(arrHor11), 2)];
[xX11t lineaX11t] = makeAtomsLine(A1, A2);
plot(xX11t, lineaX11t)





[xX12 lineaX12] = makeAtomsLine(C134, C35);
listaHor12 = crossLine2dAtom(MA, xX12, lineaX12, r, ep, 1);
arrHor12 = arrangeAtomsLine(xX12, lineaX12, MA, listaHor12, 0.001);

namefile = strcat(sHor, '12');
fileID = fopen(namefile, 'w');
for i=1:length(listaHor12)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrHor12(i,1), arrHor12(i, 2), listaHor12(i));
end
fprintf(fileID, '\n');
fclose(fileID);
horizontalesL(12, 1) = arrHor12(1,2);

for i=1:length(arrHor12) 

   nMA(i + offset, :) = [listaHor12(i) arrHor12(i, 1) arrHor12(i, 2)];

end
offset = length(nMA);
%plot(nMA(:,2), nMA(:,3), 'ro')
A1 = [arrHor12(1, 1) arrHor12(1, 2)];
A2 = [arrHor12(length(arrHor12), 1) arrHor12(length(arrHor12), 2)];
[xX12t lineaX12t] = makeAtomsLine(A1, A2);
plot(xX12t, lineaX12t)



[xX13 lineaX13] = makeAtomsLine(C132, C31);
listaHor13 = crossLine2dAtom(MA, xX13, lineaX13, r, ep, 1);
arrHor13 = arrangeAtomsLine(xX13, lineaX13, MA, listaHor13, 0.001);

namefile = strcat(sHor, '13');
fileID = fopen(namefile, 'w');
for i=1:length(listaHor13)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrHor13(i,1), arrHor13(i, 2), listaHor13(i));
end
fprintf(fileID, '\n');
fclose(fileID);
horizontalesL(13, 1) = arrHor13(1,2);

for i=1:length(arrHor13) 

   nMA(i + offset, :) = [listaHor13(i) arrHor13(i, 1) arrHor13(i, 2)];

end
offset = length(nMA);
%plot(nMA(:,2), nMA(:,3), 'ro')
A1 = [arrHor13(1, 1) arrHor13(1, 2)];
A2 = [arrHor13(length(arrHor13), 1) arrHor13(length(arrHor13), 2)];
[xX13t lineaX13t] = makeAtomsLine(A1, A2);
plot(xX13t, lineaX13t)





[xY01 lineaY01] = makeAtomsLineY(C124, C122);
listaVert01 = crossLine2dAtom(MA, xY01, lineaY01, r, ep, 2);
arrVert01 = arrangeAtomsLineY(xY01, lineaY01, MA, listaVert01, nMA);

namefile = strcat(sVert, '1');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert01)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert01(i,1), arrVert01(i, 2), listaVert01(i));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(1, 1) = arrVert01(1,1);

for i=1:length(arrVert01) 

   [ii jj] = ismember(listaVert01(i), nMA);

   nMA(jj, :) = [listaVert01(i) arrVert01(i, 1) arrVert01(i, 2)];

end
%plot(nMA(:,2), nMA(:,3), 'ro')
A1 = [arrVert01(1, 1) arrVert01(1, 2)];
A2 = [arrVert01(length(arrVert01), 1) arrVert01(length(arrVert01), 2)];
[xY01t lineaY01t] = makeAtomsLineY(A1, A2);
plot(xY01t, lineaY01t)


[xY02 lineaY02] = makeAtomsLineY(C111, C117);
listaVert02 = crossLine2dAtom(MA, xY02, lineaY02, r, ep, 2);
arrVert02 = arrangeAtomsLineY(xY02, lineaY02, MA, listaVert02, nMA);

namefile = strcat(sVert, '2');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert02)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert02(i,1), arrVert02(i, 2), listaVert02(i));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(2, 1) = arrVert02(1,1);

for i=1:length(arrVert02) 

   [ii jj] = ismember(listaVert02(i), nMA);

   nMA(jj, :) = [listaVert02(i) arrVert02(i, 1) arrVert02(i, 2)];

end
%plot(nMA(:,2), nMA(:,3), 'ro')
A1 = [arrVert02(1, 1) arrVert02(1, 2)];
A2 = [arrVert02(length(arrVert02), 1) arrVert02(length(arrVert02), 2)];
[xY02t lineaY02t] = makeAtomsLineY(A1, A2);
plot(xY02t, lineaY02t)


[xY03 lineaY03] = makeAtomsLineY(H91, C134);
listaVert03 = crossLine2dAtom(MA, xY03, lineaY03, r, ep, 2);
arrVert03 = arrangeAtomsLineY(xY03, lineaY03, MA, listaVert03, nMA);

namefile = strcat(sVert, '3');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert03)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert03(i,1), arrVert03(i, 2), listaVert03(i));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(3, 1) = arrVert03(1,1);

for i=1:length(arrVert03) 

   [ii jj] = ismember(listaVert03(i), nMA);

   nMA(jj, :) = [listaVert03(i) arrVert03(i, 1) arrVert03(i, 2)];

end
%plot(nMA(:,2), nMA(:,3), 'ro')
A1 = [arrVert03(1, 1) arrVert03(1, 2)];
A2 = [arrVert03(length(arrVert03), 1) arrVert03(length(arrVert03), 2)];
[xY03t lineaY03t] = makeAtomsLineY(A1, A2);
plot(xY03t, lineaY03t)




[xY04 lineaY04] = makeAtomsLineY(H69, C132);
listaVert04 = crossLine2dAtom(MA, xY04, lineaY04, r, ep, 2);
arrVert04 = arrangeAtomsLineY(xY04, lineaY04, MA, listaVert04, nMA);

namefile = strcat(sVert, '4');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert04)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert04(i,1), arrVert04(i, 2), listaVert04(i));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(4, 1) = arrVert04(1,1);

for i=1:length(arrVert04) 

   [ii jj] = ismember(listaVert04(i), nMA);

   nMA(jj, :) = [listaVert04(i) arrVert04(i, 1) arrVert04(i, 2)];

end
%plot(nMA(:,2), nMA(:,3), 'ro')
A1 = [arrVert04(1, 1) arrVert04(1, 2)];
A2 = [arrVert04(length(arrVert04), 1) arrVert04(length(arrVert04), 2)];
[xY04t lineaY04t] = makeAtomsLineY(A1, A2);
plot(xY04t, lineaY04t)


[xY05 lineaY05] = makeAtomsLineY(C66, C129);
listaVert05 = crossLine2dAtom(MA, xY05, lineaY05, r, ep, 2);
arrVert05 = arrangeAtomsLineY(xY05, lineaY05, MA, listaVert05, nMA);

namefile = strcat(sVert, '5');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert05)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert05(i,1), arrVert05(i, 2), listaVert05(i));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(5, 1) = arrVert05(1,1);

for i=1:length(arrVert05) 

   [ii jj] = ismember(listaVert05(i), nMA);

   nMA(jj, :) = [listaVert05(i) arrVert05(i, 1) arrVert05(i, 2)];

end
%plot(nMA(:,2), nMA(:,3), 'ro')
A1 = [arrVert05(1, 1) arrVert05(1, 2)];
A2 = [arrVert05(length(arrVert05), 1) arrVert05(length(arrVert05), 2)];
[xY05t lineaY05t] = makeAtomsLineY(A1, A2);
plot(xY05t, lineaY05t)



[xY06 lineaY06] = makeAtomsLineY(C54, C128);
listaVert06 = crossLine2dAtom(MA, xY06, lineaY06, r, ep, 2);
arrVert06 = arrangeAtomsLineY(xY06, lineaY06, MA, listaVert06, nMA);

namefile = strcat(sVert, '6');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert06)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert06(i,1), arrVert06(i, 2), listaVert06(i));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(6, 1) = arrVert06(1,1);

for i=1:length(arrVert06) 

   [ii jj] = ismember(listaVert06(i), nMA);

   nMA(jj, :) = [listaVert06(i) arrVert06(i, 1) arrVert06(i, 2)];

end
%plot(nMA(:,2), nMA(:,3), 'ro')
A1 = [arrVert06(1, 1) arrVert06(1, 2)];
A2 = [arrVert06(length(arrVert06), 1) arrVert06(length(arrVert06), 2)];
[xY06t lineaY06t] = makeAtomsLineY(A1, A2);
plot(xY06t, lineaY06t)



[xY07 lineaY07] = makeAtomsLineY(C52, C126);
listaVert07 = crossLine2dAtom(MA, xY07, lineaY07, r, ep, 2);
arrVert07 = arrangeAtomsLineY(xY07, lineaY07, MA, listaVert07, nMA);

namefile = strcat(sVert, '7');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert07)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert07(i,1), arrVert07(i, 2), listaVert07(i));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(7, 1) = arrVert07(1,1);

for i=1:length(arrVert07) 

   [ii jj] = ismember(listaVert07(i), nMA);

   nMA(jj, :) = [listaVert07(i) arrVert07(i, 1) arrVert07(i, 2)];

end
%plot(nMA(:,2), nMA(:,3), 'ro')
A1 = [arrVert07(1, 1) arrVert07(1, 2)];
A2 = [arrVert07(length(arrVert07), 1) arrVert07(length(arrVert07), 2)];
[xY07t lineaY07t] = makeAtomsLineY(A1, A2);
plot(xY07t, lineaY07t)





[xY08 lineaY08] = makeAtomsLineY(C141, C125);
listaVert08 = crossLine2dAtom(MA, xY08, lineaY08, r, ep, 2);
arrVert08 = arrangeAtomsLineY(xY08, lineaY08, MA, listaVert08, nMA);

namefile = strcat(sVert, '8');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert08)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert08(i,1), arrVert08(i, 2), listaVert08(i));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(8, 1) = arrVert08(1,1);

for i=1:length(arrVert08) 

   [ii jj] = ismember(listaVert08(i), nMA);

   nMA(jj, :) = [listaVert08(i) arrVert08(i, 1) arrVert08(i, 2)];

end
%plot(nMA(:,2), nMA(:,3), 'ro')
A1 = [arrVert08(1, 1) arrVert08(1, 2)];
A2 = [arrVert08(length(arrVert08), 1) arrVert08(length(arrVert08), 2)];
[xY08t lineaY08t] = makeAtomsLineY(A1, A2);
plot(xY08t, lineaY08t)


[xY09 lineaY09] = makeAtomsLineY(C136, C115);
listaVert09 = crossLine2dAtom(MA, xY09, lineaY09, r, ep, 2);
arrVert09 = arrangeAtomsLineY(xY09, lineaY09, MA, listaVert09, nMA);

namefile = strcat(sVert, '9');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert09)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert09(i,1), arrVert09(i, 2), listaVert09(i));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(9, 1) = arrVert09(1,1);

for i=1:length(arrVert09) 

   [ii jj] = ismember(listaVert09(i), nMA);

   nMA(jj, :) = [listaVert09(i) arrVert09(i, 1) arrVert09(i, 2)];

end
%plot(nMA(:,2), nMA(:,3), 'ro')
A1 = [arrVert09(1, 1) arrVert09(1, 2)];
A2 = [arrVert09(length(arrVert09), 1) arrVert09(length(arrVert09), 2)];
[xY09t lineaY09t] = makeAtomsLineY(A1, A2);
plot(xY09t, lineaY09t)



[xY10 lineaY10] = makeAtomsLineY(C140, C112);
listaVert10 = crossLine2dAtom(MA, xY10, lineaY10, r, ep, 2);
arrVert10 = arrangeAtomsLineY(xY10, lineaY10, MA, listaVert10, nMA);

namefile = strcat(sVert, '10');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert10)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert10(i,1), arrVert10(i, 2), listaVert10(i));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(10, 1) = arrVert10(1,1);

for i=1:length(arrVert10) 

   [ii jj] = ismember(listaVert10(i), nMA);

   nMA(jj, :) = [listaVert10(i) arrVert10(i, 1) arrVert10(i, 2)];

end
%plot(nMA(:,2), nMA(:,3), 'ro')
A1 = [arrVert10(1, 1) arrVert10(1, 2)];
A2 = [arrVert10(length(arrVert10), 1) arrVert10(length(arrVert10), 2)];
[xY10t lineaY10t] = makeAtomsLineY(A1, A2);
plot(xY10t, lineaY10t)



[xY11 lineaY11] = makeAtomsLineY(C137, C94);
listaVert11 = crossLine2dAtom(MA, xY11, lineaY11, r, ep, 2);
arrVert11 = arrangeAtomsLineY(xY11, lineaY11, MA, listaVert11, nMA);

namefile = strcat(sVert, '11');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert11)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert11(i,1), arrVert11(i, 2), listaVert11(i));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(11, 1) = arrVert11(1,1);

for i=1:length(arrVert11) 

   [ii jj] = ismember(listaVert11(i), nMA);

   nMA(jj, :) = [listaVert11(i) arrVert11(i, 1) arrVert11(i, 2)];

end
%plot(nMA(:,2), nMA(:,3), 'ro')
A1 = [arrVert11(1, 1) arrVert11(1, 2)];
A2 = [arrVert11(length(arrVert11), 1) arrVert11(length(arrVert11), 2)];
[xY11t lineaY11t] = makeAtomsLineY(A1, A2);
plot(xY11t, lineaY11t)


[xY12 lineaY12] = makeAtomsLineY(C139, C92);
listaVert12 = crossLine2dAtom(MA, xY12, lineaY12, r, ep, 2);
arrVert12 = arrangeAtomsLineY(xY12, lineaY12, MA, listaVert12, nMA);

namefile = strcat(sVert, '12');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert12)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert12(i,1), arrVert12(i, 2), listaVert12(i));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(12, 1) = arrVert12(1,1);

for i=1:length(arrVert12) 

   [ii jj] = ismember(listaVert12(i), nMA);

   nMA(jj, :) = [listaVert12(i) arrVert12(i, 1) arrVert12(i, 2)];

end
%plot(nMA(:,2), nMA(:,3), 'ro')
A1 = [arrVert12(1, 1) arrVert12(1, 2)];
A2 = [arrVert12(length(arrVert12), 1) arrVert12(length(arrVert12), 2)];
[xY12t lineaY12t] = makeAtomsLineY(A1, A2);
plot(xY12t, lineaY12t)


[xY13 lineaY13] = makeAtomsLineY(C138, C70);
listaVert13 = crossLine2dAtom(MA, xY13, lineaY13, r, ep, 2);
arrVert13 = arrangeAtomsLineY(xY13, lineaY13, MA, listaVert13, nMA);

namefile = strcat(sVert, '13');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert13)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert13(i,1), arrVert13(i, 2), listaVert13(i));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(13, 1) = arrVert13(1,1);

for i=1:length(arrVert13) 

   [ii jj] = ismember(listaVert13(i), nMA);

   nMA(jj, :) = [listaVert13(i) arrVert13(i, 1) arrVert13(i, 2)];

end
%plot(nMA(:,2), nMA(:,3), 'ro')
A1 = [arrVert13(1, 1) arrVert13(1, 2)];
A2 = [arrVert13(length(arrVert13), 1) arrVert13(length(arrVert13), 2)];
[xY13t lineaY13t] = makeAtomsLineY(A1, A2);
plot(xY13t, lineaY13t)


[xY14 lineaY14] = makeAtomsLineY(C142, C135);
listaVert14 = crossLine2dAtom(MA, xY14, lineaY14, r, ep, 2);
arrVert14 = arrangeAtomsLineY(xY14, lineaY14, MA, listaVert14, nMA);

namefile = strcat(sVert, '14');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert14)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert14(i,1), arrVert14(i, 2), listaVert14(i));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(14, 1) = arrVert14(1,1);

for i=1:length(arrVert14) 

   [ii jj] = ismember(listaVert14(i), nMA);

   nMA(jj, :) = [listaVert14(i) arrVert14(i, 1) arrVert14(i, 2)];

end
%plot(nMA(:,2), nMA(:,3), 'ro')
A1 = [arrVert14(1, 1) arrVert14(1, 2)];
A2 = [arrVert14(length(arrVert14), 1) arrVert14(length(arrVert14), 2)];
[xY14t lineaY14t] = makeAtomsLineY(A1, A2);
plot(xY14t, lineaY14t)



[xY15 lineaY15] = makeAtomsLineY(C90, C133);
listaVert15 = crossLine2dAtom(MA, xY15, lineaY15, r, ep, 2);
arrVert15 = arrangeAtomsLineY(xY15, lineaY15, MA, listaVert15, nMA);

namefile = strcat(sVert, '15');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert15)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert15(i,1), arrVert15(i, 2), listaVert15(i));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(15, 1) = arrVert15(1,1);

for i=1:length(arrVert15) 

   [ii jj] = ismember(listaVert15(i), nMA);

   nMA(jj, :) = [listaVert15(i) arrVert15(i, 1) arrVert15(i, 2)];

end
%plot(nMA(:,2), nMA(:,3), 'ro')
A1 = [arrVert15(1, 1) arrVert15(1, 2)];
A2 = [arrVert15(length(arrVert15), 1) arrVert15(length(arrVert15), 2)];
[xY15t lineaY15t] = makeAtomsLineY(A1, A2);
plot(xY15t, lineaY15t)



[xY16 lineaY16] = makeAtomsLineY(C68, C131);
listaVert16 = crossLine2dAtom(MA, xY16, lineaY16, r, ep, 2);
arrVert16 = arrangeAtomsLineY(xY16, lineaY16, MA, listaVert16, nMA);

namefile = strcat(sVert, '16');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert16)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert16(i,1), arrVert16(i, 2), listaVert16(i));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(16, 1) = arrVert16(1,1);

for i=1:length(arrVert16) 

   [ii jj] = ismember(listaVert16(i), nMA);

   nMA(jj, :) = [listaVert16(i) arrVert16(i, 1) arrVert16(i, 2)];

end
%plot(nMA(:,2), nMA(:,3), 'ro')
A1 = [arrVert16(1, 1) arrVert16(1, 2)];
A2 = [arrVert16(length(arrVert16), 1) arrVert16(length(arrVert16), 2)];
[xY16t lineaY16t] = makeAtomsLineY(A1, A2);
plot(xY16t, lineaY16t)



[xY17 lineaY17] = makeAtomsLineY(C65, C30);
listaVert17 = crossLine2dAtom(MA, xY17, lineaY17, r, ep, 2);
arrVert17 = arrangeAtomsLineY(xY17, lineaY17, MA, listaVert17, nMA);

namefile = strcat(sVert, '17');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert17)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert17(i,1), arrVert17(i, 2), listaVert17(i));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(17, 1) = arrVert17(1,1);

for i=1:length(arrVert17) 

   [ii jj] = ismember(listaVert17(i), nMA);

   nMA(jj, :) = [listaVert17(i) arrVert17(i, 1) arrVert17(i, 2)];

end
%plot(nMA(:,2), nMA(:,3), 'ro')
A1 = [arrVert17(1, 1) arrVert17(1, 2)];
A2 = [arrVert17(length(arrVert17), 1) arrVert17(length(arrVert17), 2)];
[xY17t lineaY17t] = makeAtomsLineY(A1, A2);
plot(xY17t, lineaY17t)



[xY18 lineaY18] = makeAtomsLineY(C47, C31);
listaVert18 = crossLine2dAtom(MA, xY18, lineaY18, r, ep, 2);
arrVert18 = arrangeAtomsLineY(xY18, lineaY18, MA, listaVert18, nMA);

namefile = strcat(sVert, '18');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert18)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert18(i,1), arrVert18(i, 2), listaVert18(i));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(18, 1) = arrVert18(1,1);

for i=1:length(arrVert18) 

   [ii jj] = ismember(listaVert18(i), nMA);

   nMA(jj, :) = [listaVert18(i) arrVert18(i, 1) arrVert18(i, 2)];

end
%plot(nMA(:,2), nMA(:,3), 'ro')
A1 = [arrVert18(1, 1) arrVert18(1, 2)];
A2 = [arrVert18(length(arrVert18), 1) arrVert18(length(arrVert18), 2)];
[xY18t lineaY18t] = makeAtomsLineY(A1, A2);
plot(xY18t, lineaY18t)



[xY19 lineaY19] = makeAtomsLineY(C48, C35);
listaVert19 = crossLine2dAtom(MA, xY19, lineaY19, r, ep, 2);
arrVert19 = arrangeAtomsLineY(xY19, lineaY19, MA, listaVert19, nMA);

namefile = strcat(sVert, '19');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert19)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert19(i,1), arrVert19(i, 2), listaVert19(i));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(19, 1) = arrVert19(1,1);

for i=1:length(arrVert19) 

   [ii jj] = ismember(listaVert19(i), nMA);

   nMA(jj, :) = [listaVert19(i) arrVert19(i, 1) arrVert19(i, 2)];

end
%plot(nMA(:,2), nMA(:,3), 'ro')
A1 = [arrVert19(1, 1) arrVert19(1, 2)];
A2 = [arrVert19(length(arrVert19), 1) arrVert19(length(arrVert19), 2)];
[xY19t lineaY19t] = makeAtomsLineY(A1, A2);
plot(xY19t, lineaY19t)


[xY20 lineaY20] = makeAtomsLineY(C41, C43);
listaVert20 = crossLine2dAtom(MA, xY20, lineaY20, r, ep, 2);
arrVert20 = arrangeAtomsLineY(xY20, lineaY20, MA, listaVert20, nMA);

namefile = strcat(sVert, '20');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert20)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert20(i,1), arrVert20(i, 2), listaVert20(i));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(20, 1) = arrVert20(1,1);

for i=1:length(arrVert20) 

   [ii jj] = ismember(listaVert20(i), nMA);

   nMA(jj, :) = [listaVert20(i) arrVert20(i, 1) arrVert20(i, 2)];

end
%plot(nMA(:,2), nMA(:,3), 'ro')
A1 = [arrVert20(1, 1) arrVert20(1, 2)];
A2 = [arrVert20(length(arrVert20), 1) arrVert20(length(arrVert20), 2)];
[xY20t lineaY20t] = makeAtomsLineY(A1, A2);
plot(xY20t, lineaY20t)



[xY21 lineaY21] = makeAtomsLineY(C45, C46);
listaVert21 = crossLine2dAtom(MA, xY21, lineaY21, r, ep, 2);
arrVert21 = arrangeAtomsLineY(xY21, lineaY21, MA, listaVert21, nMA);

namefile = strcat(sVert, '21');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert21)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert21(i,1), arrVert21(i, 2), listaVert21(i));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(21, 1) = arrVert21(1,1);

for i=1:length(arrVert21) 

   [ii jj] = ismember(listaVert21(i), nMA);

   nMA(jj, :) = [listaVert21(i) arrVert21(i, 1) arrVert21(i, 2)];

end
plot(nMA(:,2), nMA(:,3), 'ro')
A1 = [arrVert21(1, 1) arrVert21(1, 2)];
A2 = [arrVert21(length(arrVert21), 1) arrVert21(length(arrVert21), 2)];
[xY21t lineaY21t] = makeAtomsLineY(A1, A2);
plot(xY21t, lineaY21t)


namefile = 'nMA-BN.dat';
fileID = fopen(namefile, 'w');
for i=1:length(nMA(:,1))
   fprintf(fileID, '%f\t %f\t %f\t\n', nMA(i,1),...
   nMA(i,2), nMA(i,3));
end




xnMA = nMA(:,2);
ynMA = nMA(:,3);


xnMA = unique(xnMA);
ynMA = unique(ynMA);


xnMA2 = verticalesL;
ynMA2 = horizontalesL;
ynMA2(1) = [];





[XnMA YnMA] = meshgrid(xnMA, ynMA);
Zmalla = (XnMA - XnMA) + (YnMA - YnMA) + 1.5;
%figure('Name', 'malla')
figure('Name', 'mesh')
%surface(Xp, Yp, Zp, 'EdgeColor', 'none'), view(3)
mesh(Xp, Yp, Zp)
axis([-inf inf -inf inf -5 5])
hold on
surface(Xl', Yl', ZrC, 'EdgeColor', 'none'), view(3)
mesh(XnMA, YnMA, Zmalla)


figure('Name', 'mallanMA')
hold on
for i=1:length(xnMA)

   plot([xnMA(i) xnMA(i)], [ynMA(1) ynMA(length(ynMA))]);

end
for i=1:length(ynMA)

   plot([xnMA(1) xnMA(length(xnMA))], [ynMA(i) ynMA(i)]);

end
for i=1:length(xnMA2)

   plot([xnMA2(i) xnMA2(i)], [ynMA2(1) ynMA2(length(ynMA2))], 'r');

end
for i=1:length(ynMA2)

   plot([xnMA2(1) xnMA2(length(xnMA2))], [ynMA2(i) ynMA2(i)], 'r');

end
axis([limitXIzq limitXDer limitYAbj limitYArr]);
grid on

save('nMA-BN', 'nMA')

%A = importdata('nMA-BN');
A = nMA;
B = importdata('char-BN.dat');
C = importdata('frame1-BN-rota.dat');

fileID = fopen('mallanMA-BN.xyz','w');
fileID2 = fopen('mallanMA-BN.dat','w');
nAtoms = 142;

[A1 A2] = unique(A(:,1));

fprintf(fileID, '%s\n', int2str(nAtoms));
fprintf(fileID, '%s\t%d\n', 'frame', 1);

for i=1:(nAtoms)

   if (B{i}=='H')
      
      fprintf(fileID, '%c  %f  %f  %f\n', B{i}, C(i,1), C(i,2), 0.00);
      fprintf(fileID2, '%f  %f  %f\n', C(i,1), C(i,2), 0.00);

   
   else
      j = i - 28; % siempre y cuando los primeros 28 sean
                  % hidrogeno se han tomado tambien las mismas
                  % posiciones en la malla corregida que en la
                  % malla sin corregir, se han dejado las
                  % posiciones en z como cero ya que la primera
                  % intencion ha sido ver que haya orden y todos
                  % los lugares de los atomos esten ocupados.

      fprintf(fileID, '%c  %f  %f  %f\n', B{i}, A(A2(j),2), A(A2(j),3), 0.00);
      fprintf(fileID2, '%f  %f  %f\n', A(A2(j),2), A(A2(j),3), 0.00);
   end

end

fclose('all');

% Falta comparar en XYZ la malla corregida con la 
% no corregida
%

% obteniendo el mapeo del numero de atomos en la 
% malla 'corregida' con los numeros de atomos de la malla
% sin corregir
%
% sirve como verificacion de atomos
%
% ya que es probable que el atomo 3 no coincida con el atomo
% 3 de la malla original

%
%
%



D = importdata('mallanMA-BN.dat');

fileMap = fopen('mapnMA-BN.dat','w');
for i=1:nAtoms

   F = abs(X-D(i,1));
   G = abs(Y-D(i,2));

   H = abs(F + G);

   [I J] = min(H);

   fprintf(fileMap, '%f  %f\n', i, J);

end

fclose('all');











   




