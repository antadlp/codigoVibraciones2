close all
clear all
filename = 'zzout2.dat';
%filename = 'frame-500-000.dat';
XYZFile = importdata(filename);
format short e

MA = XYZFile;

%A = XYZFile(:,1);
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













H11 = [MA(11, 1) MA(11, 2) MA(11, 3)];
H07 = [MA(7, 1) MA(7, 2) MA(7, 3)];
H22 = [MA(22, 1) MA(22, 2) MA(22, 3)];
H04 = [MA(4, 1) MA(4, 2) MA(4, 3)];
H15 = [MA(15, 1) MA(15, 2) MA(15, 3)];
H14 = [MA(14, 1) MA(14, 2) MA(14, 3)];
H01 = [MA(1, 1) MA(1, 2) MA(1, 3)];
H05 = [MA(5, 1) MA(5, 2) MA(5, 3)];
H69 = [MA(69, 1) MA(69, 2) MA(69, 3)];
H47 = [MA(47, 1) MA(47, 2) MA(47, 3)];
H91 = [MA(91, 1) MA(91, 2) MA(91, 3)];
H48 = [MA(48, 1) MA(48, 2) MA(48, 3)];
C88 = [MA(88, 1) MA(88, 2) MA(88, 3)];
C42 = [MA(42, 1) MA(42, 2) MA(42, 3)];
C111 = [MA(111, 1) MA(111, 2) MA(111, 3)];
C43 = [MA(43, 1) MA(43, 2) MA(43, 3)];
C108 = [MA(108, 1) MA(108, 2) MA(108, 3)];
C44 = [MA(44, 1) MA(44, 2) MA(44, 3)];
C124 = [MA(124, 1) MA(124, 2) MA(124, 3)];
C45 = [MA(45, 1) MA(45, 2) MA(45, 3)];
C122 = [MA(122, 1) MA(122, 2) MA(122, 3)];
C46 = [MA(46, 1) MA(46, 2) MA(46, 3)];
C120 = [MA(120, 1) MA(120, 2) MA(120, 3)];
C40 = [MA(40, 1) MA(40, 2) MA(40, 3)];
C117 = [MA(117, 1) MA(117, 2) MA(117, 3)];
C41 = [MA(41, 1) MA(41, 2) MA(41, 3)];
C114 = [MA(114, 1) MA(114, 2) MA(114, 3)];
C34 = [MA(34, 1) MA(34, 2) MA(34, 3)];
C134 = [MA(134, 1) MA(134, 2) MA(134, 3)];
C35 = [MA(35, 1) MA(35, 2) MA(35, 3)];
C31 = [MA(31, 1) MA(31, 2) MA(31, 3)];
C132 = [MA(132, 1) MA(132, 2) MA(132, 3)];
C66 = [MA(66, 1) MA(66, 2) MA(66, 3)];
C129 = [MA(129, 1) MA(129, 2) MA(129, 3)];
C54 = [MA(54, 1) MA(54, 2) MA(54, 3)];
C128 = [MA(128, 1) MA(128, 2) MA(128, 3)];
C126 = [MA(126, 1) MA(126, 2) MA(126, 3)];
C52 = [MA(52, 1) MA(52, 2) MA(52, 3)];
C141 = [MA(141, 1) MA(141, 2) MA(141, 3)];
C125 = [MA(125, 1) MA(125, 2) MA(125, 3)];
C136 = [MA(136, 1) MA(136, 2) MA(136, 3)];
C115 = [MA(115, 1) MA(115, 2) MA(115, 3)];
C140 = [MA(140, 1) MA(140, 2) MA(140, 3)];
C112 = [MA(112, 1) MA(112, 2) MA(112, 3)];
C94 = [MA(94, 1) MA(94, 2) MA(94, 3)];
C137 = [MA(137, 1) MA(137, 2) MA(137, 3)];
C139 = [MA(139, 1) MA(139, 2) MA(139, 3)];
C92 = [MA(92, 1) MA(92, 2) MA(92, 3)];
C138 = [MA(138, 1) MA(138, 2) MA(138, 3)];
C70 = [MA(70, 1) MA(70, 2) MA(70, 3)];
C142 = [MA(142, 1) MA(142, 2) MA(142, 3)];
C135 = [MA(135, 1) MA(135, 2) MA(135, 3)];
C90 = [MA(90, 1) MA(90, 2) MA(90, 3)];
C133 = [MA(133, 1) MA(133, 2) MA(133, 3)];
C68 = [MA(68, 1) MA(68, 2) MA(68, 3)];
C131 = [MA(131, 1) MA(131, 2) MA(131, 3)];
C65 = [MA(65, 1) MA(65, 2) MA(65, 3)];
C30 = [MA(30, 1) MA(30, 2) MA(30, 3)];
C47 = [MA(47, 1) MA(47, 2) MA(47, 3)];
C31 = [MA(31, 1) MA(31, 2) MA(31, 3)];
C48 = [MA(48, 1) MA(48, 2) MA(48, 3)];
C35 = [MA(35, 1) MA(35, 2) MA(35, 3)];


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

figure('Name', 'mesh')
%surface(Xp, Yp, Zp, 'EdgeColor', 'none'), view(3)
mesh(Xp, Yp, Zp)
axis([-inf inf -inf inf -5 5])
hold on
surface(Xl', Yl', ZrC, 'EdgeColor', 'none'), view(3)




r = .3;
ep = 20;
offset = 0;
[xX01 lineaX01] = makeAtomsLine(H11, H07);
listaHor01= crossLine2dAtom(MA, xX01, lineaX01, r, ep, 1);
arrHor01 = arrangeAtomsLine(xX01, lineaX01, MA, listaHor01, 0.001);
namefile = strcat('horizontales', '01');
fileID = fopen(namefile, 'w');
for i=1:length(listaHor01)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrHor01(i,1), arrHor01(i, 2), MA(listaHor01(i), 3));
end
fprintf(fileID, '\n');
fclose(fileID);
horizontalesL(1, 1) = arrHor01(1,2);
for i=1:length(arrHor01)

   nMA(i, :) = [listaHor01(i) arrHor01(i, 1) arrHor01(i, 2)];

end
offset = length(nMA);
figure('Name', 'arrAtoms')
A1 = [arrHor01(1, 1) arrHor01(1, 2)];
A2 = [arrHor01(length(arrHor01), 1) arrHor01(length(arrHor01), 2)];
[xX01t lineaX01t] = makeAtomsLine(A1, A2);
plot(xX01t, lineaX01t)
hold on
axis([limitXIzq limitXDer limitYAbj limitYArr]);
grid on



[xX02 lineaX02] = makeAtomsLine(H69, H47);
listaHor02 = crossLine2dAtom(MA, xX02, lineaX02, r, ep, 1);
arrHor02 = arrangeAtomsLine(xX02, lineaX02, MA, listaHor02, 0.001);
namefile = strcat('horizontales', '02');
fileID = fopen(namefile, 'w');
for i=1:length(listaHor02)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrHor02(i,1), arrHor02(i, 2), MA(listaHor02(i), 3));
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



[xX03 lineaX03] = makeAtomsLine(H91, H48);
listaHor03 = crossLine2dAtom(MA, xX03, lineaX03, r, ep, 1);
arrHor03 = arrangeAtomsLine(xX03, lineaX03, MA, listaHor03, 0.001);

namefile = strcat('horizontales', '03');
fileID = fopen(namefile, 'w');
for i=1:length(listaHor03)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrHor03(i,1), arrHor03(i, 2), MA(listaHor03(i), 3));
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

namefile = strcat('horizontales', '04');
fileID = fopen(namefile, 'w');
for i=1:length(listaHor04)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrHor04(i,1), arrHor04(i, 2), MA(listaHor04(i), 3));
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

namefile = strcat('horizontales', '05');
fileID = fopen(namefile, 'w');
for i=1:length(listaHor05)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrHor05(i,1), arrHor05(i, 2), MA(listaHor05(i), 3));
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

namefile = strcat('horizontales', '06');
fileID = fopen(namefile, 'w');
for i=1:length(listaHor06)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrHor06(i,1), arrHor06(i, 2), MA(listaHor06(i), 3));
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

namefile = strcat('horizontales', '07');
fileID = fopen(namefile, 'w');
for i=1:length(listaHor07)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrHor07(i,1), arrHor07(i, 2), MA(listaHor07(i), 3));
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

namefile = strcat('horizontales', '08');
fileID = fopen(namefile, 'w');
for i=1:length(listaHor08)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrHor08(i,1), arrHor08(i, 2), MA(listaHor08(i), 3));
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

namefile = strcat('horizontales', '09');
fileID = fopen(namefile, 'w');
for i=1:length(listaHor09)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrHor09(i,1), arrHor09(i, 2), MA(listaHor09(i), 3));
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

namefile = strcat('horizontales', '10');
fileID = fopen(namefile, 'w');
for i=1:length(listaHor10)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrHor10(i,1), arrHor10(i, 2), MA(listaHor10(i), 3));
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

namefile = strcat('horizontales', '11');
fileID = fopen(namefile, 'w');
for i=1:length(listaHor11)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrHor11(i,1), arrHor11(i, 2), MA(listaHor11(i), 3));
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

namefile = strcat('horizontales', '12');
fileID = fopen(namefile, 'w');
for i=1:length(listaHor12)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrHor12(i,1), arrHor12(i, 2), MA(listaHor12(i), 3));
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

namefile = strcat('horizontales', '13');
fileID = fopen(namefile, 'w');
for i=1:length(listaHor13)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrHor13(i,1), arrHor13(i, 2), MA(listaHor13(i), 3));
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

namefile = strcat('verticales', '01');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert01)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert01(i,1), arrVert01(i, 2), MA(listaVert01(i), 3));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(1, 1) = arrVert01(1,2);

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

namefile = strcat('verticales', '02');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert02)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert02(i,1), arrVert02(i, 2), MA(listaVert02(i), 3));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(2, 1) = arrVert02(1,2);

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

namefile = strcat('verticales', '03');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert03)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert03(i,1), arrVert03(i, 2), MA(listaVert03(i), 3));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(3, 1) = arrVert03(1,2);

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

namefile = strcat('verticales', '04');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert04)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert04(i,1), arrVert04(i, 2), MA(listaVert04(i), 3));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(4, 1) = arrVert04(1,2);

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

namefile = strcat('verticales', '05');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert05)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert05(i,1), arrVert05(i, 2), MA(listaVert05(i), 3));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(5, 1) = arrVert05(1,2);

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

namefile = strcat('verticales', '06');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert06)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert06(i,1), arrVert06(i, 2), MA(listaVert06(i), 3));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(6, 1) = arrVert06(1,2);

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

namefile = strcat('verticales', '07');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert07)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert07(i,1), arrVert07(i, 2), MA(listaVert07(i), 3));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(7, 1) = arrVert07(1,2);

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

namefile = strcat('verticales', '08');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert08)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert08(i,1), arrVert08(i, 2), MA(listaVert08(i), 3));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(8, 1) = arrVert08(1,2);

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

namefile = strcat('verticales', '09');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert09)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert09(i,1), arrVert09(i, 2), MA(listaVert09(i), 3));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(9, 1) = arrVert09(1,2);

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

namefile = strcat('verticales', '10');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert10)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert10(i,1), arrVert10(i, 2), MA(listaVert10(i), 3));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(10, 1) = arrVert10(1,2);

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

namefile = strcat('verticales', '11');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert11)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert11(i,1), arrVert11(i, 2), MA(listaVert11(i), 3));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(11, 1) = arrVert11(1,2);

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

namefile = strcat('verticales', '12');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert12)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert12(i,1), arrVert12(i, 2), MA(listaVert12(i), 3));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(12, 1) = arrVert12(1,2);

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

namefile = strcat('verticales', '13');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert13)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert13(i,1), arrVert13(i, 2), MA(listaVert13(i), 3));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(13, 1) = arrVert13(1,2);

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

namefile = strcat('verticales', '14');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert14)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert14(i,1), arrVert14(i, 2), MA(listaVert14(i), 3));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(14, 1) = arrVert14(1,2);

for i=1:length(arrVert14) 

   [ii jj] = ismember(listaVert14(i), nMA);

   nMA(jj, :) = [listaVert14(i) arrVert14(i, 1) arrVert14(i, 2)];

end
%plot(nMA(:,2), nMA(:,3), 'ro')
A1 = [arrVert14(1, 1) arrVert14(1, 2)];
A2 = [arrVert14(length(arrVert14), 1) arrVert14(length(arrVert14), 2)];
[xY14t lineaY14t] = makeAtomsLineY(A1, A2);
plot(xY14t, lineaY14t)



[xY14 lineaY14] = makeAtomsLineY(C90, C133);
listaVert14 = crossLine2dAtom(MA, xY14, lineaY14, r, ep, 2);
arrVert14 = arrangeAtomsLineY(xY14, lineaY14, MA, listaVert14, nMA);

namefile = strcat('verticales', '14');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert14)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert14(i,1), arrVert14(i, 2), MA(listaVert14(i), 3));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(14, 1) = arrVert14(1,2);

for i=1:length(arrVert14) 

   [ii jj] = ismember(listaVert14(i), nMA);

   nMA(jj, :) = [listaVert14(i) arrVert14(i, 1) arrVert14(i, 2)];

end
%plot(nMA(:,2), nMA(:,3), 'ro')
A1 = [arrVert14(1, 1) arrVert14(1, 2)];
A2 = [arrVert14(length(arrVert14), 1) arrVert14(length(arrVert14), 2)];
[xY14t lineaY14t] = makeAtomsLineY(A1, A2);
plot(xY14t, lineaY14t)



[xY15 lineaY15] = makeAtomsLineY(C68, C131);
listaVert15 = crossLine2dAtom(MA, xY15, lineaY15, r, ep, 2);
arrVert15 = arrangeAtomsLineY(xY15, lineaY15, MA, listaVert15, nMA);

namefile = strcat('verticales', '15');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert15)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert15(i,1), arrVert15(i, 2), MA(listaVert15(i), 3));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(15, 1) = arrVert15(1,2);

for i=1:length(arrVert15) 

   [ii jj] = ismember(listaVert15(i), nMA);

   nMA(jj, :) = [listaVert15(i) arrVert15(i, 1) arrVert15(i, 2)];

end
%plot(nMA(:,2), nMA(:,3), 'ro')
A1 = [arrVert15(1, 1) arrVert15(1, 2)];
A2 = [arrVert15(length(arrVert15), 1) arrVert15(length(arrVert15), 2)];
[xY15t lineaY15t] = makeAtomsLineY(A1, A2);
plot(xY15t, lineaY15t)



[xY16 lineaY16] = makeAtomsLineY(C65, C30);
listaVert16 = crossLine2dAtom(MA, xY16, lineaY16, r, ep, 2);
arrVert16 = arrangeAtomsLineY(xY16, lineaY16, MA, listaVert16, nMA);

namefile = strcat('verticales', '16');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert16)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert16(i,1), arrVert16(i, 2), MA(listaVert16(i), 3));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(16, 1) = arrVert16(1,2);

for i=1:length(arrVert16) 

   [ii jj] = ismember(listaVert16(i), nMA);

   nMA(jj, :) = [listaVert16(i) arrVert16(i, 1) arrVert16(i, 2)];

end
%plot(nMA(:,2), nMA(:,3), 'ro')
A1 = [arrVert16(1, 1) arrVert16(1, 2)];
A2 = [arrVert16(length(arrVert16), 1) arrVert16(length(arrVert16), 2)];
[xY16t lineaY16t] = makeAtomsLineY(A1, A2);
plot(xY16t, lineaY16t)



[xY17 lineaY17] = makeAtomsLineY(C47, C31);
listaVert17 = crossLine2dAtom(MA, xY17, lineaY17, r, ep, 2);
arrVert17 = arrangeAtomsLineY(xY17, lineaY17, MA, listaVert17, nMA);

namefile = strcat('verticales', '17');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert17)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert17(i,1), arrVert17(i, 2), MA(listaVert17(i), 3));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(17, 1) = arrVert17(1,2);

for i=1:length(arrVert17) 

   [ii jj] = ismember(listaVert17(i), nMA);

   nMA(jj, :) = [listaVert17(i) arrVert17(i, 1) arrVert17(i, 2)];

end
%plot(nMA(:,2), nMA(:,3), 'ro')
A1 = [arrVert17(1, 1) arrVert17(1, 2)];
A2 = [arrVert17(length(arrVert17), 1) arrVert17(length(arrVert17), 2)];
[xY17t lineaY17t] = makeAtomsLineY(A1, A2);
plot(xY17t, lineaY17t)



[xY18 lineaY18] = makeAtomsLineY(C48, C35);
listaVert18 = crossLine2dAtom(MA, xY18, lineaY18, r, ep, 2);
arrVert18 = arrangeAtomsLineY(xY18, lineaY18, MA, listaVert18, nMA);

namefile = strcat('verticales', '18');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert18)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert18(i,1), arrVert18(i, 2), MA(listaVert18(i), 3));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(18, 1) = arrVert18(1,2);

for i=1:length(arrVert18) 

   [ii jj] = ismember(listaVert18(i), nMA);

   nMA(jj, :) = [listaVert18(i) arrVert18(i, 1) arrVert18(i, 2)];

end
%plot(nMA(:,2), nMA(:,3), 'ro')
A1 = [arrVert18(1, 1) arrVert18(1, 2)];
A2 = [arrVert18(length(arrVert18), 1) arrVert18(length(arrVert18), 2)];
[xY18t lineaY18t] = makeAtomsLineY(A1, A2);
plot(xY18t, lineaY18t)


[xY19 lineaY19] = makeAtomsLineY(C41, C43);
listaVert19 = crossLine2dAtom(MA, xY19, lineaY19, r, ep, 2);
arrVert19 = arrangeAtomsLineY(xY19, lineaY19, MA, listaVert19, nMA);

namefile = strcat('verticales', '19');
fileID = fopen(namefile, 'w');
for i=1:length(listaVert19)
   fprintf(fileID, '%f\t %f\t %f\t\n', ...
   arrVert19(i,1), arrVert19(i, 2), MA(listaVert19(i), 3));
end
fprintf(fileID, '\n');
fclose(fileID);
verticalesL(19, 1) = arrVert19(1,2);

for i=1:length(arrVert19) 

   [ii jj] = ismember(listaVert19(i), nMA);

   nMA(jj, :) = [listaVert19(i) arrVert19(i, 1) arrVert19(i, 2)];

end
plot(nMA(:,2), nMA(:,3), 'ro')
A1 = [arrVert19(1, 1) arrVert19(1, 2)];
A2 = [arrVert19(length(arrVert19), 1) arrVert19(length(arrVert19), 2)];
[xY19t lineaY19t] = makeAtomsLineY(A1, A2);
plot(xY19t, lineaY19t)



figure('Name', 'xy')
plot(xX01, lineaX01)
hold on
plot(xX02, lineaX02)
%plot([H11(1)], [H11(2)], 'o')
plot(X', Y', 'o')
plot(nMA(:,2), nMA(:,3), 'ro')
axis([limitXIzq limitXDer limitYAbj limitYArr]);
grid on


