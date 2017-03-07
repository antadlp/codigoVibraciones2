filename = 'frame-500-002.dat';
filename2 = 'frame-500-003.dat';
A = importdata(filename);
B = importdata(filename2);

format short e

fs = 10;


Z1 = A(:,3);
X1 = A(:, 1);
Y1 = A(:, 2);
At = B(:, 1);

X1 = floor(X1*fs)/fs;
Y1 = floor(Y1*fs)/fs;

xmax = max(X1);
xmin= min(X1);
x = xmin:1/fs:xmax;

ymax = max(Y1);
ymin= min(Y1);
y = ymin:1/fs:ymax;



n = 1;
p = 1;



fileIDX = fopen('verX2.dat', 'w');
fileIDY = fopen('verY2.dat', 'w');

Z = zeros(length(x), length(y));
for i=1:length(x)
   for j=1:length(y)

      Z(i,j) = NaN;

   end
end

for i=1:length(X1)
   
   fprintf(fileIDX, '\nVERIFICANDO %f\n', X1(i));

   xCheck = abs(x - X1(i));
   [xx I(i)] = min(xCheck);
     
   fprintf(fileIDX, '%f\t%f\n\n', X1(i), x(I(i)));
   
   yCheck = abs(y - Y1(i));
   [yy J(i)] = min(yCheck);

   fprintf(fileIDY, '\nVERIFICANDO %f\n', Y1(i));
   fprintf(fileIDY, '%f\t%f\n\n', Y1(i), y(J(i)));

   Z(I(i), J(i)) = Z1(i);
   z2(i,3) = Z(I(i), J(i));
   z2(i,1) = x(I(i));
   z2(i,2) = y(J(i));

end
fclose(fileIDX);
fclose(fileIDY);



r = 0.3;
Zr = zeros(length(x), length(y));

for i=1:length(x)
   for j=1:length(y)

      ZrC(i,j) = NaN;

   end
end



for l=1:length(X1)
   for i=1:length(x)
      for j=1:length(y)

         ecC = (x(i) - x(I(l)))^2 + (y(j) - y(J(l)))^2;

         if (ecC <= r^2)
        
            Zr(i,j,l) = sqrt(r^2 - ecC);
            ZrC(i,j)= sqrt(r^2 - ecC) + Z(I(l),J(l)); 

         else

            Zr(i,j,l) = NaN;

         end

      end
   end
end

%xl = -2*pi:0.3:2*pi;
%yl = xl;
%
%[Xl Yl] = meshgrid(xl, yl);
%
%Zl = sin(Xl) + cos(Yl);
%
%figure
%%surf(Xl, Yl, Zl)
%surface(Xl, Yl, Zl, 'EdgeColor', 'none'), view(3)
%
xl = x;
yl = y;
Zl = ZrC;

rZ = zeros(length(xl), length(yl));
for i=1:length(xl)
   for j=length(yl)
      rZ(i, j) = NaN;
   end
end

Rz = rotz(20);

for i = 1:length(xl)
   for j = 1:length(yl)

      R = Rz*[xl(i); yl(j); Zl(i, j)];

      if (Zl(i,j) ~= NaN))

         mIJ(i,j) = 


      rX(i,j) = R(1);
      rY(i,j) = R(2);
      rZ(i,j) = R(3);

   end
end

rXmin = minValueMatrix(rX);
rXmax= maxValueMatrix(rX);
rYmin = minValueMatrix(rY);
rYmax= maxValueMatrix(rY);

rXX = rXmin:0.1:rXmax;
rYY = rYmin:0.1:rYmax;

rZZ = zeros(length(rXX), length(rYY));

for i = 1:length(xl)
   for j = 1:length(yl)

      [XXMin Ixx] = min(abs(rXX - rX(i,j)));
      [YYMin Iyy] = min(abs(rYY - rY(i,j)));

      rZZ(Ixx, Iyy) = rZ(i, j);

   end
end







      




   


[X Y] = meshgrid(x, y);
%figure('Name', 'surface')
%surface(X', Y', Z, 'EdgeColor', 'none'), view(3)
[RX RY] = meshgrid(rXX, rYY);
figure('Name', 'rota')
%surf(RX', RY', rZZ);
surface(RX', RY', rZZ, 'EdgeColor', 'none'), view(3)
axis([-inf inf -inf inf -5 5])
  
figure('Name', 'mesh')
mesh(X', Y', Z*6)
hold on
surface(X', Y', ZrC, 'EdgeColor', 'none'), view(3)
axis([-inf inf -inf inf -5 5])


%figure('Name', 'atomos')
%hold on
%for l=1:length(X1)
%  
%   surface(X', Y', Zr(:,:,l), 'EdgeColor', 'none'), view(3)
%
%end
%axis([-inf inf -inf inf -5 5])
 
figure('Name', 'atomos2')
surface(X', Y', ZrC, 'EdgeColor', 'none'), view(3)
axis([-inf inf -inf inf -5 5])

for i=1:length(X1)
%   disp(Z(I(i),J(i)));
   fprintf('%i\t%f\n', i, Z(I(i),J(i)));
end




[Xs Ixs] = sort(X1);
[Ys Iys] = sort(Y1);
Xa = unique(Xs);
Ya = unique(Ys);
Za = zeros(length(Xa), length(Ya));
for i=1:length(Xa)
   for j=1:length(Ya)

      Za(i,j) = NaN;

   end
end


for i=1:length(X1)
   
   xChecka = abs(Xa - X1(i));
   [xxa Ia(i)] = min(xChecka);
     
   yChecka = abs(Ya - Y1(i));
   [yya Ja(i)] = min(yChecka);

   Za(Ia(i), Ja(i)) = Z1(i);

end

ZrCa = zeros(length(Xa), length(Ya));
for i=1:length(Xa)
   for j=1:length(Ya)

      ZrCa(i,j) = NaN;

   end
end

r=0.7;
for l=1:length(X1)
   for i=1:length(Xa)
      for j=1:length(Ya)

         ecC = (Xa(i) - Xa(Ia(l)))^2 + (Ya(j) - Ya(Ja(l)))^2;

         if (ecC <= r^2)
        
            ZrCa(i,j)= sqrt(r^2 - ecC) + Za(Ia(l),Ja(l)); 

         end

      end
   end
end

[XA YA] = meshgrid(Xa, Ya);
figure('Name', 'atomos3')
surface(XA', YA', ZrCa, 'EdgeColor', 'none'), view(3)
axis([-inf inf -inf inf -5 5])

figure('Name', 'mesh3')
ZrCa = zeros(length(Xa), length(Ya));
mesh(XA', YA', ZrCa);
hold on
surface(X', Y', ZrC, 'EdgeColor', 'none'), view(3)
axis([-inf inf -inf inf -5 5])



