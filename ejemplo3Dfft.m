fs = 100;
fr1 = 2;
fr2 = 7;
fr3 = 9;
fr4 = 1;

t = 0:1/fs:2-1/fs;                      % 10 sec sample
x = t;
y = x;
[X1 Y1] = meshgrid(x, y);

Z1 = (1.3)*sin(2*pi*2*fr1*X1)+((1.7)*sin(2*pi*2*fr2*(X1-2))) ...
+ (1.3)*sin(2*pi*fr3*Y1)+((1.7)*sin(2*pi*fr4*(Y1-2)));

[M N] = size(Z1);

figure('Name', 'Z1')
surface(X1, Y1, Z1,'EdgeColor', 'none'), view(3)

Zc1 = Z1(floor(M/3), :);
Xc1 = X1(floor(M/3), :);

figure('Name', 'Zc1')
plot(Xc1, Zc1)

Zc2 = Z1(floor(2*M/3), :);
Xc2 = X1(floor(2*M/3), :);
figure('Name', 'Zc2')
plot(Xc2, Zc2)

Zc3 = Z1(:, floor(M/3));
Yc3 = Y1(:, floor(M/3));
figure('Name', 'Zc3')
plot(Yc3, Zc3)


Zc4 = Z1(:, floor(2*M/3));
Yc4 = Y1(:, floor(2*M/3));
figure('Name', 'Zc4')
plot(Yc4, Zc4)

Z1f = fft2(Z1);
Z1sh = fftshift(Z1f);
Z1P = Z1sh.*conj(Z1sh);
%Z1P = Z1f.*conj(Z1f);
x2 = (-M/2:M/2-1)*(fs/M);
y2 = x2;

%x2 = (1:M)*(fs/M);
%y2 = x2;
[X2 Y2] = meshgrid(x2, y2);


figure('Name', 'Z1P')
maxP = maxValueMatrix(Z1P);
k = 1.00;
Z1P = Z1P/(k*maxP);
surface(X2, Y2, Z1P,'EdgeColor', 'none'), view(3)
maxP = maxValueMatrix(Z1P);

N = M;

l = 1;
for i=1:M
   for j=1:N
      if (Z1P(i,j) >= 0.2*maxP)
%         y(l, :) = [x2(i) y2(j)];
%         yx(l) = x2(i);
%         yy(l) = y2(j);
         ii(l) = i;
         jj(l) = j;
%         yx(l) = 2*pi*(i-1)/M;
%         yy(l) = 2*pi*(j-1)/N;
         yx(l) = x2(i);
         yy(l) = y2(j);
         l = l + 1;
      end
   end
end

v = unique(yx);
w = unique(yy);
numFreq = 10;

if (length(v) < numFreq)

   vv = v;

else

   for i=1:numFreq 

      vv(i) = v(length(v) - (i-1));

   end

end



if (length(w) < numFreq)

   ww = w;

else

   for i=1:numFreq 

      ww(i) = w(length(w) - (i-1));

   end

end

ii
jj

vv
ww

dummie = Z1P;
for i=1:20
   
   [aa bb cc] = maxValueMatrix2(dummie);
   Pf(i,:) = [x2(bb) y2(cc)];
   Pf2(i) = aa;
   dummie(bb, cc) = -inf;

end

Pf


















