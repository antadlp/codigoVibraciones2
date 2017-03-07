
%t = 0:1/fs:10-1/fs;                      % 10 sec sample

fs = 100;                                % Sample frequency (Hz)
fr1 = 2;
fr2 = 7;

P1 = 1/fr1;
P2 = 1/fr2;

tf = max([P1 P2]);

%t = 0:.01:10;                      % 10 sec sample
t = 0:1/fs:10-1/fs;                      % 10 sec sample

x = (1.3)*sin(2*pi*fr1*t)+((1.7)*sin(2*pi*fr2*(t-2)));          % 40 Hz component
%  + 2.5*gallery('normaldata',size(t),4); % Gaussian noise;

figure(1), plot(t, x)
axis([t(1) 12*P2 min(x) max(x)])

%
%P = gcd(P1,P1); 
%
%figure(1), plot(t, x)
%axis([t(1) P min(x) max(x)])

m = length(x);          % Window length
%n = pow2(nextpow2(m))  % Transform length
n = 512;
y = fft(x,n);           % DFT
f = (0:n-1)*(fs/n);     % Frequency range
power = y.*conj(y)/n;   % Power of the DFT

figure(2), plot(f,power)
xlabel('Frequency (Hz)')
ylabel('Power')
title('{\bf Periodogram}')

y0 = fftshift(y);          % Rearrange y values
f0 = (-n/2:n/2-1)*(fs/n);  % 0-centered frequency range
power0 = y0.*conj(y0)/n;   % 0-centered power

fft1 = fft(x);
power1 = fft1.*conj(fft1);
figure(5), plot(t, power1)
figure(6), plot(t, real(fft1))
y1 = fftshift(power1);          % Rearrange y values
t = 0:1/fs:10-1/fs;                      % 10 sec sample
tff = length(y1);
dt = abs(t(2) - t(1));

t = 0:(tff-1);
%t = (-tff/2:1/fs:tff/2);

figure(7), plot(t, y1)

figure(3), plot(f0,power0)
xlabel('Frequency (Hz)')
ylabel('Power')
title('{\bf 0-Centered Periodogram}')

t = 0:1/fs:10-1/fs;                      % 10 sec sample
fs = 100;
t0 = t(1);
t = t0:1/fs:(12*P2-(1/fs));
x1 = t;
y1 = x1;
[X1, Y1] = meshgrid(x1, y1);
Z1 = (1.3)*sin(2*pi*fr1*X1)+((1.7)*sin(2*pi*fr2*(X1-2))) ...
+ (1.3)*sin(2*pi*fr1*Y1)+((1.7)*sin(2*pi*fr2*(Y1-2)));

Z2 = fft2(Z1);
Z2s = fftshift(Z2);
Z20 = log2(1 + abs(Z2s));
%f = (0:n-1)*(fs/n);     % Frequency range
%f0 = (-n/2:n/2-1)*(fs/n);  % 0-centered frequency range

[M N] = size(Z2);
x2 = (-M/2:M/2-1)*(fs/M);
y2 = x2;
[X2 Y2] = meshgrid(x2, y2);


figure(4), surface(X1, Y1, Z1,'EdgeColor', 'none'), view(3)
figure('Name', 'power 2'), surface(X2, Y2, Z20,'EdgeColor', 'none'), view(3)

Z3 = abs(Z2s.^2);
Z3log = log(Z3 + eps);
figure('Name', 'power 2web'), surface(X2, Y2, Z3log,'EdgeColor', 'none'), view(3)

Z4 = Z2s.*conj(Z2s);
figure('Name', 'power 2conj'), surface(X2, Y2, Z4,'EdgeColor', 'none'), view(3)

Z4 = Z2s.*conj(Z2s)/M;
figure('Name', 'power 2conjM'), surface(X2, Y2, Z4,'EdgeColor', 'none'), view(3)

Z5 = Z2s.*conj(Z2s)/M;
Z5(Z5==0)=NaN;
figure('Name', 'power 2conjNaN'), surface(X2, Y2, Z5,'EdgeColor', 'none'), view(3)

Z6 = Z2s.*conj(Z2s);
maxP = maxValueMatrix(Z6);
k = 0.05;
Z6N = Z6/(k*maxP);

umbral = 1/100;
for i=1:M
   for j=1:N
      if (Z6N(i, j) <= umbral)
         Z6N(i, j) = 0;
      else
         Z6N(i, j) = Z6N(i, j);
      end
   end
end

Z6N(Z6N==0)=NaN;
figure('Name', 'power 2conjNaNN'), surface(X2, Y2, Z6N,'EdgeColor', 'none'), view(3)

Z7 = Z2s.*conj(Z2s);

umbral = .5*maxP;
for i=1:M
   for j=1:N
      if (Z7(i, j) <= umbral)
         Z7(i, j) = 0;
      else
         Z7(i, j) = Z7(i, j);
      end
   end
end

Z7(Z7==0)=NaN;
figure('ZName', 'Z7'), surface(X2, Y2, Z7,'EdgeColor', 'none'), view(3)

Z8 = Z6/(k*maxP);
%Z8(Z8<=10)=NaN;
figure('Name', 'Z8'), surface(X2, Y2, Z8,'EdgeColor', 'none'), view(3)

filename = 'frame500-001.dat';
A = importdata(filename);

Z9 = A(3, :);
[M N] = size(Z9);
x9 = (-M/2:M/2-1)*(fs/M);
y9 = x9;
[X2 Y2] = meshgrid(x2, y2);








