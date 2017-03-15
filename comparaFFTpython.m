% Compare with
%http://stackoverflow.com/questions/25735153/plotting-a-fast-fourier-transform-in-python
%Paul H
%

N = 600;

T = 1/800;
x = linspace(0.0, N*T, N);
y = sin(50*2*pi*x) + 0.5*sin(80*2*pi*x);
yf = fft(y);
xf = linspace(0, 1/(2*T), N/2);

plot(xf, 2/N*abs(yf(1:(floor(N/2)))))


for i=1:10
   disp(yf(i))
end

 


