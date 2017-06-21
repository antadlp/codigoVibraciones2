Fs = 1500;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 1500;             % Length of signal
t = (0:L-1)*T;        % Time vector

S = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t);

X = S + 2*randn(size(t));

figure(1), plot(1000*t(1:50),X(1:50))
title('Signal Corrupted with Zero-Mean Random Noise')
xlabel('t (milliseconds)')
ylabel('X(t)')


Y = fft(X);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
figure(2), plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

clear all
deltat = 0.01;
t1 = deltat:deltat:100;
y1 = sin(2*pi*(5.3)*t1) + 0.3*sin(2*pi*9*t1) + sin(2*pi*20*t1);

N = length(t1)

yf = fft(y1);
P2 = abs(yf/N);


P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);

deltat = abs(t1(10)-t1(11));
f = (1/deltat)*(0:(N/2))/N;

lf = length(f)

figure(3), plot(f, P1)
grid on
