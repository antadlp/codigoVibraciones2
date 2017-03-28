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

t1=1:500;
N = length(t1);
y1 = sin(2*pi*15*t1);

yf = fft(y1);
ypow = abs(yf/N);
ypow1 = ypow(1:N/2+1);
ypow1(2:end-1) = 2*ypow1(2:end-1)
f = (0:(N/2))/N;

A = length(f)

B = length(ypow1)


figure(3), plot(f, ypow1)





