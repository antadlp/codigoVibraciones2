function animacionTerminosFourier



%save('terB12N12a.mat', 'mfbeB12N12a')
%mfbeB12N12a(i,:) = manualFBT001(ffteB12N12a, I, 1);

load('terB12N12a.mat');
en = importdata('datosMallaPython/eB12N12a.dat');
Yp = mfbeB12N12a;


A = length(en)
a = length(Yp(1,:))

ap = abs(A-a) + 1


a = 1000;
b = 1500;

a+ap

B = length(en((a+ap):(b+ap)))
C = length(Yp(1,a:b))

Ae = a+ap-1;
Be = b+ap-1;

t = length(mfbeB12N12a(:,1));
t
x = 1:length(Yp(1,a:b));
y = real(Yp(1,a:b));

xmin = min(x);
xmax = max(x);
ymin = min(real(Yp(end,a:b)))
ymax = max(real(Yp(end,a:b)))


plot(x, en(Ae:Be), 'r')



axis([xmin xmax ymin ymax])



hold on
h = plot(x,y, 'YDataSource', 'y');


for m=2:10
    m
    axis([xmin xmax ymin ymax])
    y = real(Yp(m,a:b));
    refreshdata(h,'caller')
    drawnow; 
    pause(.7)
end

%
%y = funcionPulsoCuadrado(t, t(1), to, A);
%
%%g = GraficaVelocidadExp(t, A, A);
%%
%%plot(t, g, 'r')
%
%%hold on
%
%%plot([0 0], [0 5*A], 'r')
%
%h = plot(t, y, 'YDataSource', 'y', 'LineWidth', 2);
%
%
%
%%for m = 2:length(t)
%for m = 2:3
%	
%	%axis([t(1) t(length(t)) -1 (A + 5)])
%	axis([t(1) t(length(t)) -1 ((1.5)*A)])
%	y = funcionPulsoCuadrado(t, t(m), to, A);
%	refreshdata(h, 'caller')
%	drawnow;
%	pause(.0005)
%end
%
%figure(2)
%y = funcionPulsoCuadrado(t, -6, to, A);
%filename='funcionPulsoCuadradoFijo1.txt';
%fid = fopen(filename, 'w');
%fprintf(fid, '%f  %f\n', [t; y]);
%plot(t, y, 'b')
%axis([t(1) t(length(t)) -1 ((1.5)*A)])
%grid on
%
%

