
function GraficaPulsoCuadrado2(t, to, A)

y = funcionPulsoCuadrado(t, t(1), to, A);

%g = GraficaVelocidadExp(t, A, A);
%
%plot(t, g, 'r')

%hold on

%plot([0 0], [0 5*A], 'r')

h = plot(t, y, 'YDataSource', 'y', 'LineWidth', 2);



%for m = 2:length(t)
for m = 2:3
	
	%axis([t(1) t(length(t)) -1 (A + 5)])
	axis([t(1) t(length(t)) -1 ((1.5)*A)])
	y = funcionPulsoCuadrado(t, t(m), to, A);
	refreshdata(h, 'caller')
	drawnow;
	pause(.0005)
end

figure(2)
y = funcionPulsoCuadrado(t, -6, to, A);
filename='funcionPulsoCuadradoFijo1.txt';
fid = fopen(filename, 'w');
fprintf(fid, '%f  %f\n', [t; y]);
plot(t, y, 'b')
axis([t(1) t(length(t)) -1 ((1.5)*A)])
grid on



