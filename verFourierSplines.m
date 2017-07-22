function verFourierSplines(A)

location = '/media/antadlp/sda13/mallaPaulina/';
folder = 'datos-malla-BN/';
location = strcat(location,folder);

nMA = importdata(strcat(location,'nMA-BN.dat'));
%nMA = importdata('../../mallaPaulina/datos-malla-B2N2/nMA-B2N2.dat');
%A = importdata('zp-B2N2-exp.mat');
xnMA = nMA(:,2);
ynMA = nMA(:,3);
xnMA = unique(xnMA);
ynMA = unique(ynMA);
fs = 17;
x = min(xnMA):1/fs:max(xnMA);
x(1)
x(end)

%length(x)

N = length(x);
xf = ((0:(N/2))-1)/(2*pi*N);
%yf = fft(A(116,:,1));

movieName = strcat('movie-fft','-BN-2500.avi')

v = VideoWriter(movieName);
v.Quality = 100;
v.FrameRate = 13;
open(v);

Z = A(1,:,1);
%h = surf(Xpol,Ypol,Z, 'ZDataSource', 'Z', 'EdgeColor', 'none');
h = plot(x, Z, 'YDataSource', 'Z');
axis([-inf inf -1 1]);
hold on

frames=1:2500;
a = frames(1)
b = frames(end)
for m=a:b
    m
    Z = A(1,:,m);
    refreshdata(h,'caller')
    drawnow; 
    frame = getframe(gcf);
    writeVideo(v,frame);
    pause(.055)
end

close(v)
for i=1:2
   
   i
   yf = fft(A(1,:,i));
   ypow = abs(yf/N);
   ypow2 = ypow((N/2):N);
   plot(xf,ypow2)
   hold on

%   for i=1:10
%      [mx(i) y(i)] = max(ypow2);
%      ypow2(y(i)) = -Inf;
%   end
%

   [mx(i) I] = max(ypow2);

   I;
   disp(x(I))
   freq(i) = x(I);
   
   ss=sum(freq)/2500;

   disp(ss) 

end

[alturaI tiempoI] = max(A(1,1,1:500))
[alturaD tiempoD] = max(A(1,end,1:500))





