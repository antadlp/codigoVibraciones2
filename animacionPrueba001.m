%function AnimacionD1(t,x,p)
function animacionPrueba001


%[m n] = size(p);
%
%y = p(1,:);
%
%h = plot(x,y,'YDataSource','y','LineWidth',2);
%axis([x(1) x(length(x)) -.5 .5])

%for m=2:length(t)
%    y = p(m,:);
%    refreshdata(h,'caller')
%    drawnow; 
%    pause(.05)
%end
%
k=2;
l=3;
w = 0.5;

x = -k*pi:0.06:k*pi;
y = -l*pi:0.06:l*pi;
t = 0:100;

[X Y] = meshgrid(x,y);

for i=1:length(t)
   Zp(:,:,i) = sin(X-w*(t(i)));
end

%surface(X2, Y2, Z1Pk,'EdgeColor', 'none'), view(3)
%h = plot(x,y,'YDataSource','y','LineWidth',2);
Z = Zp(:,:,1);
h = surf(X,Y,Z, 'ZDataSource', 'Z', 'EdgeColor', 'none');

for m=2:length(t)
    Z = Zp(:,:,m);
    refreshdata(h,'caller')
    drawnow; 
    pause(.15)
end


%
%figure 
%Z = peaks;
%surf(Z)
%axis tight manual 
%ax = gca;
%ax.NextPlot = 'replaceChildren';
%loops = 40;
%F(loops) = struct('cdata',[],'colormap',[]);
%for j = 1:loops
%    X = sin(j*pi/10)*Z;
%    surf(X,Z)
%    drawnow
%    F(j) = getframe;
%end
%
%
