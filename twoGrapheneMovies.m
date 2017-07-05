function twoGrapheneMovies(I, C, frames)

A = I{1};
B = I{2};
clear I;

An = C{1};
Bn = C{2};

v = VideoWriter('two-subplot-Graphene.avi');
v.Quality = 70;
open(v);

fs = 10;

nMA = importdata(An);
xnMA = nMA(:,2);
ynMA = nMA(:,3);
xnMA = unique(xnMA);
ynMA = unique(ynMA);
x = min(xnMA):1/fs:max(xnMA);
y = min(ynMA):1/fs:max(ynMA);
[XA, YA] = meshgrid(x,y);

nMA = importdata(Bn);
xnMA = nMA(:,2);
ynMA = nMA(:,3);
xnMA = unique(xnMA);
ynMA = unique(ynMA);
x = min(xnMA):1/fs:max(xnMA);
y = min(ynMA):1/fs:max(ynMA);
[XB, YB] = meshgrid(x,y);

Zpa = A;
clear A;
Za = 1*Zpa(:,:,1);

Zpb = B;
clear B;
Zb = 1*Zpb(:,:,1);


subplot(2,1,1); ha = surf(XA, YA, Za, 'ZDataSource', 'Za', ...
'EdgeColor', 'none');

subplot(2,1,2); hb = surf(XB, YB, Zb, 'ZDataSource', 'Zb', ...
'EdgeColor', 'none');
axis([-inf inf -inf inf -1 1])
hold on

a = frames(1);
b = frames(end);

for m=a:b
   m
   Za = 1*Zpa(:,:,m);
   Zb = 1*Zpb(:,:,m);
   refreshdata(ha, 'caller')
   refreshdata(hb, 'caller')
   drawnow;
   frame = getframe(gcf);
   writeVideo(v,frame)
   pause(0.35)
end

close(v)

   



