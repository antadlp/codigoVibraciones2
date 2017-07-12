function twoGrapheneMovies(I, C, frames)

A = I{1};
B = I{2};
C = I{3};
D = I{4};
E = I{5};
F = I{6};

clear I;

An = C{1};
Bn = C{2};
Cn = C{3};
Dn = C{4};
En = C{5};
Fn = C{6};


v = VideoWriter('six-subplot-Graphene.avi');
v.Quality = 90;
v.FrameRate = 13;
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

nMA = importdata(Cn);
xnMA = nMA(:,2);
ynMA = nMA(:,3);
xnMA = unique(xnMA);
ynMA = unique(ynMA);
x = min(xnMA):1/fs:max(xnMA);
y = min(ynMA):1/fs:max(ynMA);
[XC, YC] = meshgrid(x,y);

nMA = importdata(Dn);
xnMA = nMA(:,2);
ynMA = nMA(:,3);
xnMA = unique(xnMA);
ynMA = unique(ynMA);
x = min(xnMA):1/fs:max(xnMA);
y = min(ynMA):1/fs:max(ynMA);
[XD, YD] = meshgrid(x,y);

nMA = importdata(En);
xnMA = nMA(:,2);
ynMA = nMA(:,3);
xnMA = unique(xnMA);
ynMA = unique(ynMA);
x = min(xnMA):1/fs:max(xnMA);
y = min(ynMA):1/fs:max(ynMA);
[XE, YE] = meshgrid(x,y);

nMA = importdata(Fn);
xnMA = nMA(:,2);
ynMA = nMA(:,3);
xnMA = unique(xnMA);
ynMA = unique(ynMA);
x = min(xnMA):1/fs:max(xnMA);
y = min(ynMA):1/fs:max(ynMA);
[XF, YF] = meshgrid(x,y);


Zpa = A;
clear A;
Za = 1*Zpa(:,:,1);

Zpb = B;
clear B;
Zb = 1*Zpb(:,:,1);

Zpc = C;
clear C;
Zc = 1*Zpc(:,:,1);

Zpd = D;
clear D;
Zd = 1*Zpd(:,:,1);

Zpe = E;
clear E;
Ze = 1*Zpe(:,:,1);

Zpf = F;
clear F;
Zf = 1*Zpf(:,:,1);

subplot(2,3,1); ha = surf(XA, YA, Za, 'ZDataSource', 'Za', ...
'EdgeColor', 'none');
axis([-inf inf -inf inf -0.6 0.6])
hold on

subplot(2,3,2); hb = surf(XB, YB, Zb, 'ZDataSource', 'Zb', ...
'EdgeColor', 'none');
axis([-inf inf -inf inf -0.6 0.6])

subplot(2,3,3); hc = surf(XC, YC, Zc, 'ZDataSource', 'Zc', ...
'EdgeColor', 'none');
axis([-inf inf -inf inf -0.6 0.6])
hold on

subplot(2,3,6); hd = surf(XD, YD, Zd, 'ZDataSource', 'Zd', ...
'EdgeColor', 'none');
axis([-inf inf -inf inf -0.6 0.6])

subplot(2,3,5); he = surf(XE, YE, Ze, 'ZDataSource', 'Ze', ...
'EdgeColor', 'none');
axis([-inf inf -inf inf -0.6 0.6])
hold on

subplot(2,3,4); hf = surf(XF, YF, Zf, 'ZDataSource', 'Zf', ...
'EdgeColor', 'none');
axis([-inf inf -inf inf -0.6 0.6])

fig = gcf;
set(fig, 'Position',[0, 0, 900, 900]);
hold on

a = frames(1);
b = frames(end);

for m=a:b
   m
   Za = 1*Zpa(:,:,m);
   Zb = 1*Zpb(:,:,m);
   Zc = 1*Zpc(:,:,m);
   Zd = 1*Zpd(:,:,m);
   Ze = 1*Zpe(:,:,m);
   Zf = 1*Zpf(:,:,m);
   refreshdata(ha, 'caller')
   refreshdata(hb, 'caller')
   drawnow;
   frame = getframe(gcf);
   writeVideo(v,frame)
   pause(0.0005)
end

close(v)

   



