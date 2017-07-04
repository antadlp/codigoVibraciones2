t  = 0:pi/10:2*pi;
%[X, Y, Z] = cylinder(4*cos(t));
%
%subplot(2,2,1); mesh(X); title('X');
%subplot(2,2,2); mesh(Y); title('Y');
%subplot(2,2,3); mesh(Z); title('Z');
%subplot(2,2,4); mesh(X,Y,Z); title('X,Y,Z');
%
%%Now animate a single plot
%

for i=1:100
   [Xp(:,:,i), Yp(:,:,i), Zp(:,:,i)] = cylinder(4*cos(t*i));
end

%
%figure(1), mesh(Xp(:,:,6))
%figure(2), mesh(Yp(:,:,6))
%figure(3), mesh(Zp(:,:,6))
%figure(4), mesh(Xp(:,:,6),Yp(:,:,6),Zp(:,:,6))
%
%X = Xp(:,:,1);
%h = mesh(X, 'ZDataSource', 'X');
%hold on
%
%for m=2:40
%   m
%   X = Xp(:,:,m);
%   refreshdata(h,'caller')
%   drawnow;
%   pause(0.5)
%end

%All works till here!!!
%lets see now with the subplots

%subplot(2,1,1); mesh(Xp(:,:,1)); title('X');
%subplot(2,1,2); mesh(Yp(:,:,1)); title('Y');

X = Xp(:,:,1);
Y = Yp(:,:,1);
subplot(2,1,1); hx = mesh(X, 'ZDataSource', 'X');
subplot(2,1,2); hy = mesh(Y, 'ZDataSource', 'Y');
hold on

for m=1:100
   m
   X = Xp(:,:,m);
   Y = Yp(:,:,m);
   refreshdata(hx,'caller')
   refreshdata(hy,'caller')
   drawnow;
   pause(0.35)
end

%AWSOMEEEE IT WOORKS :'D

%now see if i can save it to a movie

