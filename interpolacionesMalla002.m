function interpolacionesMalla002

close all
clear all

fs = 17;

location = '/media/antadlp/sda13/mallaPaulina/';
folder = 'datos-malla-BN/';
location = strcat(location,folder);

nMA = importdata(strcat(location,'nMA-BN.dat'));


xnMA = nMA(:,2);
ynMA = nMA(:,3);

xnMA = unique(xnMA);
ynMA = unique(ynMA);


[XNMA YNMA] = meshgrid(xnMA, ynMA);
x = min(xnMA):1/fs:max(xnMA);
y = min(ynMA):1/fs:max(ynMA);
[Xpol Ypol] = meshgrid(x,y);


%s0 = '/home/antadlp/malla/gridXY-BN/'
%s0 = '../../mallaPaulina/datos-malla-BN/gridXY-BN/'

s0 = strcat(location,'gridXY-BN/');
for al=1:2500


%   s = strcat('/home/antadlp/malla/separados-BN-zz-dat/frame', int2str(al), '.dat');
   s = strcat(location,'separados-BN-zz-dat/frame', int2str(al), '.dat');
   filename = s;
   zframe = importdata(filename);

   for i=2:13
      io = (i-1);
   %for i=2:3
   
   %   Ahr =0;
   %   s = 0;
   %   s2 = 0;
   %   tAhr = 0;
   %   toSpl = 0;
   %   clear Ahr, s, s2, tAhr, toSpl;
   
      s = strcat(s0,'horizontales-BN-', int2str(i));
      s2 = strcat('splineHor', i);
      filename = s;
      Ahr = importdata(filename);
   
      %Ahr contiene la informacion de la posicion de cada atomo sobre
      %cada linea horizontal, por ejemplo la linea equis va tener valores 
      %diferentes para la primer columna (valores de x-coord) y valores iguales
      %para la segunda columna que serian las posiciones y-coord por ser una linea
      %horizontal la tercer columna es el numero de atomo que corresponde a esas
      %coordenadas

      toSpl(:,1) = Ahr(:,1); %toSpl las que se van a interpolar
      toSpl(:,2) = zframe(Ahr(:,3),3); %se va interpolar las posiciones z, del atomo i
      %Ahr(i,3) con valor zframe(atomoi, z-values)
      
      if (Ahr(1,1) > xnMA(1))
   
         [ii jj] = min(abs(xnMA - Ahr(1,1)));
   
         for j=1:(jj-1)
   
            tAhr(j,1) = xnMA(j);
            tAhr(j,2) = 0;
   
         end
   
         ll = 1;
         for j=jj:(jj + length(Ahr(:,1)) - 1)
   
            tAhr(j,1) = Ahr(ll,1);
            tAhr(j,2) = toSpl(ll,2);
            ll = ll + 1;
         
         end
   
      end
   
      if ( Ahr(length(Ahr(:,1)), 1) < xnMA(length(xnMA)) )
   
         [ii jj] = min(abs(xnMA - Ahr(length(Ahr(:,1)),1)));
   
         rr = length(xnMA) - jj;
   
         tt = length(tAhr(:,1));
         tt = tt + 1;
         ss = 1;
         for q=tt:(tt + rr - 1)
   
            tAhr(q,1) = xnMA(jj + ss);
            tAhr(q,2) = 0.0;
            ss = ss + 1;
         end
   
      end
   
      if (exist('tAhr'))
   
         toSpl = tAhr;
   
      end
   
     

      interHor(io, :) = spline(toSpl(:,1), toSpl(:,2), x);
      %interHor(io, :) = pchip(toSpl(:,1), toSpl(:,2), xnMA);

%      if (length(toSpl(:,1)) == 14)
%
%         toSplX14(:,io) = toSpl(:,1);
%         toSplXZ14(:,io) = toSpl(:,2);
%
%      elseif (length(toSpl(:,1)) == 13)
%
%         toSplX13(:,io) = toSpl(:,1);
%         toSplXZ13(:,io) = toSpl(:,2);
%
%      elseif (length(toSpl(:,1)) == 12)
%
%         toSplX12(:,io) = toSpl(:,1);
%         toSplXZ12(:,io) = toSpl(:,2);
%
%      elseif (length(toSpl(:,1)) == 11)
%
%         toSplX11(:,io) = toSpl(:,1);
%         toSplXZ11(:,io) = toSpl(:,2);
%
%      end
%
%



%      length(toSpl(:,1))
%      fprintf('\nio es: %d\n', io);
%
%      figure('Name', 'verSpliness')
%      plot(toSpl(:,1),toSpl(:,2), 'ro')
%      hold on
%      plot(toSpl(:,1),toSpl(:,2))
%      plot(x, interHor(io,:), 'k')
%      %plot(xnMA, interHor(io,:), 'k')
%      grid on
      
      %interHor(i, :) = spline(toSpl(:,1), toSpl(:,2), xnMA);
   
%      figure('Name', 'verSpliness')
%      plot(Ahr(:,1), zframe(Ahr(:,3),3), 'ro')
%      hold on
%      grid on
%      plot(xnMA, interHor(i,:))
%      




      clear toSpl;
   
      if (exist('tAhr'))
         clear tAhr;
      end
   
   
   
   end
   
   
      
%   clear toSpl;
%   for i=1:21
%   %for i=2:3
%   
%   %   Ahr =0;
%   %   s = 0;
%   %   s2 = 0;
%   %   tAhr = 0;
%   %   toSpl = 0;
%   %   clear Ahr, s, s2, tAhr, toSpl;
%   
%      s = strcat(s0,'verticales-BN-', int2str(i));
%      s2 = strcat('splineVert', i);
%      filename = s;
%      Avt = importdata(filename);
%   
%      Avt(:,2) = floor(Avt(:,2)*10000)/10000;
%      ynMA = floor(ynMA*10000)/10000;
%   
%    
%   
%   
%      toSpl(:,1) = Avt(:,2);
%      toSpl(:,2) = zframe(Avt(:,3),3);
%   
%      if (Avt(1,2) > ynMA(1))
%   
%         [ii jj] = min(abs(ynMA - Avt(1,2)));
%   
%         for j=1:(jj-1)
%   
%            tAvt(j,1) = ynMA(j);
%            tAvt(j,2) = 0;
%   
%         end
%   
%         ll = 1;
%         for j=jj:(jj + length(Avt(:,2)) - 1)
%   
%            tAvt(j,1) = Avt(ll,2);
%            tAvt(j,2) = toSpl(ll,2);
%            ll = ll + 1;
%         
%         end
%   
%      end
%   
%      [nn mm] = size(Avt);
%   
%      if (Avt(nn, 2) < ynMA(length(ynMA)))
%   
%         [ii jj] = min(abs(ynMA - Avt(length(Avt(:,2)),2)));
%   
%         rr = length(ynMA) - jj;
%   
%         tt = length(tAvt(:,1));
%         tt = tt + 1;
%         ss = 1;
%         for q=tt:(tt + rr - 1)
%   
%            tAvt(q,1) = ynMA(jj + ss);
%            tAvt(q,2) = 0.0;
%            ss = ss + 1;
%         end
%   
%      end
%   
%      if (exist('tAvt'))
%   
%         toSpl = tAvt;
%   
%      end
%   
%   %   figure('Name', strcat('interpolVert',int2str(i)))
%   %   plot(toSpl(:,1), toSpl(:,2), 'r')
%   %   hold on
%   %   plot(toSpl(:,1), toSpl(:,2), 'ro')
%   %   grid on
%   %   plot(ynMA, 0, '-go')
%   %   plot(Avt(:,2), Avt(:,3))
%   
%   %   interHor(i, :) = spline(toSpl(:,1), toSpl(:,2), xnMA);
%      interVert(i, :) = pchip(toSpl(:,1), toSpl(:,2), ynMA);
%   %    interVert(i, :) = spline(toSpl(:,1), toSpl(:,2), y);
%   
%   %   plot(ynMA, interVert(i,:), 'ko--')
%   
%      clear toSpl;
%   
%      if (exist('tAvt'))
%         clear tAvt;
%      end
%   
%   
%   
%   end
%   
     
   %filename = s;
   %fileID = fopen(filename, 'w');

%   length(xnMA)
%   size(interHor)
%   size(interVert)
   
%   for ii=1:length(xnMA)
%      for jj=1:length(ynMA)
%   
%         %interTotal(ii,jj,al) = interHor(jj,ii);
%         %fprintf('%f %f', interVert(ii,jj), interHor(jj,ii));
%         %fprintf('\n');
%
%         if (jj==1 | jj==length(ynMA))
%            interTotal(ii,jj,al) = interHor(jj,ii);
%         else
%            interTotal(ii,jj,al) = (interVert(ii,jj) + interHor(jj,ii))/2;
%         end
%
%   
%      end
%   end

   for ii=1:length(interHor(1,:))
      for jj=1:length(interHor(:,1))
   
            interTotal(ii,jj,al) = interHor(jj,ii);
   
      end
   end

   
   [XNMA YNMA] = meshgrid(x, ynMA);
   x = min(xnMA):1/fs:max(xnMA);
   y = min(ynMA):1/fs:max(ynMA);
   [Xpol Ypol] = meshgrid(x,y);




   inter2(:,:,al) = interp2(XNMA,YNMA, interTotal(:,:,al)',Xpol,Ypol);

   al
   %s2 = '/home/antadlp/malla/mallaInter-BN-exp/frame';
   s2 = strcat(location,'/mallaInter-BN-exp/frame');
   filename2 = strcat(s2, int2str(al));
   fileID2 = fopen(filename2, 'w');

%   for cc = 1:length(inter2(:,1,1))
%
%      fprintf(fileID2, '%f  ', inter2(cc,:,al));
%      fprintf(fileID2, '\n');
%
%   end

   %fclose(fileID2);
%   for iq=1:21
%      fprintf('%f %f\n', interHor(1,iq), interVert(iq,1));
%   end

%   figure('Name', 'corte')
%   plot(x, inter2(length(y),:,1))
%   hold on
%   plot(toSplX14(:,12), toSplXZ14(:,12), 'ro')
%
%   fprintf('\ntoSpl\n');
%   for iq=1:length(toSplX14(:,12))
%      fprintf('%f %f\n',toSplX14(iq,12), toSplXZ14(iq,12))
%   end




%   [x2, y2, Pf2(j,:), Pf(:,:,j)] = analisisfft001(fs, inter2(:,:,j), Xpol, Ypol, s);

end

save('zp-BN-exp-2500', 'inter2')

%t = 1:100-exp;
%for i=1:length(t)
%
%   s2 = '/home/toshiba/out003/mallaInter';
%   filename2 = strcat(s2, int2str(i));
%   Zp(:,:,i) = importdata(filename2);
%
%end
%
%Z = -1*Zp(:,:,1);
%h = surf(Xpol,Ypol,Z, 'ZDataSource', 'Z', 'EdgeColor', 'none');
%axis([-inf inf -inf inf -1 1])
%hold on
%
%for m=2:length(t)
%    Z = -1*Zp(:,:,m);
%    refreshdata(h,'caller')
%    drawnow; 
%    pause(.15)
%end
%
%

%dummie = Pf2;
%for i=1:10
%   
%   [aa bb cc] = maxValueMatrix2(dummie);
%   PfT(i,:) = Pf(cc,:,bb);
%   dummie(bb, cc) = -inf;
%
%end
%
%PfT
%
%

   
   



%
%[XNMA YNMA] = meshgrid(xnMA, ynMA);
%figure('Name', 'interPolTotal')
%surface(XNMA, YNMA, interTotal', 'EdgeColor', 'none'), view(3)
%axis([-inf inf -inf inf -2 2])
%
%figure('Name', 'interPolTotalM')
%mesh(XNMA, YNMA, interTotal')
%axis([-inf inf -inf inf -2 2])
%
%fs = 10;
%x = min(xnMA):1/fs:max(xnMA);
%y = min(ynMA):1/fs:max(ynMA);
%[Xpol Ypol] = meshgrid(x,y);
%inter2 = interp2(XNMA,YNMA,interTotal',Xpol,Ypol);
%
%figure('Name', 'inter2total')
%surface(Xpol,Ypol,inter2,'EdgeColor', 'none'), view(3)
%axis([-inf inf -inf inf -2 2])
%

%[v, w, vv, ww, Pf] = analisisfft001(fs, inter2, Xpol, Ypol, s);

%v
%w
%vv
%ww
%Pf


%^^^
%figure('Name', 'xy')
%plot(xX03, lineaX03)
%hold on
%plot(xX02, lineaX02)
%%plot([H11(1)], [H11(2)], 'o')
%plot(X', Y', 'o')
%plot(nMA(:,2), nMA(:,3), 'ro')
%axis([limitXIzq limitXDer limitYAbj limitYArr]);
%grid on
%


