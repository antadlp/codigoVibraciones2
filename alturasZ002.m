function alturasZ002


%Frecuentemente al tener graficas con nombres propios (figure('Name','nombre'))
%las ventanas no se sobreponen sino que se abre una nueva, pudiendo tener N*M
%ventanas abiertas de graficas. Siendo N el numero de graficas por script y M
%el numero de veces que se ha abierto el script.

close all
clear all

filename = '/home/toshiba/movie-B12N12b-zz.xyz';
fileID3 = fopen(filename, 'w');

for j = 1:5000
   
   s = '/home/toshiba/separados-B12N12b-dat/'
   filename = strcat(s,'frame',int2str(j),'.dat');
   XYZFile = importdata(filename);
   format short e

   X = XYZFile(:,1);
   Y = XYZFile(:,2);
   Z = XYZFile(:,3);

   fs = 10;
   X = floor(X*fs)/fs;
   Y = floor(Y*fs)/fs;

   MA = XYZFile;

   P = [MA(11,1) MA(11,2) MA(11,3)];
   Q = [MA(7,1) MA(7,2) MA(7,3)];
   R = [MA(1,1) MA(1,2) MA(1,3)];
   PQ = P - Q;
   PR = P - R;
   normal = cross(PQ, PR);

   xmax = max(X);
   xmin= min(X);
   x = xmin:1/fs:xmax;

   ymax = max(Y);
   ymin= min(Y);
   y = ymin:1/fs:ymax;

   extraLimitA = abs(MA(11,1) - MA(11,1));
   extraLimitB = abs(MA(15,1) - MA(15,1));

   extraLimit = max([extraLimitA extraLimitB]);
   extraLimit2 = extraLimit/2;

   limitXIzq = xmin - extraLimit2;
   limitXDer = xmax + extraLimit2;

   limitYArr = ymax + extraLimit2;
   limitYAbj = ymin - extraLimit2;

   xp = limitXIzq:01:limitXDer;
   yp = limitYAbj:01:limitYArr;

   [Xp Yp] = meshgrid(xp, yp);
   [Zp N] = makePlano001(Xp, Yp, P, normal');

   nAtoms = 142;
   for i=1:nAtoms
      
      Zp = [MA(i,1) MA(i,2) MA(i,3)];
      Zd(i) = distanciaPuntoPlano(Zp, N);

   end

   atoms = importdata('char-B12N12b.dat');
   aGrid = importdata('mallanMA-B12N12b.dat');
   mp = importdata('mapnMA-B12N12b.dat');

   s = '/home/toshiba/separados-B12N12b-zz-xyz/';
   s2 = '/home/toshiba/separados-B12N12b-zz-dat/';
   filename1 = strcat(s,'frame',int2str(j),'.xyz');
   filename2 = strcat(s2,'frame',int2str(j),'.dat');
   fileID = fopen(filename1, 'w');
   fileID2 = fopen(filename2, 'w');

   fprintf(fileID, '%s\n', int2str(nAtoms));
   fprintf(fileID, '%s\t%d\n', 'frame', j);
   fprintf(fileID3, '%s\n', int2str(nAtoms));
   fprintf(fileID3, '%s\t%d\n', 'frame', j);
   for i=1:(nAtoms)

      l = mp(i,2);
      if (atoms{i}=='H')

         fprintf(fileID, '%c  %f  %f  %f\n', atoms{i}, aGrid(l,1), aGrid(l,2), 0.00);
         fprintf(fileID3, '%c  %f  %f  %f\n', atoms{i}, aGrid(l,1), aGrid(l,2), 0.00);
         fprintf(fileID2, '%f  %f  %f\n', aGrid(l,1), aGrid(l,2), 0.00);
      
      else
         
         fprintf(fileID, '%c  %f  %f  %f\n', atoms{i}, aGrid(l,1), aGrid(l,2),...
         Zd(i));
 
         fprintf(fileID3, '%c  %f  %f  %f\n', atoms{i}, aGrid(l,1), aGrid(l,2),...
         Zd(i));

         fprintf(fileID2, '%f  %f  %f\n', aGrid(l,1), aGrid(l,2), Zd(i));
      end

   end
   fclose(fileID)
   fclose(fileID2)

   j


end

   fclose('all')








   






