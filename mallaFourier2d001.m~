function mallaFourier2d001

close all
clear all

fs = 10;

nMA = importdata('nMA-B12N12b.dat');


xnMA = nMA(:,2);
ynMA = nMA(:,3);

xnMA = unique(xnMA);
ynMA = unique(ynMA);


[XNMA YNMA] = meshgrid(xnMA, ynMA);
x = min(xnMA):1/fs:max(xnMA);
y = min(ynMA):1/fs:max(ynMA);
[Xpol Ypol] = meshgrid(x,y);

load('zp-B12N12b.mat')
Zp = inter2;


s = 'a';
for j=1:length(Zp(1,1,:))
%for j=1:200
   [x2, y2, Pf2(j,:), Pf(:,:,j)] = analisisfft003(fs, Zp(:,:,j), Xpol, Ypol, s);
   j
end

sX = 'frecuenciasX-B12N12b.dat';
%filenameX = strcat(sX,int2str(i));
fileIDX = fopen(sX, 'w');

sY = 'frecuenciasY-B12N12b.dat';
%filenameY = strcat(sY,int2str(i));
fileIDY = fopen(sY, 'w');

for i=1:length(Zp(1,1,:))
%for i=1:200

   for j=1:length(Pf(:,1,1))
     
      fprintf(fileIDX, '%f\t', Pf(j,1,i));
      fprintf(fileIDY, '%f\t', Pf(j,2,i));

   end

   fprintf(fileIDX, '\n');
   fprintf(fileIDY, '\n');


end

save('Pf-B12N12b','Pf')





   
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
%fs = 10;
%x = min(xnMA):1/fs:max(xnMA);
%y = min(ynMA):1/fs:max(ynMA);
%[Xpol Ypol] = meshgrid(x,y);
%




