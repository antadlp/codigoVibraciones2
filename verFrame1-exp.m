close all
clear all

fs = 10;

nMA = importdata('nMA-B12N12a.dat');


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


%[x2, y2, Pf2(j,:), Pf(:,:,j)] = analisisfft003(fs, Zp(:,:,j), Xpol, Ypol, s);
[x2 y2 Pf2 Pf] = analisisfft001(fs, Zp(:,:,1), Xpol, Ypol, 'S');



