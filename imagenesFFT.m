% 2D FFT Demo
% 
% http://matlabgeeks.com/tips-tutorials/how-to-do-a-2-d-fourier-transform-in-matlab/
%
%

% treeA320x240.jpg

[imageA mapA] = imread('treeA320x240', 'jpg');
[imageB mapB] = imread('treeB320x240', 'jpg');

figure, imshow(imageA)
title('Image A - Greek Church')

figure, imshow(imageB)
title('Image A - Aishawarya Rai')

fftA = fft2(double(imageA));
fftB = fft2(double(imageB));


A = abs(fftshift(fftA));
A = 1000.*A;

figure, imshow(fftshift(fftA))

%, colormap gray
%title('Image A FFT2 Magnitude')
% figure, imshow(angle(fftshift(fftA)),[-pi pi]), colormap gray
% title('Image A FFT2 Phase')
% figure, imshow(abs(fftshift(fftB)),[24 100000]), colormap gray
% title('Image B FFT2 Magnitude')
% figure, imshow(angle(fftshift(fftB)),[-pi pi]), colormap gray
% title('Image B FFT2 Phase')





fftC = abs(fftA).*exp(i*angle(fftB));
fftD = abs(fftB).*exp(i*angle(fftA));

imageC = ifft2(fftC);
imageD = ifft2(fftD);

cmin = min(min(abs(imageC)));
cmax = max(max(abs(imageC)));

dmin = min(min(abs(imageD)));
dmax = max(max(abs(imageD)));


% figure, imshow(abs(imageC), [cmin cmax]), colormap gray
% title('Image C  Magnitude')
% figure, imshow(abs(imageD), [dmin dmax]), colormap gray
% title('Image D  Magnitude')

%saveas(1,'imageA2.png')
%saveas(2,'imageB2.png')
%saveas(3,'imageAfftmag.png')
%saveas(4,'imageAfftpha.png')
%saveas(5,'imageBfftmag.png')
%saveas(6,'imageBfftpha.png')
%saveas(7,'imageC.png')
%saveas(8,'imageD.png')



