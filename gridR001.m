function F = gridR001

filename = 'outRotadosGridR';
datosAll = importdata(filename);


for i=1:1000

   frames(:,i) = datosAll((((i-1)*8*23) + 1):(i*8*23),1);

end

namefile = 'FFT2-REAL.dat';
namefile2 = 'FFT2-IMAGINARIO.dat';
fileID = fopen(namefile, 'w');
fileID2 = fopen(namefile2, 'w');
for i=1:1000

   s = strcat('  FRAME_',int2str(i));

   fprintf(fileID, s);
   fprintf(fileID, '\n');

   fprintf(fileID2, s);
   fprintf(fileID2, '\n');


   matrizFrames(:,:,i) = reshape(frames(:,i), [23, 8]);
   F = fft2(matrizFrames(:,:,i));

      for j=1:23

         z1 = real(F(j,1));
         z2 = real(F(j,2));
         z3 = real(F(j,3));
         z4 = real(F(j,4));
         z5 = real(F(j,5));
         z6 = real(F(j,6));
         z7 = real(F(j,7));
         z8 = real(F(j,8));

         fprintf(fileID, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f', ...
         z1, z2, z3, z4, z5, z6, z7, z8);


         z1 = imag(F(j,1));
         z2 = imag(F(j,2));
         z3 = imag(F(j,3));
         z4 = imag(F(j,4));
         z5 = imag(F(j,5));
         z6 = imag(F(j,6));
         z7 = imag(F(j,7));
         z8 = imag(F(j,8));

         fprintf(fileID2, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f', ...
         z1, z2, z3, z4, z5, z6, z7, z8);

         fprintf(fileID, '\n');
         fprintf(fileID2, '\n');
      end
   
   fprintf(fileID, '\n');
   fprintf(fileID2, '\n');

end


%
%
%namefile = strcat('verticales', '21');
%fileID = fopen(namefile, 'w');
%for i=1:length(listaVert21)
%   fprintf(fileID, '%f\t %f\t %f\t\n', ...
%   arrVert21(i,1), arrVert21(i, 2), MA(listaVert21(i), 3));
%end
%fprintf(fileID, '\n');
%fclose(fileID);
%v
%
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



