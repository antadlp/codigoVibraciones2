A = importdata('nMA');
B = importdata('char-GP.dat');
C = importdata('frame1-mil-rota.dat');
% Empato las posiciones de la malla 'arreglada' (malla preparada para ser
% interpolada, viene de movie-ejemplo) realizando la misma transformacion del
% programa rota123.f. De esta forma tengo manera de comparar las posiciones de
% la malla arreglada con la malla a analizar. No es necesario compara todo el
% movie-analizar, unicamente basta con el primer frame y asi tener el mapeo de
% los numeros cardinales entre las dos mallas.

% Por pasos seria:
% 1. Rotar frame1 del movie a analizar con rota123.f (ojo no es con rota123-2.f)
% 2. Cargar data de nMA (#cardinal atomo con posiciones x-y arregladas) X. La
% otra seria usar el scritpt que crea horizontales y verticales para cada movie
% a analizar, es decir cada movie a analizar tendria sus horizontales y
% verticales, es decir su malla-arreglada.

% Por pasos usando cada movie su datos de horizontales y verticales, seria:
% 1. Rotar frame1 del movie a analizar con rota123.f (ojo no es con rota123-2.f)
% 2. Generar datos horizontales y verticales para el movie de interes.
% 3. Generar malla arreglada .xyz, .dat.
% 4. Generar archivo de alturas.





fileID = fopen('mallanMA.xyz','w');
fileID2 = fopen('mallanMA.dat','w');
nAtoms = 142;

[A1 A2] = unique(A(:,1));

fprintf(fileID, '%s\n', int2str(nAtoms));
fprintf(fileID, '%s\t%d\n', 'frame', 1);

for i=1:(nAtoms)

   if (B{i}=='H')
      
      fprintf(fileID, '%c  %f  %f  %f\n', B{i}, C(i,1), C(i,2), 0.00);
      fprintf(fileID2, '%f  %f  %f\n', C(i,1), C(i,2), 0.00);

      i
   
   else
      j = i - 28

      fprintf(fileID, '%c  %f  %f  %f\n', B{i}, A(A2(j),2), A(A2(j),3), 0.00);
      fprintf(fileID2, '%f  %f  %f\n', A(A2(j),2), A(A2(j),3), 0.00);
   end

end

fclose(fileID);





   

