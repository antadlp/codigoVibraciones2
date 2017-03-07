A = importdata('nMA');
B = importdata('char-GP.dat');
C = importdata('frame1-mil-rota.dat');

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





   

