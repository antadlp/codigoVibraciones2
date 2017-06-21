function imprimeMatriz001

for i=1:5
   for j=1:6
      B(i,j) = (i*6 -6) + j;
   end
end




for i=1:3
   A(:,:,i) = B;
end



for i=1:3

   s = 'matrixS';
   filename = strcat(s, int2str(i));
   fileID = fopen(filename, 'w');

   for j = 1:length(A(:,1,1))

      fprintf(fileID, '%f   ', A(j,:,i));
      fprintf(fileID, '\n');
   end
end



      





