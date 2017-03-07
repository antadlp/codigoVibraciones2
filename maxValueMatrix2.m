function [y II JJ] = maxValueMatrix(A)

[M N] = size(A);

for i=1:M
   
   [B(i) I(i)] = max(A(i,:));
   
end

[y J] = max(B);

%J columna
%I fila

II = J;
JJ = I(J);





