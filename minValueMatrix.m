function y = minValueMatrix(A)

[M N] = size(A);

for i=1:M
   B(i) = min(A(i,:));
end

y = min(B);



