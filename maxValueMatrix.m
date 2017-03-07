function y = maxValueMatrix(A)

[M N] = size(A);

for i=1:M
   B(i) = max(A(i,:));
end

y = max(B);



