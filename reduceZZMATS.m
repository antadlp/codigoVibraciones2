function B = reduceZZMATS(A,s)

m = length(A(1,1,:));

j=1
for i=1:1:m
   B(:,:,j) = A(:,:,i);
   j = j + 1;
end

save(s,'B')



