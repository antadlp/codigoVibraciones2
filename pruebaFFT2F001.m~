N = 5;
x=1:N;
y = x;

fr1 = 0.1;
fr2 = 0.5;


for r=1:N
   for s=1:N
      z(r,s) = 2*sin(2*pi*fr1*x(r)) + 3*sin(2*pi*fr2*y(s));
   end
end

zf = fft2XY(z);
zr = fft2UV(zf);



for r=1:N
   for s=1:N
      fprintf('%f\t%f\n',real(z(r,s)), real(zr(r,s)));
   end
end


for r=1:N
   for s=1:N
      fprintf('%f\t%f\n',imag(z(r,s)), imag(zr(r,s)));
   end
end


