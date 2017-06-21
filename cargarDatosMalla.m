%function AnimacionD1(t,x,p)
function cargarDatosMalla

t = 1:3100;
for i=1:length(t)

   s2 = '/home/toshiba/out003/mallaInter';
   filename2 = strcat(s2, int2str(i));
   Zp(:,:,i) = importdata(filename2);
   i

end

save('zp-GP.mat', 'Zp');

fprintf('done\n');

