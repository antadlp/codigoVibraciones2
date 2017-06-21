function zz = fft2UV(Z)

% m renglones => x
% n columnas => y

[M N] = size(Z);

MN = 1/(M*N);

for u=1:M
   for v=1:N
%      fprintf('\n\n')
%      fprintf('u: %f\n', u)
%      fprintf('v: %f\n', v)
      zz(u,v) = 0;
      for x=1:M
         for y=1:N
%            S = ['Zi: ', num2str(Z(u,v))];
%            disp(S)
            zz(u,v) = zz(u,v) + MN*Z(x,y)*(cos(2*pi*((u-1)*(x-1)/M+(v-1)*(y-1)/N))+j*sin(2*pi*((u-1)*(x-1)/M+(v-1)*(y-1)/N)));
         end
      end
   end
end








