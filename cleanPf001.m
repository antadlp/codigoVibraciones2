clear all
load('Pf-B12N12a.mat')

for i=1:1

%   Pfn0X = Pf(Pf(:,1,i)~=0)
%   Pfn0Y = Pf(Pf(:,2,i)~=0)

   Pfn0X = Pf(:,1,i);
   Pfn0Y = Pf(:,2,i);

   for j=1:length(Pfn0X)

      PfCleanX(j,1) = j;
      PfCleanX(j,2) = Pfn0X(j);

   end

   for j=1:length(Pfn0Y)

      PfCleanY(j,1) = j;
      PfCleanY(j,2) = Pfn0Y(j);

   end

   PfCX(:,1) = PfCleanX(:,1);
   PfCX(:,2) = PfCleanX(:,2);

   PfCY(:,1) = PfCleanY(:,1);
   PfCY(:,2) = PfCleanY(:,2);

end

for i=2:length(Pf(1,1,:))

%   Pfn0X = Pf(Pf(:,1,i)~=0)
%   Pfn0Y = Pf(Pf(:,2,i)~=0)

   Pfn0X = Pf(:,1,i);
   Pfn0Y = Pf(:,2,i);

   for j=1:length(Pfn0X)

      PfCleanX(j,1) = j;
      PfCleanX(j,2) = Pfn0X(j);

   end

   for j=1:length(Pfn0Y)

      PfCleanY(j,1) = j;
      PfCleanY(j,2) = Pfn0Y(j);

   end

   PfCX = cat(1, PfCX, PfCleanX);
   PfCY = cat(1, PfCY, PfCleanY);

end

save('PfCX-B12N12a', 'PfCX')
save('PfCY-B12N12a', 'PfCY')










