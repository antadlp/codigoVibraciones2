function mostrarGraficas


eGP = load('datosMallaPython/eGP.dat');
%eBN = load('datosMallaPython/eBN.dat');
%eB23N23 = load('datosMallaPython/eB23N23.dat');
%eB2N2 = load('datosMallaPython/eB2N2.dat');
%eB3N3 = load('datosMallaPython/eB3N3.dat');
%eB6N6 = load('datosMallaPython/eB6N6.dat');
eB12N12a = load('datosMallaPython/eB12N12a.dat');
%eB12N12b = load('datosMallaPython/eB12N12b.dat');
%eB12N12c = load('datosMallaPython/eB12N12c.dat');
%eB12N12d = load('datosMallaPython/eB12N12d.dat');
%eB12N12e = load('datosMallaPython/eB12N12e.dat');
%

cutA = 100;
L = length(eGP);

%ffteGP = fullFFT001(eGP(cutA:L));
%ffteBN = fullFFT001(eBN(cutA:L));
%ffteB23N23 = fullFFT001(eB23N23(cutA:L));
%ffteB2N2 = fullFFT001(eB2N2(cutA:L));
%ffteB3N3 = fullFFT001(eB3N3(cutA:L));
%ffteB6N6 = fullFFT001(eB6N6(cutA:L));
ffteB12N12a = fullFFT001(eB12N12a(cutA:L));
%ffteB12N12b = fullFFT001(eB12N12b(cutA:L));
%ffteB12N12c = fullFFT001(eB12N12c(cutA:L));
%ffteB12N12e = fullFFT001(eB12N12e(cutA:L));


%ffbeGP = fullFBT001(ffteGP);
%ffbeBN = fullFBT001(ffteBN);
%ffbeB23N23 = fullFBT001(ffteB23N23);
%ffbeB2N2 = fullFBT001(ffteB2N2);
%ffbeB3N3 = fullFBT001(ffteB3N3);
%ffbeB6N6 = fullFBT001(ffteB6N6);
ffbeB12N12a = fullFBT001(ffteB12N12a);
%ffbeB12N12b = fullFBT001(ffteB12N12b);
%ffbeB12N12c = fullFBT001(ffteB12N12c);
%ffbeB12N12e = fullFBT001(ffteB12N12e);
%

t = cutA:L;

%
%figure('Name', 'eGP')
%plot(t, real(ffbeGP), 'LineWidth', 2)
%hold on
%plot(t, eGP(cutA:L), '--r')
%

%figure('Name', 'eBN')
%plot(t, real(ffbeBN), 'LineWidth', 2)
%hold on
%plot(t, eBN(cutA:L), '--r')


%figure('Name', 'eB23N23')
%plot(t, real(ffbeB23N23), 'LineWidth', 2)
%hold on
%plot(t, eB23N23(cutA:L), '--r')
%
%
%figure('Name', 'eB2N2')
%plot(t, real(ffbeB2N2), 'LineWidth', 2)
%hold on
%plot(t, eB2N2(cutA:L), '--r')
%
%
%figure('Name', 'eB3N3')
%plot(t, real(ffbeB3N3), 'LineWidth', 2)
%hold on
%plot(t, eB3N3(cutA:L), '--r')
%
%
%figure('Name', 'eB6N6')
%plot(t, real(ffbeB6N6), 'LineWidth', 2)
%hold on
%plot(t, eB6N6(cutA:L), '--r')
%

figure('Name', 'eB12N12a')
plot(t, real(ffbeB12N12a), 'LineWidth', 2)
hold on
plot(t, eB12N12a(cutA:L), '--r')

%
%figure('Name', 'eB12N12b')
%plot(t, real(ffbeB12N12b), 'LineWidth', 2)
%hold on
%plot(t, eB12N12b(cutA:L), '--r')
%
%
%figure('Name', 'eB12N12c')
%plot(t, real(ffbeB12N12c), 'LineWidth', 2)
%hold on
%plot(t, eB12N12c(cutA:L), '--r')
%
%
%figure('Name', 'eB12N12d')
%plot(t, real(ffbeB12N12d), 'LineWidth', 2)
%hold on
%plot(t, eB12N12d(cutA:L), '--r')
%
%

[mx I] = getMxfreq(ffteB12N12a, 20);
mfbeB12N12a = manualFBT001(ffteB12N12a, I, 1);
figure('Name', 'freq eB12N12a')
plot(t, eB12N12a(cutA:L), 'LineWidth', 2)
hold on
plot(t, real(mfbeB12N12a), 'r')


cter = 1;
Lp = length(t)

j=1;
for i=1:100:Lp

   [mx I] = getMxfreq(ffteB12N12a, i);
   mfbeB12N12a(j,:) = real(manualFBT001(ffteB12N12a, I, 1));
   %figure(2), plot(t, real(mfbeB12N12a), 'r')
   j = j + 1;
   i

end

save('terB12N12a.mat', 'mfbeB12N12a')


%save('Pf-B12N12b','Pf')
   







