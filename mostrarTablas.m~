function mostrarGraficas


eGP = load('datosMallaPython/eGP.dat');
eBN = load('datosMallaPython/eBN.dat');
eB23N23 = load('datosMallaPython/eB23N23.dat');
eB2N2 = load('datosMallaPython/eB2N2.dat');
eB3N3 = load('datosMallaPython/eB3N3.dat');
eB6N6 = load('datosMallaPython/eB6N6.dat');
eB12N12a = load('datosMallaPython/eB12N12a.dat');
eB12N12b = load('datosMallaPython/eB12N12b.dat');
eB12N12c = load('datosMallaPython/eB12N12c.dat');
eB12N12d = load('datosMallaPython/eB12N12d.dat');
eB12N12e = load('datosMallaPython/eB12N12e.dat');


cutA = 100;
L = length(eGP);

ffteGP = fft(eGP(cutA:L));
ffteBN = fft(eBN(cutA:L));
ffteB23N23 = fft(eB23N23(cutA:L));
ffteB2N2 = fft(eB2N2(cutA:L));
ffteB3N3 = fft(eB3N3(cutA:L));
ffteB6N6 = fft(eB6N6(cutA:L));
ffteB12N12a = fft(eB12N12a(cutA:L));
ffteB12N12b = fft(eB12N12b(cutA:L));
ffteB12N12c = fft(eB12N12c(cutA:L));
ffteB12N12d = fft(eB12N12d(cutA:L));
ffteB12N12e = fft(eB12N12e(cutA:L));

numMallas = 11;


%[mx y] = getMxfreq(yf, n)

for i=1:7
   
   [mx IeGP] = getMxFreq(ffteGP, 10);
   [mx IeBN] = getMxFreq(ffteBN, 10);
   [mx IeB23N23] = getMxFreq(ffteB23N23, 10);
   [mx IeB2N2] = getMxFreq(ffteB2N2, 10);
   [mx IeB3N3] = getMxFreq(ffteB3N3, 10);
   [mx IeB6N6] = getMxFreq(ffteB6N6, 10);
   [mx IeB12N12a] = getMxFreq(ffteB12N12a, 10);





  







