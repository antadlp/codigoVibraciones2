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


cutA = 101;
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

   
   [mx(1,:) IeGP] = getMxfreq(ffteGP, 10);
   [mx(2,:) IeBN] = getMxfreq(ffteBN, 10);
   [mx(3,:) IeB23N23] = getMxfreq(ffteB23N23, 10);
   [mx(4,:) IeB2N2] = getMxfreq(ffteB2N2, 10);
   [mx(5,:) IeB3N3] = getMxfreq(ffteB3N3, 10);
   [mx(6,:) IeB6N6] = getMxfreq(ffteB6N6, 10);
   [mx(7,:) IeB12N12a] = getMxfreq(ffteB12N12a, 10);
   [mx(8,:) IeB12N12b] = getMxfreq(ffteB12N12b, 10);
   [mx(9,:) IeB12N12c] = getMxfreq(ffteB12N12c, 10);
   [mx(10,:) IeB12N12d] = getMxfreq(ffteB12N12d, 10);
   [mx(11,:) IeB12N12e] = getMxfreq(ffteB12N12e, 10);

   IeGP = (IeGP-1)/5000;
   





filename='tablaLatexFreq.tex';
fid = fopen(filename, 'w');
%fprintf(fid,'algoalgo\\');
%ffprintf(fid,'%\documentclass[a4paper,10pt,twoside]{report}');
fprintf(fid,'\n');
fprintf(fid,'\\documentclass[a4paper, landscape]{article}%');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{geometry}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{lipsum}');
fprintf(fid,'\n');
fprintf(fid,'\n');
%fprintf(fid,'\\newgeometry{');
%fprintf(fid,'\n');
%fprintf(fid,'    top=1in,');
%fprintf(fid,'\n');
%fprintf(fid,'    bottom=1in,');
%fprintf(fid,'\n');
%fprintf(fid,'    outer=1in,');
%fprintf(fid,'\n');
%fprintf(fid,'    inner=1.5in,');
%fprintf(fid,'\n');
%fprintf(fid,'}');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage[utf8]{inputenc}');
fprintf(fid,'\n');
fprintf(fid,'%\\usepackage[Lenny]{fncychap}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage {amsmath}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage {amssymb}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage [usenames]{color}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage {graphicx}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{float}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage[rftl]{floatflt}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage[T1]{fontenc}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{units}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage[spanish]{babel}');
fprintf(fid,'\n');
fprintf(fid,'\\numberwithin{equation}{section}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{booktabs}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{nopageno}');
fprintf(fid,'\n');
fprintf(fid,'%\\usepackage{txfonts}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{graphics}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{tabularx}');
fprintf(fid,'\n');
fprintf(fid,'%\\usepackage{bibunits}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{bibentry}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{rotating}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{textcomp}');
fprintf(fid,'\n');
fprintf(fid,'\\definecolor{ColorLinks}{RGB}{46,139,87}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{titlesec, blindtext, color}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{url}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage[pdftex]{hyperref}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{makeidx}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{caption}');
fprintf(fid,'\n');
fprintf(fid,'\\makeindex');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{appendix}');
fprintf(fid,'\n');
fprintf(fid,'%\\usepackage{sectsty}');
fprintf(fid,'\n');
fprintf(fid,'%\\usepackage[super,square]{natbib}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{footnote}');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\\hypersetup{');
fprintf(fid,'\n');
fprintf(fid,'colorlinks,%');
fprintf(fid,'\n');
fprintf(fid,'citecolor=black,%');
fprintf(fid,'\n');
fprintf(fid,'filecolor=black,%');
fprintf(fid,'\n');
fprintf(fid,'linkcolor=black,%');
fprintf(fid,'\n');
fprintf(fid,'urlcolor=black');
fprintf(fid,'\n');
fprintf(fid,'}');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\\definecolor{RojoTitulo}{RGB}{219,32,32}');
fprintf(fid,'\n');
fprintf(fid,'\\definecolor{AzulSeccion}{RGB}{42,76,189}');
fprintf(fid,'\n');
fprintf(fid,'\\definecolor{VerdeSubseccion}{RGB}{46,139,87}');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{titlesec}');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\\titleformat*{\\section}{\\normalfont\\Large\\bfseries\\color{blue}}');
fprintf(fid,'\n');
fprintf(fid,'\\titleformat*{\\subsection}{\\normalfont\\large\\bfseries\\color{VerdeSubseccion}}');
fprintf(fid,'\n');
fprintf(fid,'\\titleformat*{\\subsubsection}{\\normalfont\\normalsize\\bfseries\\color{black}}');
fprintf(fid,'\n');
fprintf(fid,'\\definecolor{MyBlack}     {cmyk} {0.0,0.0,0.0,1.0}');
fprintf(fid,'\n');
fprintf(fid,'\\definecolor{MyBlue}      {cmyk} {1.0,0.3,0.0,0.6}');
fprintf(fid,'\n');
fprintf(fid,'\\definecolor{MyGreen}     {cmyk} {0.8,0.0,0.4,0.5}');
fprintf(fid,'\n');
fprintf(fid,'\\definecolor{MyLightGray} {cmyk} {0.0,0.0,0.0,0.1}');
fprintf(fid,'\n');
fprintf(fid,'\\definecolor{MyRed}       {cmyk} {0.0,0.6,1.0,0.0}');
fprintf(fid,'\n');
fprintf(fid,'\\definecolor{MyPink}      {cmyk} {0.0,0.8,0.0,0.0}');
fprintf(fid,'\n');
fprintf(fid,'\\definecolor{MyWhite}     {cmyk} {0.0,0.0,0.0,0.0}');
fprintf(fid,'\n');
fprintf(fid,'\\definecolor{MyYellow}    {cmyk} {0.0,0.0,0.8,0.2}');
fprintf(fid,'\n');
fprintf(fid,'\\definecolor{Mytry}       {cmyk} {0.3,0.1,0.3,0.0}');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\\renewcommand\\spanishtablename{\\scriptsize{TABLA}}');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\\begin{document}');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\n');


fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\\begin{table}[!htb]');
fprintf(fid,'\n');
fprintf(fid,'\\captionsetup{width=0.7\\textwidth, font={singlespacing}}');
fprintf(fid,'\n');
fprintf(fid,'\\begin{center}');
fprintf(fid,'\n');

% GP BN 

fprintf(fid,'\\begin{tabular} {c c c c c c c c c c c c c c');
fprintf(fid, ' c c c c c c c}');
fprintf(fid,'\n');
fprintf(fid,'\\multicolumn {11} {l}');
fprintf(fid, ' {{\\bf Tabla .} {Frequencias'); 
fprintf(fid,' de Mallas de Grafeno}}');
fprintf(fid,'\\\\');
fprintf(fid,'\n');
fprintf(fid,'\\\\');
fprintf(fid,'\n');
fprintf(fid,'\\hline \\hline');
fprintf(fid,'\n');
fprintf(fid,'GP && BN && B23N23 && B3N2 && B3N3 && B6N6'); 
fprintf(fid,' && B12N12a && B12N12b && B12N12c && B12N12d && B12N12e \\\\');
fprintf(fid,'\n');
%fprintf(fid,'$nm^3$ && KHartree && MBar && KBar\\\\');
%fprintf(fid,'\n');
fprintf(fid,'\\cline{1-1} \\cline{3-3} \\cline{5-5}');
fprintf(fid,'\\cline{7-7} \\cline{9-9} \\cline{11-11}');
fprintf(fid,'\\cline{13-13} \\cline{15-15} \\cline{17-17}');
fprintf(fid,'\\cline{19-19} \\cline{21-21}');
fprintf(fid,'\n');

for i=1:10
   fprintf(fid,'%0.3d && %d && %d && %d && %d && %d && %d && %d && %d && %d && %d\\\\', ...
   [IeGP(i)'; IeBN(i)'; IeB23N23(i)'; IeB3N3(i)'; IeB6N6(i)';...
   IeB12N12a(i)'; IeB12N12b(i)'; IeB12N12c(i)'; IeB12N12c(i)';...
   IeB12N12d(i)'; IeB12N12e(i)']);
   fprintf(fid,'\n');
end

fprintf(fid,'\\hline');
fprintf(fid,'\n');
fprintf(fid,'\\end{tabular}');
fprintf(fid,'\n');
fprintf(fid,'\\caption{\\footnotesize{Mallas');
%fprintf(fid,' energía y presión');
%fprintf(fid,' de las moleculas bajo confinamiento');
%fprintf(fid,' a $T= 350\\ K$.$^\\ast$}}');
%fprintf(fid,'\n');
fprintf(fid,'\\label{dosdos}}}');
fprintf(fid,'\n');
fprintf(fid,'\\end{center}');
fprintf(fid,'\n');
fprintf(fid,'\\end{table}');
fprintf(fid,'\n');
fprintf(fid,'\n');

fprintf(fid,'\\end{document}');
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



filename='tablaLatexFreq2.tex';
fid = fopen(filename, 'w');
%fprintf(fid,'algoalgo\\');
%ffprintf(fid,'%\documentclass[a4paper,10pt,twoside]{report}');
fprintf(fid,'\n');
fprintf(fid,'\\documentclass[a4paper, landscape]{article}%');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{geometry}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{lipsum}');
fprintf(fid,'\n');
fprintf(fid,'\n');
%fprintf(fid,'\\newgeometry{');
%fprintf(fid,'\n');
%fprintf(fid,'    top=1in,');
%fprintf(fid,'\n');
%fprintf(fid,'    bottom=1in,');
%fprintf(fid,'\n');
%fprintf(fid,'    outer=1in,');
%fprintf(fid,'\n');
%fprintf(fid,'    inner=1.5in,');
%fprintf(fid,'\n');
%fprintf(fid,'}');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage[utf8]{inputenc}');
fprintf(fid,'\n');
fprintf(fid,'%\\usepackage[Lenny]{fncychap}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage {amsmath}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage {amssymb}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage [usenames]{color}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage {graphicx}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{float}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage[rftl]{floatflt}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage[T1]{fontenc}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{units}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage[spanish]{babel}');
fprintf(fid,'\n');
fprintf(fid,'\\numberwithin{equation}{section}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{booktabs}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{nopageno}');
fprintf(fid,'\n');
fprintf(fid,'%\\usepackage{txfonts}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{graphics}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{tabularx}');
fprintf(fid,'\n');
fprintf(fid,'%\\usepackage{bibunits}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{bibentry}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{rotating}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{textcomp}');
fprintf(fid,'\n');
fprintf(fid,'\\definecolor{ColorLinks}{RGB}{46,139,87}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{titlesec, blindtext, color}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{url}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage[pdftex]{hyperref}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{makeidx}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{caption}');
fprintf(fid,'\n');
fprintf(fid,'\\makeindex');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{appendix}');
fprintf(fid,'\n');
fprintf(fid,'%\\usepackage{sectsty}');
fprintf(fid,'\n');
fprintf(fid,'%\\usepackage[super,square]{natbib}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{footnote}');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\\hypersetup{');
fprintf(fid,'\n');
fprintf(fid,'colorlinks,%');
fprintf(fid,'\n');
fprintf(fid,'citecolor=black,%');
fprintf(fid,'\n');
fprintf(fid,'filecolor=black,%');
fprintf(fid,'\n');
fprintf(fid,'linkcolor=black,%');
fprintf(fid,'\n');
fprintf(fid,'urlcolor=black');
fprintf(fid,'\n');
fprintf(fid,'}');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\\definecolor{RojoTitulo}{RGB}{219,32,32}');
fprintf(fid,'\n');
fprintf(fid,'\\definecolor{AzulSeccion}{RGB}{42,76,189}');
fprintf(fid,'\n');
fprintf(fid,'\\definecolor{VerdeSubseccion}{RGB}{46,139,87}');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{titlesec}');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\\titleformat*{\\section}{\\normalfont\\Large\\bfseries\\color{blue}}');
fprintf(fid,'\n');
fprintf(fid,'\\titleformat*{\\subsection}{\\normalfont\\large\\bfseries\\color{VerdeSubseccion}}');
fprintf(fid,'\n');
fprintf(fid,'\\titleformat*{\\subsubsection}{\\normalfont\\normalsize\\bfseries\\color{black}}');
fprintf(fid,'\n');
fprintf(fid,'\\definecolor{MyBlack}     {cmyk} {0.0,0.0,0.0,1.0}');
fprintf(fid,'\n');
fprintf(fid,'\\definecolor{MyBlue}      {cmyk} {1.0,0.3,0.0,0.6}');
fprintf(fid,'\n');
fprintf(fid,'\\definecolor{MyGreen}     {cmyk} {0.8,0.0,0.4,0.5}');
fprintf(fid,'\n');
fprintf(fid,'\\definecolor{MyLightGray} {cmyk} {0.0,0.0,0.0,0.1}');
fprintf(fid,'\n');
fprintf(fid,'\\definecolor{MyRed}       {cmyk} {0.0,0.6,1.0,0.0}');
fprintf(fid,'\n');
fprintf(fid,'\\definecolor{MyPink}      {cmyk} {0.0,0.8,0.0,0.0}');
fprintf(fid,'\n');
fprintf(fid,'\\definecolor{MyWhite}     {cmyk} {0.0,0.0,0.0,0.0}');
fprintf(fid,'\n');
fprintf(fid,'\\definecolor{MyYellow}    {cmyk} {0.0,0.0,0.8,0.2}');
fprintf(fid,'\n');
fprintf(fid,'\\definecolor{Mytry}       {cmyk} {0.3,0.1,0.3,0.0}');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\\renewcommand\\spanishtablename{\\scriptsize{TABLA}}');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\\begin{document}');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\n');


fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\\begin{table}[!htb]');
fprintf(fid,'\n');
fprintf(fid,'\\captionsetup{width=0.7\\textwidth, font={singlespacing}}');
fprintf(fid,'\n');
fprintf(fid,'\\begin{center}');
fprintf(fid,'\n');

% GP BN 

fprintf(fid,'\\begin{tabular} {c c c c c c c c c c c}');
%fprintf(fid, ' c c c c c c c}');
fprintf(fid,'\n');
fprintf(fid,'\\multicolumn {11} {l}');
fprintf(fid, ' {{\\bf Tabla .} {Frequencias'); 
fprintf(fid,' de Mallas de Grafeno}}');
fprintf(fid,'\\\\');
fprintf(fid,'\n');
fprintf(fid,'\\\\');
fprintf(fid,'\n');
fprintf(fid,'\\hline \\hline');
fprintf(fid,'\n');
fprintf(fid,'GP && B2N2 && B3N3 && B6N6 && B23N23 && BN\\\\'); 
%fprintf(fid,' &&  &&  &&  &&  &&  \\\\');
fprintf(fid,'\n');
%fprintf(fid,'$nm^3$ && KHartree && MBar && KBar\\\\');
%fprintf(fid,'\n');
fprintf(fid,'\\cline{1-1} \\cline{3-3} \\cline{5-5}');
fprintf(fid,'\\cline{7-7} \\cline{9-9} \\cline{11-11}');
%fprintf(fid,'\\cline{13-13} \\cline{15-15} \\cline{17-17}');
%fprintf(fid,'\\cline{19-19} \\cline{21-21}');
fprintf(fid,'\n');

for i=1:10
   fprintf(fid,'%0.3d && %d && %d && %d && %d && %d \\\\', ...
   [IeGP(i)'; IeB2N2(i)'; IeB3N3(i)'; IeB6N6(i)'; IeB23N23(i)';...
   IeBN(i)']);
   fprintf(fid,'\n');
end

fprintf(fid,'\\hline');
fprintf(fid,'\n');
fprintf(fid,'\\end{tabular}');
fprintf(fid,'\n');
fprintf(fid,'\\caption{\\footnotesize{Mallas');
%fprintf(fid,' energía y presión');
%fprintf(fid,' de las moleculas bajo confinamiento');
%fprintf(fid,' a $T= 350\\ K$.$^\\ast$}}');
%fprintf(fid,'\n');
fprintf(fid,'\\label{dosdos}}}');
fprintf(fid,'\n');
fprintf(fid,'\\end{center}');
fprintf(fid,'\n');
fprintf(fid,'\\end{table}');
fprintf(fid,'\n');
fprintf(fid,'\n');


%
%
% tabla abcde
%
%

fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\\begin{table}[!htb]');
fprintf(fid,'\n');
fprintf(fid,'\\captionsetup{width=0.7\\textwidth, font={singlespacing}}');
fprintf(fid,'\n');
fprintf(fid,'\\begin{center}');
fprintf(fid,'\n');
fprintf(fid,'\\begin{tabular} {c c c c c c c c c}');
%fprintf(fid, ' c c c c c c c}');
fprintf(fid,'\n');
fprintf(fid,'\\multicolumn {9} {l}');
fprintf(fid, ' {{\\bf Tabla .} {Frequencias'); 
fprintf(fid,' de Mallas de Grafeno}}');
fprintf(fid,'\\\\');
fprintf(fid,'\n');
fprintf(fid,'\\\\');
fprintf(fid,'\n');
fprintf(fid,'\\hline \\hline');
fprintf(fid,'\n');
%fprintf(fid,'GP && BN && B23N23 && B3N2 && B3N3 && B6N6'); 
fprintf(fid,'B12N12a && B12N12b && B12N12c && B12N12d && B12N12e \\\\');
fprintf(fid,'\n');
%fprintf(fid,'$nm^3$ && KHartree && MBar && KBar\\\\');
%fprintf(fid,'\n');
fprintf(fid,'\\cline{1-1} \\cline{3-3} \\cline{5-5}');
fprintf(fid,'\\cline{7-7} \\cline{9-9}');
%\\cline{11-11}');
%fprintf(fid,'\\cline{13-13} \\cline{15-15} \\cline{17-17}');
%fprintf(fid,'\\cline{19-19} \\cline{21-21}');
fprintf(fid,'\n');

for i=1:10
   fprintf(fid,'%d && %d && %d && %d && %d \\\\', ...
   [IeB12N12a(i)'; IeB12N12b(i)'; IeB12N12c(i)'; IeB12N12d(i)';...
   IeB12N12e(i)']);   
   fprintf(fid,'\n');
end

fprintf(fid,'\\hline');
fprintf(fid,'\n');
fprintf(fid,'\\end{tabular}');
fprintf(fid,'\n');
fprintf(fid,'\\caption{\\footnotesize{Mallas');
%fprintf(fid,' energía y presión');
%fprintf(fid,' de las moleculas bajo confinamiento');
%fprintf(fid,' a $T= 350\\ K$.$^\\ast$}}');
%fprintf(fid,'\n');
fprintf(fid,'\\label{dosdos}}}');
fprintf(fid,'\n');
fprintf(fid,'\\end{center}');
fprintf(fid,'\n');
fprintf(fid,'\\end{table}');
fprintf(fid,'\n');
fprintf(fid,'\n');



fprintf(fid,'\\end{document}');
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




filename='tablaLatexFreq3.tex';
fid = fopen(filename, 'w');
%fprintf(fid,'algoalgo\\');
%ffprintf(fid,'%\documentclass[a4paper,10pt,twoside]{report}');
fprintf(fid,'\n');
fprintf(fid,'\\documentclass[a4paper, landscape]{article}%');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{geometry}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{lipsum}');
fprintf(fid,'\n');
fprintf(fid,'\n');
%fprintf(fid,'\\newgeometry{');
%fprintf(fid,'\n');
%fprintf(fid,'    top=1in,');
%fprintf(fid,'\n');
%fprintf(fid,'    bottom=1in,');
%fprintf(fid,'\n');
%fprintf(fid,'    outer=1in,');
%fprintf(fid,'\n');
%fprintf(fid,'    inner=1.5in,');
%fprintf(fid,'\n');
%fprintf(fid,'}');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage[utf8]{inputenc}');
fprintf(fid,'\n');
fprintf(fid,'%\\usepackage[Lenny]{fncychap}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage {amsmath}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage {amssymb}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage [usenames]{color}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage {graphicx}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{float}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage[rftl]{floatflt}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage[T1]{fontenc}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{units}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage[spanish]{babel}');
fprintf(fid,'\n');
fprintf(fid,'\\numberwithin{equation}{section}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{booktabs}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{nopageno}');
fprintf(fid,'\n');
fprintf(fid,'%\\usepackage{txfonts}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{graphics}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{tabularx}');
fprintf(fid,'\n');
fprintf(fid,'%\\usepackage{bibunits}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{bibentry}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{rotating}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{textcomp}');
fprintf(fid,'\n');
fprintf(fid,'\\definecolor{ColorLinks}{RGB}{46,139,87}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{titlesec, blindtext, color}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{url}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage[pdftex]{hyperref}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{makeidx}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{caption}');
fprintf(fid,'\n');
fprintf(fid,'\\makeindex');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{appendix}');
fprintf(fid,'\n');
fprintf(fid,'%\\usepackage{sectsty}');
fprintf(fid,'\n');
fprintf(fid,'%\\usepackage[super,square]{natbib}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{footnote}');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\\hypersetup{');
fprintf(fid,'\n');
fprintf(fid,'colorlinks,%');
fprintf(fid,'\n');
fprintf(fid,'citecolor=black,%');
fprintf(fid,'\n');
fprintf(fid,'filecolor=black,%');
fprintf(fid,'\n');
fprintf(fid,'linkcolor=black,%');
fprintf(fid,'\n');
fprintf(fid,'urlcolor=black');
fprintf(fid,'\n');
fprintf(fid,'}');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\\definecolor{RojoTitulo}{RGB}{219,32,32}');
fprintf(fid,'\n');
fprintf(fid,'\\definecolor{AzulSeccion}{RGB}{42,76,189}');
fprintf(fid,'\n');
fprintf(fid,'\\definecolor{VerdeSubseccion}{RGB}{46,139,87}');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{titlesec}');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\\titleformat*{\\section}{\\normalfont\\Large\\bfseries\\color{blue}}');
fprintf(fid,'\n');
fprintf(fid,'\\titleformat*{\\subsection}{\\normalfont\\large\\bfseries\\color{VerdeSubseccion}}');
fprintf(fid,'\n');
fprintf(fid,'\\titleformat*{\\subsubsection}{\\normalfont\\normalsize\\bfseries\\color{black}}');
fprintf(fid,'\n');
fprintf(fid,'\\definecolor{MyBlack}     {cmyk} {0.0,0.0,0.0,1.0}');
fprintf(fid,'\n');
fprintf(fid,'\\definecolor{MyBlue}      {cmyk} {1.0,0.3,0.0,0.6}');
fprintf(fid,'\n');
fprintf(fid,'\\definecolor{MyGreen}     {cmyk} {0.8,0.0,0.4,0.5}');
fprintf(fid,'\n');
fprintf(fid,'\\definecolor{MyLightGray} {cmyk} {0.0,0.0,0.0,0.1}');
fprintf(fid,'\n');
fprintf(fid,'\\definecolor{MyRed}       {cmyk} {0.0,0.6,1.0,0.0}');
fprintf(fid,'\n');
fprintf(fid,'\\definecolor{MyPink}      {cmyk} {0.0,0.8,0.0,0.0}');
fprintf(fid,'\n');
fprintf(fid,'\\definecolor{MyWhite}     {cmyk} {0.0,0.0,0.0,0.0}');
fprintf(fid,'\n');
fprintf(fid,'\\definecolor{MyYellow}    {cmyk} {0.0,0.0,0.8,0.2}');
fprintf(fid,'\n');
fprintf(fid,'\\definecolor{Mytry}       {cmyk} {0.3,0.1,0.3,0.0}');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\\renewcommand\\spanishtablename{\\scriptsize{TABLA}}');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\\begin{document}');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\n');


fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\\begin{table}[!htb]');
fprintf(fid,'\n');
fprintf(fid,'\\captionsetup{width=0.7\\textwidth, font={singlespacing}}');
fprintf(fid,'\n');
fprintf(fid,'\\begin{center}');
fprintf(fid,'\n');
fprintf(fid,'\\begin{tabular} {c c c c c c c c c c c c c c');
fprintf(fid, ' c c c c c c c}');
fprintf(fid,'\n');
fprintf(fid,'\\multicolumn {11} {l}');
fprintf(fid, ' {{\\bf Tabla .} {Frequencias'); 
fprintf(fid,' de Mallas de Grafeno}}');
fprintf(fid,'\\\\');
fprintf(fid,'\n');
fprintf(fid,'\\\\');
fprintf(fid,'\n');
fprintf(fid,'\\hline \\hline');
fprintf(fid,'\n');
fprintf(fid,'GP && BN && B23N23 && B3N2 && B3N3 && B6N6'); 
fprintf(fid,' && B12N12a && B12N12b && B12N12c && B12N12d && B12N12e \\\\');
fprintf(fid,'\n');
%fprintf(fid,'$nm^3$ && KHartree && MBar && KBar\\\\');
%fprintf(fid,'\n');
fprintf(fid,'\\cline{1-1} \\cline{3-3} \\cline{5-5}');
fprintf(fid,'\\cline{7-7} \\cline{9-9} \\cline{11-11}');
fprintf(fid,'\\cline{13-13} \\cline{15-15} \\cline{17-17}');
fprintf(fid,'\\cline{19-19} \\cline{21-21}');
fprintf(fid,'\n');

for i=1:10
   fprintf(fid,'%d && %d && %d && %d && %d && %d && %d && %d && %d && %d && %d\\\\', ...
   [IeGP(i)'; IeBN(i)'; IeB23N23(i)'; IeB3N3(i)'; IeB6N6(i)';...
   IeB12N12a(i)'; IeB12N12b(i)'; IeB12N12c(i)'; IeB12N12c(i)';...
   IeB12N12d(i)'; IeB12N12e(i)']);
   fprintf(fid,'\n');
end

fprintf(fid,'\\hline');
fprintf(fid,'\n');
fprintf(fid,'\\end{tabular}');
fprintf(fid,'\n');
fprintf(fid,'\\caption{\\footnotesize{Mallas');
%fprintf(fid,' energía y presión');
%fprintf(fid,' de las moleculas bajo confinamiento');
%fprintf(fid,' a $T= 350\\ K$.$^\\ast$}}');
%fprintf(fid,'\n');
fprintf(fid,'\\label{dosdos}}}');
fprintf(fid,'\n');
fprintf(fid,'\\end{center}');
fprintf(fid,'\n');
fprintf(fid,'\\end{table}');
fprintf(fid,'\n');
fprintf(fid,'\n');

fprintf(fid,'\\end{document}');



