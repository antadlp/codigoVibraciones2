function [MatrizDatos P cv] = procesoDatosAll001(opcion)

makeGnuplot = opcion(1);
makeLatexTable = opcion(2);

%%%%%%
%Entrada de Datos
%%%%%%

%Datos Generales
%beta1 = lsDInp(1);
%iPasoEnergia = lsDInp(2);
%ShowPlotAvRadius = lsDInp(3);
%ShowPlotMajorRadius = lsDInp(4);
%ShowPlotMajorRadiusConBetas = lsDInp(5);
%ShowPlotTemp = lsDInp(6);
%ShowPlotRmsTProm = lsDInp(7);
%ShowPlotRmsTrefAPartir = lsDInp(8);
%rmsBegin = lsDInp(9);
%Tref = lsDInp(10);
%ShowPlotRmsTRef = lsDInp(11);
%AbsDes = lsDInp(12);
%Tmax = lsDInp(13);
%Tmin = lsDInp(14);
%ShowPlotEspectro = lsDInp(15);
%CadaCuantos = lsDInp(16);
%ShowPlotPotEne = lsDInp(17);
%filtradoPorDelta = lsDInp(18);
%AbsDes2 = lsDInp(19);
%
%
% Ej = [beta1 iPasoEnergia ShowPlotAvRadius ShowPlotMajorRadius
%    ShowPlotMajorRadiusConBetas ShowPlotTemp ShowPlotRmsTProm 
%    ShowPlotRmsTrefAPartir rmsBegin Tref ShowPlotRmsTRef AbsDes Tmax Tmin
%    ShowPlotEspectro CadaCuantos ShowPlotPotEne filtradoPorDelta AbsDes2
%
%    ...]
%    
%Ej = [(1)8.0 (2)3000 (3)0 (4)1 (5)0 (6)1 (7)0 (8)0 (9)3000 (10)350 (11)1
%(12)35 (13)400 (14)200 (15)0 (16)10 (17)1 (18)1 (19)20]

%Ej = [8.0 3000 0 1 0 1 0 0 3000 350 1 35 400 200 0 10 1 1 20]


%
%lsDatos = [pasostotales, iPasoEnergia, numPasos, VolPromNm3, fTprom, ...
%fpromEne, PresDin, PresDinBar, FPromInt, FPromExt];
%

%Datos


%inpB6 = [6.0 100 0 0 0 1 0 0 100 350 0 35 400 200 0 10 0 1 100];
%B6 = procesoDatosOne001('cBeta7_25', inpB6);
%

Tmin = 300;
Tmax = 450;
AbsDes = 10;

inpB70 = [7.0 1000 0 0 0 0 0 0 3100 350 0 AbsDes Tmax Tmin 0 1 0 1 35];
B70 = procesoDatosOne003('Beta7_0', inpB70);
%3100

inpB725 = [7.25 1000 0 0 0 0 0 0 1000 350 0 AbsDes Tmax Tmin 0 1 0 1 35];
B725 = procesoDatosOne003('Beta7_25', inpB725);

inpB75 = [7.5 1000 0 0 0 0 0 0 3100 350 0 AbsDes Tmax Tmin 0 1 0 1 35];
B75 = procesoDatosOne003('Beta7_5', inpB75);

inpB775 = [7.75 1000 0 0 0 0 0 0 2700 350 0 AbsDes Tmax Tmin 0 1 0 1 35];
B775 = procesoDatosOne003('Beta7_75', inpB775);
%2700

inpB80 = [8.0 1000 0 0 0 0 0 0 4000 350 0 AbsDes Tmax Tmin 0 1 0 1 35];
B80 = procesoDatosOne003('Beta8_0', inpB80);
%4000
%
%inpB825 = [8.25 1000 0 0 0 1 0 0 3000 350 0 AbsDes Tmax Tmin 1 1 0 1 35];
%B825 = procesoDatosOne003('Beta8_25', inpB825);
%
%
%inpB85 = [8.5 1000 0 0 0 0 0 0 3000 350 0 35 400 200 0 10 0 1 3];
%B85 = procesoDatosOne001('Beta8_5', inpB85);
%

%MatrizDatos = [B70; B75; B775; B80; B825; B85];
MatrizDatos = [B70; B725; B75; B775; B80];

%VolPromNm3 = [B70(4) B775(4) B80(4) B825(4) B85(4)];
%fpromEne = [B70(6) B775(6) B80(6) B825(6) B85(6)];
%

%VolPromNm3 = [B70(4) B75(4) B775(4) B80(4) B825(4) B85(4)];
%fpromEne = [B70(6) B75(6) B775(6) B80(6) B825(6) B85(6)];
%

FPromInt(1) = B70(9);
FPromInt(2) = B725(9);
FPromInt(3) = B75(9);
FPromInt(4) = B775(9);
FPromInt(5) = B80(9);

FPromExt(1) = B70(10);
FPromExt(2) = B725(10);
FPromExt(3) = B75(10);
FPromExt(4) = B775(10);
FPromExt(5) = B80(10);


VolPromNm3 = [B70(4) B725(4) B75(4) B775(4) B80(4)];
fpromEne = [B70(6) B725(6) B75(6) B775(6) B80(6)];

tarE = 1206;

fpromEne = fpromEne + tarE;

[fitresult, gof] = createFitMatlabPrueba001(VolPromNm3, fpromEne);


cv = coeffvalues(fitresult);


a = cv(1);
b = cv(2);
c = cv(3);
d = cv(4);

V = VolPromNm3;
Vxi = min(V);
Vxf = max(V);

Vx = Vxi:.0001:Vxf;


P = -exp(d*V.^(2/3)).*((d*b./(3*V.^(1/3))) +  (d*c./(3*V.^(1/3)))  + ...
(c./(3*V.^(2/3))));

curvaFit = a+b*exp(d*(Vx.^(1/3)))+c*(Vx.^(1/3)).*exp(d*(Vx.^(1/3)));

%plot(Vx, curvaFit);
%
%[status,cmdout] = system('cd')
%dir2 ='/home/assusadmin/toshibaAResurrected/beta2/allpuntosDatos/scriptPruebaGnuplot.sh';
%[status,cmdout] = system(dir2)
%
Eh = 4.35974417*1e-18;
Nm3 = 1e-27;

P = P*Eh/Nm3;

P = P/(1e5);

Pdin = MatrizDatos(:,8);

Psw = P/(1e6);
Pdinsw = Pdin/(1e3);
fpromEnesw = fpromEne/(1e3);

filename='tablaDatosLatex.dat';
fid = fopen(filename, 'w');
fprintf(fid, '%f %f %.3f %.3f\n', [VolPromNm3; fpromEnesw; Psw; Pdinsw']);
fclose(fid);

filename='curvaFittingVvsE.dat';
fid = fopen(filename, 'w');
fprintf(fid, '%f %f\n', [Vx; curvaFit]);
fclose(fid);

filename='tablaDatosEcEstadoVvsEn.dat';
fid = fopen(filename, 'w');
fprintf(fid, '%f  %f\n', [VolPromNm3; fpromEne]);
fclose(fid);

filename='tablaLatex.tex';
fid = fopen(filename, 'w');
%fprintf(fid,'algoalgo\\');
%ffprintf(fid,'%\documentclass[a4paper,10pt,twoside]{report}');
fprintf(fid,'\n');
fprintf(fid,'\\documentclass[oneside]{book}%');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{geometry}');
fprintf(fid,'\n');
fprintf(fid,'\\usepackage{lipsum}');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\\newgeometry{');
fprintf(fid,'\n');
fprintf(fid,'    top=1in,');
fprintf(fid,'\n');
fprintf(fid,'    bottom=1in,');
fprintf(fid,'\n');
fprintf(fid,'    outer=1in,');
fprintf(fid,'\n');
fprintf(fid,'    inner=1.5in,');
fprintf(fid,'\n');
fprintf(fid,'}');
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
fprintf(fid,'\\begin{table}[!htb]');
fprintf(fid,'\n');
fprintf(fid,'\\captionsetup{width=0.7\\textwidth, font={singlespacing}}');
fprintf(fid,'\n');
fprintf(fid,'\\begin{center}');
fprintf(fid,'\n');
fprintf(fid,'\\begin{tabular} {c c c c c c c}');
fprintf(fid,'\n');
fprintf(fid,'\\multicolumn {7} {l} {{\\bf Tabla \\ref{unodos}.} {Volumen'); 
fprintf(fid,' promedio, energía y presión de las moleculas}}');
%fprintf(fid,'confinamiento a $T=350 K$}}');
fprintf(fid,'\\\\');
fprintf(fid,'\n');
fprintf(fid,'\\multicolumn {7} {l} {bajo confinamiento a $T= 350\\ K$.$^\\ast$}');
fprintf(fid,'\\\\');
fprintf(fid,'\n');
fprintf(fid,'\\hline \\hline');
fprintf(fid,'\n');
fprintf(fid,'Volumen && energía && Presión estática && Presión dinámica\\\\');
fprintf(fid,'\n');
fprintf(fid,'$nm^3$ && K Hartree && MBar && KBar\\\\');
fprintf(fid,'\n');
fprintf(fid,'\\cline{1-1} \\cline{3-3} \\cline{5-5} \\cline{7-7} ');
fprintf(fid,'\n');
fprintf(fid,'1.033742 && -1.035234 && 2.0  && 3.7 \\\\ ');
fprintf(fid,'\n');
fprintf(fid,'0.911486 && -1.026682 && 5.2  && 4.1 \\\\ ');
fprintf(fid,'\n');
fprintf(fid,'0.830515 && -1.012948 && 10.1 && 4.3 \\\\ ');
fprintf(fid,'\n');
fprintf(fid,'0.762473 && -0.990814 && 18.6 && 4.7 \\\\ ');
fprintf(fid,'\n');
fprintf(fid,'0.706373 && -0.960962 && 29.5 && 4.9 \\\\ ');
fprintf(fid,'\n');
fprintf(fid,'\\hline');
fprintf(fid,'\n');
fprintf(fid,'\\end{tabular}');
fprintf(fid,'\n');
fprintf(fid,'\\caption{\\footnotesize{Volumen promedio,');
fprintf(fid,' energía y presión');
fprintf(fid,' de las moleculas bajo confinamiento');
fprintf(fid,' a $T= 350\\ K$.$^\\ast$}}');
fprintf(fid,'\n');
fprintf(fid,'\\label{unodos}');
fprintf(fid,'\n');
fprintf(fid,'\\end{center}');
fprintf(fid,'\n');
fprintf(fid,'\\end{table}');
fprintf(fid,'\n');
fprintf(fid,'\n');


%fprintf(fid, '%f %f %.3f %.3f\n', [VolPromNm3; fpromEnesw; Psw; Pdinsw']);

fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\\begin{table}[!htb]');
fprintf(fid,'\n');
fprintf(fid,'\\captionsetup{width=0.7\\textwidth, font={singlespacing}}');
fprintf(fid,'\n');
fprintf(fid,'\\begin{center}');
fprintf(fid,'\n');
fprintf(fid,'\\begin{tabular} {c c c c c c c}');
fprintf(fid,'\n');
fprintf(fid,'\\multicolumn {7} {l} {{\\bf Tabla \\ref{dosdos}.} {Volumen'); 
fprintf(fid,' promedio, energía y presión de las moleculas}}');
%fprintf(fid,'confinamiento a $T=350 K$}}');
fprintf(fid,'\\\\');
fprintf(fid,'\n');
fprintf(fid,'\\multicolumn {7} {l} {bajo confinamiento a $T= 350\\ K$.$^\\ast$}');
fprintf(fid,'\\\\');
fprintf(fid,'\n');
fprintf(fid,'\\hline \\hline');
fprintf(fid,'\n');
fprintf(fid,'Volumen && energía && Presión estática && Presión dinámica\\\\');
fprintf(fid,'\n');
fprintf(fid,'$nm^3$ && KHartree && MBar && KBar\\\\');
fprintf(fid,'\n');
fprintf(fid,'\\cline{1-1} \\cline{3-3} \\cline{5-5} \\cline{7-7} ');
fprintf(fid,'\n');
%fprintf(fid,'%f && %f && %.3f && %.3f \\\\', [VolPromNm3(7); fpromEnesw(7); Psw(7); Pdinsw(7)']);
%fprintf(fid,'\n');
%fprintf(fid,'%f && %f && %.3f && %.3f \\\\', [VolPromNm3(6); fpromEnesw(6); Psw(6); Pdinsw(6)']);
%fprintf(fid,'\n');
fprintf(fid,'%f && %f && %.3f && %.3f \\\\', [VolPromNm3(5); fpromEnesw(5); Psw(5); Pdinsw(5)']);
fprintf(fid,'\n');
fprintf(fid,'%f && %f && %.3f && %.3f \\\\', [VolPromNm3(4); fpromEnesw(4); Psw(4); Pdinsw(4)']);
fprintf(fid,'\n');
fprintf(fid,'%f && %f && %.3f && %.3f \\\\', [VolPromNm3(3); fpromEnesw(3); Psw(3); Pdinsw(3)']);
fprintf(fid,'\n');
fprintf(fid,'%f && %f && %.3f && %.3f \\\\', [VolPromNm3(2); fpromEnesw(2); Psw(2); Pdinsw(2)']);
fprintf(fid,'\n');
fprintf(fid,'%f && %f && %.3f && %.3f \\\\', [VolPromNm3(1); fpromEnesw(1); Psw(1); Pdinsw(1)']);
fprintf(fid,'\n');
fprintf(fid,'\\hline');
fprintf(fid,'\n');
fprintf(fid,'\\end{tabular}');
fprintf(fid,'\n');
fprintf(fid,'\\caption{\\footnotesize{Volumen promedio,');
fprintf(fid,' energía y presión');
fprintf(fid,' de las moleculas bajo confinamiento');
fprintf(fid,' a $T= 350\\ K$.$^\\ast$}}');
fprintf(fid,'\n');
fprintf(fid,'\\label{dosdos}');
fprintf(fid,'\n');
fprintf(fid,'\\end{center}');
fprintf(fid,'\n');
fprintf(fid,'\\end{table}');
fprintf(fid,'\n');
fprintf(fid,'\n');


fprintf(fid,'\\begin{figure}[h]');
fprintf(fid,'\n');
fprintf(fid,'\\captionsetup{width=1\\textwidth, font={singlespacing}}');
fprintf(fid,'\n');
fprintf(fid,'\\centering');
fprintf(fid,'\n');
fprintf(fid,'\\resizebox{1.0\\textwidth}{!}{\\input{curvaEcEstado001}}');
fprintf(fid,'\n');
fprintf(fid,'\\label{fig:figure2}');
fprintf(fid,'\n');
fprintf(fid,'\\caption{\\scriptsize{caption}}');
fprintf(fid,'\n');
fprintf(fid,'\\end{figure}');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\n');

fprintf(fid,'\n');
fprintf(fid,'\\end{document}');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\n');
%fprintf(fid,'2 2 2 2');
fclose(fid);



%FORCES FORCES FORCES FORCES

NumGrupoBarras = 5;
sepEjeY =0.50;
anchoBarra = 1;
sepBarJuntas = 0.00;
sepBarGrupo = 0.5;

for i=1:NumGrupoBarras
	
	xBarInt(i) = sepEjeY + (anchoBarra/2)*(i) + sepBarJuntas*(i-1) ...
	+ (1.5)*anchoBarra*(i-1) + sepBarGrupo*(i-1);
end

for i=1:NumGrupoBarras
	
	xBarExt(i) = sepEjeY + anchoBarra*(i) + sepBarJuntas*(i) ...
	+ (anchoBarra/2)*(i) + (anchoBarra/2)*(i-1) + sepBarGrupo*(i-1);

end

hold off
figure('Name', 'barras forces')
for i=1:NumGrupoBarras

	%Lineas Izquierdas barras, primero bar izq 2da linea bar der
	plot([(xBarInt(i)-(anchoBarra/2)) (xBarInt(i)-(anchoBarra/2))], [0 FPromInt(i)], 'b')
	hold on
	plot([(xBarExt(i)-(anchoBarra/2)) (xBarExt(i)-(anchoBarra/2))], [0 FPromExt(i)], 'r')
	
	%Lineas Derechas barras, primero bar izq 2da linea bar der
	plot([(xBarInt(i)+(anchoBarra/2)) (xBarInt(i)+(anchoBarra/2))], [0 FPromInt(i)], 'b')
	plot([(xBarExt(i)+(anchoBarra/2)) (xBarExt(i)+(anchoBarra/2))], [0 FPromExt(i)], 'r')

	%tapaderas, primero bar izq luego bar der
	plot([(xBarInt(i)-(anchoBarra/2)) (xBarInt(i)+(anchoBarra/2))], [FPromInt(i) FPromInt(i)], 'b')
	plot([(xBarExt(i)-(anchoBarra/2)) (xBarExt(i)+(anchoBarra/2))], [FPromExt(i) FPromExt(i)], 'r')

end


ff = max(FPromExt);
ff2 = max(FPromInt);

maxff = max([ff ff2]);

xBarInt(2)
(anchoBarra/2)

xAx1 = xBarInt(3) - (anchoBarra/2) - 0.97*sepEjeY;
xAx2 = xBarExt(5) + anchoBarra;
yAx1 = 0;
yAx2 = maxff;

grid on
axis([xAx1 xAx2 yAx1 yAx2])
hold off


filename='dataBarsExt001.dat';
fid = fopen(filename, 'w');
fprintf(fid, '%f  %f\n', [xBarExt; FPromExt]);
fclose(fid);

filename='dataBarsInt001.dat';
fid = fopen(filename, 'w');
fprintf(fid, '%f  %f\n', [xBarInt; FPromInt]);
fclose(fid);


filename='graficaBarras001.gnu';
fid = fopen(filename, 'w');
fprintf(fid,'reset');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'# epslatex');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'#set terminal epslatex size 8.2cm,5cm color colortext 10');
fprintf(fid,'\n');
fprintf(fid,'#set terminal epslatex font "Gill Sans,9" linewidth 4 rounded fontscale 1.0 ');
fprintf(fid,'\n');
fprintf(fid,'set terminal epslatex linewidth 4 rounded fontscale 1.0 ');
fprintf(fid,'\n');
fprintf(fid,'set output "graficaBarras001.tex"');
fprintf(fid,'\n');
fprintf(fid,'set size ratio 0.625');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'set encoding iso_8859_1');
fprintf(fid,'\n');
fprintf(fid,'#set border linewidth 2');
fprintf(fid,'\n');
fprintf(fid,'#set style line 5 linecolor rgb "#dd181f" linetype 1 linewidth 1');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'#"#0060ad" azul');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'set style line 80 lt rgb "#808080"');
fprintf(fid,'\n');
fprintf(fid,'set style line 81 lt 0 # dashed');
fprintf(fid,'\n');
fprintf(fid,'set style line 81 lt rgb "#808080" # grey');
fprintf(fid,'\n');
fprintf(fid,'set grid back linestyle 81');
fprintf(fid,'\n');
fprintf(fid,'set border 3 back linestyle 80 # Remove border on top and right. These ');
fprintf(fid,'\n');
fprintf(fid,'#');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'set xtics nomirror ');
fprintf(fid,'\n');
fprintf(fid,'set ytics nomirror ');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'set style line 1 lt rgb "#A00000" lw 2 pt 1');
fprintf(fid,'\n');
fprintf(fid,'set style line 2 lt rgb "#00A000" lw 2 pt 6');
fprintf(fid,'\n');
fprintf(fid,'set style line 3 lt rgb "#5060D0" lw 2 pt 2');
fprintf(fid,'\n');
fprintf(fid,'set style line 4 lt rgb "#F25900" lw 2 pt 9');
fprintf(fid,'\n');
fprintf(fid,'set style line 5 linecolor rgb "#dd181f" linetype 1 linewidth 1');
fprintf(fid,'\n');
fprintf(fid,'set style line 6 linecolor rgb "#5060D0" linetype 1 linewidth 0.3');
fprintf(fid,'\n');
fprintf(fid,'set style line 7 linecolor rgb "#F25900" linetype 1 linewidth 0.3 ');
fprintf(fid,'\n');
fprintf(fid,'set style line 8 linecolor rgb "#00A000"  linetype 1 linewidth 1');
fprintf(fid,'\n');
fprintf(fid,'set style line 9 linecolor rgb "#A00000"  linetype 1 linewidth 1');
fprintf(fid,'\n');
fprintf(fid,'set style line 10 linecolor rgb "#20B2AA"  linetype 1 linewidth 0.3');
fprintf(fid,'\n');
fprintf(fid,'set style line 11 linecolor rgb "#4B0082"  linetype 1 linewidth 1');
fprintf(fid,'\n');
fprintf(fid,'set style line 12 linecolor rgb "#F25900" pt 63 ');
fprintf(fid,'\n');
fprintf(fid,'set style line 14 linecolor rgb "#008EA7"  linetype 1 linewidth 0.3');
fprintf(fid,'\n')
fprintf(fid,'set style line 15 linecolor rgb "#204253"  linetype 1 linewidth 0.3');
fprintf(fid,'\n');

fprintf(fid,'#set ylabel ');
fprintf(fid,'\n');
fprintf(fid,'#set xlabel');
fprintf(fid,'\n');
fprintf(fid,'set xrange [%f:%f]',xAx1,xAx2);
fprintf(fid,'\n');
fprintf(fid,'set yrange [%f:%f]',yAx1,yAx2);
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'set format "$%%g$"');
fprintf(fid,'\n');
fprintf(fid,'#set key top right');
fprintf(fid,'\n');
fprintf(fid,'set key off');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\n');

dbarras = (xBarExt(1)-xBarInt(1))/2;
xA(1) = xBarInt(1) + dbarras;
xA(2) = xBarInt(2) + dbarras;
xA(3) = xBarInt(3) + dbarras;
xA(4) = xBarInt(4) + dbarras;
xA(5) = xBarInt(5) + dbarras;


%fprintf(fid,'set ytics ("%.3f" %.3f, "%.3f" %.3f, "%.3f" %.3f, "%.3f" %.3f, "%.3f" %.3f, "%.3f" %.3f) nomirror', FPromInt(3),FPromInt(3)...
%,FPromExt(3),FPromExt(3),FPromInt(4),FPromInt(4),FPromExt(4),FPromExt(4),FPromInt(5),FPromInt(5),FPromExt(5),FPromExt(5));



fprintf(fid,'\n');
fprintf(fid,'set xtics ("ro=7.0" %.3f, "ro=7.25" %.3f, "ro=7.5" %.3f, "ro=7.75" %.3f, "ro=8.0" %.3f) nomirror', xA(1), xA(2), xA(3), xA(4), xA(5));
fprintf(fid,'\n');
fprintf(fid,'set boxwidth %f', anchoBarra);
fprintf(fid,'\n');
fprintf(fid,'set style fill solid ');
fprintf(fid,'\n');
fprintf(fid,'#set style fill transparent solid 0.2 noborder');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'plot "dataBarsExt001.dat" using 1:2 with boxes ls 14,\\');
fprintf(fid,'\n');
fprintf(fid,'     "dataBarsInt001.dat" using 1:2 with boxes ls 15');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'#plot "dataExampleBars.dat" using 1:2 ls 10 with lines title "Particulas Confinadas" #, \\');
fprintf(fid,'\n');
fprintf(fid,'#"tempCage.dat" using 1:2 ls 7 with lines title "Particulas Confinadas"');
fprintf(fid,'\n');
%dataBarsExt001.dat
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\n');


%fprintf(fid,'%.2f,t lw .3 lc rgb "%s" , \\',t(i), hexcd1(i,:));
%ps2pdf -dEPSCrop graficaBarrasPrueba001.eps graficaBarrasPrueba001.pdf
fclose(fid);

%[status,cmdout] = system('ps2pdf -dEPSCrop graficaBarras001.eps graficaBarras001.pdf');

