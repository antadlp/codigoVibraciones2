function animacionTerminosFourier



%save('terB12N12a.mat', 'mfbeB12N12a')
%mfbeB12N12a(i,:) = manualFBT001(ffteB12N12a, I, 1);

load('terB12N12a.mat');

t = length(mfbeB12N12a(:,1));
x = 1:length(mfbeB12N12a(1,:));
Yp = mfbeB12N12a;
Y = Yp(1,:);
x
Y
h = plot(x,Y, 'YDataSource', 'Y')
%axis([-inf inf -inf inf -1 1])
hold on

for m=2:length(t)
    m
    Y = Yp(m,:);
    refreshdata(h,'caller')
    drawnow; 
    pause(.035)
end


