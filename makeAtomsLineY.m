function [x y] = makeAtomsLineY(P1, P2)

ymin = min([P1(2) P2(2)]);
ymax = max([P1(2) P2(2)]);

y = ymin:.01:ymax;

invm = (P2(1) - P1(1))/(P2(2) - P1(2));

x = ((y - P1(2))*invm)+P1(1); 


