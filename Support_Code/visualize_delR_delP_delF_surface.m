% This is a script to loop thru different pairs of (R,P) values and
% calculate F-measure for each and delF for each pair.  


numgrid = 20;

r = linspace(0,1,numgrid);
p = linspace(0,1,numgrid);

% method 1 (Recall, Precision, F-measure)
R1 = repmat( r, numgrid, 1);
P1 = repmat( p, numgrid, 1)';
F1 = 2*P1.*R1./(P1+R1+((P1+R1)==0));

% method 2 (Recall, Precision, F-measure)
R2 = repmat( r, numgrid, 1);
P2 = repmat( p, numgrid, 1)';
F2 = 2*P2.*R2./(P2+R2+((P2+R2)==0));

% what is delta P (P1-P2) and delta R (R1-R2)?
delP = P1-P2';
delR = R1-R2';
delF = F1-F2';



