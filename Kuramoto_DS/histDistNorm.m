function [H] = histDistNorm(DistsPW,titleStr)


[xi,yi] = hist(DistsPW(1,:),100);
[xo,yo] = hist(DistsPW(2,:),100);
H=figure; hold on
bar(yo,xo,'b')
bar(yi,xi,'r')
axis([0 1 0 max([xi,xo])])
legend('across-Cluster Distance','in-Cluster Distance')
title(titleStr)

titleStr
[min([yi,yo]), max([yi,yo])]