function [THopt,Iopt] = optimize_threshold2(img,truth,N,iter)

% Syntax: THopt = optimize_threshold2(img,truth,N,iter)
%
% Find the Max Information and the Optimal Threshold that leads to it.
% This method is not just a line search, but an iterative series of line
% searches that zooms in on max value each time and redoes line search on
% progressively smaller intervals.
%
% Find optimum from min to max using N points in between
% Shrink length of interval by 1/2 centered around new optimum. (Maybe increase N by 10x)
% Iterate a number of times defined by iter ..

if ~exist('N','var')
    N=10;
end
if ~exist('N','var')
    iter = 10; 
end

j=0; 
converge = 1e-2;
THopt_old = 1e5; THopt=1;
info = zeros(1,N);

A = min(img(:));
B = max(img(:));

while(abs(THopt-THopt_old)>converge && j<iter)
    j=j+1;

    Thr = linspace(A,B,N);
    L = abs(max(Thr)-min(Thr));

    for i = 1:N
        seg = img>Thr(i);
        info(i) = calc_infoB(seg,truth);
    end

    ind = find(info==max(info));
    %
    THopt_old = THopt;
    THopt = Thr(ind);
    THopt = THopt(1);

%     figure, plot(Thr,info),xlabel('Threshold'), ylabel('Information')
    
    A=max(THopt-L/4, A);
    B=min(THopt+L/4, B);
    N=2*N;

end

Iopt = info(ind);
Iopt = Iopt(1);

% keyboard