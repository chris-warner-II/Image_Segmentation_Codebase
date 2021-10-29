function [Fx, Fy] = gradientB(X,circ)

% This function will return the magnitude or horizontal (Fx) and vertical
% (Fy) gradients in the input X.  If X is a linear variable - that is
% circ=0, then this runs similar to matlab's gradient function.  It is
% different because it returns magnitude of gradient as if you are
% exporting abs(Fx) and abs(Fy).
%
% If circ=1, the function make sure that any horizontal or vertical gradients
% are less than pi before averaging pairs of side-by-side points together
% to get out Fx and Fy.
%
% Written by Chris Warner 9/15 to mimic gradient with some differences.


% This created a simple test image that I used to construct this function.
if(0)
    X(1,:) = [ 0 0 0 0 0 0 0 0 0 ]; 
    X(2,:) = [ 0 0 0 0 0 0 0 0 0 ];
    X(3,:) = [ 0 0 0 2 2 2 0 0 0 ];
    X(4,:) = [ 0 0 0 1 2 2 0 0 0 ];
    X(5,:) = [ 0 0 0 2 2 2 0 0 0 ];
    X(6,:) = [ 0 0 0 0 0 0 0 0 0 ];
    X(7,:) = [ 0 0 0 0 0 0 0 0 0 ];

    figure, imagesc(X)
end


y = abs(diff(X));
x = abs(diff(X')');


% have to break in here and make sure no jumps are larger than pi. Because
% if we found a jump > pi, we should have gone around the other way.
if(circ)
    
    yn(:,:,1) = y;             % size of vertical gradient given the way we did travel around circle.
    yn(:,:,2) = abs(2*pi - y); % size of vertical gradient if we travel otherway around circle.
    yd = min(yn,[],3);         % vertical gradient given shorter way around circle
    
    xn(:,:,1) = x;             % size of horizontal gradient given the way we did travel around circle.
    xn(:,:,2) = abs(2*pi - x); % size of horizontal gradient if we travel otherway around circle.
    xd = min(xn,[],3);         % horizontal gradient given shorter way around circle
    
else
    
    yd = y;
    xd = x;                    % do nothing if not circular.
    
end


% Now we do same as gradient function. Pad zeros on either side to shift
% gradients from between pixels to on pixels and then average together the
% diff on either side of a pixel to be its value at that pixel.
xdb(:,:,1) = [ zeros(size(X,1),1), xd];
xdb(:,:,2) = [xd, zeros(size(X,1),1) ];
Fx = mean(xdb,3);


ydb(:,:,1) = [ zeros(1,size(X,2)); yd];
ydb(:,:,2) = [yd; zeros(1,size(X,2)) ];
Fy = mean(ydb,3);



if(0)
    disp('compare [Fx,Fy] to [x2,y2] output from gradient function')
    [x2,y2] = gradient(X);
end
