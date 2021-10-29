function [ DEV, ICE ] = power_method( Q, ximg, yimg, EvecTrue, initdist, method, project, thresh, breakout, dispsho, DynMov, plotErr, HistMov)

% syntax: DEV = power_method( Q, ximg, yimg, EvecTrue, initdist, method, project, dispsho, DynMov, plotErr, HistMov)
%
% Output:
% DEV = Dominant Eigenvector (one with largest positive eigenvalue) (2nd largest if projecting)
% ICE = [Iteration,Change,Error]. Concatination of Vectors of [Iteration #, Iterative Change, Error from EvecTrue]
%
% Input:
% [ximg,yimg] = size of [x,y] dimensions (in pixels) of image.  To reshape Eigenvector.
% EvecTrue = passing in true eigenvector (from 'eig' function call) for error calculation
% initdist = initial distribution or vector to input to PM algorithm.
%            (random, delta function, uniform) - maybe add Gaussian?
% project = {0,1} a flag to project out all 1's DEV for Negative Graph
%           Laplacian to find next (i.e., 2nd) Largest Eigenvector.
% dispsho = {0,1} a flag to plot dynamics of PM algorithm (eigenvector vs. iteration #)
% saveMov = {0,1} a flag to save plots of PM dynamics to a movie.
% plotErr = {0,1} a flag to plot iteration vs. error between PM DEV and true Evec vs Evec change (convergence)
% HistMov = {0,1} a flag to plot a Pixel-wise Histogram of diff between PM & True Evecs.  Should be peaked around zero.
%
% TODO: Plot a curve (hopefully monotonically decreasing) of iteration vs. error between PM DEV and true Eigenvector

%% Initializations and flags on how to run & display
if strcmp(initdist,'random')
    x = rand(size(Q,1),1); % initialize eigenvector to uniformly chosen random distrib.
end

if strcmp(initdist,'delta')
    x = zeros(size(Q,1),1); % initialize eigenvector with randomly placed delta function
    x( round(1200*rand(1)) ) = 1;
end

if strcmp(initdist,'uniform')
    x = ones(size(Q,1),1); % initialize eigenvector with uniform distribution across nodes.
    % is this a useful distribution ?
end

if strcmp(initdist,'gaussian')
    % A Gaussian distribution for initialization. Is this a useful distribution ?
end

if strcmp(initdist,'zeros')
    x = zeros(size(Q,1),1); % initialize eigenvector with uniform distribution across nodes.
    % is this a useful distribution ?
end

%
if ~exist('thresh','var')
    thresh = 1e-7; % 
end
if ~exist('breakout','var')
    breakout = 200000;
end
%
change = 1;
i=1; j=0; % loop/iteration variables below.
%
shownum = 50; % show a plot every X iterations.

% (B). Open a file to save a movie of Power Method Dynamics
if(DynMov && dispsho)
    vDyn = VideoWriter('PM_Dynamics','Uncompressed AVI');
    vDyn.FrameRate = 10;
    open(vDyn);
end

% (D). Open a file to save a movie of Histogram of Pixel-wise Error Distributions
if(HistMov && dispsho)
    vHist = VideoWriter('Err_Hist','Uncompressed AVI');
    vHist.FrameRate = 10;
    open(vHist);
end

% (C). Initialize a 3 x N vector where N = breakout/shownum
    if(plotErr) 
        N = breakout/shownum;
        ICE = zeros(3,N);
    end

% Project out 1's vector
if project
    v = ones(size(x));                   % the 1's vector
    P = eye(numel(x)) - (v*v')/(v'*v);   % projection matrix orthogonal to v   
    Q = Q*P;
    %Px = P*x;        % check: Px vector projected orthogonal to 1's vector
end


%% Use iterative Power Method to solve for Dominant Eigenvector
while ( (change > thresh) && (i < breakout) )
    
    % The Guts of the Power Method
    y = Q*x;
    y = y/max(abs(y));
    i = i+1; 
    change = mean(abs(x-y)); % for breakout and (C) plot below
    x=y;

    % (A). To visualize dynamics of power method
    if( dispsho && ( i==1 || (mod(i,shownum)==0) ) )
        i
        j = j+1;
        DEV = reshape(y, ximg, yimg);
        if(  mean(mean(abs(EvecTrue - DEV))) < mean(mean(abs(EvecTrue + DEV)))  )
            diffEvecs = DEV - EvecTrue;
            ErrTrue = mean(mean(abs(diffEvecs)));
        else
            diffEvecs = DEV + EvecTrue;
            ErrTrue = mean(mean(abs(diffEvecs)));
        end
        
        % Plot Dominant Eigenvector for current iteration of Power Method
        hDyn = figure; imagesc(DEV); % NOTE: We are no longer taking log of eigenvector.
        axis off 
        colorbar('FontSize',16,'FontWeight','Bold'); %,[-1,1]  
        title(['Power Method: ',method,': Iter #',num2str(i),' - \Delta = ',num2str(change,'%0.1e'),' - Err = ',num2str(ErrTrue,'%0.1g')],'FontSize',20,'FontWeight','Bold')

        % Plot Histogram of distribution of differences
        [num,bin] = hist(abs(diffEvecs(:)),1000);
        hHist = figure; 
        plot( bin, num )
        axis([0 1 0 max(num)])
        title('Pixel-wise Error between True & PM Evec','FontSize',20,'FontWeight','Bold')
        xlabel('\Delta Diff in Pixel Values','FontSize',16,'FontWeight','Bold')
        ylabel(['# Pixels (out of ',num2str(ximg*yimg),') with that diff value'],'FontSize',16,'FontWeight','Bold')
        %thr = mean(q) + std(q);
        %disp(['Mean + 1-STD of histogram of difference between Power Method & Analytical Eigenvector: ', num2str(thr)])
        
        % (B). Save each frame of the PM Dynamics movie file
        if(DynMov) % note: must have dispshow == 1 also to get here.
            currFrame = getframe(hDyn);
            writeVideo(vDyn,currFrame);
        end
        
        % (D). Save each frame of the Pixelwise Error Histogram movie file
        if(HistMov) % note: must have dispshow == 1 also to get here.
            currFrame = getframe(hHist);
            writeVideo(vHist,currFrame);
        end
        
        % (C). Record iteration # vs. change vs. error between PM Evec & True one.
        if(plotErr)
        	ICE(1,j) = i;
            ICE(2,j) = change;
            ICE(3,j) = ErrTrue;
        end
        %
        if(mod(i,10*shownum)==0) % clean up figures periodically
            close all
        end
        
    end

end

% (B). Close the Evec Dynamics movie file
if(DynMov && dispsho)
    close(vDyn);
end

% (D). Close the Histogram movie file
if(HistMov && dispsho)
    close(vHist);
end

% (C). Plot Iteration # vs. change vs. error between PM Evec & True one. 
if(plotErr && dispsho)
    %legend('(\Delta) Iterative Change of Evec','(\epsilon) Error from True Evec','FontSize',16,'FontWeight','Bold')
    figure, subplot(211), semilogy(ICE(1,:),ICE(2,:));
    title(['Power Method: ',method,': Convergence Dynamics'],'FontSize',20,'FontWeight','Bold')
    ylabel('Log of Iterative Change','FontSize',16,'FontWeight','Bold')
    %
    subplot(212), plot(ICE(1,:),ICE(3,:));
    xlabel('Iteration #','FontSize',16,'FontWeight','Bold')
    ylabel('Error','FontSize',16,'FontWeight','Bold')

end

disp(['Power Method Brokeout at Iteration #',num2str(i),' - Change = ',num2str(change)])

DEV = reshape(y, yimg, ximg);

if(~exist('ICE','var'))
    ICE = 0;
end

% set zeros in img_seg to small number (not zero) because log(zero) = -Inf
% DEV(DEV==0) = min(min(DEV(DEV~=0)));
% DO I REALLY NEED THIS?  (I think so because of Log of Dominant
% Eigenvalue)  Do I really need to do that?
