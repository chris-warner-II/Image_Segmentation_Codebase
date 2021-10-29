function pp= vizPhaseLock(theta,tau,location,frameRate,dec,t_start,t_fin)

% Syntax: pp = vizPhaseLock(theta,tau,location,dec,t_start,t_fin)
%
% This function takes in THETA (an array of the phase of all oscillators
% for all time points) and plots a movie of the temporal evolution of
% phase. Oscillators will be placed in location given by topography
% specified in LOCATION.  Time may be decimated or downsampled by DEC.
% 
% Note: do I want to decimate or downsample?  I can run into problems with
% aliasing probably.
%
% THIS DUDE NEEDS TO BE UPDATED IF IM GONNA USE IT AGAIN.

T = size(theta,1);
%
if ~exist('t_start','var')
    t_start = 1;
end
%
if ~exist('t_fin','var')
    t_fin = T;
end
%
if ~exist('dec','var')
    dec = 1;
end
%
if ~exist('frameRate','var')
    frameRate = 10;
end

% Circular Colormap that wraps around for phase so 0 = 2*pi
circ_hsv = colormap(hsv);
circ_hsv = [circ_hsv; [1,0,0] ];


%% Preallocate arrays for data and movie structure.
writerObj = VideoWriter('Kuramoto_Phase_Viz.avi');
writerObj.FrameRate = frameRate;
open(writerObj);



%% 
for t=t_start:dec:t_fin % Loop through timesteps

    % Create movie of polar phase plot of oscillators
    pp = figure;
    
    x = reshape(theta(t,:),size(location)); % This is not completely general.
    imagesc(x), caxis([0, 2*pi]), colorbar
    colormap(circ_hsv); % because hsv is kinda circular with pink and red almost same.
    title(['Oscillator Phase (t = ',num2str(t*tau),' A.U.)'],'FontSize',20,'FontWeight','Bold');
    set(gca,'XTick',[1:size(location,2)],'YTick',[1:size(location,1)],'FontSize',16,'FontWeight','Bold')
    xlabel(['Cluster #'],'FontSize',16,'FontWeight','Bold')
    ylabel(['Osc #'],'FontSize',16,'FontWeight','Bold')

    % Maybe this would be compelling with a 3d stem or scatter plot !!
    % YEA, use SCATTER3 with IMAGE UNDERLAYED!!


    mov = getframe(pp);
    writeVideo(writerObj,mov);
%     if(mod(t,20)==1)
%        keyboard
       close(pp) % to keep figures from overloading memory of compu.
%     end

end



%% Create AVI movie of phase oscillators around circle
close(writerObj);

