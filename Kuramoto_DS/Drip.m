function [freq, rStats, fezStats] = Drip(w, s, c, T, reach, dynMov)

% syntax: [x] = Drip(w, s, T, x0);
%
% This function models the discrete-time, discrete-space "Dripping
% Handrail" Dynamical System Model.  It will simulate an oscillatory 
% handrail system with a stable periodic limit cycle.
%
% Inputs: w = The natural frequency of drop in isolation
%         s = rate of each droplet fall (in isolation) - resonant frequency of uncoupled osc.
%         T = number of timesteps to consider.  Length of x vector.
%
% Outputs: x = a vector of time evolution of the single drop system.  It is
%              oscillatory, with water collecting until it reaches 1 and 
%              then dropping off to be 0 and start all over again.
%
%           Author: Chris Warner 5/2013 



%% Foreplay.
N = numel(w);  % Number of coupled oscillators

x = zeros(T,N); % preallocate matrix to hold phases of oscillators at each timestep
theta = zeros(T,N); % preallocate matrix to hold phases of oscillators at each timestep
fwd = zeros(N); % preallocate matrix
shft = zeros(2*N-1,N); % preallocate matrix

% vector of initial conditions (can be uniform, completely random, have some Gaussian distribution, ...)
x(1,:) = rand(1,N);    % initial conditions of oscillator phases
% x(1,:) = zeros(1,N); % to start all oscillators at zero phase
% x(1,:) = rand.*ones(1,N); % to start all oscillators at same random phase of oscillation.
theta(1,:) = 2*pi.*x(1,:);

if(N==1)
   reach = 0;
   rch = 1;
else
   rch = reach; % rch is divided by below.
end

% Preallocate arrays for data and movie structure.
if(dynMov)
    writerObj = VideoWriter(['DripMov.avi']); % _w',p,'_sp_c',p,'_r',p,'
    writerObj.FrameRate = 2;
    open(writerObj);
end


%% Evolve state of drop by linear equation with nonlinear (mod1) constraint.
for t = 1:(T-1)
    
    % create shifted versions of the x vector in a matrix to calculate coupling easily
    % Note: this construction assumes wrap-around boundaries of 1D osc. chain.
    if(reach)
        fwd = circulant(x(t,:),1);
        shft = [fwd(2:end,:);fwd(2:end,:)];
        nebrs = shft(N-reach:(N-1)+reach,:)./(2*reach);
    else
        nebrs = 0;
    end
    
    % This is the Dripping Handrail Circle Map Lattice Equation !!
    x(t+1,:) = mod( x(t,:) + (w') + (s').*x(t,:) + (c').*sum(nebrs,1) , 1 );
    % convert x to a phase (0,1) -> (0,2pi) so we can judge coherence radius and make polar movie.
    theta(t+1,:) = 2*pi.*x(t+1,:);
    
    % Create movie of polar phase plot of oscillators
    if(dynMov)
        pp = figure; hold on,
        title(['Oscillator Phase (t = ',num2str(t),')'],'FontSize',20,'FontWeight','Bold');
        p = polar(theta(t,:),ones(1,N),'r.'); % rho 2nd 1 was N
        p = polar(theta(t,2),ones(1,1),'b.');
%         p = polar(theta(t,3),ones(1,1),'g.');
        axis([-1 1 -1 1]);
        delete(findall(ancestor(p,'figure'),'HandleVisibility','off','type','line','-or','type','text')); 
        set(gca,'nextplot','replacechildren');
        mov = getframe(pp);
        writeVideo(writerObj,mov);
        close % to keep figures from overloading memory

    end

end

% Calculate average phases and coherence radius of oscillators
phi = mean(theta');
phiStd = std(theta');
rx = mean(cos(theta')); % THESE MEANS DO MEAN FIELD SHIT.  DO I REALLY WANT THIS?
ry = mean(sin(theta'));
r = sqrt(rx.^2 + ry.^2);
rStats = [mean(r), std(r)]; % average Coherence Radius

%% Calculate Period of each oscillator using Fourier Transform
N1 = 1.2*T;
X1 = abs(fft(x,N1));
F1=[0: N1-1]/N1;

for i=1:N
    pk = find( X1(:,i) == max(X1( round(0.01*N1):end - round(0.01*N1), i ) ) );
    freq(i) = F1(pk(1));
end

% figure, plot(F1,X1,'-x'),title(num2str(freq))
% pause

phse = mod(diff(x,1,2),1); % This diff move only works for 2 osc.  Think about more.  Pairwise phase-diffs?
fezStats = [mean(phse),std(phse)];
% keyboard


%% Plot radius of convergence and average phase with STD errorbars
if(0)
    h=figure;
    subplot(311), plot(theta); %errorbar(phi,phiStd)
    title('Average phase of all Oscillators','FontSize',20,'FontWeight','Bold')
    subplot(312), plot(r,'r')
    title('Average Coherence Radius of Oscillators','FontSize',20,'FontWeight','Bold')
    subplot(313), 
    %
    bar([0.1:0.1:1],y);
    title('Histogram of Convergence Radius','FontSize',20,'FontWeight','Bold')
    xlabel(['[w = ',num2str(w'),'] [s = ',num2str(s'),'] [c = ',num2str(c'),']'],'FontSize',18,'FontWeight','Bold')
    set(h,'Position',get(0,'ScreenSize')); % make right size for my laptop
    
    period_N
    keyboard
end


%% Plot the 1D map of single drip function f.
if(0) % NOT WORKING, BUT MEH WHO CARES...
    a = linspace(0,1,100);
    f = mod( w + s.*a , 1 );
    figure, hold on
    plot(a,f,'g','LineWidth',2), 
    line([0 1], [0 1])
    title ('The Map Function','FontSize',20,'FontWeight','Bold')      
end


%% Is there something interesting in here?  It looks perty.
if(0)
    polar(theta,ones(size(theta)))
end


%% Movies of Polar Plot of Oscillator Dynamics
% close Circle Dynamics movie file
if(dynMov)
    close(writerObj);
end