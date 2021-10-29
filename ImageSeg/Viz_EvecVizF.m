% Script to plot the Eigenvector Visualization function EvecVizF for
% different values of slope.  I also want to look at shifting the location
% of the steepest part of the curve away from zero (to the mean or
% something).

x = linspace(-1,1,1000);
slope = [1 1e-1 1e-2 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8 1e-9 1e-10 1e-11 1e-12 1e-13 1e-14 1e-15];
slope_leg{1} = 'null';

y = zeros(numel(slope),numel(x));
for i=1:numel(slope)
    y(i,:) = real(EvecVizF(x,slope(i)));
    slope_leg{i+1} = num2str(slope(i)); 
end

figure, hold on

plot(x,x,'k--')
plot(x,y);

legend(slope_leg,'Location','SouthEast')

xlabel('Original Eigenvector Value')
ylabel('Transformed Eigenvector Value')


% THINK ABOUT HOW TO SHIFT STEEP PART OF THIS FUNCTION FROM 0 (MAYBE TO
% MEAN PIXEL VALUE OR SOMETHING)
