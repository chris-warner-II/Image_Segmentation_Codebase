% This script will plot r on the x-axis from 0 to 1.  It will plot the 3
% different potential measures of circular standard deviation on the y-axis
% to see how they compare to one another.

N=100; % hard coded for now.  We are just inputting r value.

r = [0:0.01:1];

s = sqrt(2*(1-r));      % Circular STD from Variance : bounds = [0,1]
s0 = sqrt(-2*log(r));   % Circular STD : bounds = [0, inf]


% finally STD defined as 1/sqrt(kappa) : bounds = [0,inf]
for i = 1:numel(r)
    
    i
    
    R = r(i);
    
    if R < 0.53
      kappa = 2*R + R^3 + 5*R^5/6;
    elseif R>=0.53 && R<0.85
      kappa = -.4 + 1.39*R + 0.43/(1-R);
    else
      kappa = 1/(R^3 - 4*R^2 + 3*R);
    end

    if N<15 && N>1
      if kappa < 2
        kappa = max(kappa-2*(N*kappa)^-1,0);    
      else
        kappa = (N-1)^3*kappa/(N^3+N);
      end
    end

    s1(i) = 1./sqrt(kappa);
    
end



figure, hold on,
plot(r,s,'b','LineWidth',2)
plot(r,s0,'g','LineWidth',2)
plot(r,s1,'r','LineWidth',2)

title('Three Ways to Compute Circular Standard Deviation','FontSize',20,'FontWeight','Bold')

xlabel('resultant vector length (r)','FontSize',18,'FontWeight','Bold')
ylabel('''STD''','FontSize',18,'FontWeight','Bold')

legend({'$\sqrt{2(1-r)}$','$\sqrt{-2ln(r)}$','$1/\sqrt{\kappa}$'},'Interpreter','Latex')

set(gca,'FontSize',16,'FontWeight','Bold')

keyboard