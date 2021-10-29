%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% In this simple example, we assume there is a mixture of 2 dimensional
% gaussian variables, where mu and covariance are unknown. But we know the
% covariances are diagonal and isotropic. Therefore what we don't know 
% are mu and the scalar factor sigma. We use dirichlet process to model
% the problem and do clustering. We assume it has a conjugate prior for
% mu and sigma. Since the likelihood is gaussian, the conjugate prior 
% should be normal-gamma distribution. More specifically, sigma has gamma
% distribution and mu has multivariate student-t distribution. The purpose
% of using dirichlet process is that we do not want to specify the number
% of components in the mixture, but instead give a prior over 1 to
% infinite. Then Gibbs sampler is used to draw sample from posterior
% distribution based on the observation. See [2] for detail about the
% algorithm.
% 
% most implementation is encapsulated in DirichMix class (DirichMix.m)
%
% distributable under GPL
% written by Zhiyuan Weng, Nov 26 2011
%
%
% Reference:
% [1] D Fink, "A Compendium of Conjugate Priors"
% [2] R Neal, "Markov Chain Sampling Methods for Dirichlet Process Mixture Models"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear
% clc







figure, plot(phaseAtClk')

x = cos( phaseAtClk(:,55) );
xp = [x,x];


dirich = DirichMix; % construct an object of the class
dirich.InputData(xp(:,1:2));
dirich.DoIteration(10); % 100 iterations



% subplot(1,2,2)

figure,
dirich.PlotData
title('clustering results');
axis([0, 2*pi, 0, 2*pi]);
axis square;

