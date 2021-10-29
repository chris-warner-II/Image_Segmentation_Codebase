function [paramEsts,paramCIs,acov,stderr] = truncNormStats_MLEest(x,xTruncL,xTruncR)




% Script to play around with finding the parameters (mu, sig) for a normal
% distribution on a finite interval using MLE.

if ~exist('xTruncL','var')
    xTruncL=0;
end

if ~exist('xTruncR','var')
    xTruncR=1;
end


mu = mean(x); % + (rand-0.5).*mean(x);  % these will be wrong, but maybe ok starting point
sigma = std(x); % + (rand-0.5).*std(x);


pdf_truncnorm = @(x,mu,sigma) normpdf(x,mu,sigma) ./ (normcdf(xTruncR,mu,sigma)-normcdf(xTruncL,mu,sigma));


start = [mu,sigma]; % iuno. Starting guess. At empirical mean & Std.



options = statset('MaxIter',1000, 'MaxFunEvals',1000, 'FunValCheck','off'); %

[paramEsts,paramCIs] = mle(x, 'pdf',pdf_truncnorm, 'start',start, 'lowerbound',[-Inf 0], 'options', options); % 'lowerbound' was just 'lower'


% OK, I MIGHT HAVE TO CHANGE THIS TO A STUDENT-T DISTRIBUTION FOR ALL.  IT
% WILL BE MORE CORRECT FOR CASES WITH FEW SAMPLES.  IT WILL RELAX TO THE
% GAUSSIAN DISTRIBUTION FOR LARGE NUMBER OF SAMPLES.  
%
% THE CONFIDENCE INTERVALS OF STD SHOULD BE NON-NEGATIVE. AND SHOULD ALSO
% NOT BE SYMMETRIC.



acov = mlecov(paramEsts, x, 'pdf',pdf_truncnorm);
stderr = sqrt(diag(acov));


