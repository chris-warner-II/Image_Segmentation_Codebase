% This is the simplest ass test of the Graph Laplacian (Negative GL &
% Normalized GL too).  And it seems to fail in my eyes.  Do you agree? Does
% this surprise you?

%% Make a very simple input image.
% img = rand(2);

% img = [0.2370    0.0198; 0.0076    0.3131]
% img = [0.1 0.8; 0.1 0.9]
% img = [0.1 0.8 0.1; 0.7 0.9 0.8; 0.2 0.75 0.1] % 5th largest Evec looks best.
% img = [0.8 0.2 0.1; 0.7 0.9 0.4; 0.8 0.75 0.7] % Xth largest Evec looks best.
img = [1 0 0; 1 1 0; 1 1 1] % Binary.

ximg = size(img,1);
yimg = size(img,2);


%% Define Weights Matrix
sigma = 0.0001; %0.1
sigdist = 1;
threshdist = 0.55;
W = calc_weights(img,sigma,sigdist,threshdist);

%% Calculate Graph Laplacian Matrix
% normlze = 0;
% neg = 1;
% Q = compute_Laplacian(W,normlze,neg,img);

W = full(W)                       % weights
D = diag(sum(W))                  % vertex incidences
L  = D - W                        % Laplacian
LS = L + 10.*eye(size(L));        % Shift Evals

% % Normalized Graph Laplacian
% Di = eye(size(W));          % ones on the diagonal
% degW = sum(W,1);            % degree of incidence into each vertex
% ODi = -W./sqrt(degW'*degW);  % off-diagonal elements for NGL
% for i = 1:numel(degW)
%    ODi(i,i) = 0; 
% end
% Ln = Di + ODi
% LnS = Ln - 5.*eye(size(Ln));  % Shift Evals

%% Eigensystem Calculations
disp('Graph Laplacian')
[Evec,Eval]=eig(L);
diag(Eval)'
Evec
Eval1 = Eval(1,1)
Evec1 = reshape(Evec(:,1),ximg,yimg)
Eval2 = Eval(2,2)
Evec2 = reshape(Evec(:,2),ximg,yimg)

% disp('Shifted Graph Laplacian')
% [EvecS,EvalS]=eig(LS);
% diag(EvalS)'
% EvecS
% EvalS1 = EvalS(1,1)
% EvecS1 = reshape(EvecS(:,1),ximg,yimg)
% EvalS2 = EvalS(2,2)
% EvecS2 = reshape(EvecS(:,2),ximg,yimg)

% disp('Normalized Graph Laplacian')
% [EvecN,EvalN]=eig(Ln);
% diag(EvalN)'
% EvecN
% Eval1N = EvalN(1,1)
% Evec1N = reshape(EvecN(:,1),ximg,yimg)
% Eval2N = EvalN(2,2)
% Evec2N = reshape(EvecN(:,2),ximg,yimg)
% 
% disp('Shifted Normalized Graph Laplacian')
% [EvecSN,EvalSN]=eig(LnS);
% diag(EvalSN)'
% EvecSN
% Eval1SN = EvalSN(1,1)
% Evec1SN = reshape(EvecSN(:,1),ximg,yimg)
% Eval2SN = EvalSN(2,2)
% Evec2SN = reshape(EvecSN(:,2),ximg,yimg)

figure, imagesc(img)
