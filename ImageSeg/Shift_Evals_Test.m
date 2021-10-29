X = rand(9);               % Some Matrix
Xs = X - 5.*eye(size(X));  % Shift Evals

[Evec,Eval]=eig(X);        % Eigensystems
[EvecS,EvalS]=eig(Xs);

diag(Eval)'
Evec

diag(EvalS)'
EvecS