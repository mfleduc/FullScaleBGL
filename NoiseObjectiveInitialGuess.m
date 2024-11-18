function val = NoiseObjectiveInitialGuess( params,covmodel,Phi,Y, varargin)
% Objective function for determining noise parameters in the basis
% graphical lasso noise model.
%   more detail to come
if ~isempty(varargin)
    precisionflag = varargin{1};
else
    precisionflag=0;
end
val=0;
n = size(Y,2) ;
% l = size(Phi,1);
% % del = 1+tol ;
Sigma = covmodel(params(2:end));
L = chol(params(1)*(Phi*Phi') + Sigma);
trpart = norm(L\Y,'fro')^2/n;
ldpart = 2*sum(log(diag(L)));
val = ldpart+trpart;
end


 
