function val = NoiseObjectiveFirst( params,covmodel,Phi,Y, varargin)
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
l = size(Phi,2);
% del = 1+tol ;
Sigma = covmodel(params(2:end));
L = chol(Sigma , 'lower') ; 
if ~precisionflag
    ldsigma = 2*sum(log(diag(L))); 
    LiY = L\Y ; 
    LiPhi = L\Phi ;
    trS = norm(LiY,'fro')^2*1/n;
elseif precisionflag
    ldsigma = -2*sum(log(diag(L)));
    LiY = L'*Y;
    LiPhi = L'*Phi;
    trS = norm(LiY,'fro')^2*1/n;
end
QpLiPhisq = params(1)*eye(l)+LiPhi'*LiPhi ; 
LQ = chol(QpLiPhisq,'lower');
ldQ = 2*sum(log(diag(LQ))) ;
YSiPhi = (LiY)'*(LiPhi);
bigmatrix = 1/n*norm(LQ\YSiPhi','fro')^2;
trbig = sum(diag(bigmatrix));
val = ldsigma+ldQ+trS-trbig - l*log(params(1));
end
