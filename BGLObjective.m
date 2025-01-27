function val = BGLObjective( logtausq,Q,Phi,Y, varargin)
% Objective function for determining noise parameters in the basis
% graphical lasso noise model.
%   more detail to come
if ~isempty(varargin)
    precisionflag = varargin{1};
else
    precisionflag=0;
end
val=0;
tausq = 10^logtausq;
n = size(Y,2) ;
% del = 1+tol ;
ldsigma = log(tausq)*n;

 LiY = Y/sqrt(tausq) ; 
 LiPhi = Phi/sqrt(tausq) ;
 trS = norm(LiY,'fro')^2*1/n;

QpLiPhisq = Q+LiPhi'*LiPhi ; 
LQ = chol(QpLiPhisq,'lower');
ldQ = 2*sum(log(diag(LQ))) ;
YSiPhi = (LiY)'*(LiPhi);
bigmatrix = 1/n*norm(LQ\YSiPhi','fro')^2;
trbig = sum(diag(bigmatrix));
val = ldsigma+ldQ+trS-trbig;
end
