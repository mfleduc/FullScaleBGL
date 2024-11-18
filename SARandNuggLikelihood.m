function val = SARandNuggLikelihood(B,tausq,Y)
% Fits a covariance model of the form (B^TB)^-1+tau^2 I
% where B is a spatial autoregrression matrix. See the LatticeKrig paper
% (Nychka et al 2015, a Gaussian Process Model for Large Spatial Datasets)
% for info
n = size(Y,2);
ldB = logdet(B);
L = chol(speye(size(Y,1))+tausq*(B'*B),'lower');
ldItB = sum(log(diag(L)));
trSQi = 1/n*norm(B*Y,'fro')^2;
trbig = tausq/n*norm(L\((B'*B)*Y),'fro')^2;

val = -2*ldB+ldItB+trSQi - trbig;

end