function out = FitWendlandCorr(x, d, Y, stdprof)
%% 
n = size(Y,2);
Sigma = WendlandCov(x, d, 'corr', stdprof);
L = chol(Sigma,'lower');
ldpart = 2*sum(log(diag(L)));
LiY = L\Y;
trPart =  1/n*norm(LiY,'fro')^2;
out = ldpart+trPart;
end