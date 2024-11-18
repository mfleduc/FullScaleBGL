function out = FitMOW(x, d, Y ,stprof)
%% normalize data to zero mean, unit variance
n = size(Y,2);
if length(x)==5
Sigma = MixOfWendlandCov(x, d,'cov',stprof);
elseif length(x)==4
params = [x, 1-x(3)];
Sigma = MixOfWendlandCov(params, d,'cov',stprof);
end
L = chol(Sigma,'lower');
ldpart = 2*sum(log(diag(L)));
LiY = L\Y;
trPart =  1/n*norm(LiY,'fro')^2;
out = ldpart+trPart;
end