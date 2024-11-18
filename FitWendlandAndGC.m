function [out] = FitWendlandAndGC(x,d)
%FITWENDLANDANDGC Summary of this function goes here
Y = Y./stdprof;
n = size(Y,2);
Sigma = gentemp([x(1:2),-0.1], 2*sind(d/2)) + WendlandCov(x(3:end), d, 'cov');
L = chol(Sigma,'lower');
ldpart = 2*sum(log(diag(L)));
LiY = L\Y;
trPart =  1/n*norm(LiY,'fro')^2;
out = ldpart+trPart;
end

