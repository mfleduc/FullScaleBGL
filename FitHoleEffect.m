function out = FitHoleEffect(x, d, Y, stdprof)
%%
Y = Y./stdprof;
n = size(Y,2);
Sigma = gentemp([x,-0.1], d);
L = chol(Sigma,'lower');
ldpart = 2*sum(log(diag(L)));
LiY = L\Y;
trPart =  1/n*norm(LiY,'fro')^2;
out = ldpart+trPart;
end