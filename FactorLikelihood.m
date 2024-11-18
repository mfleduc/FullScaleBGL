function val = FactorLikelihood(B, Y, precisionflag, tausq)
% Calculates likelihood that Y ~ N(0, B^TB) or
% equivalently, if precisionflag = 1, N(0, (B^TB)^{-1})
ldpart = 2*logdet(B);
n = size(Y,2);
if precisionflag
    ldpart = -ldpart;
    trS = 1/n*norm( B*Y,'fro' )^2;
else
    trS = 1/n*norm(B'\Y,'fro')^2;
end
val = ldpart+trS;
end