function Sigma = NuggetEffectCovariance(logtausq,npts)
%% Nugget effect covariance model
%logtausq = log10(tau^2), tau^2 the nugget variance
Sigma = 10^(logtausq(1))*speye(npts);
end