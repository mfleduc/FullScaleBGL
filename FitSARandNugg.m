function val = FitSARandNugg( x, d, n, Y , stdprof)
% Fitting parameterized SAR and nugget covariance model
%
[~,B] = MarkovRandomField(x,d,n,1, stdprof);
val = SARandNuggLikelihood(B,x(1), Y);
end