function [Q,B] = MarkovRandomField(x,d,n,stdprof)
%% Note: This outputs Sigma^-1, so make sure you use it with the correct likelihood function!!!
% keyboard
x(3)=10^x(3);
p = size(d,1);
B = speye(p)*(1+x(2)^2)+d/n;
sdevs =spdiags(1./stdprof,0,p,p);
B = sqrt(x(3))*B*sdevs;
Q = B'*B;
% Q=1;
% Q = Q-Q*((speye(size(d,1))/x(1)+Q)\Q);
end
% fnToMin = @(x)NoiseObjective(x, covmodel, Q ,Ahat, yhat,1) ; 