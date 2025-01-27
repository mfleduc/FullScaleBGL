function [paramvec,ps] = FitLatticeKrig(data,PTP,Phi,mesh)
%% Fits a LatticeKrig model to the given data
%% We will do cross-validation outside this function.
lvls = length(mesh);
npts=501;
fn = @(x)LatticeKrigLikelihood(x,mesh,Phi,PTP,data);
optparams=struct;
optparams.algorithm = 'simanneal';
% nparams = 3+length(mesh);
ranges = 10.^linspace( -2,2 , npts);
for ii = 1:length(mesh)
    ranges = [ranges;10.^linspace( -2,2,npts )];
end
ranges = [ranges; 10.^linspace(-3.9,0.5,npts) ;10.^linspace(-10,1,npts)];
% x0 = [250, 250, 250, 250, 250, 250];
x0 = [];
[~,ps,paramvec] = SimulatedAnnealing(fn,ranges,5000,0.9,1e-5,x0);



end

