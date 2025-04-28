%% Full scale basis graphical lasso
function fullScaleBGLWendland(lambdas,lvls)
% clear variables; close all;
fprintf('Running Full-Scale Basis Graphical Lasso with %d levels of resolution and lambdas %.4f to %.4f\n'...
    ,lvls,lambdas(1),lambdas(end))
fprintf('*****WENDLAND COVARIANCE*****\n');
%% 
optparams = struct; 
optparams.algorithm='simanneal';
optparams.ranges = [[-10,1];%log10(tau^2)
                    [1, 20];%theta1
                    [0.01,3]];%alpha2
% optparams.ranges = [-10,1];
optparams.maxiters = 500; %number of simulated annealing steps
optparams.npts = 501; %Number of lattice points
optparams.x0 = [-4,10,0.5];%Initial guess
optparams.pchange = 0.9;

fprintf('********************Fitting residual model**********************\n');

[spx,spy,spz] = sph2cart(pi/180*LON(:), pi/180*LAT(:),1);
spc=[spx,spy,spz];
d = real(acosd(spc*spc')) ; %d = d - diag(diag(d));
%imaginary values on the diagonal sometimes        
covmodel = @(x)WendlandCov(x,d,'cov');
fnname = func2str(covmodel{1});
fnname = fnname(5:end);
ndx = strfind(fnname, '('); 
fnname = fnname(1:ndx(1)-1);
[val,obj,np] = CalculateInitialGuess(covmodel,input.Psi,input.data, optparams);
fprintf('Done\n');

for lambda = input.lambda       
    alpha = np(1);
    [Qest] = CalculateQ(input.data, input.Psi, lambda, input.tol,covmodel,np(2:end), 1 );
    save(fullfile(input.svDir,sprintf('FSBGL results %d levels lambda %.4f covmodel %s.mat',input.lvls,lambda , fnname)),'alpha', 'Qest', 'np', 'optparams', 'covmodel','lambda','-v7.3')
end
end
 

