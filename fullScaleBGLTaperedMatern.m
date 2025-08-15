%% Full scale basis graphical lasso
function fullScaleBGLTaperedMatern(input)
fprintf('Running Full-Scale Basis Graphical Lasso with %d levels of resolution and lambdas %.4f to %.4f\n'...
    ,input.lvls,input.lambda(1),input.lambda(end))
%% Parameters of high res component
fprintf('*******TAPERED MATERN********\n');
optparams = struct; 
optparams.algorithm='simanneal';
optparams.ranges = [[-3,3];%log10(tau^2)
                    [0.001,10];%sigma^2
                       [0.001,2];%Taper range
                       [0.1,1];%Smoothness
                       [0.001,4]];%Matern range
% optparams.ranges = [-10,1];
optparams.maxiters = 5000; %number of simulated annealing steps
optparams.npts = 501; %Number of lattice points
optparams.x0 = [0, 1, 1, 0.5, 1];%Initial guess
optparams.pchange = 0.95;
fprintf('********************Fitting residual model**********************\n');

% [spx,spy,spz] = sph2cart(pi/180*input.LON(:), pi/180*input.LAT(:),1);
% spc=[spx,spy,spz];
% d = real(acosd(spc*spc')) ; %d = d - diag(diag(d));
%imaginary values on the diagonal sometimes 
d = sqrt(abs(input.LAT(:)-input.LAT(:)').^2+abs(input.LON(:)-input.LON(:)').^2);
covmodel = @(x)TaperedMatern( d,x);
fnname = func2str(covmodel);
fnname = fnname(5:end);
ndx = strfind(fnname, '('); 
fnname = fnname(1:ndx(1)-1);
[val,obj,np] = CalculateInitialGuess(covmodel,input.Psi,input.data, optparams);
% np = [-1.108,-0.816,1032,27.55,0.6346,7.5214];
fprintf('Done\n');

for lambda = input.lambda       
    alpha = np(1);
    [Qest] = CalculateQ(input.data, input.Psi, lambda, input.tol,covmodel,np(2:end), 1 );
    save(fullfile(input.svDir,sprintf('FSBGL results %d levels lambda %.4f covmodel %s.mat',input.lvls,lambda , fnname)),'input','alpha', 'Qest', 'np', 'optparams', 'covmodel','lambda','-v7.3')
end
end
 

