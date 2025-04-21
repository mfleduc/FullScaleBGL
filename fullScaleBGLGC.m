function fullScaleBGLGC(input)
%% Full scale basis graphical lasso
% clear variables; close all;
fprintf('Running Full-Scale Basis Graphical Lasso with %d levels of resolution and lambdas %.4f to %.4f\n'...
    ,input.lvls,input.lambda(1),input.lambda(end))
fprintf('*******GASPARI-COHN COVARIANCE********\n');
%% Example optimization parameters
optparams = struct; 
optparams.algorithm='simanneal';
optparams.ranges = [[-3,3];[0.001,1500];[.0001,0.5]];
optparams.maxiters = 5000; %number of simulated annealing steps
optparams.npts = 501; %Number of lattice points
optparams.x0 = [0,750,0.25];%Initial guess
optparams.pchange = 0.9;
%%
fprintf('********************Fitting residual model**********************\n');

[spx,spy,spz] = sph2cart(pi/180*input.LON(:), pi/180*input.LAT(:),1);
spc=[spx,spy,spz];
d = real(acosd(spc*spc')) ; %d = d - diag(diag(d));
%imaginary values on the diagonal sometimes        
d = 2*sind(d/2); %chordal distance required for gaspari cohn correlation fn
covmodel = @(x)gentemp([x(2:3),-0.1],d,x(1));%+gentemp([x(3:4),-0.1],d);
% covmodel = @(x)gentemp([x(1:2),-0.1], 2*sind(d/2)) + WendlandCov(x(3:end), d, 'cov');
[val,obj,np ] = CalculateInitialGuess(covmodel{1}, A ,testdata, optparams);
fprintf('Done\n');
alpha=np(1);
fnname = func2str(covmodel{1});
fnname = fnname(5:end);
ndx = strfind(fnname, '('); 
fnname = fnname(1:ndx(1)-1);
disp(toc)
for l1 = 1:length(lambdas) %% cross-validate for each lambda
    lambda = lambdas(l1);
    Qest = cell(1,1);
    % np = [];
    fprintf('lambda = %.6f, calculating Q|Sigma\n',lambda);
    %% We want to fit a Wendland covariance
    
    [Qest{1}] = CalculateQ(testdata, A, lambda, ...
        2e-2,covmodel{1},np, strcmpi(fnname,'MarkovRandomField') );
    likelihood = NoiseObjective( np(2:end), covmodel{1},...
        Qest{1},A,testdata, strcmpi(fnname,'MarkovRandomField'));
    svdir = 'FSBGL results/temperature/GC/';
    mkdir(svdir);
     save(sprintf([svdir, ...
         'FSBGL results %d levels lambda %.4f covmodel %s.mat'],lvls,lambda,fnname ),...
    'Qest', 'np', 'optparams', 'covmodel','lambda','alpha','likelihood','-v7.3')
     disp(toc)
end

end

