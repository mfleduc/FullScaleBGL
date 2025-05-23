function fullScaleBGLMOW(input)
fprintf('Running Full-Scale Basis Graphical Lasso with %d levels of resolution and lambdas %.4f to %.4f\n'...
    ,input.lvls,input.lambda(1),input.lambda(end))
fprintf('*******MMIXTURE OF WENDLANDS********\n');
%% Parameters for determining residual model parameters
% Not needed for the BGL
optparams = struct; 
optparams.algorithm='simanneal';
optparams.ranges = [[-3,3];%log10(tau^2)
                    [1, 20];%theta1
                    [0.01,600];%alpha1
                    [1,20];%theta2
                    [0.01,600]]; %alpha2
optparams.maxiters = 5000; %number of simulated annealing steps
optparams.npts = 501; %Number of lattice points
optparams.x0 = [0, 10, 300, 10, 300];%Initial guess  7.65, 301.205,
optparams.pchange = 0.9;
%% Distance measures
fprintf('********************Fitting residual model**********************\n');
[spx,spy,spz] = sph2cart(pi/180*input.LON(:), pi/180*input.LAT(:),1);
spc=[spx,spy,spz];
d = real(acosd(spc*spc')) ; %d = d - diag(diag(d));
%imaginary values on the diagonal sometimes        
covmodel = @(x)MixOfWendlandCov(x,d,'cov');      
%% 
fnname = func2str(covmodel);
fnname = fnname(5:end);
ndx = strfind(fnname, '(');
fnname = fnname(1:ndx-1);
% separation = rand([1,64800]);
% ndcs = separation < 0.25;dtmp = d(ndcs,ndcs);
[val,obj,np] = CalculateInitialGuess(covmodel,input.Psi,input.data, optparams);
fprintf('Done\n');
for lambda = input.lambda       
    alpha = np(1);
    [Qest] = CalculateQ(input.data, input.Psi, lambda, input.tol,covmodel,np(2:end), 1 );
    save(fullfile(input.svDir,sprintf('FSBGL results %d levels lambda %.4f covmodel %s.mat',input.lvls,lambda , fnname)),'alpha', 'Qest', 'np', 'optparams', 'covmodel','lambda','-v7.3')
end
end

