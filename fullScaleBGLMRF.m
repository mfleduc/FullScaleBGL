%% Full scale basis graphical lasso
% clear variables; close all;
function fullScaleBGLMRF(input, optparams)
fprintf('Running Full-Scale Basis Graphical Lasso with %d levels of resolution and lambdas %.4f to %.4f\n'...
    ,input.lvls,input.lambda(1),input.lambda(end))
fprintf('*******MARKOV RANDOM FIELD********\n');
%% Parameters for determining residual model p[arameters
% Not needed for the BGL
optparams = struct; 
optparams.algorithm= 'simanneal';
optparams.ranges = [...
                    [0.0001,10];%kappa
                    [-8,1]];%1/alpha for MRF
% optparams.ranges = [-10,1];
optparams.maxiters = 1000; %number of simulated annealing steps
optparams.npts = 501; %Number of lattice points
optparams.x0 = [4, -4];%Initial guess
optparams.pchange = 0.9;
%% Split into time steps
fprintf('********************Fitting residual model**********************\n');
LAT=input.LAT;LON=input.LON;
vals=[];rows=[];cols=[];
for ii = 1:length(LAT)
    ndcs = find((abs(LAT-LAT(ii))==1 & (LON-LON(ii))==0 )|(abs(LON-LON(ii))==1 & (LAT-LAT(ii))==0 )...
        |(abs(LON-LON(ii))==359 & (LAT-LAT(ii))==0 ));
    nums = min(length(ndcs),4);
    vals(end+1:end+nums)=-1;
    rows(end+1:end+nums)=ii;
    cols(end+1:end+nums) = ndcs(1:nums);
end
d = sparse(rows,cols,vals,length(LAT),length(LAT));
covmodel = @(x)MarkovRandomField([0,x],d,4,1);
p = size(input.Psi,1);
% covmodel = @(x)NuggetEffectCovariance(x,p);
fnname = func2str(covmodel);
fnname = fnname(5:end);
ndx = strfind(fnname, '(');
fnname = fnname(1:ndx-1);
% separation = rand([1,64800]);
% ndcs = separation < 0.25;dtmp = d(ndcs,ndcs);

%% 
[val,obj,np] = CalculateInitialGuess(covmodel,input.Psi,input.data, optparams);
fprintf('Done\n');
for lambda = input.lambda       
    alpha = np(1);
    [Qest] = CalculateQ(input.data, input.Psi, lambda, input.tol,covmodel,np(2:end), 1 );
    save(fullfile(input.svDir,sprintf('FSBGL results %d levels lambda %.4f covmodel %s.mat',input.lvls,lambda , fnname)),'alpha', 'Qest', 'np', 'optparams', 'covmodel','lambda','-v7.3')
end
end