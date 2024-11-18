%% Full scale basis graphical lasso
% clear variables; close all;
function fullScaleBGLMRF(lambdas, lvl)
fprintf('Running Full-ScaleBasis Graphical Lasso with %d levels of resolution and lambdas %.4f to %.4f\n'...
    ,lvl,lambdas(1),lambdas(end))
fprintf('*******MARKOV RANDOM FIELD********');
nneeds = [12,60,252,1020]+1;
rng(6516156)
cd('/glade/work/mleduc/GP emulation/needlet fbgl');
addpath( genpath('/glade/work/mleduc/Ionosphere_inversion/NeedMat/') );
addpath( genpath('/glade/work/mleduc/Ionosphere_inversion/Spherical-Harmonic-Transform/') );
addpath( genpath('/glade/work/mleduc/Ionosphere_inversion/spherepts/') );
addpath( genpath('/glade/work/mleduc/BigQUIC_release/'))
addpath( genpath('/glade/work/mleduc/MultivariateBasisGraphicalLasso'))
addpath( genpath('/glade/work/mleduc/FMGL-0/'))
addpath( genpath('/glade/work/mleduc/GP emulation/needlet fbgl/NeedletBGL'));
addpath( genpath('/glade/work/mleduc/S2-Sampling-Toolbox'));

load(sprintf('needlets degree resolution j=5.mat'), 'A');
load('mean of beta algorithm residual transform alt final.mat' ,'coeffs', ...
    'a00','c00', 'testdata','finalfit','LAT', 'LON' );
% lambdas = [0.7,0.4,0.3];
%% Split into time steps
timendx = 1;
testdata = reshape(testdata, [],25,6);
data = testdata - mean(testdata, 2);
data = data(:,:,1)./std(data(:,:,1),[],2);
vals=[];rows=[];cols=[];
for ii = 1:length(LAT)
    ndcs = find((abs(LAT-LAT(ii))==1 & (LON-LON(ii))==0 )|(abs(LON-LON(ii))==1 & (LAT-LAT(ii))==0 )...
        |(abs(LON-LON(ii))==359 & (LAT-LAT(ii))==0 ));
    nums = min(length(ndcs),4);
    vals(end+1:end+nums)=-1;
    rows(end+1:end+nums)=ii;
    cols(end+1:end+nums) = ndcs(1:nums);
end
d = sparse(rows,cols,vals,64800,64800);
covmodel = @(x)MarkovRandomField([0,x],d,4,1);
p = size(A,1);
% covmodel = @(x)NuggetEffectCovariance(x,p);
fnname = func2str(covmodel);
fnname = fnname(5:end);
ndx = strfind(fnname, '(');
fnname = fnname(1:ndx-1);
% separation = rand([1,64800]);
% ndcs = separation < 0.25;dtmp = d(ndcs,ndcs);

%% 
optparams = struct; 
optparams.algorithm= 'simanneal';
optparams.ranges = [[0.001,5];%kappa
                    [0.001, 3]];%1/alpha
% optparams.ranges = [-10,1];
optparams.maxiters = 1000; %number of simulated annealing steps
optparams.npts = 501; %Number of lattice points
optparams.x0 = [0.5,1.5];%Initial guess
optparams.pchange = 0.9;

A = [ones(64800,1)/2/sqrt(pi), A]; Ahat = A(:,1:nneeds(lvl));
for lambda = lambdas
    
    [val,obj,np] = CalculateInitialGuess(covmodel, Ahat,data, optparams);
    alpha = np(1);
    [Qest] = CalculateQ(data, Ahat, lambda, 1e-2,covmodel,np(2:end), 1 );
    save(sprintf('FSBGL Nugget first/MRF/FSBGL results %d levels time %d lambda %.4f covmodel %s.mat',lvl,timendx,lambda , fnname),'alpha', 'Qest', 'np', 'optparams', 'timendx','covmodel','lambda','-v7.3')
end
end