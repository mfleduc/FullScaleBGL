%% Full scale basis graphical lasso
% clear variables; close all;
function fullScaleBGL(lambdas,lvls)
% clear variables; close all;
fprintf('Running Full-Scale Basis Graphical Lasso with %d levels of resolution and lambdas %.4f to %.4f\n'...
    ,lvls,lambdas(1),lambdas(end))
fprintf('MRF covariance model\n');
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
lambdas = [1,0.5,0.25,0.1,0.05,0.03,0.025];
%% Split into time steps
timendx = 1;
testdata = reshape(testdata, [],25,6);
data = testdata - mean(testdata, 2);
data = data(:,:,1)./std(data(:,:,1),[],2);
%% We want to fit a Markov Random Field to the residuals
% 4 nearest neighbors. Need to maybe think about this?
p = size(testdata, 1);
d = speye(size(testdata,1));
for ii = 1:p-1
    these = ii+1:p;
    latdiff = abs(LAT(ii)-LAT(ii+1:p));
    londiff = abs(LON(ii)-LON(ii+1:p));
    ndcs = find((latdiff ==1 & londiff==0)|(latdiff==0 &(londiff==1|londiff==359)));
    d(ii,these(ndcs)) = -1;
    d(these(ndcs),ii) = -1;
end
d = d - speye(p);
n = 4;
covmodel = @(x)MarkovRandomField([0,x],d,4,1);
% covmodel = @(x)NuggetEffectCovariance(x,p);
fnname = func2str(covmodel);
fnname = fnname(5:end);
ndx = strfind(fnname, '('); 
fnname = fnname(1:ndx-1);
% separation = rand([1,64800]);
% ndcs = separation < 0.25;dtmp = d(ndcs,ndcs);

%% 
optparams = struct; 
optparams.algorithm='simanneal';
optparams.ranges = [[0.001,5];%kappa
                    [0.001, 3]];%1/alpha
% optparams.ranges = [-10,1];
optparams.maxiters = 500; %number of simulated annealing steps
optparams.npts = 501; %Number of lattice points
optparams.x0 = [0.25,1];%Initial guess
optparams.pchange = 0.9;
tausq = 0.0051;
A = [ones(64800,1)/2/sqrt(pi), A]; Ahat = A(:,1:nneeds(lvls));
for lambda = lambdas
    [Qest, np] = CalculateQAndSigma(data, Ahat, lambda, 1e-2,covmodel,optparams, strcmpi(fnname,'MarkovRandomField') ,0.0051);
    save(sprintf('FSBGL Nugget first/FSBGL results %d levels time %d lambda %.4f covmodel %s nugg first.mat',lvls,timendx,lambda , fnname), 'Qest', 'np', 'optparams', 'timendx','covmodel','lambda','tausq','-v7.3')
end
end