function FullScaleApproximations(pctvar)
%% For performing full-scale approximations in the style of Stein(2008) and Sang et al (2012)
% Write the process as a sum of a KL expansion and a residual process
% pctvar = fraction of variance explained by low rank process
cd('/glade/work/mleduc/GP emulation/needlet fbgl');
addpath( genpath('/glade/work/mleduc/Ionosphere_inversion/NeedMat/') );
addpath( genpath('/glade/work/mleduc/Ionosphere_inversion/Spherical-Harmonic-Transform/') );
addpath( genpath('/glade/work/mleduc/Ionosphere_inversion/spherepts/') );
addpath( genpath('/glade/work/mleduc/BigQUIC_release/'))
addpath( genpath('/glade/work/mleduc/MultivariateBasisGraphicalLasso'))
addpath( genpath('/glade/work/mleduc/FMGL-0/'))
addpath( genpath('/glade/work/mleduc/GP emulation/needlet fbgl/NeedletBGL'));
addpath( genpath('/glade/work/mleduc/S2-Sampling-Toolbox'));

load('mean of beta algorithm residual transform alt final.mat' , ...
     'testdata','LAT', 'LON' );


% data = testdata ;
data = (testdata(:,1:25)-mean(testdata(:,1:25),2))./std(testdata(:,1:25),[],2); %Setting to marginal standard normal
m = size(data,2) ;
%% Perform PCA, project onto selected basis, determine residuals
[U,S] = svd(data(:,1:25)/sqrt(m),'econ');
totalvar = sum(diag(S));
cumulativeVar = cumsum(diag(S))/totalvar;
ndx = find(cumulativeVar > pctvar);
Phi = U(:,1:ndx(1)) ;
Shat = S(1:ndx(1),1:ndx(1));
% coeffs = data'*Phi;
% residual = data - Phi*coeffs';
%% Now fit a model to the residuals
%optimization parameters
optparams = struct; 
optparams.algorithm='simanneal';
optparams.ranges = [[-10,1];%log10(tau^2)
                    [1, 20];%theta1
                    [0.01,1];%alpha1
                    [1,20];%theta2
                    [0.01,1];%alpha2
                    [1,20];%theta3
                    [0.01,1]];%alpha3
% optparams.ranges = [-10,1];
optparams.maxiters = 500; %number of simulated annealing steps
optparams.npts = 501; %Number of lattice points
optparams.x0 = [-4,13,0.15,7,0.15,3,0.15];%I1itial guess
optparams.pchange = 0.9;
%
ncvs = 5;
nptseach = 8000;
ndcs = randperm( size(data,1) );
covmodel = {};
np = [];
for nn = 1:ncvs
    trainset = ndcs((nn-1)*nptseach+1:nn*nptseach);
    LATd = LAT(trainset);
    LONd = LON(trainset);
    [spx,spy,spz] = sph2cart(pi/180*LONd, pi/180*LATd,1);
    spc=[spx,spy,spz];
    d = real(acosd(spc*spc'));
    covmodel{nn} = @(x)MixOfWendlandCov(x,d,'cov');
    % np(nn,:) = CalculateSigma(covmodel{nn},inv(Shat.^2),Phi(trainset,:),data(trainset,1:25),optparams,0);
    
end
save(sprintf('FSBGL Nugget first/Full scale/Full scale results 3Wendlands pctvar=%.4f.mat',pctvar),'pctvar','optparams','covmodel','np','ndcs','Phi','Shat')


end