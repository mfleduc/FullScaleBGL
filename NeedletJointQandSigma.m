clear variables; close all ;
%% Needlet Basis Graphical Lasso
cd('/glade/work/mleduc/GP emulation/needlet fbgl');
addpath( genpath('/glade/work/mleduc/Ionosphere_inversion/NeedMat/') );
addpath( genpath('/glade/work/mleduc/Ionosphere_inversion/Spherical-Harmonic-Transform/') );
addpath( genpath('/glade/work/mleduc/Ionosphere_inversion/spherepts/') );
addpath( genpath('/glade/work/mleduc/BigQUIC_release/'));
addpath( genpath('/glade/work/mleduc/MultivariateBasisGraphicalLasso'));
addpath( genpath('/glade/work/mleduc/GP emulation/needlet fbgl/NeedletBGL'));
% 
load('needlets degree resolution j=5.mat', 'A') ; 
load('mean of beta algorithm residual transform alt final.mat', 'coeffs', ...
    'a00','c00', 'testdata','finalfit','LAT', 'LON', 'lambda' );
%% The data we feed in is the test data - mean model, which in our case is given by 
% A*\hat{\beta} where \hat{\beta} is the sample mean of the coefficients. 
% We also have to fit a profile to the standard deviations. We can sort
% that out. 
%% Parameters for optimization
optparams = struct; 
optparams.algorithm='simanneal';
optparams.ranges = [[-6,-3];... %log10(tausq)
                    [0.5,10];... %angular decorr distance
                    [0.0001, 1]];%fraction of variance explained by the noise model
optparams.maxiters = 200; %number of simulated annealing steps
optparams.npts = 501; %Number of lattice points
optparams.x0 = [1e-4,5,0.8];%Initial guess
optparams.pchange = 0.9;
%%
ndcs = find(rand([1,64800])<0.03);
Ahat = A(ndcs, 1:252);
% noisepart = testdata - finalfit; % Use this to inform the standard deviation profiles. 
% stdevprofile = std(noisepart, [],2);
testdata = reshape(testdata, size(testdata,1), size(coeffs, 2), size(coeffs, 3));
Qest = {};
noiseparams = [];
d = zeros(length(ndcs));
for ii = 1:length(ndcs)
    d(ii:length(ndcs),ii) = distance([LAT(ndcs(ii)),LON(ndcs(ii))],...
        [LAT(ndcs(ii:length(ndcs))), LON(ndcs(ii:length(ndcs)))]);
end
d = d+d';
for ii = 1:size(testdata, 3)
    Y = testdata(ndcs,:,ii) - Ahat(:,1:252)*mean(coeffs(1:252,:,ii),2) - mean(a00(:,ii)+c00(:,ii))/2/sqrt(pi);
    stdevprofile = std(Y,[],2) ;
    Y = (Y-mean(Y,2))  ;
    covmodel = @(x)WendlandCov(x, d, 'cov',stdevprofile);
    % covmodel = 'wendland';
    fprintf('Lambda = %.4f \n', 5);
    fprintf('Time step %d\n',ii);
    % keyard
    [Qest{ii}, noiseparams(:,ii)] = CalculateQAndSigma(Y, Ahat(:,1:252), 15, 0.05,covmodel, optparams) ;
    
end