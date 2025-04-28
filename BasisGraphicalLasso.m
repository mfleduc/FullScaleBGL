%% Full scale basis graphical lasso
% clear variables; close all;
function BasisGraphicalLasso(input, varargin)
fprintf('Running Basis Graphical Lasso with %d levels of resolution and lambdas %.4f to %.4f\n'...
    ,input.lvls,input.lambda(1),input.lambda(end))
% nneeds = [12,60,252,1020]+1;
% rng(6516156)
% cd('/glade/work/mleduc/GP emulation/needlet fbgl');
% addpath( genpath('/glade/work/mleduc/Ionosphere_inversion/NeedMat/') );
% addpath( genpath('/glade/work/mleduc/Ionosphere_inversion/Spherical-Harmonic-Transform/') );
% addpath( genpath('/glade/work/mleduc/Ionosphere_inversion/spherepts/') );
% addpath( genpath('/glade/work/mleduc/BigQUIC_release/'))
% addpath( genpath('/glade/work/mleduc/MultivariateBasisGraphicalLasso'))
% addpath( genpath('/glade/work/mleduc/FMGL-0/'))
% addpath( genpath('/glade/work/mleduc/GP emulation/needlet fbgl/NeedletBGL'));
% addpath( genpath('/glade/work/mleduc/S2-Sampling-Toolbox'));
% 
% load(sprintf('needlets degree resolution j=5.mat'), 'A');
% load('mean of beta algorithm residual transform alt final.mat'  , ...
%      'testdata' );
% lambdas = [0.7,0.4,0.3];
%% Split into time steps
% timendx = 1;
% testdata = reshape(testdata, [],25,6);
% data = testdata - mean(testdata, 2);
% data = data(:,:,1)./std(data(:,:,1),[],2);
p = size(input.Psi,1);
covmodel = @(x)NuggetEffectCovariance(x,p);
%% 
 % 
[alpha,np] = FitNugget(input.Psi, input.data);
for lambda = input.lambda
    fprintf('*************Lambda = %.6f************\n', lambda)
    [Qest] = CalculateQ(input.data, input.Psi, input.lambda, input.tol,covmodel,np, 0 );
    save(fullfile(input.svDir,sprintf('BGL results %d levels lambda %.4f.mat',input.lvls,lambda )),'Qest', 'np', 'alpha', 'input','covmodel','lambda','-v7.3')
end
end