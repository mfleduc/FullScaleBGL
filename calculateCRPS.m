%% Parameter cross validation when necessaryGCMglade
% function crossValidation()
clear variables;close all;
rng(65161660)
cd('/glade/work/mleduc/GP emulation/needlet fbgl');
addpath( genpath('/glade/work/mleduc/Ionosphere_inversion/NeedMat/') );
addpath( genpath('/glade/work/mleduc/Ionosphere_inversion/Spherical-Harmonic-Transform/') );
addpath( genpath('/glade/work/mleduc/Ionosphere_inversion/spherepts/') );
addpath( genpath('/glade/work/mleduc/BigQUIC_release/'))
addpath( genpath('/glade/work/mleduc/MultivariateBasisGraphicalLasso'))
addpath( genpath('/glade/work/mleduc/FMGL-0/'))
addpath( genpath('/glade/work/mleduc/GP emulation/needlet fbgl/NeedletBGL'));
addpath( genpath('/glade/work/mleduc/S2-Sampling-Toolbox'));
%

load(sprintf('needlets degree resolution j=5.mat'), 'A');
A = [ones(size(A,1),1)/2/sqrt(pi),A];
ndcs = randperm(size(A,1),8000);
load('mean of beta algorithm residual transform alt final.mat' , ...
    'testdata','finalfit','LAT', 'LON' ,'coeffs','a00','c00');
data = (testdata(:,1:25)-mean(testdata(:,1:25),2))./std(testdata(:,1:25),[],2);
data = data(ndcs,:);

Ahat = A(ndcs,:);
[XH,YH] = sph2hammer(pi/180*(LON(ndcs)-180),pi/180*LAT(ndcs));
%
[spx,spy,spz] = sph2cart(pi/180*LON(ndcs), pi/180*LAT(ndcs),1);
spc=[spx,spy,spz];
d = 2*sind(0.5*real(acosd(spc*spc') ));
aicresults = load('CV results mix of GC 3 levels nugg first.mat');
for kk = 1:length(aicresults.data)
    Qtmp = aicresults.data{kk}.Qest{aicresults.ndcs(kk)};
    noiseparams =  aicresults.data{kk}.np( aicresults.ndcs(kk),: );
    Sigma = gentemp( [noiseparams(1:2),-0.1],d ) + gentemp( [noiseparams(3:4),-0.1],d );
    lowrankcov = Ahat(:,1:253)*(Qtmp\Ahat(:,1:253)');
    pwstdevs = sqrt(diag(lowrankcov)+diag(Sigma));
    CRPS(kk) = 
end