%% Parameter cross validation when necessaryGCMglade
function crossValidation()
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
A = [ones(size(A,1),1)/2/sqrt(pi),A];
load('mean of beta algorithm residual transform alt final.mat' , ...
    'testdata','finalfit','LAT', 'LON' ,'coeffs','a00','c00');

% files = {'FSBGL results 3 levels time 1 lambda 1.0000 covmodel 3 MixOfWendlandCov.mat',...
%     'FSBGL results 3 levels time 1 lambda 0.5000 covmodel 3 MixOfWendlandCov.mat'};
% files = {'FSBGL results 4 levels time 1 lambda 0.1000 covmodel 3 MixOfWendlandCov.mat',...
%     'FSBGL results 4 levels time 1 lambda 0.2500 covmodel 3 MixOfWendlandCov.mat',...
%     'FSBGL results 4 levels time 1 lambda 0.5000 covmodel 3 MixOfWendlandCov.mat',...
%     'FSBGL results 4 levels time 1 lambda 1.0000 covmodel 3 MixOfWendlandCov.mat'} ;
files = dir('FSBGL Nugget first/1Wendland/FSBGL results 3 levels time 1 lambda *.mat');
% files = dir('FSBGL Nugget first/FSBGL results 3 levels time 1 lambda * covmodel WendlandCov nugg first.mat');
files = fullfile({files(1:end).folder},{files(1:end).name});
likelihoods = [];
data=cell(1,length(files));
% data = load(files{1}, 'np','ndcs','nptseach');
for kk = 1:length(files)
    data{kk} = load(files{kk},'Qest','lambda','np','ndcs','nptseach');
    testset = data{kk}.ndcs(5*data{kk}.nptseach+1:end);
    LATd = LAT(testset);LONd = LON(testset);
    datapts = (testdata(:,1:25) - mean(testdata(:,1:25),2))./std(testdata(:,1:25),[],2) ;
    datapts = datapts(testset,:);
    [spx,spy,spz] = sph2cart(pi/180*LONd, pi/180*LATd,1);
    spc=[spx,spy,spz];
    d = real(acosd(spc*spc') );% 
    % d = 2*sind(d/2) ;
    % covmodel = @(x)gentemp([x(1:2),-0.1],d); %+gentemp([x(3:4),-0.1],d)  ;
    Ahat = A(testset,1:size(data{kk}.Qest{1},1));
    if length(data{kk}.np(1,:))==4
        npndcs = 2:4;
    else
        npndcs = 1:3;
    end
    covmodel = @(x)WendlandCov(x,d,'cov') ;
    for jj = 1:length(data{kk}.Qest)

        likelihoods(kk,jj) = NoiseObjective(data{kk}.np(jj,npndcs), covmodel, ...
            data{kk}.Qest{jj},Ahat , datapts);
        likelihoods(kk,jj) = likelihoods(kk,jj) - logdet(data{kk}.Qest{jj});
    end
    lambdas(kk) = data{kk}.lambda;
end
[~,ndcs] = min(likelihoods,[],2);

for kk = 1:length(data)
    nll(kk) = likelihoods(kk,ndcs(kk));
    ks(kk) = nnz(triu(data{kk}.Qest{ndcs(kk)})) + length(data{kk}.np(ndcs(kk),2:end)) ;
end
aic = (2*ks+2*nll);

save('FSBGL Nugget first/1Wendland/CV results 1Wendland 3 levels.mat', 'ndcs', 'aic','data','likelihoods','ks','lambdas');
% end