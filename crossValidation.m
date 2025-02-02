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
% files = dir('FSBGL Nugget first/3Wendlands/FSBGL results 3 levels time 1 lambda *.mat');
% files = dir('FSBGL Nugget first/GC/FSBGL results 3 levels time 1 lambda *.mat');
files = dir('FSBGL Nugget first/PGL/PGL results 3 levels * MixOfWendland.mat');
files = fullfile({files(1:end).folder},{files(1:end).name});
pgl = true;
likelihoods = [];ks=[];
data=cell(1,length(files));
if ~pgl
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
        % covmodel = @(x)gentemp([x(1:2),-0.1],d)+gentemp([x(3:4),-0.1],d)  ;
        Ahat = A(testset,1:size(data{kk}.Qest{1},1));
        % if length(data{kk}.np(1,:))==4
            npndcs = 2:8;
        % else
        %     npndcs = 1:3;
        % % end
        % npndcs = 1:4;
        covmodel = @(x)MixOfWendlandCov(x,d,'cov') ;
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
        ks(kk) = nnz(triu(data{kk}.Qest{ndcs(kk)})) + length(data{kk}.np(ndcs(kk),npndcs)) ;
    end
    aic = (2*ks+2*nll);
    
    save('FSBGL Nugget first/Mix GC/CV results Mix GC 3 levels.mat', 'ndcs', 'aic','data','likelihoods','ks','lambdas');
    % end
else
    lambdas = [];
    % likelihoods()
    testset = randperm(64800,24800);
    Ahat = A(testset,1:253);
    LATd = LAT(testset);LONd = LON(testset);
    datapts = (testdata(:,1:25) - mean(testdata(:,1:25),2))./std(testdata(:,1:25),[],2) ;
    datapts = datapts(testset,:);
    [spx,spy,spz] = sph2cart(pi/180*LONd, pi/180*LATd,1);
    spc=[spx,spy,spz];
    d = real(acosd(spc*spc') );% 
    covmodel = @(x)MixOfWendlandCov(x,d,'cov');
    for kk = 1:length(files)
        data{kk} = load(files{kk});
        % keyboard
        lambdas(end+1:end+length(data{kk}.lambdas)) = data{kk}.lambdas;
        for jj = 1:length(data{kk}.lambdas)
            likelihoods = [likelihoods;zeros(1,5)];
            ks(end+1) = nnz(triu(data{kk}.Qs{jj})) + size(data{kk}.np,3);
            for ii = 1:5
                likelihoods(end,ii) = NoiseObjective(data{kk}.np(jj,ii,:), covmodel, ...
                data{kk}.Qs{jj},Ahat , datapts);
                likelihoods(end,ii) = likelihoods(end,ii) - logdet(data{kk}.Qs{jj});
            end
        end
    end
    nll = likelihoods(:,1);
    datapts = (testdata(:,1:25) - mean(testdata(:,1:25),2))./std(testdata(:,1:25),[],2) ;
    for ii = 2:5
        ndcs = randperm(64800,24800);
        % testset = randperm(64800,24800);
        Ahat = A(ndcs,1:253);
        LATd = LAT(ndcs);LONd = LON(ndcs);
        thisdata = datapts(ndcs,:);
        [spx,spy,spz] = sph2cart(pi/180*LONd, pi/180*LATd,1);
        spc=[spx,spy,spz];
        d = real(acosd(spc*spc') );% 
        covmodel = @(x)MixOfWendlandCov(x,d,'cov');
        for kk = 1:length(files)
            nll(ii,kk) = NoiseObjective(data{kk}.np(1,1,:), covmodel, ...
                data{kk}.Qs{1},Ahat , thisdata);
            nll(ii,kk) = nll(ii,kk)-logdet(data{kk}.Qs{1});
        end
    end
end