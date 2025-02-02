%% Calculating the AIC for model determination
function CalculateAIC(files,cvresults,varargin)
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
normdata = (testdata(:,1:25,1) - mean(testdata(:,1:25,1),2))./std(testdata(:,1:25,1),[],2);
if ~isempty(varargin)
    pglflag = varargin{1};
else
    pglflag = 0;
end

%% First part: Testing FS-BGL 
% files = {'FSBGL results 4 levels time 1 lambda 0.5000 covmodel gentemp.mat',...
%     'FSBGL results 4 levels time 1 lambda 0.2500 covmodel gentemp.mat',...
%     'FSBGL results 4 levels time 1 lambda 1.0000 covmodel gentemp.mat',...
%     'FSBGL results 4 levels time 1 lambda 0.1000 covmodel gentemp.mat'};
%     'FSBGL results 3 levels time 1 lambda 0.0500 covmodel MarkovRandomField.mat',...
%     'FSBGL results 3 levels time 1 lambda 0.0300 covmodel MarkovRandomField.mat',...
%     'FSBGL results 3 levels time 1 lambda 0.0250 covmodel MarkovRandomField.mat'};
% 
% files = {'FSBGL Nugget results/FSBGL results 4 levels time 1 lambda 1.0000 covmodel NuggetEffectCovariance.mat';...
%     'FSBGL Nugget results/FSBGL results 4 levels time 1 lambda 0.5000 covmodel NuggetEffectCovariance.mat';...
%     'FSBGL Nugget results/FSBGL results 4 levels time 1 lambda 0.7000 covmodel NuggetEffectCovariance.mat'};
%     'FSBGL Nugget results/FSBGL results 3 levels time 1 lambda 0.5000 covmodel NuggetEffectCovariance.mat'; ...
%     'FSBGL Nugget results/FSBGL results 3 levels time 1 lambda 0.0250 covmodel NuggetEffectCovariance.mat'; ...
%     'FSBGL Nugget results/FSBGL results 3 levels time 1 lambda 0.0500 covmodel NuggetEffectCovariance.mat'};
% files = dir('FSBGL Nugget first/FSBGL results 3 levels time 1 lambda * covmodel mix gentemp*.mat');

% files = dir('/glade/work/mleduc/GP emulation/needlet fbgl/FSBGL results/nugget cov/FSBGL results 3*.mat');
% files = dir('FSBGL results 4 levels time 1 lambda *covmodel 3 MixOfWendlandCov.mat');
% files = fullfile({files(1:end).folder},{files(1:end).name});
if ~pglflag
    lambdas = zeros(1,length(files));
    Qs = cell(1,length(files));
    noiseparams = [];
    % covmodel = @(x)NuggetEffectCovariance(x,64800);
    nloglike = zeros(1,length(files))';
    ks = zeros(1,length(files))';
    if exist('cvresults','var')==1 && ~isempty(cvresults)
        cvres = load( cvresults );
        for ii = 1:length(files)
            results = load(files{ii});
            whichone = cvres.ndcs(ii);
            lambdas(ii) = results.lambda;
            for jj = 1:10
                % ndcs = 1:64800 ;
                ndcs = randperm(64800,8000);
                LATd = LAT(ndcs);LONd = LON(ndcs);
                [spx,spy,spz] = sph2cart(pi/180*LONd, pi/180*LATd,1);
                spc=[spx,spy,spz];
                d = real(acosd(spc*spc') );
                covmodel = eval(func2str(results.covmodel));
                ks(ii,jj) = size(results.Qest{whichone},1) + 0.5*(nnz(results.Qest{whichone})...
                    -size(results.Qest{whichone},1))...
                    +length(results.np(whichone,:));
                nloglike(ii,jj) = NoiseObjective(results.np(whichone,:), covmodel,results.Qest{whichone},...
                    A(ndcs,1:size(results.Qest{whichone},1)), normdata(ndcs,:),0);
                nloglike(ii,jj) = nloglike(ii,jj) - logdet(results.Qest{whichone});
            end
        end
    
    else   
        for ii = 1:length(files)
            results = load(files{ii});
            % covmodel = @(x)NuggetEffectCovariance(x,24800);
            lambdas(ii) = results.lambda;
            
            for jj=1 
                % ndcs = randperm(64800,24800) ;
                ndcs = 1:64800;
                nloglike(ii,jj) = NoiseObjective(results.np(2:end), results.covmodel,results.Qest,...
                    A(ndcs,1:size(results.Qest,1)), normdata(ndcs,:),0,ndcs);
                nloglike(ii,jj) = nloglike(ii,jj) - logdet(results.Qest);
                ks(ii,jj) = nnz(triu(results.Qest))+length(results.np);
            end
        end
    end
elseif pglflag
    keyboard
    lambdas = [];
    for ii = 1:length(files)
        results = load(files{ii});
        
    end
end
AIC = 2*(ks+nloglike);
figure;semilogx(lambdas, mean(AIC,2), 'k.');grid on;
xlabel('\lambda')
ylabel('AIC')
% title('AIC, FS-BGL with mixture of Gaspari-Cohn, 3 levels')
[folder,~,~] = fileparts(files{1});
save(fullfile(folder,'AIC results.mat'),'nloglike','ks','AIC','lambdas' );
%% Testing PGL reults
% files = {'3 levels of resolution mean of beta algorithm alt quic time step 1.mat'};
% data = load(files{1});
% nll = zeros(1,length(data.lambda));
% ks = zeros(1,length(data.lambda));
% for ii = 1:length(data.lambda)
%     ks(ii) = nnz(data.Qest{ii})/2 + size(data.Qest{ii},1)/2;
%     L = chol( data.Qest{ii} );
%     ldpart = -2*sum(log(diag(L)));
%     thesecoeffs = [a00(:,1)'+c00(:,1)';data.coeffs(:,:,1)];
%     trpart = 1/size(thesecoeffs,2)*norm( L\(thesecoeffs-mean(thesecoeffs,2)),'fro' );
%     nll(ii) = ldpart+trpart;
% end
% AIC = 2*(ks+nll);
% figure;
% semilogx(data.lambda, AIC,'k.')
end