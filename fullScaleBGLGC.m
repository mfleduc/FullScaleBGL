function fullScaleBGLGC(lambdas,lvls)
%% Full scale basis graphical lasso
% clear variables; close all;
fprintf('Running Full-Scale Basis Graphical Lasso with %d levels of resolution and lambdas %.4f to %.4f\n'...
    ,lvls,lambdas(1),lambdas(end))
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
load('mean of beta algorithm residual transform alt final.mat' , ...
     'testdata','LAT', 'LON' );
% lambdas = [1,0.5,0.25,0.1,0.05,0.03,0.025];
%% Split into time steps
timendx = 1;
testdata = reshape(testdata, [],25,6); 
data = testdata - mean(testdata, 2);
data = data(:,:,timendx)./std(data(:,:,timendx),[],2);
%% 
optparams = struct; 
optparams.algorithm='simanneal';
% optparams.ranges = [[0.001,1];... %alpha
%     [0.001,1];...  % 0.5*chordal decorrelation distance (c parameter)
%     [-9,-1];... %log10(nugget)
%     [0.001,20];...% Wendland decorrelation distance]
%     [0.001,1]]; %alpha for Wendland
% optparams.ranges = [-10,1];
optparams.ranges = [[0.001,1];[.0001,0.5]];
optparams.maxiters = 500; %number of simulated annealing steps
optparams.npts = 501; %Number of lattice points
optparams.x0 = [0.25,0.07];%Initial guess
optparams.pchange = 0.9;
%%
A = [ones(64800,1)/2/sqrt(pi), A]; 
ncvs = 5;
nptseach = 8000;
np = nan(ncvs, size(optparams.ranges,1)+1);
ndcs = randperm(size(A,1));
covmodel = cell(1,ncvs);
for nn = 1:ncvs
        fprintf('********************Fitting residual model, Training Set %d **********************\n',nn);
        thesendcs = ndcs((nn-1)*nptseach+1:nn*nptseach);
        Ahat = A(thesendcs,1:nneeds(lvls));
        thisdata = data(thesendcs,:);
        % [~,tausq] = FitNugget(Ahat,thisdata);
        LATd=LAT(thesendcs);LONd=LON(thesendcs);
        [spx,spy,spz] = sph2cart(pi/180*LONd, pi/180*LATd,1);
        spc=[spx,spy,spz];
        d = real(acosd(spc*spc')) ; %d = d - diag(diag(d));
 %imaginary values on the diagonal sometimes        
        d = 2*sind(d/2); %chordal distance required for gaspari cohn correlation fn
        covmodel{nn} = @(x)gentemp([x(1:2),-0.1],d);%+gentemp([x(3:4),-0.1],d);
        % covmodel = @(x)gentemp([x(1:2),-0.1], 2*sind(d/2)) + WendlandCov(x(3:end), d, 'cov');
        
        fprintf('Fitting residual parameters...');
        [val,obj,np(nn,:)] = CalculateInitialGuess(covmodel{nn}, Ahat,thisdata, optparams);
        fprintf('Done\n');
end
fnname = func2str(covmodel{1});
fnname = fnname(5:end);
ndx = strfind(fnname, '('); 
fnname = fnname(1:ndx(1)-1);
for l1 = 1:length(lambdas) %% cross-validate for each lambda
    lambda = lambdas(l1);
    Qest = cell(ncvs,1);
    % np = [];
    fprintf('lambda = %.6f, calculating Q|Sigma\n');
    %% We want to fit a Wendland covariance
    for nn = 1:ncvs
        fprintf('***************Training set %d***************\n');
        thesendcs = ndcs((nn-1)*nptseach+1:nn*nptseach);
        Ahat = A(thesendcs,1:nneeds(lvls));    
        thisdata = data(thesendcs,:);
        [Qest{nn}] = CalculateQ(thisdata, Ahat, lambda, ...
            1e-2,covmodel{nn},np(nn,2:end), strcmpi(fnname,'MarkovRandomField') );
    end
    
         save(sprintf('FSBGL Nugget first/GC/lower tolerance/FSBGL results %d levels time %d lambda %.4f covmodel gentemp.mat',lvls,timendx,lambda ),...
        'Qest', 'np', 'optparams', 'timendx','covmodel','lambda','ndcs','nptseach','-v7.3')
end

end

