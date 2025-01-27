function ProjectedGraphicalLasso(lambdas,lvls,varargin)
%% Runs the projected graphical lasso using the needlet transform
%Inputs:    
%   data: struct with fields
%         Y: dataset, sptial location x realizations
%         LON: longitude of datapoints, radians
%         LAT: latitude of datapoints, 0 = north, radians


%% Input parser
rng(34563)
cd('/glade/work/mleduc/GP emulation/needlet fbgl');
addpath('/glade/work/mleduc/GP emulation/needlet fbgl/NeedletBGL');
validresidmodels = {'nugget', 'MRF','Wendland','MixOfWendland'};
ip = inputParser();
    % addRequired(ip, 'data',@(x)(isstruct(x)));%&&isfield(data,'Y')&&isfield(data,'LAT')&&isfield(data,'LON')));
    addRequired(ip, 'lambdas',@(lambdas)(isnumeric(lambdas)&&all(isreal(lambdas))&&(min(lambdas)>0) ));
    addRequired(ip, 'lvls',@(lvls)(isscalar(lvls)&&isreal(lvls)&&(lvls==floor(lvls))));
    addOptional(ip, 'residualmodel','MixOfWendland',@(x)any(validatestring(x,validresidmodels)));
    addOptional(ip, 'crossvalidate',false,@(x)(x==true ||x==false))
ip.parse(lambdas,lvls,varargin{:});
residualmodel = ip.Results.residualmodel;
%% Adding paths to necessary matlab code
fprintf('Running the Projected Graphical Lasso with %d levels of resolution and lambda = %.6f to %.6f \n',lvls,lambdas(1),lambdas(end));
addpath( genpath('/glade/work/mleduc/Ionosphere_inversion/NeedMat/') );
addpath( genpath('/glade/work/mleduc/Ionosphere_inversion/Spherical-Harmonic-Transform/') );
addpath( genpath('/glade/work/mleduc/Ionosphere_inversion/spherepts/') );
addpath( genpath('/glade/work/mleduc/MultivariateBasisGraphicalLasso')); %QUIC is in here
%% 

load('pgl data.mat','data');
nrealizations = size(data.Y,2);
nneeds = 1+[12,60,252,1020];
%% Calculate the transform
try
    load('weights 1 deg resolution.mat');
    alm = spharmonic_tran_irr(data.LAT,data.LON, data.Y,2^(lvls),weights);
catch excep
    fprintf([excep.message,'\n']);
    alm = spharmonic_tran_irr(data.LAT,data.LON, data.Y,2^(lvls));
end

a00 = squeeze(alm(1,2^(lvls)+1,:));
beta{1} = spneedlet_tran(alm(:,:,1), 2^(lvls), 2) ;
b=[];
for ii = 1:lvls
    b = [b;beta{1}{ii}];
end
for jj = 2:size(alm,3)
    bp=[];
    beta{jj} = spneedlet_tran(alm(:,:,jj), 2^(lvls), 2) ;
    for ii = 1:(lvls)
        bp = [bp;beta{jj}{ii}];
    end
    b = [b,bp];
end
b = [a00';b];%Adding the first shpericla harmonic coefficients
%The needlets are not a frame without Y_00 = sqrt(1/4pi)
% mnb = mean(b,2);
% bdemeaned = b-mnb;
% stdevb = std(b,[],2);
% bnorm = bdemeaned./stdevb;
Qs = cell( 1,length(lambdas) );
Sb = b*b'/size(b,2);
for ll = 1:length(lambdas)
    penalty = lambdas(ll)*(ones(size(Sb)) - eye(size(Sb,1)));
    Qs{ll} = QUIC(Sb,penalty, 1e-6);
end
%% Now: Fit the residual model
try
    load(sprintf('needlets degree resolution j=5.mat'), 'A');
catch
    A = get_A(2, 0, lvls-1, data.LAT, data.LON, 1e3);
end
A = [ones(size(A,1),1)/2/sqrt(pi),A];
fittedFld = A(:,1:nneeds(lvls))*b;
residuals = data.Y - fittedFld;
switch residualmodel
    case  'nugget'
        npts = size(data.Y,1);
        likelihood = @(x)(npts*log(x.^2)+1./(nrealizations*x.^2)*norm(residuals,'fro')^2) ;
        np = fminbnd(likelihood,0,100);
    case 'MixOfWendland'
        np=[];
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
        optparams.x0 = [-4,12,0.1,5,0.3,2,0.4];%I1itial guess
        optparams.pchange = 0.9;
        ncvs=5;nptseach=8000;
        covmodel = {};
        ndcs = randperm(size(data.Y,1));
        for nn = 1:ncvs
            fprintf('***************Training set %d***************\n',nn);
            thesendcs = ndcs((nn-1)*nptseach+1:nn*nptseach);
            Ahat = A(thesendcs,1:nneeds(lvls));    
            thisdata = data.Y(thesendcs,:);
            
        % [~,tausq] = FitNugget(Ahat,thisdata);
            LATd=data.LAT(thesendcs);LONd=data.LON(thesendcs);
            [spx,spy,spz] = sph2cart(LONd,pi/2-LATd,1);
            spc=[spx,spy,spz];
            d = real(acosd(spc*spc')) ; %d = d - diag(di
            covmodel{nn} = @(x)MixOfWendlandCov(x,d,'cov');  
            for ll = 1:length(lambdas)
                fprintf('Lambda = %.4f \n',lambdas(ll));
                [np(ll,nn,:)] = CalculateSigma(covmodel{nn},Qs{ll}, Ahat,data.Y(thesendcs,:) ,optparams,0);
                             
            end
        end
    otherwise
        error('Not implemented')
end
svname = sprintf('FSBGL Nugget first/PGL/PGL results %d levels lambda %.6f to %.6f residuals %s.mat',lvls,lambdas(1),lambdas(end), residualmodel);
save(svname, 'lambdas', 'data', 'Qs', 'residualmodel','np','alm','b');
end