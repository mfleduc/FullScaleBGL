function ResidualReanalysis()
%% Checking whether the residual models change significantly if you recalculate the parameters knowing Q
% From Mitch: They don't change. I am doubtful.
cd('/glade/work/mleduc/GP emulation/needlet fbgl');
addpath('/glade/work/mleduc/GP emulation/needlet fbgl/NeedletBGL');
files = dir('FSBGL Nugget first/3Wendlands/FSBGL results 3 levels*.mat');
files = fullfile(files(1).folder,{files(1:end).name});
cvres = load('FSBGL Nugget first/3Wendlands/CV results 3Wendlands 3 levels.mat');

load(sprintf('needlets degree resolution j=5.mat'), 'A');
A = [ones(64800,1)/2/sqrt(pi), A]; 
load('mean of beta algorithm residual transform alt final.mat'  ...
    ,'testdata','LAT', 'LON' );
dataset = (testdata(:,1:25) - mean(testdata(:,1:25),2))./std(testdata(:,1:25),[],2);

%load each file
for ff = 1:length(files)
    prevresults = load(files{ff},'Qest','ndcs','optparams','np','covmodel','nptseach','ndcs','lambda');
    lambdandx = find(prevresults.lambda == cvres.lambdas);
    Q = prevresults.Qest{cvres.ndcs(lambdandx)};
    np = prevresults.np(lambdandx,:);prevresults.optparams.x0 = np(2:end);
    trainset = cvres.ndcs(lambdandx)*prevresults.nptseach+1:(cvres.ndcs(lambdandx)+1)*prevresults.nptseach ;
    traindata = dataset(trainset,:);
    trainbasis = A(trainset,1:size(Q,1));
    covmodel = prevresults.covmodel{cvres.ndcs(lambdandx)};
    updatednp = CalculateSigma( covmodel, Q, trainbasis, traindata, prevresults.optparams, 0) ;

end


end