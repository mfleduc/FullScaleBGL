function CompareToDifferentDataset( )
%COMPARETODIFFERENTDATASET Summary of this function goes here
%   Detailed explanation goes here
rng(654016)
load('mean of beta algorithm residual transform alt final.mat' , 'LAT','LON','testdata');
load(sprintf('needlets degree resolution j=5.mat'), 'A');
tmp=load('/glade/work/mleduc/GP emulation/needlet fbgl/FSBGL Nugget first/nugget cov/BGL results 3 levels lambda 0.0500.mat');
% load('/glade/work/mleduc/GP emulation/needlet fbgl/FSBGL Nugget first/Full scale/Full scale results 3Wendlands pctvar=0.9000.mat')
% tmp = load('FSBGL Nugget first/Mix GC/FSBGL results 3 levels time 1 lambda 0.5000 covmodel mix gentemp.mat') ;
% tmp = load('FSBGL Nugget first/PGL/PGL results 3 levels lambda 0.005000 to 0.005000 residuals MixOfWendland.mat') ;
Phi = [ones(size(A,1),1)/2/sqrt(pi),A(:,1:252)];
% load('FSBGL Nugget first/LatticeKrig/LatticeKrig results basis_.mat');
% load('FSBGL Nugget first/LatticeKrig/LatticeKrig results_Q.mat')
Q = tmp.Qest;
% Q = tmp.Qs{1};
% Q = inv(Shat);
timesteps = 1:5:125;
outsidelikelihood = [];
for tt = 1:length(timesteps)
    cvdata = testdata(:,timesteps(tt):24+timesteps(tt))...
        - mean(testdata(:,timesteps(tt):24+timesteps(tt)),2);
    ncvdata = cvdata./std(cvdata,[],2);
    
    ncvs = 1;
    nptseach = 24800;
    testndcs = randperm( size(cvdata,1) ,ncvs*nptseach);
    testndcs = reshape(testndcs,ncvs,nptseach);
    % load(sprintf('needlets degree resolution j=5.mat'), 'A');
    % load('/glade/work/mleduc/GP emulation/needlet fbgl/FSBGL Nugget first/Full scale/Full scale results 3Wendlands pctvar=0.9000.mat');
    
    for nn = 1:ncvs
        LATd = LAT(testndcs(nn,:));LONd = LON(testndcs(nn,:));
        [spx,spy,spz] = sph2cart(pi/180*LONd, pi/180*LATd,1);
        spc=[spx,spy,spz];
        d = real(acosd(spc*spc'));
        %d = 2*sind(d/2);
        %covmodel = @(x)gentemp([x(1:2),-0.1],d)+gentemp([x(3:4),-0.1],d)  ;
        % covmodel = @(x)MixOfWendlandCov(x,d,'cov');
        covmodel = @(x)NuggetEffectCovariance(x,24800);
        outsidelikelihood(tt,nn) = NoiseObjective( (tmp.np),covmodel,Q ,Phi(testndcs(nn,:),:),ncvdata(testndcs(nn,:),:),0) ;
        outsidelikelihood(tt,nn) = outsidelikelihood(tt,nn) -logdet(Q ) ;
    end
end
% save('FSBGL Nugget first/Mix GC/CV results Mix GC 3 levels.mat','outsidelikelihood','-append')
end

